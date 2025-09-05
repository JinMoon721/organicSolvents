#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <limits>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <cctype>
#include <unordered_map>

// ------------------------ Basic data structs ------------------------
struct Atom {
  float x, y, z;
};

struct Box {
  double Lx{0}, Ly{0}, Lz{0}; // orthorhombic lengths
};

struct Frame {
  std::vector<Atom> atoms;
  Box box;
};

struct Angle { 
  float zloc;  // z location of terminal carbon
  float ccCos; // angle from +z axis of terminal cc bond
  float chCos; // angle from +z axis of the mean terminal ch bonds
};

struct Targets {
  int t; // terminal 
  int c; // carbon next to terminal
  int h1; // hydrogen1 at terminal
  int h2;
  int h3;
  Targets operator+(int k) const { return {t+k, c+k, h1+k, h2+k, h3+k}; }
  Targets& operator+=(int k) { t+=k; c+=k; h1+=k; h2+=k; h3+=k; return *this; }
};

// ------------------------ Molecule Target Selection -------------
enum class Key { DES, DBS, DEO, DEA, Unknown };
Key toKey(std::string s) {
  for (char& ch : s ) ch = std::toupper(static_cast<unsigned char>(ch));
  if (s == "DES") return Key::DES;
  if (s == "DBS") return Key::DBS;
  if (s == "DEO") return Key::DEO;
  if (s == "DEA") return Key::DEA;
  return Key::Unknown;
}

Targets findTargets(const std::string& s, bool left) {
  Targets init{};
  switch ( toKey(s)) {
    case Key::DES:
      init = left ? Targets{1, 2, 19, 20, 21} : Targets{18, 17, 44, 43, 42}; break;
    case Key::DBS:
      init = left ? Targets{1, 2, 23, 24, 25} : Targets{22, 21, 56, 55, 54}; break;
    case Key::DEO:
      init = left ? Targets{1, 2, 11, 12, 13} : Targets{10, 9, 20, 19, 18}; break;
    case Key::DEA:
      init = left ? Targets{1, 2, 15, 16, 17} : Targets{14, 13, 32, 31, 30}; break;
    default:
      throw std::runtime_error("Unknown key");
  }

  init += -1;
  return init;
}
// ----------------------- START READ LAMMPS INPUT ---------------------
// One row from the Atoms section (fields you asked for)
struct AtomRow {
  int index = 0;   // LAMMPS atom id (1-based)
  int mol   = 0;   // molecule id (0 if missing)
  int type  = 0;   // atom type
  double charge = 0.0; // charge (0 if missing)
};

// Parsed result container
struct LammpsData {
  std::vector<AtomRow> atoms;                 // all Atoms rows we read
  std::unordered_map<int, double> masses;     // type -> mass (from Masses)
  std::string atomsStyleHint;                 // e.g., "full", "atomic", "charge", if seen
};

// --------------------- Small utilities ---------------------

// Trim leading/trailing whitespace in-place.
static inline void trimInplace(std::string& s) {
  size_t i = 0;
  while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
  size_t j = s.size();
  while (j > i && std::isspace(static_cast<unsigned char>(s[j-1]))) --j;
  s.assign(s, i, j - i);
}

// Remove everything after a '#' (LAMMPS inline comment).
static inline void stripComment(std::string& s) {
  size_t pos = s.find('#');
  if (pos != std::string::npos) s.erase(pos);
}

// Case-insensitive prefix check.
static inline bool startsWithCi(const std::string& s, const std::string& prefix) {
  if (s.size() < prefix.size()) return false;
  for (size_t i = 0; i < prefix.size(); ++i) {
    if (std::tolower(static_cast<unsigned char>(s[i])) !=
        std::tolower(static_cast<unsigned char>(prefix[i]))) return false;
  }
  return true;
}

// Strict integer parse (entire token must be an int).
static inline bool parseIntStrict(const std::string& tok, int& out) {
  const char* p = tok.c_str();
  char* end = nullptr;
  long v = std::strtol(p, &end, 10);
  if (end == p || *end != '\0') return false;
  out = static_cast<int>(v);
  return true;
}

// Strict double parse (entire token must be a number).
static inline bool parseDoubleStrict(const std::string& tok, double& out) {
  const char* p = tok.c_str();
  char* end = nullptr;
  double v = std::strtod(p, &end);
  if (end == p || *end != '\0') return false;
  out = v;
  return true;
}

// Tokenize a whitespace string into vector<string>.
static inline std::vector<std::string> splitWs(const std::string& s) {
  std::vector<std::string> out;
  std::istringstream iss(s);
  std::string w;
  while (iss >> w) out.push_back(w);
  return out;
}

// Heuristic: does a cleaned line look like a known section header (besides Atoms/Masses)?
static inline bool isOtherSectionHeader(const std::string& s) {
  static const char* names[] = {
    "Velocities",
    "Bonds", "Angles", "Dihedrals", "Impropers",
    "Pair Coeffs", "Bond Coeffs", "Angle Coeffs",
    "Dihedral Coeffs", "Improper Coeffs", "PairIJ Coeffs",
    "Atoms", "Masses"  // we will handle these explicitly elsewhere
  };
  for (const char* n : names) {
    if (startsWithCi(s, n)) return true;
  }
  // Fallback: alphabetic line that isn't data (be conservative)
  if (!s.empty() && std::isalpha(static_cast<unsigned char>(s[0]))) return true;
  return false;
}

// Try to extract a style hint from a header line like "Atoms # full"
static inline std::string extractStyleHint(const std::string& headerLine) {
  size_t hash = headerLine.find('#');
  if (hash == std::string::npos) return std::string();
  std::string tail = headerLine.substr(hash + 1);
  trimInplace(tail);
  // Tail might be like: "full" or "charge" or "molecular"
  // Strip trailing words like "atoms" if present (rare)
  return tail;
}

// --------------------- Atoms parsing ---------------------

// Parse one Atoms data line to (index, mol, type, charge) using style rules.
// - If style is known ("full", "molecular", "charge", "atomic"), follow that.
// - If style is empty, infer from token pattern; missing fields default to 0 / 0.0.
static inline bool parseAtomsRow(const std::vector<std::string>& tok,
                                 const std::string& styleHint,
                                 AtomRow& out)
{
  // Defensive defaults
  AtomRow row;

  // Known styles first
  if (!styleHint.empty()) {
    // Atoms # full : id mol type q x y z ...
    if (startsWithCi(styleHint, "full")) {
      if (tok.size() < 4) return false;
      if (!parseIntStrict(tok[0], row.index)) return false;
      if (!parseIntStrict(tok[1], row.mol))   return false;
      if (!parseIntStrict(tok[2], row.type))  return false;
      if (!parseDoubleStrict(tok[3], row.charge)) return false;
      out = row; return true;
    }
    // Atoms # molecular : id mol type x y z ...
    if (startsWithCi(styleHint, "molecular")) {
      if (tok.size() < 3) return false;
      if (!parseIntStrict(tok[0], row.index)) return false;
      if (!parseIntStrict(tok[1], row.mol))   return false;
      if (!parseIntStrict(tok[2], row.type))  return false;
      row.charge = 0.0; out = row; return true;
    }
    // Atoms # charge : id type q x y z ...
    if (startsWithCi(styleHint, "charge")) {
      if (tok.size() < 3) return false;
      if (!parseIntStrict(tok[0], row.index)) return false;
      if (!parseIntStrict(tok[1], row.type))  return false;
      if (!parseDoubleStrict(tok[2], row.charge)) return false;
      row.mol = 0; out = row; return true;
    }
    // Atoms # atomic : id type x y z ...
    if (startsWithCi(styleHint, "atomic")) {
      if (tok.size() < 2) return false;
      if (!parseIntStrict(tok[0], row.index)) return false;
      if (!parseIntStrict(tok[1], row.type))  return false;
      row.mol = 0; row.charge = 0.0; out = row; return true;
    }
    // Unknown hint → fall through to inference
  }

  // No hint (or unrecognized): infer:
  // Try "full"-like first: id, mol, type, charge, ...
  bool ok = false;
  if (tok.size() >= 4) {
    AtomRow t{};
    ok = parseIntStrict(tok[0], t.index)
      && parseIntStrict(tok[1], t.mol)
      && parseIntStrict(tok[2], t.type)
      && parseDoubleStrict(tok[3], t.charge);
    if (ok) { out = t; return true; }
  }
  // Try "charge"-like: id, type, charge, ...
  if (tok.size() >= 3) {
    AtomRow t{};
    ok = parseIntStrict(tok[0], t.index)
      && parseIntStrict(tok[1], t.type)
      && parseDoubleStrict(tok[2], t.charge);
    if (ok) { t.mol = 0; out = t; return true; }
  }
  // Try "atomic"-like: id, type, ...
  if (tok.size() >= 2) {
    AtomRow t{};
    ok = parseIntStrict(tok[0], t.index)
      && parseIntStrict(tok[1], t.type);
    if (ok) { t.mol = 0; t.charge = 0.0; out = t; return true; }
  }

  return false; // couldn’t parse
}

// --------------------- File reader ---------------------

// Parse a LAMMPS data file, collecting Atoms (id,mol,type,charge) and Masses (type->mass).
LammpsData readLammpsData(const std::string& path) {
  std::ifstream fin(path);
  if (!fin) throw std::runtime_error("Cannot open file: " + path);

  LammpsData out;
  enum State { None, InAtoms, InMasses };
  State state = None;

  std::string line;
  while (std::getline(fin, line)) {
    // Keep a copy for header detection; but strip comments for data parsing.
    std::string raw = line;
    stripComment(line);
    trimInplace(line);
    if (line.empty()) continue;

    // Section enters
    if (startsWithCi(raw, "Atoms")) {
      state = InAtoms;
      out.atomsStyleHint = extractStyleHint(raw); // may be empty
      continue;
    }
    if (startsWithCi(raw, "Masses")) {
      state = InMasses;
      continue;
    }

    // If a new (other) section begins, leave current section.
    if (isOtherSectionHeader(raw)) {
      state = None;
      continue;
    }

    // Parse content by section
    if (state == InAtoms) {
      // Tokenize current Atoms data row
      std::vector<std::string> tok = splitWs(line);
      if (tok.empty()) continue;

      AtomRow row;
      if (parseAtomsRow(tok, out.atomsStyleHint, row)) {
        out.atoms.push_back(row);
      } else {
        // Silent skip or log a warning:
        // std::cerr << "Warn: could not parse Atoms line: " << raw << "\n";
      }
      continue;
    }

    if (state == InMasses) {
      // Expect: "type mass" (possibly more tokens we ignore)
      std::vector<std::string> tok = splitWs(line);
      if (tok.size() < 2) continue;
      int type = 0; double mass = 0.0;
      if (parseIntStrict(tok[0], type) && parseDoubleStrict(tok[1], mass)) {
        out.masses[type] = mass; // last one wins if duplicated
      }
      continue;
    }

    // Otherwise (state == None): outside sections → ignore
  }

  return out;
}


double massFromType(const std::vector<std::pair<int, double>>& type2mass, int type, double def = -1.0) {
  auto it = std::lower_bound(type2mass.begin(), type2mass.end(), type, [](const auto& p, int key){return p.first < key; });
  return (it != type2mass.end() && it->first == type) ? it->second : def;
}

float wrapPBC(float x, float L) { // wrap x into [0, L)
  float y = std::fmod(x, L);
  if (y<0) y+=L;
  return y;
}

float comZ(const std::vector<Atom>& atoms, const std::vector<int>& atomTypes, const std::vector<std::pair<int, double>>& type2mass){
  if (atoms.size() != atomTypes.size() ) {
    throw std::runtime_error("atoms/types size mismatch");
  }
  long double M = 0.0L, sz = 0.0L;
  for (size_t i=0; i<atoms.size(); i++) {
    const double m = massFromType(type2mass, atomTypes[i]);
    sz += m*atoms[i].z;
    M += m;
  }
  return (M == 0.0L) ? 0.0 : float(sz / M);
}

float comZPBC(const std::vector<Atom>& atoms, const std::vector<int>& atomTypes, const std::vector<std::pair<int, double>>& type2mass, float Lz){
  if (atoms.size() != atomTypes.size() ) {
    throw std::runtime_error("atoms/types size mismatch");
  }
  if(!(Lz > 0.0) ) {
    throw std::runtime_error("Lz must be positive");
  }
  long double C = 0.0L, S = 0.0L, M=0.0L;
  const long double twoPi = 2.0L * acosl(-1.0L);
  const long double k = twoPi / Lz;

  for (size_t i =0; i< atoms.size(); i++) {
    const double m = massFromType(type2mass, atomTypes[i]);
    const long double theta = k* wrapPBC(atoms[i].z, Lz);
    C += m * std::cos(theta);
    S += m * std::sin(theta);
    M += m;
  }
  if (M == 0.0L) return 0.0;

  long double ang = std::atan2(S, C);
  if (ang < 0) ang += twoPi;
  return float (ang * (Lz / twoPi));
}


void writeLammpsDump(std::ostream& out, long timestep, const std::vector<Atom>& atoms, const std::vector<int> atomTypes, const Box& box) {
  out << "ITEM: TIMESTEP\n" << timestep << "\n";
  out << "ITEM: NUMBER OF ATOMS\n" << atoms.size() << "\n";
  out << "ITEM: BOX BOUNDS pp pp pp\n";
  out << std::setprecision(10) << std::fixed;
  float bl=0;
  out << bl << " " << box.Lx << "\n";
  out << bl << " " << box.Ly << "\n";
  out << bl << " " << box.Lz << "\n";

  out << "ITEM: ATOMS id type x y z\n";
  out << std::setprecision(10) << std::fixed;
  for( int atom=0; atom< atoms.size(); atom++ ) {
    out << atom+1 << " " << atomTypes[atom] << " " << atoms[atom].x << " " <<  atoms[atom].y << " " << atoms[atom].z << "\n"; 
  }
}




// ------------------------ START READ DSD FILE -------------------------------
static inline uint32_t bswap32(uint32_t x) {
  return ((x & 0xFF000000u) >> 24) |
         ((x & 0x00FF0000u) >> 8)  |
         ((x & 0x0000FF00u) << 8)  |
         ((x & 0x000000FFu) << 24);
}
static inline uint64_t bswap64(uint64_t x) {
  return ((x & 0xFF00000000000000ull) >> 56) |
         ((x & 0x00FF000000000000ull) >> 40) |
         ((x & 0x0000FF0000000000ull) >> 24) |
         ((x & 0x000000FF00000000ull) >> 8)  |
         ((x & 0x00000000FF000000ull) << 8)  |
         ((x & 0x0000000000FF0000ull) << 24) |
         ((x & 0x000000000000FF00ull) << 40) |
         ((x & 0x00000000000000FFull) << 56);
}

template <typename T>
static inline void byteswap_inplace(T& v) {
  if (sizeof(T) == 4) {
    uint32_t t; std::memcpy(&t, &v, 4); t = bswap32(t); std::memcpy(&v, &t, 4);
  } else if (sizeof(T) == 8) {
    uint64_t t; std::memcpy(&t, &v, 8); t = bswap64(t); std::memcpy(&v, &t, 8);
  } else {
    // unsupported size: no-op
  }
}

template <typename T>
static inline void maybe_swap_buffer(T* ptr, size_t n, bool need_swap) {
  if (!need_swap) return;
  for (size_t i = 0; i < n; ++i) byteswap_inplace(ptr[i]);
}

// Read one Fortran unformatted record: [int32 len] [payload] [int32 len]
static std::vector<char> read_fortran_record(std::ifstream& f, bool& need_swap) {
  uint32_t len1 = 0;
  if (!f.read(reinterpret_cast<char*>(&len1), 4)) throw std::runtime_error("Unexpected EOF (record start).");

  uint32_t len = len1;
  if (len != 84 && len != 164 && len > (1u<<26)) { // heuristic to detect swapped length
    uint32_t s = bswap32(len1);
    if (s < (1u<<24)) { need_swap = !need_swap; len = s; }
  }

  std::vector<char> payload(len);
  if (!f.read(payload.data(), len)) throw std::runtime_error("Unexpected EOF (record payload).");

  uint32_t len2 = 0;
  if (!f.read(reinterpret_cast<char*>(&len2), 4)) throw std::runtime_error("Unexpected EOF (record end).");
  if (need_swap) { len1 = bswap32(len1); len2 = bswap32(len2); }
  if (len1 != len2) throw std::runtime_error("Record length mismatch.");
  return payload;
}

// ------------------------ DCD Reader ------------------------
class DCDReader {
public:
  explicit DCDReader(const std::string& path) : ifs_(path, std::ios::binary) {
    if (!ifs_) throw std::runtime_error("Cannot open file: " + path);
    parse_header();
  }

  int natoms() const { return natoms_; }
  int nframes() const { return nset_; }
  bool has_unitcell() const { return has_unitcell_; }
  bool has_fixed_atoms() const { return has_fixed_atoms_; }

  std::vector<Frame> read_all() {
    std::vector<Frame> out;
    out.reserve(nset_ > 0 ? nset_ : 64);

    int frame_index = 0;
    while (ifs_.peek() != std::char_traits<char>::eof()) {
      Frame fr;
      read_one_frame(fr, frame_index);

      // Fixed-atoms: fill unchanged atoms from previous full coords
      if (has_fixed_atoms_) {
        if (frame_index == 0) {
          prev_full_coords_ = fr.atoms; // seed
        } else {
          for (int i = 0; i < natoms_; ++i) {
            if (!is_free_mask_[i]) fr.atoms[i] = prev_full_coords_[i];
          }
          prev_full_coords_ = fr.atoms; // update seed
        }
      }

      out.push_back(std::move(fr));
      ++frame_index;
    }
    return out;
  }

private:
  std::ifstream ifs_;
  bool need_swap_ = false;

  int natoms_ = 0;
  int nset_   = 0;
  bool has_unitcell_ = false;

  // Fixed-atoms bookkeeping (tentative, verified)
  bool has_fixed_atoms_ = false;
  int nfixed_ = 0;
  int nfreat_ = 0;
  std::vector<int> free_idx_;          // 0-based indices of free atoms
  std::vector<uint8_t> is_free_mask_;  // 0/1 mask per atom
  std::vector<Atom> prev_full_coords_; // last full coordinates for fixed-atom fill

  void parse_header() {
    // ---- Record 1: "CORD" + 20 ints ----
    auto rec1 = read_fortran_record(ifs_, need_swap_);
    if (rec1.size() < 4 + 20*4) throw std::runtime_error("DCD header too short.");

    char magic[5] = {0,0,0,0,0};
    std::memcpy(magic, rec1.data(), 4);
    if (std::strncmp(magic, "CORD", 4) != 0) {
      throw std::runtime_error("Not a DCD file (missing 'CORD').");
    }

    const int nctrl = 20;
    const int32_t* icntrl = reinterpret_cast<const int32_t*>(rec1.data() + 4);
    std::vector<int32_t> ctrl(nctrl);
    for (int i = 0; i < nctrl; ++i) {
      int32_t v = icntrl[i];
      if (need_swap_) byteswap_inplace(v);
      ctrl[i] = v;
    }

    nset_ = ctrl[0];
    const int iflag = ctrl[7];
    int nfixed_tentative = ctrl[9];

    has_unitcell_ = (iflag & 0x04) != 0;

    // ---- Record 2: title block ----
    auto rec2 = read_fortran_record(ifs_, need_swap_);
    if (rec2.size() < 4) throw std::runtime_error("Corrupt title block.");

    // ---- Record 3: number of atoms ----
    auto rec3 = read_fortran_record(ifs_, need_swap_);
    if (rec3.size() < 4) throw std::runtime_error("Corrupt atom count block.");
    int32_t n = 0;
    std::memcpy(&n, rec3.data(), 4);
    if (need_swap_) byteswap_inplace(n);
    if (n <= 0) throw std::runtime_error("Invalid atom count in DCD.");
    natoms_ = n;

    // ---- (Tentative) fixed-atoms detection, with verification ----
    has_fixed_atoms_ = false;
    nfixed_ = 0;
    nfreat_ = 0;
    free_idx_.clear();
    is_free_mask_.clear();

    if (nfixed_tentative > 0 && nfixed_tentative < natoms_) {
      // Save position, try to read IFREAT
      std::streampos pos = ifs_.tellg();
      int nfreat_try = natoms_ - nfixed_tentative;
      bool ok = false;

      try {
        auto rec_free = read_fortran_record(ifs_, need_swap_);
        if (rec_free.size() == static_cast<size_t>(nfreat_try * 4)) {
          free_idx_.resize(nfreat_try);
          const int32_t* raw = reinterpret_cast<const int32_t*>(rec_free.data());
          ok = true;
          for (int i = 0; i < nfreat_try; ++i) {
            int32_t idx = raw[i];
            if (need_swap_) byteswap_inplace(idx);
            if (idx <= 0 || idx > natoms_) { ok = false; break; }
            free_idx_[i] = idx - 1; // 1-based -> 0-based
          }
        }
      } catch (...) {
        ok = false;
      }

      if (ok) {
        has_fixed_atoms_ = true;
        nfixed_ = nfixed_tentative;
        nfreat_ = natoms_ - nfixed_;
        is_free_mask_.assign(natoms_, 0);
        for (int k = 0; k < nfreat_; ++k) is_free_mask_[ free_idx_[k] ] = 1;
      } else {
        // Rewind and treat as non-fixed (header field was misleading)
        ifs_.clear();
        ifs_.seekg(pos);
      }
    }
  }

  void read_unit_cell(Box& box) {
    // Some writers use 6 doubles, others 6 floats; accept either.
    auto rec = read_fortran_record(ifs_, need_swap_);
    if (rec.size() == 6*sizeof(double)) {
      double cell[6];
      std::memcpy(cell, rec.data(), 6*sizeof(double));
      maybe_swap_buffer(cell, 6, need_swap_);
      box.Lx = cell[0]; // A
      box.Ly = cell[2]; // B
      box.Lz = cell[4]; // C
    } else if (rec.size() == 6*sizeof(float)) {
      float cellf[6];
      std::memcpy(cellf, rec.data(), 6*sizeof(float));
      maybe_swap_buffer(cellf, 6, need_swap_);
      box.Lx = cellf[0];
      box.Ly = cellf[2];
      box.Lz = cellf[4];
    } else if (rec.size() == 0) {
      // Empty cell block (shouldn't happen): keep zeros
    } else {
      throw std::runtime_error("Unexpected unit cell record size.");
    }
  }

  void read_xyz_full(std::vector<Atom>& atoms) {
    auto recx = read_fortran_record(ifs_, need_swap_);
    auto recy = read_fortran_record(ifs_, need_swap_);
    auto recz = read_fortran_record(ifs_, need_swap_);
    
    const size_t want_f = static_cast<size_t>(natoms_) * sizeof(float);
    const size_t want_d = static_cast<size_t>(natoms_) * sizeof(double);
    
    if (recx.size() == want_f && recy.size() == want_f && recz.size() == want_f) {
      // float payloads
      std::vector<float> buf(natoms_);
      std::memcpy(buf.data(), recx.data(), recx.size());
      maybe_swap_buffer(buf.data(), natoms_, need_swap_);
      for (int i = 0; i < natoms_; ++i) atoms[i].x = buf[i];
    
      std::memcpy(buf.data(), recy.data(), recy.size());
      maybe_swap_buffer(buf.data(), natoms_, need_swap_);
      for (int i = 0; i < natoms_; ++i) atoms[i].y = buf[i];
    
      std::memcpy(buf.data(), recz.data(), recz.size());
      maybe_swap_buffer(buf.data(), natoms_, need_swap_);
      for (int i = 0; i < natoms_; ++i) atoms[i].z = buf[i];
    
    } else if (recx.size() == want_d && recy.size() == want_d && recz.size() == want_d) {
      // double payloads
      std::vector<double> bufD(natoms_);
      std::memcpy(bufD.data(), recx.data(), recx.size());
      maybe_swap_buffer(bufD.data(), natoms_, need_swap_);
      for (int i = 0; i < natoms_; ++i) atoms[i].x = static_cast<float>(bufD[i]);
    
      std::memcpy(bufD.data(), recy.data(), recy.size());
      maybe_swap_buffer(bufD.data(), natoms_, need_swap_);
      for (int i = 0; i < natoms_; ++i) atoms[i].y = static_cast<float>(bufD[i]);
    
      std::memcpy(bufD.data(), recz.data(), recz.size());
      maybe_swap_buffer(bufD.data(), natoms_, need_swap_);
      for (int i = 0; i < natoms_; ++i) atoms[i].z = static_cast<float>(bufD[i]);
    
    } else {
      throw std::runtime_error("XYZ array size mismatch (full): got {" +
                               std::to_string(recx.size()) + "," +
                               std::to_string(recy.size()) + "," +
                               std::to_string(recz.size()) + "} bytes");
    }
  }
  void read_xyz_free(std::vector<Atom>& atoms) {
    auto recx = read_fortran_record(ifs_, need_swap_);
    auto recy = read_fortran_record(ifs_, need_swap_);
    auto recz = read_fortran_record(ifs_, need_swap_);
    
    const size_t want_f = static_cast<size_t>(nfreat_) * sizeof(float);
    const size_t want_d = static_cast<size_t>(nfreat_) * sizeof(double);
    
    if (recx.size() == want_f && recy.size() == want_f && recz.size() == want_f) {
      std::vector<float> buf(nfreat_);
    
      std::memcpy(buf.data(), recx.data(), recx.size());
      maybe_swap_buffer(buf.data(), nfreat_, need_swap_);
      for (int i = 0; i < nfreat_; ++i) atoms[ free_idx_[i] ].x = buf[i];
    
      std::memcpy(buf.data(), recy.data(), recy.size());
      maybe_swap_buffer(buf.data(), nfreat_, need_swap_);
      for (int i = 0; i < nfreat_; ++i) atoms[ free_idx_[i] ].y = buf[i];
    
      std::memcpy(buf.data(), recz.data(), recz.size());
      maybe_swap_buffer(buf.data(), nfreat_, need_swap_);
      for (int i = 0; i < nfreat_; ++i) atoms[ free_idx_[i] ].z = buf[i];
    
    } else if (recx.size() == want_d && recy.size() == want_d && recz.size() == want_d) {
      std::vector<double> bufD(nfreat_);
    
      std::memcpy(bufD.data(), recx.data(), recx.size());
      maybe_swap_buffer(bufD.data(), nfreat_, need_swap_);
      for (int i = 0; i < nfreat_; ++i) atoms[ free_idx_[i] ].x = static_cast<float>(bufD[i]);
    
      std::memcpy(bufD.data(), recy.data(), recy.size());
      maybe_swap_buffer(bufD.data(), nfreat_, need_swap_);
      for (int i = 0; i < nfreat_; ++i) atoms[ free_idx_[i] ].y = static_cast<float>(bufD[i]);
    
      std::memcpy(bufD.data(), recz.data(), recz.size());
      maybe_swap_buffer(bufD.data(), nfreat_, need_swap_);
      for (int i = 0; i < nfreat_; ++i) atoms[ free_idx_[i] ].z = static_cast<float>(bufD[i]);
    
    } else {
      throw std::runtime_error("XYZ array size mismatch (free): got {" +
                               std::to_string(recx.size()) + "," +
                               std::to_string(recy.size()) + "," +
                               std::to_string(recz.size()) + "} bytes");
    }
  }

  static inline void parse_cell_record(const std::vector<char>& rec, bool need_swap, Box& box) {
    double v[6] = {0,0,0,0,0,0};
    
    if (rec.size() == 6*sizeof(double)) {
      std::memcpy(v, rec.data(), 6*sizeof(double));
      maybe_swap_buffer(v, 6, need_swap);
    } else if (rec.size() == 6*sizeof(float)) {
      float vf[6];
      std::memcpy(vf, rec.data(), 6*sizeof(float));
      maybe_swap_buffer(vf, 6, need_swap);
      for (int i = 0; i < 6; ++i) v[i] = vf[i];
    } else {
      throw std::runtime_error("Unexpected unit cell record size.");
    }
    
    // Default CHARMM/NAMD/LAMMPS convention: (A, gamma, B, beta, C, alpha)
    double Lx = v[0];
    double Ly = v[2];
    double Lz = v[4];
    
    // Heuristics for quirky writers:
    // If C came out 0 or tiny, try the last slot (some paths put C at index 5),
    // or any entry that must be a length (e.g., > 180, which cannot be an angle).
    auto is_angle_like = [](double a) {
      // many writers use 0 or ~90 for angles; treat [0..180] as angle-like
      return (a >= 0.0 && a <= 180.0);
    };
    
    if (Lz <= 1e-9) {
      if (v[5] > 1e-9 && (!is_angle_like(v[5]) || v[5] > 180.0)) {
        Lz = v[5];                       // take the last slot if it looks like a length
      } else {
        // As a last resort, scan all positions for a plausible length not used yet
        // Prefer values > 180 (cannot be angle), else largest positive not equal to ~90.
        int candidates[6] = {0,1,2,3,4,5};
        double best = 0.0;
        for (int k : candidates) {
          if (k == 0 || k == 2 || k == 4) continue; // already assigned
          double x = v[k];
          if (x > best && (!is_angle_like(x) || x > 180.0)) best = x;
        }
        if (best > 0.0) Lz = best;
      }
    }
    
    box.Lx = Lx;
    box.Ly = Ly;
    box.Lz = Lz;
  }


  void read_one_frame(Frame& fr, int frame_index) {
    fr.atoms.resize(natoms_);
    fr.box = Box{}; // default zeros
    
    // Read the first record of the frame. It might be:
    //   (a) unit cell (6 floats or 6 doubles), or
    //   (b) X array
    auto rec0 = read_fortran_record(ifs_, need_swap_);
    
    const size_t wantX_f = static_cast<size_t>(natoms_) * sizeof(float);
    const size_t wantX_d = static_cast<size_t>(natoms_) * sizeof(double);
    
    bool rec0_is_cell = (rec0.size() == 6*sizeof(double)) || (rec0.size() == 6*sizeof(float));
    
    std::vector<char> recx, recy, recz;
    
    if (rec0_is_cell) {
      // We discovered a per-frame unit cell even if the header didn't say so
      parse_cell_record(rec0, need_swap_, fr.box);
    
      // Now read the three coord arrays
      recx = read_fortran_record(ifs_, need_swap_);
      recy = read_fortran_record(ifs_, need_swap_);
      recz = read_fortran_record(ifs_, need_swap_);
    } else {
      // No cell block here; rec0 is actually X
      recx = std::move(rec0);
      recy = read_fortran_record(ifs_, need_swap_);
      recz = read_fortran_record(ifs_, need_swap_);
    }
    
    // Accept either float or double coordinate payloads and convert to float
    if (recx.size() == wantX_f && recy.size() == wantX_f && recz.size() == wantX_f) {
      std::vector<float> buf(natoms_);
    
      std::memcpy(buf.data(), recx.data(), recx.size());
      maybe_swap_buffer(buf.data(), natoms_, need_swap_);
      for (int i = 0; i < natoms_; ++i) fr.atoms[i].x = buf[i];
    
      std::memcpy(buf.data(), recy.data(), recy.size());
      maybe_swap_buffer(buf.data(), natoms_, need_swap_);
      for (int i = 0; i < natoms_; ++i) fr.atoms[i].y = buf[i];
    
      std::memcpy(buf.data(), recz.data(), recz.size());
      maybe_swap_buffer(buf.data(), natoms_, need_swap_);
      for (int i = 0; i < natoms_; ++i) fr.atoms[i].z = buf[i];
    
    } else if (recx.size() == wantX_d && recy.size() == wantX_d && recz.size() == wantX_d) {
      std::vector<double> bufD(natoms_);
    
      std::memcpy(bufD.data(), recx.data(), recx.size());
      maybe_swap_buffer(bufD.data(), natoms_, need_swap_);
      for (int i = 0; i < natoms_; ++i) fr.atoms[i].x = static_cast<float>(bufD[i]);
    
      std::memcpy(bufD.data(), recy.data(), recy.size());
      maybe_swap_buffer(bufD.data(), natoms_, need_swap_);
      for (int i = 0; i < natoms_; ++i) fr.atoms[i].y = static_cast<float>(bufD[i]);
    
      std::memcpy(bufD.data(), recz.data(), recz.size());
      maybe_swap_buffer(bufD.data(), natoms_, need_swap_);
      for (int i = 0; i < natoms_; ++i) fr.atoms[i].z = static_cast<float>(bufD[i]);
    
    } else {
      throw std::runtime_error(
        "XYZ array size mismatch (full): got {" +
        std::to_string(recx.size()) + "," +
        std::to_string(recy.size()) + "," +
        std::to_string(recz.size()) + "} bytes");
    }
  }
};
//----------------------- END READ DCD -------------------------------

//------------------------ START STATISTICS ---------------------------
struct Histogram {
  std::vector<float> edges;
  std::vector<float> counts;
  float binWidth(std::size_t i=0) const {
    return edges.size() > 1? (edges[i+1]-edges[i]) : 0.0;
  }
};

Histogram makeHistogram(const std::vector<float>& data, size_t nbins, float minEdge, float maxEdge){
  Histogram h;
  h.edges.resize(nbins + 1);
  h.counts.assign(nbins, 0.0);

  if (nbins == 0) return h;
  if (maxEdge == minEdge) maxEdge = minEdge + 1;

  const float width = (maxEdge - minEdge)/ static_cast<float>(nbins);
  for( size_t i =0; i<= nbins; i++ ) h.edges[i] = minEdge + i* width;

  for(size_t k=0; k<data.size(); k++) {
    float x = data[k];
    long idx = static_cast<long>(std::floor( (x-minEdge) / width));
    if (idx < 0 ) idx =0;
    if (idx > static_cast<long>(nbins)) idx = static_cast<long>(nbins)-1;
    h.counts[static_cast<size_t>(idx)] += 1.0;
  }
  return h;
}

Histogram makeHistogramAuto(const std::vector<float>& data, size_t nbins) {
  if ( data.empty() || nbins == 0) return {};
  auto mm = std::minmax_element(data.begin(), data.end());
  float mn = *mm.first;
  float mx = *mm.second;
  if ( mx == mn ) {mn -= 0.5; mx += 0.5;}

  return makeHistogram(data, nbins, mn, mx);
}
std::vector<float> normalizePDF(const Histogram& h) {
  float N = std::accumulate(h.counts.begin(), h.counts.end(), 0.0);
  std::vector<float> pdf(h.counts.size(), 0.0);
  if ( N > 0 && h.edges.size() >=2) {
    for (size_t i =0; i< h.counts.size(); i++) {
      float bw = h.edges[i+1] - h.edges[i];
      if (bw > 0 ) pdf[i] = h.counts[i] / (N*bw);
    }
  }
  return pdf;
}

void printHistogram( const Histogram& h, std::ostream& os = std::cout) {
  std::vector<float> pdf = normalizePDF(h);
  os.setf(std::ios::fixed);
  os << std::setprecision(10);
  int nr = h.edges.size();
  int half = nr/2;
  for( size_t i =half; i< nr; i++ ) {
    float x = 0.5 * (h.edges[i] + h.edges[i+1]);
    os << x << "\t" << (pdf[i] + pdf[2*half-i-1])/2.0 << "\n";
  }
}

std::vector<float> histPDFatEdges(const Histogram& h) {
  size_t nb = h.counts.size();
  std::vector<float> y; 
  y.assign(nb+1, 0.0);

  float N = std::accumulate(h.counts.begin(), h.counts.end(), 0.0);
  if ( N<=0.0 || h.edges.size() != nb+1 ) return y;

  for (size_t i =0; i< nb ; i++){
    float bw = h.edges[i+1] - h.edges[i];
    y[i] = (bw > 0.0) ? (h.counts[i] / N / bw) : 0.0;
  }
  y[nb] = y[nb - (nb > 0 ? 1: 0)];
  return y;
}

float linInterp(float x0, float y0, float x1, float y1, float x) {
  float dx = x1-x0;
  if (dx==0.0) return 0.5*(y0+y1);
  float t= (x-x0)/dx;
  return y0 + t* (y1-y0);
}

float integratePDF(const Histogram& h, float xi, float xf) {
  if (h.edges.size() < 2 || h.counts.size() +1 != h.edges.size()) return 0.0;
  if (xi > xf ) std::swap(xi, xf);

  const std::vector<float>& x = h.edges;
  std::vector<float> y = histPDFatEdges(h);

  float a = x.front();
  float b = x.back();
  if ( xf <= a || xi >= b) return 0.0;
  xi = std::max(xi, a);
  xf = std::min(xf, b);

  std::vector<float>::const_iterator itL = std::upper_bound(x.begin(), x.end(), xi);
  std::vector<float>::const_iterator itR = std::upper_bound(x.begin(), x.end(), xf);
  size_t j = (itL == x.begin() ? 0 : static_cast<size_t>((itL - x.begin()) -1));
  size_t k = (itR == x.begin() ? 0 : static_cast<size_t>((itR - x.begin()) -1));
  if( j >= x.size() -1) j = x.size() -2;
  if( k >= x.size() -1) k= x.size() -2;
  float y_xi = linInterp(x[j], y[j], x[j+1], y[j+1], xi);
  float y_xf = linInterp(x[k], y[k], x[k+1], y[k+1], xf);
  if ( j== k) { 
    return 0.5 * (y_xi + y_xf) * (xf - xi);
  }

  float area = 0.0;
  area += 0.5 * (y_xi + y[j+1]) * (x[j+1] - xi);

  for ( size_t m = j+1; m< k; m++ ) {
    area += 0.5 * (y[m]+y[m+1]) * (x[m+1] - x[m]);
  }
  area += 0.5 * (y[k] + y_xf) * (xf - x[k]);
  return area;
}


int binIndexFromEdges(float x, const std::vector<float>& edges){
  if( x< edges.front() || x > edges.back() ) return -1;
  if ( x == edges.back()) return int(edges.size())-2;
  std::vector<float>::const_iterator it = std::upper_bound(edges.begin(), edges.end(), x);
  int idx = int (it - edges.begin()) -1;
  if (idx < 0 || idx >= int(edges.size()) -1 ) return -1;
  return idx;
}

struct momentStats {
  std::vector<float> mean, se;
};

struct Hist2D {
  std::vector<float> zedges, aedges;
  std::vector<float> counts;
  int nz = 0, na = 0;
  float totalSamples = 0.0;

  float& at(int iz, int ia) { return counts[size_t(iz) * size_t(na) + size_t(ia)];}
  float  at(int iz, int ia) const { return counts[size_t(iz) * size_t(na) + size_t(ia)];}

  float zwidth(int iz ) const { return zedges[iz+1] - zedges[iz];}
  float awidth(int ia ) const { return aedges[ia+1] - aedges[ia];}

  float zcenter(int iz ) const { return 0.5*(zedges[iz+1] + zedges[iz]);}
  float acenter(int ia ) const { return 0.5*(aedges[ia+1] + aedges[ia]);}

  void normalizeMass() {
    if ( totalSamples <= 0.0) return;
    for (size_t i=0; i< counts.size(); i++) counts[i] /= totalSamples;
  }

  void normalizeDensity() {
    if (totalSamples <= 0.0) return;
    for (int iz =0 ; iz < nz; iz++) {
      for (int ia=0; ia < na; ia++ ) {
        float area = zwidth(iz) * awidth(ia);
        if (area > 0) {
          at(iz, ia) /= (totalSamples * area);
        }
      }
    }
  }

  void printHist(const std::string& path) const {
    FILE* f = std::fopen(path.c_str(), "w");
    if (!f) throw std::runtime_error("Fail to open file: " + path);
    std::fprintf(f, "z, cos(theta)\n");
    for( int iz=0; iz< nz; iz++ ) {
      for (int ia=0; ia<na; ia++ ) {
        std::fprintf(f, "%.5f\t%.5f\t%.5f\n", zcenter(iz), acenter(ia), at(iz, ia));
      }
    }
    std::fclose(f);
  }

  std::vector<float> aSliceDensity(float zq) const {
    int iz = binIndexFromEdges(zq, zedges);
    if ( iz< 0) throw std::runtime_error("z value out of range");
    std::vector<float> pdf (na, 0.0);
    float rowSum = 0.0;
    for (int ia = 0; ia< na; ia++) rowSum += at(iz,ia);
    if (rowSum > 0.0) {
      for( int ia=0; ia < na ; ia++ ) {
        float da = awidth(ia);
        if (da > 0.0) pdf[ia] = at(iz, ia) / rowSum / da;
      }
    }
    return pdf;
  }

  void printZsliceHist(float zq, std::ostream& os = std::cout,int precision=10 ) const {
    int iz = binIndexFromEdges(zq, zedges);
    if( iz<0) throw std::runtime_error("z value is out of range");
    std::vector<float> v = aSliceDensity(zq);
    os.setf(std::ios::fixed);
    os << std::setprecision(precision);
    os << "# z-bin center = " << zcenter(iz) << " (query z = " << zq << ")\n";
    os << "# angle center, pdf\n";
    for( size_t i =0; i< na; i++ ) {
      os << acenter(i) << "\t" << v[i] << "\n";
    }
  }

  std::vector<float> kthMomentByZ(int k) const {
    std::vector<float> out(size_t(nz), std::numeric_limits<double>::quiet_NaN());
    for( int iz = 0; iz < nz ; iz ++ ) {
      long double num = 0.0L, den = 0.0L;
      for( int ia = 0; ia < na; ia++){
        const double ac = acenter(ia);
        long double w = at(iz, ia) * awidth(ia);
        if ( w <= 0.0L) continue;
        num += std::pow(ac, k) * w;
        den += w;
      }
      if (den > 0.0L) out[size_t(iz)] = float(num/den);
    }
    return out;
  }
  momentStats kthMomentStatsByZ(int k, const std::vector<float>& raw) const {
    if (k<0) throw std::runtime_error("k must be >= 0");
    if (zedges.size() != nz +1 || aedges.size() != na+1 || raw.size() != nz*na ) {
      throw std::runtime_error("size mismatch");
    }
    momentStats out;
    out.mean.assign(nz, std::numeric_limits<float>::quiet_NaN());
    out.se.assign(nz, std::numeric_limits<float>::quiet_NaN());

    for(int iz =0; iz < nz; iz++) {
      long double N = 0.0L, sumY=0.0L, sumY2=0.0L;
      for( int ia=0; ia<na; ia++) {
        const size_t idx = size_t(iz) * size_t(na) + size_t(ia);
        const long double n = raw[idx];
        if ( n <= 0.0L) continue;
        const long double ac = acenter(ia);
        const long double y = std::pow(ac, k);
        sumY += n*y;
        sumY2+= n*y*y;
        N += n;
      }
      if ( N >=1.0L) {
        const long double mean = sumY / N;
        out.mean[size_t(iz)] = float(mean);
        if ( N > 1.0L) {
          const long double varY = (sumY2 - N*mean*mean) / (N-1.0L);
          out.se[size_t(iz)] = (varY > 0.0L ) ? float(std::sqrt(varY/N)) : 0.0;
        } else {
          out.se[size_t(iz)] = std::numeric_limits<float>::infinity();
        }
      }
    }
    return out;
  }

  // print first and third moments
  /*
  void printMoments(std::ostream& os = std::cout, int precision=10 ) const {
    auto m1 = kthMomentByZ(1);
    auto m3 = kthMomentByZ(3);
    os.setf(std::ios::fixed);
    os << std::setprecision(precision);
    os << "# z center, first moment, third moment\n";
    int st = static_cast<int>(nz/2);
    for( size_t i =st; i< nz; i++ ) {
      os << zcenter(i) << "\t" << (m1[i] - m1[2*st-i-1])/2.0 << "\t" << (m3[i] - m3[2*st-i-1])/2.0 << "\n";
    }
  }
  */
  void printMoments(std::ostream& os, momentStats& stat1, momentStats& stat3, int precision=10 ) const {
    os.setf(std::ios::fixed);
    os << std::setprecision(precision);
    os << "# z center, first moment, standard error, third moment, standard error\n";
    int st = static_cast<int>(nz/2);
    for( size_t i =st; i< nz; i++ ) {
      os << zcenter(i) << "\t" 
        << (stat1.mean[i] - stat1.mean[2*st-i-1])/2.0 << "\t" << (stat1.se[i] + stat1.se[2*st-i-1])/2.0 << "\t" 
        << (stat3.mean[i] - stat3.mean[2*st-i-1])/2.0 << "\t" << (stat3.se[i] + stat3.se[2*st-i-1])/2.0 << "\n";
    }
  }
};



std::pair<float, float> minmaxComponent(const std::vector<Angle>& a, bool use_z) {
  float lo = std::numeric_limits<float>::infinity();
  float hi = -std::numeric_limits<float>::infinity();
  for( size_t i =0; i<a.size(); i++ ) {
    float v = use_z ? float(a[i].zloc) : float(a[i].chCos);
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  if (!std::isfinite(lo) || !std::isfinite(hi)) { lo = 0.0; hi = 1.0;}
  if (hi==lo) hi = lo + 1;
  return std::pair<float, float> (lo, hi);
}

Hist2D makeHist2D(const std::vector<Angle>& data, int nz, int na) {
  if (nz <= 0 || na <= 0) throw std::runtime_error("nz/na must be positive\n");
  float zmin, zmax;
  auto zmm = minmaxComponent(data, true);
  zmin = zmm.first; zmax = zmm.second;

  float amin, amax;
  auto amm = minmaxComponent(data, false);
  amin = amm.first; amax = amm.second;


  Hist2D H;
  H.nz = nz; H.na = na;
  H.zedges.resize(size_t(nz)+1);
  H.aedges.resize(size_t(na)+1);
  H.counts.assign(size_t(nz)*size_t(na), 0.0);
  float dz = (zmax-zmin)/float(nz);
  float da = (amax-amin)/float(na);
  for(int i = 0; i<=nz; i++) H.zedges[i] = zmin + dz*float(i);
  for(int j = 0; j<=na; j++) H.aedges[j] = amin + da*float(j);
  return H;
}

float wrapMinImage( float x, float L ) {
  return x - L * std::floor((x+0.5*L)/L);
}

Hist2D makeHist2DSymZ(const std::vector<Angle>& data, float center, float Lz, float zBinWidth, int na) {
  if (!(Lz > 0.0) || !(zBinWidth > 0.0) ||  na <= 0) {
    throw std::runtime_error("nz/na must be positive\n");
  }

  // find max absolute displacement from center
  float zmaxDisp = 0.0;
  for (const auto& a : data) {
    float d = wrapMinImage(wrapPBC(a.zloc, Lz) - center, Lz);
    float ab = std::fabs(d);
    if (ab > zmaxDisp) zmaxDisp = ab;
  }
  zmaxDisp = std::min(zmaxDisp, 0.5f * Lz);

  // choose bin count per side
  int nSide = std::max(1, int(std::ceil(zmaxDisp / zBinWidth)));
  float R = nSide * zBinWidth;
  int nz = nSide*2;

  std::cout << "max z : " << zmaxDisp << " and nz : " << nz << "\n";

  float amin, amax;
  auto amm = minmaxComponent(data, false);
  amin = amm.first; amax = amm.second;

  Hist2D H;
  H.nz = nz; H.na = na;
  H.zedges.resize(size_t(nz)+1);
  H.aedges.resize(size_t(na)+1);
  H.counts.assign(size_t(nz)*size_t(na), 0.0);
  float dz = zBinWidth;
  float da = (amax-amin)/float(na);
  for(int i = 0; i<=nz; i++) H.zedges[i] = center-R + dz*float(i);
  for(int j = 0; j<=na; j++) H.aedges[j] = amin + da*float(j);
  return H;
}

void hist2DAccumulate(Hist2D& H, const std::vector<Angle>& data, float weight = 1.0) {
  for (size_t k=0; k< data.size(); k++ ) {
    float z = float(data[k].zloc) ;
    float a = float(data[k].chCos);
    if (!std::isfinite(z) || !std::isfinite(a)) continue;
    int iz = binIndexFromEdges(z, H.zedges);
    int ia = binIndexFromEdges(a, H.aedges);
    if (iz >= 0 && ia >= 0) {
      H.at(iz, ia) += weight;
      H.totalSamples += weight;
    }
  }
}

void hist2DAdd(Hist2D& H, float z, float a, float weight=1.0) {
  if(!std::isfinite(z) || !std::isfinite(a)) return;
  int iz = binIndexFromEdges(z, H.zedges);
  int ia = binIndexFromEdges(a, H.aedges);
  if (iz >= 0 && ia >= 0) {
    H.at(iz, ia) += weight;
    H.totalSamples += weight;
  }
}


// ---- Make an empty histogram with existing edges (no auto-ranging) ----
Hist2D makeEmptyHistWithEdges(const std::vector<float>& zEdges,
                              const std::vector<float>& aEdges) {
  if (zEdges.size() < 2 || aEdges.size() < 2) {
    throw std::runtime_error("makeEmptyHistWithEdges: edges too short");
  }
  Hist2D H;
  H.nz = int(zEdges.size()) - 1;
  H.na = int(aEdges.size()) - 1;
  H.zedges = zEdges;
  H.aedges = aEdges;
  H.counts.assign(size_t(H.nz) * size_t(H.na), 0.0f);
  H.totalSamples = 0.0f;
  return H;
}

// ---- raw-counts kth moment of 'a' per z-row (no density; uses bin centers) ----
static std::vector<float> kthMomentByZFromCounts(const Hist2D& H, int k) {
  std::vector<float> out(size_t(H.nz), std::numeric_limits<float>::quiet_NaN());
  for (int iz = 0; iz < H.nz; ++iz) {
    long double N = 0.0L, sumY = 0.0L;
    for (int ia = 0; ia < H.na; ++ia) {
      const long double n  = H.at(iz, ia);
      if (n <= 0.0L) continue;
      const long double ac = H.acenter(ia);
      const long double y  = std::pow(ac, k);
      sumY += n * y;
      N    += n;
    }
    if (N > 0.0L) out[size_t(iz)] = float(sumY / N);
  }
  return out;
}

// ---- Block-average: mean & SE over blocks, per z-bin (uses raw counts) ----
struct BlockMomentStats {
  std::vector<float> mean;  // size = nz
  std::vector<float> se;    // size = nz
};

BlockMomentStats blockAverageMomentByZ(
    const std::vector<std::vector<Angle>>& frames, // frames[t] = angles for that frame
    int k,
    int B,
    const std::vector<float>& zEdges,
    const std::vector<float>& aEdges)
{
  if (B < 2) throw std::runtime_error("blockAverageMomentByZ: B must be >= 2");
  const int F = int(frames.size());
  if (F < B) throw std::runtime_error("blockAverageMomentByZ: F < B");
  // Prepare per-block histograms with identical edges
  std::vector<Hist2D> HB;
  HB.reserve(B);
  for (int b = 0; b < B; ++b) HB.push_back(makeEmptyHistWithEdges(zEdges, aEdges));

  // Even split frames into B blocks: block b => [start[b], start[b+1])
  std::vector<int> start(B + 1, 0);
  for (int b = 0; b <= B; ++b) start[b] = (b * F) / B;

  // Accumulate raw counts per block
  for (int b = 0; b < B; ++b) {
    for (int f = start[b]; f < start[b + 1]; ++f) {
      hist2DAccumulate(HB[b], frames[size_t(f)], /*weight=*/1.0f);
    }
  }

  const int nz = int(zEdges.size()) - 1;
  BlockMomentStats out;
  out.mean.assign(size_t(nz), std::numeric_limits<float>::quiet_NaN());
  out.se  .assign(size_t(nz), std::numeric_limits<float>::quiet_NaN());

  std::vector<float> mb(size_t(B), std::numeric_limits<float>::quiet_NaN());
  std::vector<char>  has(size_t(B), 0);

  // For each z-row: compute block means, then mean & SE across blocks that have data
  for (int iz = 0; iz < nz; ++iz) {
    int M = 0;
    for (int b = 0; b < B; ++b) {
      // row count in this block
      long double Nrow = 0.0L;
      for (int ia = 0; ia < HB[b].na; ++ia) Nrow += HB[b].at(iz, ia);
      if (Nrow > 0.0L) {
        // kth moment for this row from this block
        long double N = 0.0L, sumY = 0.0L;
        for (int ia = 0; ia < HB[b].na; ++ia) {
          const long double n  = HB[b].at(iz, ia);
          if (n <= 0.0L) continue;
          const long double ac = HB[b].acenter(ia);
          const long double y  = std::pow(ac, k);
          sumY += n * y; N += n;
        }
        mb[size_t(b)]  = (N > 0.0L) ? float(sumY / N) : std::numeric_limits<float>::quiet_NaN();
        has[size_t(b)] = 1;
        ++M;
      } else {
        has[size_t(b)] = 0;
      }
    }

    if (M >= 1) {
      long double mu = 0.0L;
      for (int b = 0; b < B; ++b) if (has[size_t(b)]) mu += mb[size_t(b)];
      mu /= M;
      out.mean[size_t(iz)] = float(mu);

      if (M >= 2) {
        long double s2 = 0.0L;
        for (int b = 0; b < B; ++b) if (has[size_t(b)]) {
          const long double d = mb[size_t(b)] - mu;
          s2 += d * d;
        }
        out.se[size_t(iz)] = float(std::sqrt(s2 / (M * (M - 1))));
      } // else leave NaN
    }
  }
  return out;
}




// -------------- END STATISTICS ---------------

// -------------- START BASIC -------------------
float applyPBC(float x, float box){
  float hbox = box/2.0;
  float wrapped = fmod(x + hbox, box);
  if (wrapped < 0) wrapped += box;
  return wrapped - hbox;
}

float distance(const Atom& a, const Atom& b, const Box& box) {
  float dx, dy, dz, rsq;
  dx=applyPBC( a.x - b.x, box.Lx);
  dy=applyPBC( a.y - b.y, box.Ly);
  dz=applyPBC( a.z - b.z, box.Lz);
  rsq = dx*dx +dy*dy +dz*dz;
  return std::sqrt(rsq);
}

Angle computeTerminalAngle( Targets& target, const std::vector<Atom>& atoms, const Box& box) {
  Angle result;
  // compute C-C vector angle
  float dz;
  dz=applyPBC(atoms[target.t].z - atoms[target.c].z, box.Lz);
  result.ccCos = dz/distance(atoms[target.t], atoms[target.c], box);
  //std::cout << "target : " << target.t << " carbon : " << target.c << " dz " << dz << " result " << result.ccCos << "\n";

  // compute averaged C-H vector angle, assume CH bond length fixed
  float dx1, dx2, dx3, dy1, dy2, dy3, dz1, dz2, dz3, dx, dy;
  dx1 = applyPBC(atoms[target.h1].x - atoms[target.t].x, box.Lx);
  dx2 = applyPBC(atoms[target.h2].x - atoms[target.t].x, box.Lx);
  dx3 = applyPBC(atoms[target.h3].x - atoms[target.t].x, box.Lx);
  dx = dx1+dx2+dx3;

  dy1 = applyPBC(atoms[target.h1].y - atoms[target.t].y, box.Ly);
  dy2 = applyPBC(atoms[target.h2].y - atoms[target.t].y, box.Ly);
  dy3 = applyPBC(atoms[target.h3].y - atoms[target.t].y, box.Ly);
  dy = dy1 + dy2 + dy3;

  dz1 = applyPBC(atoms[target.h1].z - atoms[target.t].z, box.Lz);
  dz2 = applyPBC(atoms[target.h2].z - atoms[target.t].z, box.Lz);
  dz3 = applyPBC(atoms[target.h3].z - atoms[target.t].z, box.Lz);
  dz = dz1 + dz2 + dz3;

  float rsq = dx*dx + dy*dy + dz*dz;
  result.chCos = dz/std::sqrt(rsq);

  result.zloc = atoms[target.t].z;
  
  return result;
}




// ------------------------ MAIN ------------------------
int main(int argc, char** argv) {
  if (argc != 6 ) {
    std::cerr << "Error: Num Args not Match\n" ;
    std::cerr << "Usage : " << argv[0] << "dumpDir targetMol numMols timestep(ps) eqtime(ns)\n";
    return 1;
  }
  std::string dump = std::string(argv[1])  + argv[2] + ".dcd";
  int nummols = std::stoi(argv[3]);
  float timestep = std::stof(argv[4]);
  float eqtime = std::stof(argv[5]);

  // Read DCD file
  std::cout << "---- Reading dcd file at " << dump << " ----\n";
  DCDReader reader(dump);
  int numatoms = reader.natoms();
  std::cout << "number of atoms: " << numatoms << "\n";
  std::cout << "number of mols : " << nummols << "\n";
  std::cout << "number of frames (original): " << reader.nframes() << "\n";
  std::cout << "Initial " << eqtime << " ns are used for equilibration\n";
  int eqsnap=static_cast<int>(1000*eqtime/timestep);
  auto frames = reader.read_all();
  size_t numsnap = frames.size();
  std::cout << "number of frames (after EQ): " << numsnap - eqsnap << "\n\n";
  int oneMolecule = numatoms / nummols;


  // Read Lammps input file to get type, mass, mol etc
  std::string lammpsInput = std::string("../../../") + argv[2] +"/system.data";
  LammpsData data = readLammpsData(lammpsInput);

  std::cout << "---- Reading Lammps input at " << lammpsInput << " ----\n";
  std::cout << "Atoms style : " 
            << (data.atomsStyleHint.empty() ? "<none>" : data.atomsStyleHint) << "\n";
  std::cout << "Read atoms: " << data.atoms.size() << "\n";

  //generate type vector
  std::vector<int> atomTypes;
  if (!data.atoms.empty()) {
    for (size_t i = 0; i < data.atoms.size(); ++i) {
      const AtomRow& a = data.atoms[i];
      atomTypes.push_back(a.type);
    }
  }

  // Output video file
  std::ofstream ofs("video.lammpstrj");
  if (!ofs) {
    std::cerr << "Error. Cannot open the file\n";
    return 1;
  }

  std::cout << "Read masses: " << data.masses.size() << " types\n";
  std::vector<std::pair<int,double>> type2mass;
  if (!data.masses.empty()) {
    // print a map : type to mass 
    type2mass.reserve(data.masses.size());
    for (const auto& kv : data.masses) type2mass.emplace_back(kv.first, kv.second);
    std::sort(type2mass.begin(), type2mass.end(), [](const std::pair<int,double>& a, const std::pair<int,double>& b){ return a.first < b.first;});
    for (size_t i = 0; i < type2mass.size(); ++i) {
      std::cout << "  type " << type2mass[i].first << " mass " << type2mass[i].second << "\n";
    }
  }

  auto initL = findTargets(argv[2], true); // left
  auto initR = findTargets(argv[2], false); // right

  numsnap = numsnap-eqsnap;
  std::vector<Angle> angles; 
  angles.reserve(nummols*2*numsnap);
  for( size_t time=eqsnap; time<numsnap+eqsnap; time++ ) {
    auto frame = frames[time].atoms;
    const auto box = frames[time].box;
    // first, find pbc center of mass, shift it to the center of box
    float zc = comZPBC(frame, atomTypes, type2mass, box.Lz); 
    for ( auto& atom : frame ) atom.z = wrapPBC(atom.z - zc + box.Lz/2.0, box.Lz);
    // second, find geometric center of mass, shift box 
    float zc2 = comZ(frame, atomTypes, type2mass);
    for ( auto& atom : frame ) atom.z = wrapPBC(atom.z - zc2 + box.Lz/2.0, box.Lz);
    //std::cout << "Time " << time << " z center was " << zc << " now goes to " << zc2 << " \n";

    // find true center of mass, and shift again
    for( int mol=0; mol < nummols; mol++ ) {
      auto indexL = initL + static_cast<int>(oneMolecule*mol);
      auto indexR = initR + static_cast<int>(oneMolecule*mol);
      Angle left  = computeTerminalAngle(indexL, frame, box);
      Angle right = computeTerminalAngle(indexR, frame, box);
      angles.emplace_back(left);
      angles.emplace_back(right);
    }

    if (time % 100 == 0 ) {
      writeLammpsDump(ofs, time, frame, atomTypes,  box);
    }
  }
  ofs.close();

  /// makeHist2D(std::vector<Angle> data, int nz, int na)
  //auto H = makeHist2D(angles, 30, 50);

  /// makeHist2DSymZ(std::vector<Angle> data, float center, float Lz, float zBinWidth, int na)
  auto H =  makeHist2DSymZ(angles, 100, 200, 2, 50);
  hist2DAccumulate(H, angles);
  std::vector<float> raw = H.counts; 
  H.normalizeDensity();
  H.printHist("2dhist.dat");
  std::ofstream output("slice.dat");
  H.printZsliceHist(120, output);

  // generate terminal carbon density profile
  std::vector<float> rho(angles.size());
  std::transform(angles.begin(), angles.end(), rho.begin(), [](const Angle& a){ return static_cast<float>(a.zloc);});
  auto density = makeHistogramAuto(rho, 200);
  std::ofstream out("densityProfile.dat");
  printHistogram(density, out);

  auto stats1 = H.kthMomentStatsByZ(/*k=*/1, raw);
  auto stats3 = H.kthMomentStatsByZ(/*k=*/3, raw);
  std::ofstream output2("moments.dat");
  H.printMoments(output2, stats1, stats3);


  // try block averaging
  bool blockAverage=true;
  if(blockAverage) {
    std::vector<std::vector<Angle>> traj(numsnap, std::vector<Angle>(nummols*2));
    for( size_t time =0; time < numsnap; time++ ) {
      for ( int mol =0; mol< nummols*2; mol++ ) {
        traj[time][mol] = angles[time * (nummols*2) + mol];
      }
    }
    int B =20;
    auto stats1 = blockAverageMomentByZ(traj, 1, B, H.zedges, H.aedges);
    auto stats3 = blockAverageMomentByZ(traj, 3, B, H.zedges, H.aedges);
    std::ofstream output3("BAmoments.dat");
    H.printMoments(output3, *(momentStats*)&stats1, *(momentStats*)&stats3, 8);
  }






  return 0;
}

