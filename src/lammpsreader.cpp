#include <ioTraj/lammpsreader.hpp>

#include <vector>
#include <string>
#include <fstream>
#include <sstream> // istringstream
#include <cctype> // isspace, tolower
#include <cstdlib> // strtod
#include <iomanip>


namespace ioTraj {

static inline void trimInplace(std::string& s) {
  size_t i = 0;
  while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
  size_t j = s.size();
  while (j > i && std::isspace(static_cast<unsigned char>(s[j-1]))) --j;
  s.assign(s, i, j - i);
}

static inline void stripComment(std::string& s) {
  size_t pos = s.find('#');
  if (pos != std::string::npos) s.erase(pos);
}

static inline bool startsWithCi(const std::string& s, const std::string& prefix) {
  if (s.size() < prefix.size()) return false;
  for (size_t i = 0; i < prefix.size(); ++i) {
    if (std::tolower(static_cast<unsigned char>(s[i])) !=
        std::tolower(static_cast<unsigned char>(prefix[i]))) return false;
  }
  return true;
}

static inline bool parseIntStrict(const std::string& tok, int& out) {
  const char* p = tok.c_str();
  char* end = nullptr;
  long v = std::strtol(p, &end, 10);
  if (end == p || *end != '\0') return false;
  out = static_cast<int>(v);
  return true;
}

static inline bool parseDoubleStrict(const std::string& tok, double& out) {
  const char* p = tok.c_str();
  char* end = nullptr;
  double v = std::strtod(p, &end);
  if (end == p || *end != '\0') return false;
  out = v;
  return true;
}

static inline std::vector<std::string> splitWs(const std::string& s) {
  std::vector<std::string> out;
  std::istringstream iss(s);
  std::string w;
  while (iss >> w) out.push_back(w);
  return out;
}

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

static inline std::string extractStyleHint(const std::string& headerLine) {
  size_t hash = headerLine.find('#');
  if (hash == std::string::npos) return std::string();
  std::string tail = headerLine.substr(hash + 1);
  trimInplace(tail);
  // Tail might be like: "full" or "charge" or "molecular"
  // Strip trailing words like "atoms" if present (rare)
  return tail;
}

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
void writeLammpsDump(std::ostream& out, long timestep, const std::vector<Atom>& atoms, const std::vector<int> atomTypes, const Box& box) {
  out << "ITEM: TIMESTEP\n" << timestep << "\n";
  out << "ITEM: NUMBER OF ATOMS\n" << atoms.size() << "\n";
  out << "ITEM: BOX BOUNDS pp pp pp\n";
  out << std::setprecision(10) << std::fixed;
  double bl=0;
  out << bl << " " << box.x << "\n";
  out << bl << " " << box.y << "\n";
  out << bl << " " << box.z << "\n";

  out << "ITEM: ATOMS id type x y z\n";
  out << std::setprecision(10) << std::fixed;
  for( int atom=0; atom< atoms.size(); atom++ ) {
    out << atom+1 << " " << atomTypes[atom] << " " << atoms[atom].x << " " <<  atoms[atom].y << " " << atoms[atom].z << "\n"; 
  }
}


}
