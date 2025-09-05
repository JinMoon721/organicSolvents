#include <ioTraj/atom.hpp>
#include <ioTraj/dcdreader.hpp>

#include <vector>
#include <cstring> // memcpy
#include <stdexcept> // runtime_error
#include <fstream> // ios
#include <string>
#include <cstring> // strncmp


namespace ioTraj {

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

// ---- class DCDReader
DCDReader::DCDReader(const std::string& path) : ifs_(path, std::ios::binary) {
    if (!ifs_) throw std::runtime_error("Cannot open file: " + path);
    parse_header();
}
int DCDReader::natoms() const { return natoms_; }
int DCDReader::nframes() const { return nset_; }
bool DCDReader::has_unitcell() const { return has_unitcell_; }
bool DCDReader::has_fixed_atoms() const { return has_fixed_atoms_; }

std::vector<Frame> DCDReader::read_all() {
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

void DCDReader::parse_header() {
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

void DCDReader::read_unit_cell(Box& box) {
  // Some writers use 6 doubles, others 6 floats; accept either.
  auto rec = read_fortran_record(ifs_, need_swap_);
  if (rec.size() == 6*sizeof(double)) {
    double cell[6];
    std::memcpy(cell, rec.data(), 6*sizeof(double));
    maybe_swap_buffer(cell, 6, need_swap_);
    box.x = cell[0]; // A
    box.y = cell[2]; // B
    box.z = cell[4]; // C
  } else if (rec.size() == 6*sizeof(float)) {
    float cellf[6];
    std::memcpy(cellf, rec.data(), 6*sizeof(float));
    maybe_swap_buffer(cellf, 6, need_swap_);
    box.x = cellf[0];
    box.y = cellf[2];
    box.z = cellf[4];
  } else if (rec.size() == 0) {
    // Empty cell block (shouldn't happen): keep zeros
  } else {
    throw std::runtime_error("Unexpected unit cell record size.");
  }
}

void DCDReader::read_xyz_full(std::vector<Atom>& atoms) {
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

void DCDReader::read_xyz_free(std::vector<Atom>& atoms) {
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

inline void DCDReader::parse_cell_record(const std::vector<char>& rec, bool need_swap, Box& box) {
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
  
  box.x = Lx;
  box.y = Ly;
  box.z = Lz;
}


void DCDReader::read_one_frame(Frame& fr, int frame_index) {
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
  
  // Accept either float or double coordinate payloads and convert to double
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
    for (int i = 0; i < natoms_; ++i) fr.atoms[i].x = static_cast<double>(bufD[i]);
  
    std::memcpy(bufD.data(), recy.data(), recy.size());
    maybe_swap_buffer(bufD.data(), natoms_, need_swap_);
    for (int i = 0; i < natoms_; ++i) fr.atoms[i].y = static_cast<double>(bufD[i]);
  
    std::memcpy(bufD.data(), recz.data(), recz.size());
    maybe_swap_buffer(bufD.data(), natoms_, need_swap_);
    for (int i = 0; i < natoms_; ++i) fr.atoms[i].z = static_cast<double>(bufD[i]);
  
  } else {
    throw std::runtime_error(
      "XYZ array size mismatch (full): got {" +
      std::to_string(recx.size()) + "," +
      std::to_string(recy.size()) + "," +
      std::to_string(recz.size()) + "} bytes");
  }
}

}
