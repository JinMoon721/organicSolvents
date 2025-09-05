#pragma once
#include <vector>
#include <string>
#include <ioTraj/atom.hpp>
#include <fstream> // ifs

namespace ioTraj {

static inline uint32_t bswap32(uint32_t x);
static inline uint64_t bswap64(uint64_t x);

template <typename T>
static inline void byteswap_inplace(T& v);

template <typename T>
static inline void maybe_swap_buffer(T* ptr, size_t n, bool need_swap);

// Read one Fortran unformatted record: [int32 len] [payload] [int32 len]
static std::vector<char> read_fortran_record(std::ifstream& f, bool& need_swap);


class DCDReader {
public:
  explicit DCDReader(const std::string& path);

  int natoms() const;
  int nframes() const;
  bool has_unitcell() const;
  bool has_fixed_atoms() const;

  std::vector<Frame> read_all();

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

  //public methods
  void parse_header();
  void read_unit_cell(Box& box); 
  void read_xyz_full(std::vector<Atom>& atoms);
  void read_xyz_free(std::vector<Atom>& atoms);

  static inline void parse_cell_record(const std::vector<char>& rec, bool need_swap, Box& box);
  void read_one_frame(Frame& fr, int frame_index);
};


} 
