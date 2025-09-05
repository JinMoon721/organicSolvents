#pragma once

#include <vector>
#include <unordered_map>
#include <string>

namespace ioTraj {

struct Atom {
  double x, y, z;
};

struct Box {
  double x{0}, y{0}, z{0}; // orthorhombic lengths
};

struct Frame {
  std::vector<Atom> atoms;
  Box box;
};

struct AtomRow {
  int index = 0;   // LAMMPS atom id (1-based)
  int mol   = 0;   // molecule id (0 if missing)
  int type  = 0;   // atom type
  double charge = 0.0; // charge (0 if missing)
};

struct LammpsData {
  std::vector<AtomRow> atoms;                 // all Atoms rows we read
  std::unordered_map<int, double> masses;     // type -> mass (from Masses)
  std::string atomsStyleHint;                 // e.g., "full", "atomic", "charge", if seen
};

}
