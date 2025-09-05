#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <ioTraj/atom.hpp>

namespace ioTraj {

// Trim leading/trailing whitespace in-place.
static inline void trimInplace(std::string& s);

// Remove everything after a '#' (LAMMPS inline comment).
static inline void stripComment(std::string& s);

// Case-insensitive prefix check.
static inline bool startsWithCi(const std::string& s, const std::string& prefix);

// Strict integer parse (entire token must be an int).
static inline bool parseIntStrict(const std::string& tok, int& out);

// Strict double parse (entire token must be a number).
static inline bool parseDoubleStrict(const std::string& tok, double& out);

// Tokenize a whitespace string into vector<string>.
static inline std::vector<std::string> splitWs(const std::string& s);

// Heuristic: does a cleaned line look like a known section header (besides Atoms/Masses)?
static inline bool isOtherSectionHeader(const std::string& s);

// Try to extract a style hint from a header line like "Atoms # full"
static inline std::string extractStyleHint(const std::string& headerLine);

// Parse one Atoms data line to (index, mol, type, charge) using style rules.
// - If style is known ("full", "molecular", "charge", "atomic"), follow that.
// - If style is empty, infer from token pattern; missing fields default to 0 / 0.0.
static inline bool parseAtomsRow(const std::vector<std::string>& tok,
                                 const std::string& styleHint,
                                 AtomRow& out);


// Parse a LAMMPS data file, collecting Atoms (id,mol,type,charge) and Masses (type->mass).
LammpsData readLammpsData(const std::string& path);

void writeLammpsDump(std::ostream& out, long timestep, const std::vector<Atom>& atoms, const std::vector<int> atomTypes, const Box& box);

}
