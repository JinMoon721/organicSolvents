#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>

#include <ioTraj/atom.hpp>
using namespace ioTraj;

// ----------------- namespace angles -------------------------------
namespace angles {


struct Angle { 
  double zloc;  // z location of terminal carbon
  double ccCos; // angle from +z axis of terminal cc bond
  double chCos; // angle from +z axis of the mean terminal ch bonds
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

enum class Key { DES, DBS, DEO, DEA, Unknown };
Key toKey(std::string s); 

Targets findTargets(const std::string& s, bool left); 
double massFromType(const std::vector<std::pair<int, double>>& type2mass, int type, double def = -1.0);
double comZ(const std::vector<Atom>& atoms, const std::vector<int>& atomTypes, const std::vector<std::pair<int, double>>& type2mass);
double comZPBC(const std::vector<Atom>& atoms, const std::vector<int>& atomTypes, const std::vector<std::pair<int, double>>& type2mass, double Lz);

Angle computeTerminalAngle( Targets& target, const std::vector<Atom>& atoms, const Box& box);

double applyPBC(double x, double box);
double distance(const Atom& a, const Atom& b, const Box& box);
double wrapPBC(double x, double L);

} // -----------namespace angles


