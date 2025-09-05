#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <ioTraj/atom.hpp>


namespace stats {

struct Histogram {
  std::vector<double> edges;
  std::vector<double> counts;
  double binWidth(std::size_t i=0) const {
    return edges.size() > 1? (edges[i+1]-edges[i]) : 0.0;
  }
};

Histogram makeHistogram(const std::vector<double>& data, size_t nbins, double minEdge, double maxEdge);

Histogram makeHistogramAuto(const std::vector<double>& data, size_t nbins);
std::vector<double> normalizePDF(const Histogram& h);

void printHistogram( const Histogram& h, std::ostream& os = std::cout);

std::vector<double> histPDFatEdges(const Histogram& h);

double linInterp(double x0, double y0, double x1, double y1, double x);

double integratePDF(const Histogram& h, double xi, double xf);

int binIndexFromEdges(double x, const std::vector<double>& edges);

struct momentStats {
  std::vector<double> mean, se;
};

struct Hist2D {
  std::vector<double> zedges, aedges;
  std::vector<double> counts;
  int nz = 0, na = 0;
  double totalSamples = 0.0;

  double& at(int iz, int ia) { return counts[size_t(iz) * size_t(na) + size_t(ia)];}
  double  at(int iz, int ia) const { return counts[size_t(iz) * size_t(na) + size_t(ia)];}

  double zwidth(int iz ) const { return zedges[iz+1] - zedges[iz];}
  double awidth(int ia ) const { return aedges[ia+1] - aedges[ia];}

  double zcenter(int iz ) const { return 0.5*(zedges[iz+1] + zedges[iz]);}
  double acenter(int ia ) const { return 0.5*(aedges[ia+1] + aedges[ia]);}

  void normalizeMass();
  void normalizeDensity();
  void printHist(const std::string& path) const;
  std::vector<double> aSliceDensity(double zq) const;
  void printZsliceHist(double zq, std::ostream& os = std::cout,int precision=10 ) const;
  std::vector<double> kthMomentByZ(int k) const;
  momentStats kthMomentStatsByZ(int k, const std::vector<double>& raw) const;
  void printMoments(std::ostream& os, momentStats& stat1, momentStats& stat3, int precision=10 ) const;
};



double swrapPBC(double x, double L);
std::pair<double, double> minmaxComponent(const std::vector<double>& a);
Hist2D makeHist2D(const std::vector<double>& dataz, const std::vector<double>& dataa, int nz, int na);

double wrapMinImage( double x, double L );

Hist2D makeHist2DSymZ(const std::vector<double>& dataz, const std::vector<double>& dataa, double center, double Lz, double zBinWidth, int na);

void hist2DAccumulate(Hist2D& H, const std::vector<double>& dataz, const std::vector<double>& dataa, double weight = 1.0);

void hist2DAdd(Hist2D& H, double z, double a, double weight=1.0);

// ---- Make an empty histogram with existing edges (no auto-ranging) ----
Hist2D makeEmptyHistWithEdges(const std::vector<double>& zEdges,
                              const std::vector<double>& aEdges);

// ---- raw-counts kth moment of 'a' per z-row (no density; uses bin centers) ----
static std::vector<double> kthMomentByZFromCounts(const Hist2D& H, int k);

struct BlockMomentStats {
  std::vector<double> mean;  // size = nz
  std::vector<double> se;    // size = nz
};

BlockMomentStats blockAverageMomentByZ(
    const std::vector<std::vector<double>>& framesz, // frames[t] = angles for that frame
    const std::vector<std::vector<double>>& framesa, // frames[t] = angles for that frame
    int k,
    int B,
    const std::vector<double>& zEdges,
    const std::vector<double>& aEdges);


}
