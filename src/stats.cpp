#include <ioTraj/atom.hpp>
#include <stats/stats.hpp>
#include <vector>
#include <cmath> 
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <stdexcept>
#include <limits>


namespace stats{

Histogram makeHistogram(const std::vector<double>& data, size_t nbins, double minEdge, double maxEdge){
  Histogram h;
  h.edges.resize(nbins + 1);
  h.counts.assign(nbins, 0.0);

  if (nbins == 0) return h;
  if (maxEdge == minEdge) maxEdge = minEdge + 1;

  const double width = (maxEdge - minEdge)/ static_cast<double>(nbins);
  for( size_t i =0; i<= nbins; i++ ) h.edges[i] = minEdge + i* width;

  for(size_t k=0; k<data.size(); k++) {
    double x = data[k];
    long idx = static_cast<long>(std::floor( (x-minEdge) / width));
    if (idx < 0 ) idx =0;
    if (idx > static_cast<long>(nbins)) idx = static_cast<long>(nbins)-1;
    h.counts[static_cast<size_t>(idx)] += 1.0;
  }
  return h;
}

Histogram makeHistogramAuto(const std::vector<double>& data, size_t nbins) {
  if ( data.empty() || nbins == 0) return {};
  auto mm = std::minmax_element(data.begin(), data.end());
  double mn = *mm.first;
  double mx = *mm.second;
  if ( mx == mn ) {mn -= 0.5; mx += 0.5;}

  return makeHistogram(data, nbins, mn, mx);
}
std::vector<double> normalizePDF(const Histogram& h) {
  double N = std::accumulate(h.counts.begin(), h.counts.end(), 0.0);
  std::vector<double> pdf(h.counts.size(), 0.0);
  if ( N > 0 && h.edges.size() >=2) {
    for (size_t i =0; i< h.counts.size(); i++) {
      double bw = h.edges[i+1] - h.edges[i];
      if (bw > 0 ) pdf[i] = h.counts[i] / (N*bw);
    }
  }
  return pdf;
}

void printHistogram( const Histogram& h, std::ostream& os ) {
  std::vector<double> pdf = normalizePDF(h);
  os.setf(std::ios::fixed);
  os << std::setprecision(10);
  int nr = h.edges.size();
  int half = nr/2;
  for( size_t i =half; i< nr; i++ ) {
    double x = 0.5 * (h.edges[i] + h.edges[i+1]);
    os << x << "\t" << (pdf[i] + pdf[2*half-i-1])/2.0 << "\n";
  }
}

std::vector<double> histPDFatEdges(const Histogram& h) {
  size_t nb = h.counts.size();
  std::vector<double> y; 
  y.assign(nb+1, 0.0);

  double N = std::accumulate(h.counts.begin(), h.counts.end(), 0.0);
  if ( N<=0.0 || h.edges.size() != nb+1 ) return y;

  for (size_t i =0; i< nb ; i++){
    double bw = h.edges[i+1] - h.edges[i];
    y[i] = (bw > 0.0) ? (h.counts[i] / N / bw) : 0.0;
  }
  y[nb] = y[nb - (nb > 0 ? 1: 0)];
  return y;
}

double linInterp(double x0, double y0, double x1, double y1, double x) {
  double dx = x1-x0;
  if (dx==0.0) return 0.5*(y0+y1);
  double t= (x-x0)/dx;
  return y0 + t* (y1-y0);
}

double integratePDF(const Histogram& h, double xi, double xf) {
  if (h.edges.size() < 2 || h.counts.size() +1 != h.edges.size()) return 0.0;
  if (xi > xf ) std::swap(xi, xf);

  const std::vector<double>& x = h.edges;
  std::vector<double> y = histPDFatEdges(h);

  double a = x.front();
  double b = x.back();
  if ( xf <= a || xi >= b) return 0.0;
  xi = std::max(xi, a);
  xf = std::min(xf, b);

  std::vector<double>::const_iterator itL = std::upper_bound(x.begin(), x.end(), xi);
  std::vector<double>::const_iterator itR = std::upper_bound(x.begin(), x.end(), xf);
  size_t j = (itL == x.begin() ? 0 : static_cast<size_t>((itL - x.begin()) -1));
  size_t k = (itR == x.begin() ? 0 : static_cast<size_t>((itR - x.begin()) -1));
  if( j >= x.size() -1) j = x.size() -2;
  if( k >= x.size() -1) k= x.size() -2;
  double y_xi = linInterp(x[j], y[j], x[j+1], y[j+1], xi);
  double y_xf = linInterp(x[k], y[k], x[k+1], y[k+1], xf);
  if ( j== k) { 
    return 0.5 * (y_xi + y_xf) * (xf - xi);
  }

  double area = 0.0;
  area += 0.5 * (y_xi + y[j+1]) * (x[j+1] - xi);

  for ( size_t m = j+1; m< k; m++ ) {
    area += 0.5 * (y[m]+y[m+1]) * (x[m+1] - x[m]);
  }
  area += 0.5 * (y[k] + y_xf) * (xf - x[k]);
  return area;
}


int binIndexFromEdges(double x, const std::vector<double>& edges){
  if( x< edges.front() || x > edges.back() ) return -1;
  if ( x == edges.back()) return int(edges.size())-2;
  std::vector<double>::const_iterator it = std::upper_bound(edges.begin(), edges.end(), x);
  int idx = int (it - edges.begin()) -1;
  if (idx < 0 || idx >= int(edges.size()) -1 ) return -1;
  return idx;
}


void Hist2D::normalizeMass() {
  if ( totalSamples <= 0.0) return;
  for (size_t i=0; i< counts.size(); i++) counts[i] /= totalSamples;
}

void Hist2D::normalizeDensity() {
  if (totalSamples <= 0.0) return;
  for (int iz =0 ; iz < nz; iz++) {
    for (int ia=0; ia < na; ia++ ) {
      double area = zwidth(iz) * awidth(ia);
      if (area > 0) {
        at(iz, ia) /= (totalSamples * area);
      }
    }
  }
}

void Hist2D::printHist(const std::string& path) const {
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

std::vector<double> Hist2D::aSliceDensity(double zq) const {
  int iz = binIndexFromEdges(zq, zedges);
  if ( iz< 0) throw std::runtime_error("z value out of range");
  std::vector<double> pdf (na, 0.0);
  double rowSum = 0.0;
  for (int ia = 0; ia< na; ia++) rowSum += at(iz,ia);
  if (rowSum > 0.0) {
    for( int ia=0; ia < na ; ia++ ) {
      double da = awidth(ia);
      if (da > 0.0) pdf[ia] = at(iz, ia) / rowSum / da;
    }
  }
  return pdf;
}

void Hist2D::printZsliceHist(double zq, std::ostream& os ,int precision ) const {
  int iz = binIndexFromEdges(zq, zedges);
  if( iz<0) throw std::runtime_error("z value is out of range");
  std::vector<double> v = aSliceDensity(zq);
  os.setf(std::ios::fixed);
  os << std::setprecision(precision);
  os << "# z-bin center = " << zcenter(iz) << " (query z = " << zq << ")\n";
  os << "# angle center, pdf\n";
  for( size_t i =0; i< na; i++ ) {
    os << acenter(i) << "\t" << v[i] << "\n";
  }
}

std::vector<double> Hist2D::kthMomentByZ(int k) const {
  std::vector<double> out(size_t(nz), std::numeric_limits<double>::quiet_NaN());
  for( int iz = 0; iz < nz ; iz ++ ) {
    long double num = 0.0L, den = 0.0L;
    for( int ia = 0; ia < na; ia++){
      const double ac = acenter(ia);
      long double w = at(iz, ia) * awidth(ia);
      if ( w <= 0.0L) continue;
      num += std::pow(ac, k) * w;
      den += w;
    }
    if (den > 0.0L) out[size_t(iz)] = double(num/den);
  }
  return out;
}

momentStats Hist2D::kthMomentStatsByZ(int k, const std::vector<double>& raw) const {
  if (k<0) throw std::runtime_error("k must be >= 0");
  if (zedges.size() != nz +1 || aedges.size() != na+1 || raw.size() != nz*na ) {
    throw std::runtime_error("size mismatch");
  }
  momentStats out;
  out.mean.assign(nz, std::numeric_limits<double>::quiet_NaN());
  out.se.assign(nz, std::numeric_limits<double>::quiet_NaN());

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
      out.mean[size_t(iz)] = double(mean);
      if ( N > 1.0L) {
        const long double varY = (sumY2 - N*mean*mean) / (N-1.0L);
        out.se[size_t(iz)] = (varY > 0.0L ) ? double(std::sqrt(varY/N)) : 0.0;
      } else {
        out.se[size_t(iz)] = std::numeric_limits<double>::infinity();
      }
    }
  }
  return out;
}

void Hist2D::printMoments(std::ostream& os, momentStats& stat1, momentStats& stat3, int precision ) const {
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




std::pair<double, double> minmaxComponent(const std::vector<double>& a) {
  double lo = std::numeric_limits<double>::infinity();
  double hi = -std::numeric_limits<double>::infinity();
  for( size_t i =0; i<a.size(); i++ ) {
    double v=a[i];
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  if (!std::isfinite(lo) || !std::isfinite(hi)) { lo = 0.0; hi = 1.0;}
  if (hi==lo) hi = lo + 1;
  return std::pair<double, double> (lo, hi);
}

Hist2D makeHist2D(const std::vector<double>& dataz, const std::vector<double>& dataa, int nz, int na) {
  if (nz <= 0 || na <= 0) throw std::runtime_error("nz/na must be positive\n");
  double zmin, zmax;
  auto zmm = minmaxComponent(dataz);
  zmin = zmm.first; zmax = zmm.second;

  double amin, amax;
  auto amm = minmaxComponent(dataa);
  amin = amm.first; amax = amm.second;


  Hist2D H;
  H.nz = nz; H.na = na;
  H.zedges.resize(size_t(nz)+1);
  H.aedges.resize(size_t(na)+1);
  H.counts.assign(size_t(nz)*size_t(na), 0.0);
  double dz = (zmax-zmin)/double(nz);
  double da = (amax-amin)/double(na);
  for(int i = 0; i<=nz; i++) H.zedges[i] = zmin + dz*double(i);
  for(int j = 0; j<=na; j++) H.aedges[j] = amin + da*double(j);
  return H;
}

double wrapMinImage( double x, double L ) {
  return x - L * std::floor((x+0.5*L)/L);
}

double swrapPBC(double x, double L) { // wrap x into [0, L)
  double y = std::fmod(x, L);
  if (y<0) y+=L;
  return y;
}

Hist2D makeHist2DSymZ(const std::vector<double>& dataz, const std::vector<double>& dataa, double center, double Lz, double zBinWidth, int na) {
  if (!(Lz > 0.0) || !(zBinWidth > 0.0) ||  na <= 0) {
    throw std::runtime_error("nz/na must be positive\n");
  }

  // find max absolute displacement from center
  double zmaxDisp = 0.0;
  for (const auto& a : dataz) {
    double d = wrapMinImage(swrapPBC(a, Lz) - center, Lz);
    double ab = std::fabs(d);
    if (ab > zmaxDisp) zmaxDisp = ab;
  }
  zmaxDisp = std::min(zmaxDisp, 0.5 * Lz);

  // choose bin count per side
  int nSide = std::max(1, int(std::ceil(zmaxDisp / zBinWidth)));
  double R = nSide * zBinWidth;
  int nz = nSide*2;

  std::cout << "max z : " << zmaxDisp << " and nz : " << nz << "\n";

  double amin, amax;
  auto amm = minmaxComponent(dataa);
  amin = amm.first; amax = amm.second;

  Hist2D H;
  H.nz = nz; H.na = na;
  H.zedges.resize(size_t(nz)+1);
  H.aedges.resize(size_t(na)+1);
  H.counts.assign(size_t(nz)*size_t(na), 0.0);
  double dz = zBinWidth;
  double da = (amax-amin)/double(na);
  for(int i = 0; i<=nz; i++) H.zedges[i] = center-R + dz*double(i);
  for(int j = 0; j<=na; j++) H.aedges[j] = amin + da*double(j);
  return H;
}

void hist2DAccumulate(Hist2D& H, const std::vector<double>& dataz, const std::vector<double>& dataa, double weight ) {
  for (size_t k=0; k< dataz.size(); k++ ) {
    double z = double(dataz[k]) ;
    double a = double(dataa[k]);
    if (!std::isfinite(z) || !std::isfinite(a)) continue;
    int iz = binIndexFromEdges(z, H.zedges);
    int ia = binIndexFromEdges(a, H.aedges);
    if (iz >= 0 && ia >= 0) {
      H.at(iz, ia) += weight;
      H.totalSamples += weight;
    }
  }
}

void hist2DAdd(Hist2D& H, double z, double a, double weight) {
  if(!std::isfinite(z) || !std::isfinite(a)) return;
  int iz = binIndexFromEdges(z, H.zedges);
  int ia = binIndexFromEdges(a, H.aedges);
  if (iz >= 0 && ia >= 0) {
    H.at(iz, ia) += weight;
    H.totalSamples += weight;
  }
}


// ---- Make an empty histogram with existing edges (no auto-ranging) ----
Hist2D makeEmptyHistWithEdges(const std::vector<double>& zEdges,
                              const std::vector<double>& aEdges) {
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
static std::vector<double> kthMomentByZFromCounts(const Hist2D& H, int k) {
  std::vector<double> out(size_t(H.nz), std::numeric_limits<double>::quiet_NaN());
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
    if (N > 0.0L) out[size_t(iz)] = double(sumY / N);
  }
  return out;
}


BlockMomentStats blockAverageMomentByZ(
    const std::vector<std::vector<double>>& framesz, // framesz[t] = z for that frame
    const std::vector<std::vector<double>>& framesa, // framesz[t] = angle for that frame
    int k,
    int B,
    const std::vector<double>& zEdges,
    const std::vector<double>& aEdges)
{
  if (B < 2) throw std::runtime_error("blockAverageMomentByZ: B must be >= 2");
  const int F = int(framesz.size());
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
      hist2DAccumulate(HB[b], framesz[size_t(f)], framesa[size_t(f)], /*weight=*/1.0f);
    }
  }

  const int nz = int(zEdges.size()) - 1;
  BlockMomentStats out;
  out.mean.assign(size_t(nz), std::numeric_limits<double>::quiet_NaN());
  out.se  .assign(size_t(nz), std::numeric_limits<double>::quiet_NaN());

  std::vector<double> mb(size_t(B), std::numeric_limits<double>::quiet_NaN());
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
        mb[size_t(b)]  = (N > 0.0L) ? double(sumY / N) : std::numeric_limits<double>::quiet_NaN();
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
      out.mean[size_t(iz)] = double(mu);

      if (M >= 2) {
        long double s2 = 0.0L;
        for (int b = 0; b < B; ++b) if (has[size_t(b)]) {
          const long double d = mb[size_t(b)] - mu;
          s2 += d * d;
        }
        out.se[size_t(iz)] = double(std::sqrt(s2 / (M * (M - 1))));
      } // else leave NaN
    }
  }
  return out;
}
}
