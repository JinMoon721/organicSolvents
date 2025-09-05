#include <ioTraj/atom.hpp>
#include <angles/angles.hpp>
#include <ioTraj/dcdreader.hpp>
#include <ioTraj/lammpsreader.hpp>
#include <stats/stats.hpp>

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace angles;
using namespace ioTraj;
using namespace stats;

int main(int argc, char** argv) {
  if (argc != 6 ) {
    std::cerr << "Error: Num Args not Match\n" ;
    std::cerr << "Usage : " << argv[0] << "dumpDir targetMol numMols timestep(ps) eqtime(ns)\n";
    return 1;
  }
  std::string dump = std::string(argv[1])  + argv[2] + ".dcd";
  int nummols = std::stoi(argv[3]);
  double timestep = std::stof(argv[4]);
  double eqtime = std::stof(argv[5]);

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
  std::string lammpsInput = std::string("../../") + argv[2] +"/system.data";
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
    double zc = comZPBC(frame, atomTypes, type2mass, box.z); 
    for ( auto& atom : frame ) atom.z = wrapPBC(atom.z - zc + box.z/2.0, box.z);
    // second, find geometric center of mass, shift box 
    double zc2 = comZ(frame, atomTypes, type2mass);
    for ( auto& atom : frame ) atom.z = wrapPBC(atom.z - zc2 + box.z/2.0, box.z);
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

  /// makeHist2DSymZ(std::vector<Angle> data, double center, double Lz, double zBinWidth, int na)

  std::vector<double> termZ, termAngle;
  termZ.reserve(angles.size());
  termAngle.reserve(angles.size());

  std::transform(angles.begin(), angles.end(), std::back_inserter(termZ), [](const Angle& a) {return a.zloc;});
  std::transform(angles.begin(), angles.end(), std::back_inserter(termAngle), [](const Angle& a) {return a.ccCos;});

  auto H =  makeHist2DSymZ(termZ, termAngle, 100, 200, 2, 50);
  hist2DAccumulate(H, termZ, termAngle);
  std::vector<double> raw = H.counts; 
  H.normalizeDensity();
  H.printHist("2dhist.dat");
  std::ofstream output("slice.dat");
  H.printZsliceHist(120, output);

  // generate terminal carbon density profile
  std::vector<double> rho(angles.size());
  std::transform(angles.begin(), angles.end(), rho.begin(), [](const Angle& a){ return static_cast<double>(a.zloc);});
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
    std::vector<std::vector<double>> trajZ(numsnap, std::vector<double>(nummols*2));
    std::vector<std::vector<double>> trajA(numsnap, std::vector<double>(nummols*2));
    for( size_t time =0; time < numsnap; time++ ) {
      for ( int mol =0; mol< nummols*2; mol++ ) {
        trajZ[time][mol] = angles[time * (nummols*2) + mol].zloc;
        trajA[time][mol] = angles[time * (nummols*2) + mol].ccCos;
      }
    }
    int B =20;
    auto stats1 = blockAverageMomentByZ(trajZ, trajA,  1, B, H.zedges, H.aedges);
    auto stats3 = blockAverageMomentByZ(trajZ, trajA,  3, B, H.zedges, H.aedges);
    std::ofstream output3("BAmoments.dat");
    H.printMoments(output3, *(momentStats*)&stats1, *(momentStats*)&stats3, 8);
  }

  return 0;
}

