//#define __OPEN_MP__
#define __DEBUG1__
//#define __DEBUG2__
//#define __DEBUG3__
//#define __DEBUG4__

#include <iostream>
#include <fstream>
#include "IISPH_solver.hpp"

#ifdef __OPEN_MP__
#include <omp.h>
#endif

string fileOutput = "liquidPointCloud";

//Simulation parameters
const Real solvNu = 0.08;
const Real solvH = 0.5;
const Real solvDensity = 3e3;
const Vec3f solvG = Vec3f(0, 0, -9.8);
const Real solvInitP = 0.5;
const Real solvOmega = 0.3;
const Real solvPressureError = 0.99;

int nbConsecutiveSteps = 5;
const Real dt = 0.01f;
long unsigned int timesteps = 500;
const Vec3f gridRes = Vec3f(20, 20, 25);
const Vec3f initShift = Vec3f(0, 0, 10);
const Vec3f initBlock = Vec3f(16, 16, 8);

IisphSolver solver(dt, solvNu, solvH, solvDensity, solvG, solvInitP, solvOmega, solvPressureError);


int main(int argc, char **argv) {
#ifdef __OPEN_MP__
    omp_set_num_threads(6);
#endif

  try {
    solver.initScene(gridRes, initShift, initBlock);
  } catch (length_error& e) {
    cout << e.what() << endl;
    return 0;
  }

  ofstream file;
  std::stringstream fpath;
  int fileNum = 1;
  while(true) {
    fpath.str(string());
    fpath << "./" << fileOutput << "_" << std::setw(3) << std::setfill('0') << fileNum++ << ".txt";
    if (!ifstream(fpath.str()).is_open()) {
      break;
    }
  }

  const int nbFluidPart = initBlock.x*initBlock.y*initBlock.z*8;
  const int nbWallPart = solver.wallParticleCount();
  vector<Vec3f> partPos((timesteps + 1)*nbFluidPart, Vec3f(0));

  for (int i = 0; i < nbFluidPart; ++i) {
    partPos[i] = solver.position(i + nbWallPart);
    }

  for (long unsigned int t = 1; t <= timesteps; ++t) {
    #ifdef __DEBUG1__
    cout << "Step number " << t << "\n";
    #endif
    for(int i=0; i<nbConsecutiveSteps; ++i) solver.update();

    for (int i = 0; i < nbFluidPart; ++i) {
      partPos[t*nbFluidPart + i] = solver.position(i + nbWallPart);
    }
  }

  cout << "End of the simulation, saving data..." << endl;

  file.open(fpath.str());
  file << gridRes.x << " " << gridRes.y << " " << gridRes.z << "\n";
  file << dt * nbConsecutiveSteps << " " << timesteps + 1 << "\n";

  // Next part print wall particles. It is currently removed because unused.
  /*file << nbWallPart << "\n";
  for (int i = 0; i < nbWallPart; ++i) {
  file << solver.position(i).x << " " << solver.position(i).y << " " << solver.position(i).z << "\n";
  }*/ 

  file << nbFluidPart << "\n";

  for (auto part : partPos) {
    file << part.x << " " << part.y << " " << part.z << "\n";
  }

  cout << " > Quit" << endl;
  file.close();

  return 0;
}