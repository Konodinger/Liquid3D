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

using namespace std;

string fileOutput = "liquidPointCloud";
bool solverStop = false;


//Simulation parameters
const Real solvNu = 0.08;
const Real solvH = 0.5;
const Real solvDensity = 3e3;
const Vec3f solvG = Vec3f(0, 0, -9.8);
const Real solvInitP = 0.5  ;
const Real solvOmega = 0.3;
const Real solvPressureError = 0.99;

int nbConsecutiveSteps = 5;
const Real dt = 0.005f;
long unsigned int timesteps = 500;
const Vec3f gridRes = Vec3f(32, 32, 32);
const Vec3f initShift = Vec3f(0, 0, 10);
const Vec3f initBlock = Vec3f(16, 8, 8);

IisphSolver solver(dt, solvNu, solvH, solvDensity, solvG, solvInitP, solvOmega, solvPressureError);


int main(int argc, char **argv)
{
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

  file.open(fpath.str());
  file << gridRes.x << " " << gridRes.y << " " << gridRes.z << "\n";
  file << dt * nbConsecutiveSteps << " " << timesteps + 1 << "\n";
  
  int nbWallPart = solver.wallParticleCount();

  // Next part print wall particles. It is currently removed because unused.
  /*file << nbWallPart << "\n";
  for (int i = 0; i < nbWallPart; ++i) {
  file << solver.position(i).x << " " << solver.position(i).y << " " << solver.position(i).z << "\n";
  }*/ 

  int nbPart = solver.particleCount();
  file << nbPart - nbWallPart << "\n";


  for (int i = nbWallPart; i < nbPart; ++i) {
    file << solver.position(i).x << " " << solver.position(i).y << " " << solver.position(i).z << "\n";
  }

  for (long unsigned int t = 0; t < timesteps; ++t) {
    #ifdef __DEBUG1__
    cout << "Step number " << t + 1 << "\n";
    #endif
    if (!solverStop){
      for(int i=0; i<nbConsecutiveSteps; ++i) solver.update();

      for (int i = nbWallPart; i < nbPart; ++i) {
        file << solver.position(i).x << " " << solver.position(i).y << " " << solver.position(i).z << "\n";
      }

    }
  }
  cout << " > Quit" << endl;
  file.close();
  return 0;
}