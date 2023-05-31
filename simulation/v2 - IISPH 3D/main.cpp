//#define __OPEN_MP__
#define __DEBUG1__
#define __DEBUG2__
#define __DEBUG3__
//#define __DEBUG4__

#include <stdio.h>
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
const Vec3f solvG = Vec3f(0, -9.8, 0);
const Real solvInitP = 0.5  ;
const Real solvOmega = 0.3;
const Real solvPressureError = 0.99;

int nbConsecutiveSteps = 5;
const Real dt = 0.005f;
long unsigned int timesteps = 200;
const int res_x = 48;
const int res_y = 32;
const int res_z = 48;
const int f_length = 8;
const int f_height = 8;
const int f_width = 8;

IisphSolver solver(dt, solvNu, solvH, solvDensity, solvG, solvInitP, solvOmega, solvPressureError);


int main(int argc, char **argv)
{
#ifdef __OPEN_MP__
  omp_set_num_threads(6);
#endif

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
  file << res_x << " " << res_y << " " << res_z << "\n";
  file << dt * nbConsecutiveSteps << " " << timesteps + 1 << "\n";


  solver.initScene(res_x, res_y, res_z, f_length, f_height, f_width);
  int nbWallPart = solver.wallParticleCount();

  // Next part print wall particles. It is currently removed because unused.
  /*file << nbWallPart << "\n";
  for (int i = 0; i < nbWallPart; ++i) {
  file << solver.position(i).x << " " << solver.position(i).y << " " << solver.position(i).z << "\n";
  }*/ 

  int nbPart = solver.particleCount();
  file << nbPart - nbWallPart << "\n";


  long unsigned int t = 0;
  for (int i = nbWallPart; i < nbPart; ++i) {
    file << solver.position(i).x << " " << solver.position(i).y << " " << solver.position(i).z << "\n";
  }

  while(t < timesteps) {
    #ifdef __DEBUG1__
    cout << "Step number " << t + 1 << "\n";
    #endif
    if (!solverStop){
      for(int i=0; i<nbConsecutiveSteps; ++i) solver.update();

      for (int i = nbWallPart; i < nbPart; ++i) {
        file << solver.position(i).x << " " << solver.position(i).y << " " << solver.position(i).z << "\n";
      }

      t++;
    }
  }
  cout << " > Quit" << endl;
  file.close();
  return 0;
}