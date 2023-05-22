//#define __OPEN_MP__
#define __DEBUG1__
//#define __DEBUG2__
//#define __DEBUG3__

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "IISPH_solver.hpp"

#ifdef __OPEN_MP__
#include <omp.h>
#endif

using namespace std;

string fileOutput = "liquidPointCloud";
bool gSolverStop = false;


//Simulation parameters

int nbConsecutiveSteps = 50;
float dt = 0.02f;
long unsigned int timesteps = 20;
const int res_x = 48;
const int res_y = 32;
const int res_z = 48;
const int f_length = 8;
const int f_height = 8;
const int f_width = 8;

SphSolver gSolver(0.08, 0.5, 2.5e3, Vec3f(0, -9.8, 0), 0.01, 7.0, 0.5, 0.5, 0.99);


int main(int argc, char **argv)
{
#ifdef __OPEN_MP__
  omp_set_num_threads(6);
#endif

  ofstream file;
  file.open("./" + fileOutput + ".txt");
  file << res_x << " " << res_y << " " << res_z << "\n";
  file << dt * nbConsecutiveSteps << " " << timesteps + 1 << "\n";


  gSolver.initScene(res_x, res_y, res_z, f_length, f_height, f_width);
  int nbWallPart = gSolver.wallParticleCount();

  // Next part print wall particles. It is currently removed because unused.
  /*file << nbWallPart << "\n";
  for (int i = 0; i < nbWallPart; ++i) {
  file << gSolver.position(i).x << " " << gSolver.position(i).y << " " << gSolver.position(i).z << "\n";
  }*/ 

  int nbPart = gSolver.particleCount();
  file << nbPart - nbWallPart << "\n";


  long unsigned int t = 0;
  for (int i = nbWallPart; i < nbPart; ++i) {
    file << gSolver.position(i).x << " " << gSolver.position(i).y << " " << gSolver.position(i).z << "\n";
  }

  while(t < timesteps) {
    #ifdef __DEBUG1__
    cout << "Step number " << t + 1 << "\n";
    #endif
    if (!gSolverStop){
      for(int i=0; i<nbConsecutiveSteps; ++i) gSolver.update();

      for (int i = nbWallPart; i < nbPart; ++i) {
        file << gSolver.position(i).x << " " << gSolver.position(i).y << " " << gSolver.position(i).z << "\n";
      }

      t++;
    }
  }
  cout << " > Quit" << endl;
  file.close();
  return 0;
}