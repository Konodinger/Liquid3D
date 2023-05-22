#include <stdio.h>
#include <iostream>
#include <fstream>
#include "IISPH_solver.hpp"

using namespace std;

string fileOutput = "liquidPointCloud";
bool gSolverStop = false;


//Simulation parameters

int nbConsecutiveSteps = 5;
float dt = 0.2f;
long unsigned int timesteps = 100;
const int res_x = 48;
const int res_y = 32;
const int f_height = 8;
const int f_width = 8;

SphSolver gSolver(0.08, 0.5, 1e3, Vec2f(0, -9.8), 0.01, 7.0, 0.5, 0.5, 1.);


int main(int argc, char **argv)
{

  //omp_set_num_threads( 12 );
  ofstream file;
  file.open("./" + fileOutput + ".txt");
  file << res_x << " " << res_y << " 1\n";
  file << dt * nbConsecutiveSteps << " " << timesteps << "\n";


  gSolver.initScene(res_x, res_y, f_height, f_width);
  int nbWallPart = gSolver.wallParticleCount();
  file << nbWallPart << "\n";
  for (int i = 0; i < nbWallPart; ++i) {
    file << gSolver.position(i).x << " " << gSolver.position(i).y << " 0.\n";
  }

  int nbPart = gSolver.particleCount();
  file << nbPart - nbWallPart << "\n";

  long unsigned int t = 0;
  while(t < timesteps) {
    if (!gSolverStop){
      for(int i=0; i<nbConsecutiveSteps; ++i) gSolver.update();

      for (int i = nbWallPart; i < nbPart; ++i) {
        file << gSolver.position(i).x << " " << gSolver.position(i).y << " 0.\n";
      }

      t++;
    }
  }
  cout << " > Quit" << endl;
  file.close();
  return 0;
}