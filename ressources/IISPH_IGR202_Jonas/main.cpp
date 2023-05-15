#include <stdio.h>
#include <iostream>
#include "IISPH_solver.cpp"

using namespace std;

float timer = 0.f;
bool gSolverStop = false;

SphSolver gSolver(0.08, 0.5, 1e3, Vec2f(0, -9.8), 0.01, 7.0, 0.5, 0.5, 1.);


int main(int argc, char **argv)
{

  //omp_set_num_threads( 12 );

  gSolver.initScene(48, 32, 16, 16);
  while(true) {
    if (!gSolverStop){
      for(int i=0; i<10; ++i) gSolver.update();
      timer += 0.1f;
    }
  }
  cout << " > Quit" << endl;
  return 0;
}