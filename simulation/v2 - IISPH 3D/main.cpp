#define __DEBUG1__
//#define __DEBUG2__
//#define __DEBUG3__
//#define __DEBUG4__
//#define __DEBUG5__ //For CFL conditions.

#include <iostream>
#include <fstream>
#include "IISPH_solver.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

string fileOutput = "liquidPointCloud";

//Simulation parameters
const Real solvNu = 0.02;
const Real solvH = 0.5;
const Real solvDensity = 3e3;
const Vec3f solvG = Vec3f(0, 0, -9.8);
const Real solvInitP = 0.5;
const Real solvOmega = 0.3;
const Real solvPressureError = 0.99;

const int nbConsecutiveSteps = 5;
const Real dt = 0.008f;
const Real fileDt = nbConsecutiveSteps * dt;
const Real cflNumber = 0.4f;

long unsigned int timesteps = 100;
const Vec3f gridRes = Vec3f(20, 20, 25);
const Vec3f initShift = Vec3f(2, 2, 10);
const Vec3f initBlock = Vec3f(12, 12, 8);

IisphSolver solver(dt, solvNu, solvH, solvDensity, solvG, solvInitP, solvOmega, solvPressureError);


int main(int argc, char **argv) {

    // if openmp is enabled, we print the number of threads used
#ifdef _OPENMP
    cout << "OpenMP enabled, using " << omp_get_max_threads() << " threads" << endl;
#else
    cout << "OpenMP not enabled, using only 1 thread" << endl;
#endif

    try {
        solver.initScene(gridRes, initShift, initBlock);
    } catch (length_error &e) {
        cout << e.what() << endl;
        return 0;
    }

    ofstream file;
    stringstream fpath;
    int fileNum = 1;
    do {
        fpath.str(string());
        fpath << "./" << fileOutput << "_" << setw(3) << setfill('0') << fileNum++ << ".txt";
    } while (ifstream(fpath.str()).is_open());

    const unsigned long int nbFluidPart = initBlock.x * initBlock.y * initBlock.z * 8;
    const int nbWallPart = solver.wallParticleCount();
    vector<Vec3f> partPos((timesteps + 1) * nbFluidPart, Vec3f(0));

    for (unsigned long int i = 0; i < nbFluidPart; ++i) {
        partPos[i] = solver.position(i + nbWallPart);
    }
    
    for (unsigned int t = 1; t < timesteps; ++t) {
#ifdef __DEBUG1__
        cout << "Step number " << t << endl;
#endif
        Real timeStepDt = dt;
        Real timeElapsed = 0.;
        while (timeElapsed < fileDt) {
            solver.update(timeStepDt);

            timeElapsed += timeStepDt;
            Real cflCriterion = cflNumber * solver.getKernel().supportRadius()/solver.maxVelocity().length();
            timeStepDt = min(min(dt, cflCriterion), fileDt - timeElapsed);

        }

        for (unsigned long int i = 0; i < nbFluidPart; ++i) {
            partPos[t * nbFluidPart + i] = solver.position(i + nbWallPart);
        }
    }

    cout << "End of the simulation, saving data..." << endl;

    file.open(fpath.str());
    file << gridRes.x << " " << gridRes.y << " " << gridRes.z << "\n";
    file << dt * nbConsecutiveSteps << " " << timesteps + 1 << "\n";

    file << nbFluidPart << "\n";

    for (auto part: partPos) {
        file << part.x << " " << part.y << " " << part.z << "\n";
    }

    cout << " > Quit" << endl;
    file.close();

    return 0;
}