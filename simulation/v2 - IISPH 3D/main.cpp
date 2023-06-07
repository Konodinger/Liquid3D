//#define __OPEN_MP__
#define __DEBUG1__
//#define __DEBUG2__
//#define __DEBUG3__
//#define __DEBUG4__
//#define __BENCHMARK__

#include <iostream>
#include <fstream>
#include "IISPH_solver.hpp"

#ifdef _OPENMP

#include <omp.h>

#endif

std::string fileOutput = "liquidPointCloud";

//Simulation parameters
const Real solvNu = 0.08;
const Real solvH = 0.5;
const Real solvDensity = 3e3;
const Vec3f solvG = Vec3f(0, 0, -9.8);
const Real solvInitP = 0.5;
const Real solvOmega = 0.3;
const Real solvPressureError = 0.99;

const int nbConsecutiveSteps = 5;
const Real dt = 0.01f;
long unsigned int timesteps = 500;
const Vec3f gridRes = Vec3f(20, 20, 25);
const Vec3f initShift = Vec3f(0, 0, 10);
const Vec3f initBlock = Vec3f(16, 16, 8);

IisphSolver solver(dt, solvNu, solvH, solvDensity, solvG, solvInitP, solvOmega, solvPressureError);


int main(int argc, char **argv) {

    // if openmp is enabled, we print the number of threads used
#ifdef _OPENMP
    std::cout << "OpenMP enabled, using " << omp_get_max_threads() << " threads" << std::endl;
#else
    std::cout << "OpenMP not enabled, using only 1 thread" << std::endl;
#endif

    try {
        solver.initScene(gridRes, initShift, initBlock);
    } catch (length_error &e) {
        std::cout << e.what() << std::endl;
        return 0;
    }

    ofstream file;
    std::stringstream fpath;
    int fileNum = 1;
    while (true) {
        fpath.str(std::string());
        fpath << "./" << fileOutput << "_" << std::setw(3) << std::setfill('0') << fileNum++ << ".txt";
        if (!ifstream(fpath.str()).is_open()) {
            break;
        }
    }

    const unsigned long int nbFluidPart = initBlock.x * initBlock.y * initBlock.z * 8;
    const int nbWallPart = solver.wallParticleCount();
    vector<Vec3f> partPos((timesteps + 1) * nbFluidPart, Vec3f(0));

    for (unsigned long int i = 0; i < nbFluidPart; ++i) {
        partPos[i] = solver.position(i + nbWallPart);
    }

    for (long unsigned int t = 1; t <= timesteps; ++t) {
#ifdef __DEBUG1__
        std::cout << "Step number " << t << std::endl;
#endif
        for (unsigned int i = 0; i < nbConsecutiveSteps; ++i) solver.update();

        for (unsigned long int i = 0; i < nbFluidPart; ++i) {
            partPos[t * nbFluidPart + i] = solver.position(i + nbWallPart);
        }
    }

    std::cout << "End of the simulation, saving data..." << std::endl;

    file.open(fpath.str());
    file << gridRes.x << " " << gridRes.y << " " << gridRes.z << "\n";
    file << dt * nbConsecutiveSteps << " " << timesteps + 1 << "\n";

    // Next part print wall particles. It is currently removed because unused.
    /*file << nbWallPart << "\n";
    for (int i = 0; i < nbWallPart; ++i) {
    file << solver.position(i).x << " " << solver.position(i).y << " " << solver.position(i).z << "\n";
    }*/

    file << nbFluidPart << "\n";

    for (auto part: partPos) {
        file << part.x << " " << part.y << " " << part.z << "\n";
    }

    std::cout << " > Quit" << std::endl;
    file.close();

    return 0;
}