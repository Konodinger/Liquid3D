#define __DEBUG1__
//#define __DEBUG2__
//#define __DEBUG3__
//#define __DEBUG4__
//#define __DEBUG5__  //For CFL conditions.
//#define __PRINT_VELOCITY__

#include <iostream>
#include <fstream>
#include "IISPH_solver.hpp"
#include "SFB_generator.hpp"
#include "utils.hpp"

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

//Foam parameters
const int minBubbleNeighbor = 20;
const int minFoamNeighbor = 6;
const Real kb = 0.3f;
const Real kd = 0.7f;

const int nbConsecutiveSteps = 5;
const Real dt = 0.008f;
const Real fileDt = nbConsecutiveSteps * dt;
const Real cflNumber = 0.4f;

#define DEFAULT_NB_TIMESTEPS 50
#define DEFAULT_INIT_TYPE InitType::SPHERE;

#define DEFAULT_RESOLUTION 1.0f
Vec3f gridRes = Vec3f(20, 20, 25);

IisphSolver solver(dt, solvNu, solvH, solvDensity, solvG, solvInitP, solvOmega, solvPressureError);

sfbSimulation sfbSim(&solver, solvH, dt, minBubbleNeighbor, minFoamNeighbor);


int main(int argc, char **argv) {
    if (cmdOptionExists(argv, argv + argc, "--help") || cmdOptionExists(argv, argv + argc, "-h")) {
        cout << "Usage: ./main [options]" << endl;
        cout << "Options:" << endl;
        cout << "\t--help, -h\t Display this help message" << endl;
        cout << "\t--output <value>, -o <value>\t Set the output file name (default: liquidPointCloud)" << endl;
        cout << "\t--timesteps <value>, -t <value>\t Set the number of timesteps (default: " << DEFAULT_NB_TIMESTEPS
             << ")" << endl;
        cout << "\t--init <value>, -i <value>\t Set the initial particle configuration (default: sphere) (torus, block, sphere)" << endl;
        cout << "\t--resolution <value>, -r <value>\t Set the scale of the scene (default: 1.0)" << endl;

        return 0;
    }

    unsigned long int timesteps = DEFAULT_NB_TIMESTEPS;
    char *timestepsCmdArg = getCmdOption(argv, argv + argc, "--timesteps");
    char *timestepsCmdArgShort = getCmdOption(argv, argv + argc, "-t");
    if (timestepsCmdArg) timesteps = atoi(timestepsCmdArg);
    else if (timestepsCmdArgShort) timesteps = atoi(timestepsCmdArgShort);

    InitType initType = DEFAULT_INIT_TYPE;
    char *initTypeCmdArg = getCmdOption(argv, argv + argc, "--init");
    char *initTypeCmdArgShort = getCmdOption(argv, argv + argc, "-i");
    std::string initTypeStr;
    if (initTypeCmdArg) initTypeStr = initTypeCmdArg;
    else if (initTypeCmdArgShort) initTypeStr = initTypeCmdArgShort;
    if (!initTypeStr.empty()) {
        // sphere
        if (initTypeStr == "sphere") {
            initType = InitType::SPHERE;
        } else if (initTypeStr == "torus") {
            initType = InitType::TORUS;
        } else if (initTypeStr == "block") {
            initType = InitType::BLOCK;
        } else {
            cout << "Unknown initial scene type: " << initTypeStr << endl;
            return 0;
        }
    }

    stringstream fpath;
    ofstream file;

    char *outputCmdArg = getCmdOption(argv, argv + argc, "--output");
    char *outputCmdArgShort = getCmdOption(argv, argv + argc, "-o");
    if (outputCmdArg) {
        fileOutput = outputCmdArg;
        fpath << "./" << fileOutput << ".txt";
    } else if (outputCmdArgShort) {
        fileOutput = outputCmdArgShort;
        fpath << "./" << fileOutput << ".txt";
    } else {
        int fileNum = 1;
        do {
            fpath.str(string());
            fpath << "./" << fileOutput << "_" << setw(3) << setfill('0') << fileNum++ << ".txt";
        } while (ifstream(fpath.str()).is_open());
    }

    char *resolutionCmdArg = getCmdOption(argv, argv + argc, "--resolution");
    char *resolutionCmdArgShort = getCmdOption(argv, argv + argc, "-r");
    Real resolution;
    if (resolutionCmdArg) resolution = atof(resolutionCmdArg);
    else if (resolutionCmdArgShort) resolution = atof(resolutionCmdArgShort);
    else resolution = DEFAULT_RESOLUTION;

    gridRes *= resolution;

    // if openmp is enabled, we print the number of threads used
#ifdef _OPENMP
    cout << "OpenMP enabled, using " << omp_get_max_threads() << " threads" << endl;
#else
    cout << "OpenMP not enabled, using only 1 thread" << endl;
#endif

    try {
        solver.initScene(gridRes, initType);
    } catch (length_error &e) {
        cout << e.what() << endl;
        return 0;
    }

    const unsigned long int nbFluidPart = solver.fluidParticleCount();
    const int nbWallPart = solver.wallParticleCount();
    vector<Vec3f> partPos((timesteps + 1) * nbFluidPart, Vec3f(0));
    vector<Vec3f> partVel((timesteps + 1) * nbFluidPart, Vec3f(0));

    for (unsigned long int i = 0; i < nbFluidPart; ++i) {
        partPos[i] = solver.position(i + nbWallPart);
        partVel[i] = solver.velocity(i + nbWallPart);
    }

    for (long unsigned int t = 1; t <= timesteps; ++t) {
#ifdef __DEBUG1__
        cout << "Step number " << t << endl;
#endif
        Real timeStepDt = dt;
        Real timeElapsed = 0.;
        while (timeElapsed < fileDt) {
            solver.update(timeStepDt);
            sfbSim.sfbStep(timeStepDt);

            timeElapsed += timeStepDt;
            Real cflCriterion = cflNumber * solver.getKernel().supportRadius() / solver.maxVelocity().length();
            timeStepDt = min(min(dt, cflCriterion), fileDt - timeElapsed);

        }

        for (unsigned long int i = 0; i < nbFluidPart; ++i) {
            partPos[t * nbFluidPart + i] = solver.position(i + nbWallPart);
            partVel[t * nbFluidPart + i] = solver.velocity(i + nbWallPart);
        }
    }

    cout << "End of the simulation, saving data..." << endl;

    file.open(fpath.str());
    file << gridRes.x << " " << gridRes.y << " " << gridRes.z << "\n";
    file << dt * nbConsecutiveSteps << " " << timesteps + 1 << "\n";

    file << nbFluidPart << "\n";

    for (unsigned int i = 0; i < (timesteps + 1) * nbFluidPart; ++i) {
        file << partPos[i].x << " " << partPos[i].y << " " << partPos[i].z;
#ifdef __PRINT_VELOCITY__
        << " " << partVel[i].x << " " << partVel[i].y << " " << partVel[i].z;
#endif
        file << "\n";

    }

    cout << " > Quit" << endl;
    file.close();

    return 0;
}