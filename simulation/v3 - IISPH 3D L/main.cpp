#define __DEBUG1__
//#define __DEBUG2__
//#define __DEBUG3__
//#define __DEBUG4__
//#define __DEBUG5__  //For CFL conditions.
//#define __DEBUG6__ //For debugging diffuse particle formation.
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

string fileOutputFluid = "liquidPointCloud";
string fileOutputDiffuse = "diffusePointCloud";

//Simulation parameters
#define DEFAULT_VISCOSITY 0.02f
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

#define DEFAULT_USE_LCONFIG false

#define DEFAULT_RESOLUTION 1.0f
Vec3f gridRes = Vec3f(20, 20, 25);

int main(int argc, char **argv) {
    if (cmdOptionExists(argv, argv + argc, "--help") || cmdOptionExists(argv, argv + argc, "-h")) {
        cout << "Usage: ./main [options]" << endl;
        cout << "Options:" << endl;
        cout << "\t--help, -h\t Display this help message" << endl;
        cout << "\t--output <value>, -o <value>\t Set the output file name (default: liquidPointCloud)" << endl;
        cout << "\t--timesteps <value>, -t <value>\t Set the number of timesteps (default: " << DEFAULT_NB_TIMESTEPS
             << ")" << endl;
        cout
                << "\t--init <value>, -i <value>\t Set the initial particle configuration (default: sphere) (torus, block, sphere)"
                << endl;
        cout << "\t--resolution <value>, -r <value>\t Set the resolution of the scene (default: 1.0)" << endl;
        cout << "\t--noFoam\t Disable foam generation" << endl;
        cout << "\t--viscosity <value>, -v <value>\t Set the viscosity of the fluid (default: " << DEFAULT_VISCOSITY
             << ")" << endl;
        cout << "\t--useLConfig\t Use the LConfig file to initialize the scene (default: false)" << endl;

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

    stringstream fpathFluid, fpathSpray, fpathFoam, fpathBubble;
    ofstream file;

    char *outputCmdArg = getCmdOption(argv, argv + argc, "--output");
    char *outputCmdArgShort = getCmdOption(argv, argv + argc, "-o");
    if (outputCmdArg) {
        fileOutputFluid = outputCmdArg;
        fpathFluid << "./" << fileOutputFluid << ".txt";
    } else if (outputCmdArgShort) {
        fileOutputFluid = outputCmdArgShort;
        fpathFluid << "./" << fileOutputFluid << ".txt";
    } else {
        int fileNum = 0;
        do {
            fpathFluid.str(string());
            fpathFluid << "./" << fileOutputFluid << "_" << setw(3) << setfill('0') << ++fileNum << ".txt";
            fpathSpray.str(string());
            fpathSpray << "./" << fileOutputFluid << "Spray_" << setw(3) << setfill('0') << fileNum << ".txt";
            fpathFoam.str(string());
            fpathFoam << "./" << fileOutputFluid << "Foam_" << setw(3) << setfill('0') << fileNum << ".txt";
            fpathBubble.str(string());
            fpathBubble << "./" << fileOutputFluid << "Bubble_" << setw(3) << setfill('0') << fileNum << ".txt";
        } while (ifstream(fpathFluid.str()).is_open());
    }

    char *resolutionCmdArg = getCmdOption(argv, argv + argc, "--resolution");
    char *resolutionCmdArgShort = getCmdOption(argv, argv + argc, "-r");
    Real resolution;
    if (resolutionCmdArg) resolution = atof(resolutionCmdArg);
    else if (resolutionCmdArgShort) resolution = atof(resolutionCmdArgShort);
    else resolution = DEFAULT_RESOLUTION;

    gridRes *= resolution;

    bool foamEnabled = !cmdOptionExists(argv, argv + argc, "--noFoam");

    char *viscosityCmdArg = getCmdOption(argv, argv + argc, "--viscosity");
    char *viscosityCmdArgShort = getCmdOption(argv, argv + argc, "-v");
    Real viscosity;
    if (viscosityCmdArg) viscosity = atof(viscosityCmdArg);
    else if (viscosityCmdArgShort) viscosity = atof(viscosityCmdArgShort);
    else viscosity = DEFAULT_VISCOSITY;

    // if openmp is enabled, we print the number of threads used
#ifdef _OPENMP
    cout << "OpenMP enabled, using " << omp_get_max_threads() << " threads" << endl;
#else
    cout << "OpenMP not enabled, using only 1 thread" << endl;
#endif

    IisphSolver solver(dt, viscosity, solvH, solvDensity, solvG, solvInitP, solvOmega, solvPressureError);
    sfbSimulation sfbSim(&solver, solvH, dt, minBubbleNeighbor, minFoamNeighbor);
    solver.scaleGarvity(resolution);

    bool useLConfig = cmdOptionExists(argv, argv + argc, "--useLConfig");

    try {
        solver.initScene(gridRes, initType, useLConfig);
        if (foamEnabled) sfbSim.initScene();
    } catch (length_error &e) {
        cout << e.what() << endl;
        return 0;
    }

    const unsigned long int nbFluidPart = solver.fluidParticleCount();
    const int nbWallPart = solver.wallParticleCount();
    vector<Vec3f> partPos((timesteps + 1) * nbFluidPart, Vec3f(0));
    vector<Vec3f> partVel((timesteps + 1) * nbFluidPart, Vec3f(0));
    vector<vector<Vec3f>> sprayPos(timesteps + 1, vector<Vec3f>(0));
    vector<vector<Vec3f>> foamPos(timesteps + 1, vector<Vec3f>(0));
    vector<vector<Vec3f>> bubblePos(timesteps + 1, vector<Vec3f>(0));

    const list<sfb> *diffuseList = sfbSim.diffuseList();

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
            if (foamEnabled) sfbSim.sfbStep(timeStepDt);

            timeElapsed += timeStepDt;
            Real cflCriterion = cflNumber * solver.getKernel().supportRadius() / solver.maxVelocity().length();
            timeStepDt = min(min(dt, cflCriterion), fileDt - timeElapsed);

        }

        for (unsigned long int i = 0; i < nbFluidPart; ++i) {
            partPos[t * nbFluidPart + i] = solver.position(i + nbWallPart);
            partVel[t * nbFluidPart + i] = solver.velocity(i + nbWallPart);
        }

        if (foamEnabled) {
            vector<Vec3f> sprayPosStep, foamPosStep, bubblePosStep;
            for (const sfb &diffuse: (*diffuseList)) {
                switch (diffuse.nature) {
                    case sfbType::SPRAY :
                        sprayPosStep.push_back(diffuse.position);
                        break;
                    case sfbType::FOAM :
                        foamPosStep.push_back(diffuse.position);
                        break;
                    case sfbType::BUBBLE :
                        bubblePosStep.push_back(diffuse.position);
                        break;
                }
            }

            sprayPos[t] = sprayPosStep;
            foamPos[t] = foamPosStep;
            sprayPos[t] = sprayPosStep;
        }
    }

    cout << "End of the simulation, saving data to " << fpathFluid.str() << endl;

    //Fluid part.
    file.open(fpathFluid.str());
    file << gridRes.x << " " << gridRes.y << " " << gridRes.z << "\n";
    file << dt * nbConsecutiveSteps << " " << timesteps + 1 << "\n";

    for (unsigned int i = 0; i < (timesteps + 1); ++i) {
        file << nbFluidPart << "\n";
        for (unsigned long int j = 0; j < nbFluidPart; ++j) {
            unsigned long int index = i * nbFluidPart + j;
            file << partPos[index].x << " " << partPos[index].y << " " << partPos[index].z;
#ifdef __PRINT_VELOCITY__
            file << " " << partVel[index].x << " " << partVel[index].y << " " << partVel[index].z;
#endif
            file << "\n";
        }

    }
    file.close();

    //Diffuse part.
    if (foamEnabled) {
        file.open(fpathSpray.str());

        file << gridRes.x << " " << gridRes.y << " " << gridRes.z << "\n";
        file << dt * nbConsecutiveSteps << " " << timesteps + 1 << "\n";

        for (unsigned int i = 0; i <= timesteps; ++i) {
            file << sprayPos[i].size() << "\n";
            for (unsigned long long int j = 0; j < sprayPos[i].size(); ++j) {
                file << sprayPos[i][j].x << " " << sprayPos[i][j].y << " " << sprayPos[i][j].z << "\n";
            }
        }

        file.close();

        file.open(fpathFoam.str());

        file << gridRes.x << " " << gridRes.y << " " << gridRes.z << "\n";
        file << dt * nbConsecutiveSteps << " " << timesteps + 1 << "\n";

        for (unsigned int i = 0; i <= timesteps; ++i) {
            file << foamPos[i].size() << "\n";
            for (unsigned long long int j = 0; j < foamPos[i].size(); ++j) {
                file << foamPos[i][j].x << " " << foamPos[i][j].y << " " << foamPos[i][j].z << "\n";
            }
        }

        file.close();

        file.open(fpathBubble.str());

        file << gridRes.x << " " << gridRes.y << " " << gridRes.z << "\n";
        file << dt * nbConsecutiveSteps << " " << timesteps + 1 << "\n";

        for (unsigned int i = 0; i <= timesteps; ++i) {
            file << bubblePos[i].size() << "\n";
            for (unsigned long long int j = 0; j < bubblePos[i].size(); ++j) {
                file << bubblePos[i][j].x << " " << bubblePos[i][j].y << " " << bubblePos[i][j].z << "\n";
            }
        }

        file.close();
    }


    cout << " > Quit" << endl;
    return 0;
}