#include <openvdb/openvdb.h>
#include <openvdb/Types.h>

#include <iostream>
#include <fstream>

#include "rasterize.hpp"
#include "pointGrid.hpp"

#define DEFAULT_PARTICLE_RADIUS 0.4f
#define DEFAULT_VOXEL_SIZE 0.1f
#define DEFAULT_HALF_WIDTH 3.0f

#define DEFAULT_OUTPUT_FILE_NAME "fluid"

int main(int argc, char **argv) {
    if (cmdOptionExists(argv, argv + argc, "-h")) {
        std::cout << "Usage: " << argv[0] << " -f <fileName>\n" << std::endl;

        std::cout << "Options:" << std::endl;
        std::cout << "\t-f <fileName>\tSpecify the file to read" << std::endl;
        std::cout << "\t-o <fileName>\tSpecify the output file name. Default: " << DEFAULT_OUTPUT_FILE_NAME
                  << std::endl;
        std::cout << "\t-h\tPrint this help message" << std::endl;
        std::cout << "\t-p\tGenerate point grids as well" << std::endl;
        std::cout << "\t--obj\tGenerate .obj files as well" << std::endl;
        std::cout << "\t--particleRadius <size>\tSpecify the radius of the particles. Default: "
                  << DEFAULT_PARTICLE_RADIUS << std::endl;
        std::cout << "\t--voxelSize <size>\tSpecify the size of the voxels. Default: " << DEFAULT_VOXEL_SIZE
                  << std::endl;
        std::cout << "\t--halfWidth <size>\tSpecify the half width of the narrow band the levelset. Default: "
                  << DEFAULT_HALF_WIDTH << std::endl;
        std::cout << std::endl;

        return 0;
    }

    char *fileName = getCmdOption(argv, argv + argc, "-f");

    if (fileName == nullptr) {
        std::cout << "Usage: " << argv[0] << " -f <fileName>\n";
        return 1;
    }

    char *outFileName = getCmdOption(argv, argv + argc, "-o");
    if (outFileName == nullptr) {
        outFileName = DEFAULT_OUTPUT_FILE_NAME;
    }

    bool shouldGenerateGrids = cmdOptionExists(argv, argv + argc, "-p");

    bool shouldGenerateObjFiles = cmdOptionExists(argv, argv + argc, "--obj");

    char *cmdParticleRadius = getCmdOption(argv, argv + argc, "--particleRadius");
    float particleRadius = cmdParticleRadius == nullptr ? DEFAULT_PARTICLE_RADIUS : std::stof(cmdParticleRadius);

    char *cmdVoxelSize = getCmdOption(argv, argv + argc, "--voxelSize");
    float voxelSize = cmdVoxelSize == nullptr ? DEFAULT_VOXEL_SIZE : std::stof(cmdVoxelSize);

    char *cmdHalfWidth = getCmdOption(argv, argv + argc, "--halfWidth");
    float halfWidth = cmdHalfWidth == nullptr ? DEFAULT_HALF_WIDTH : std::stof(cmdHalfWidth);

    // read particlesPositions from file particlesPositions.txt
    // each line contains 3 float numbers separated by space
    std::ifstream infile(fileName);
    if (!infile) {
        std::cout << "Error opening " << fileName << " file" << std::endl;
        return 1;
    }

    int sizeX, sizeY, sizeZ;
    infile >> sizeX >> sizeY >> sizeZ;
    std::cout << "Reading grid size: " << sizeX << " " << sizeY << " " << sizeZ << std::endl;

    float dt;
    int nbTimeSteps;
    infile >> dt >> nbTimeSteps;
    std::cout << "Reading time step: " << dt << std::endl;
    std::cout << "Reading number of time steps: " << nbTimeSteps << std::endl;

    std::vector<std::vector<openvdb::Vec3R>> particlesPositions;
    float x, y, z;
    for (int i = 0; i < nbTimeSteps; i++) {
        int nbParticles;
        infile >> nbParticles;
        std::cout << i << "th iteration: reading " << nbParticles << " particles: " << std::endl;
        particlesPositions.emplace_back();
        for (unsigned long int j = 0; j < nbParticles; j++) {
            infile >> x >> y >> z;
            particlesPositions[i].emplace_back(x, y, z);
        }
        if (nbParticles > 0)
            std::cout << "First particle of " << i << "th step: " << particlesPositions[i][0] << std::endl;

        assert(particlesPositions[i].size() == nbParticles &&
               "Number of particlesPositions read is not the same as the number of particlesPositions");

    }

    infile.close();

    assert(particlesPositions.size() == nbTimeSteps &&
           "Number of particlesPositions read is not the same as the number of particlesPositions");

    // Initialize the OpenVDB library.  This must be called at least
    // once per program and may safely be called multiple times.
    openvdb::initialize();

    for (int i = 0; i < nbTimeSteps; i++) {
        std::cout << "\n----------- Time step " << i << " -----------\n" << std::endl;

        // creates a SDF and a mesh from the particles (output in the results folder)
        rasterizeParticles(particlesPositions[i], outFileName, particleRadius, voxelSize, halfWidth, i, dt,
                           shouldGenerateObjFiles);

        // creates a point grid from the particles (output in the results folder) not compatible with Blender
        if (shouldGenerateGrids) createPointGrid(particlesPositions[i], outFileName, i);
    }

    std::cout << "\n----------- Done -----------\n" << std::endl;

    return 0;
}