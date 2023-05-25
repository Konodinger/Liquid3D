#include <openvdb/openvdb.h>
#include <openvdb/Types.h>

#include <iostream>
#include <fstream>

#include "rasterize.hpp"
#include "pointGrid.hpp"

int main(int argc, char **argv) {
    // usage: ./openvdbBridge <particlesPositions.txt>

    if (cmdOptionExists(argv, argv + argc, "-h")) {
        std::cout << "Usage: " << argv[0] << " -f <fileName>\n" << std::endl;

        std::cout << "Options:" << std::endl;
        std::cout << "\t-f <fileName>\tSpecify the file to read" << std::endl;
        std::cout << "\t-h\tPrint this help message" << std::endl;
        std::cout << "\t-p\tGenerate point grids as well" << std::endl;
        std::cout << "\t--obj\tGenerate .obj files as well" << std::endl;
        std::cout << std::endl;

        return 0;
    }

    char *fileName = getCmdOption(argv, argv + argc, "-f");

    if (fileName == nullptr) {
        std::cout << "Usage: " << argv[0] << " -f <fileName>\n";
        return 1;
    }

    bool shouldGenerateGrids = cmdOptionExists(argv, argv + argc, "-p");

    bool shouldGenerateObjFiles = cmdOptionExists(argv, argv + argc, "--obj");

    // read particlesPositions from file particlesPositions.txt
    // each line contains 3 float numbers separated by space
    std::ifstream infile(fileName);
    if (!infile) {
        std::cout << "Error opening" << fileName << "file" << std::endl;
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

    int nbParticles;
    infile >> nbParticles;
    std::cout << "Reading number of particles: " << nbParticles << std::endl;

    std::vector<std::vector<openvdb::Vec3R>> particlesPositions;
    float x, y, z;
    for (int i = 0; i < nbTimeSteps; i++) {
        particlesPositions.emplace_back();
        for (int j = 0; j < nbParticles; j++) {
            infile >> x >> y >> z;
            particlesPositions[i].emplace_back(x, y, z);
        }
        assert(particlesPositions[i].size() == nbParticles &&
               "Number of particlesPositions read is not the same as the number of particlesPositions");

    }

    infile.close();

    assert(particlesPositions.size() == nbTimeSteps &&
           "Number of particlesPositions read is not the same as the number of particlesPositions");

    std::cout << "First particle of first step: " << particlesPositions[0][0] << std::endl;

    // Initialize the OpenVDB library.  This must be called at least
    // once per program and may safely be called multiple times.
    openvdb::initialize();

    for (int i = 0; i < nbTimeSteps; i++) {
        std::cout << "\n----------- Time step " << i << " -----------\n" << std::endl;

        // creates a SDF and a mesh from the particles (output in the results folder)
        rasterizeParticles(particlesPositions[i], "fluid", i, shouldGenerateObjFiles);

        // creates a point grid from the particles (output in the results folder) not compatible with Blender
        if (shouldGenerateGrids) createPointGrid(particlesPositions[i], "fluid", i);
    }

    std::cout << "\n----------- Done -----------\n" << std::endl;

    return 0;
}