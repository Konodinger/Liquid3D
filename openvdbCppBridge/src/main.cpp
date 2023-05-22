#include <openvdb/openvdb.h>
#include <openvdb/Types.h>

#include <iostream>
#include <fstream>

#include "rasterize.hpp"
#include "pointGrid.hpp"

int main(int argc, char **argv) {
    // usage: ./openvdbBridge <particlesPositions.txt>
    std::string fileName;
    if (argc > 1) {
        fileName = argv[1];
    } else {
        fileName = "./particles.txt";
    }

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

    int nbBoundaryParticles;
    infile >> nbBoundaryParticles;
    std::cout << "Reading number of boundary particles: " << nbBoundaryParticles << std::endl;

    float x, y, z;
    for (int i = 0; i < nbBoundaryParticles; i++) {
        infile >> x >> y >> z;
        // do nothing with boundary particles
    }

    int nbParticles;
    infile >> nbParticles;
    std::cout << "Reading number of particles: " << nbParticles << std::endl;

    std::vector<std::vector<openvdb::Vec3R>> particlesPositions;
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
        rasterizeParticles(particlesPositions[i], "fluid", i);

        // creates a point grid from the particles (output in the results folder) not compatible with Blender
        createPointGrid(particlesPositions[i], "fluid", i);
    }

    std::cout << "\n----------- Done -----------\n" << std::endl;

    return 0;
}