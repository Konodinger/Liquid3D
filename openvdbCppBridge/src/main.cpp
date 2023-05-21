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

    std::vector<openvdb::Vec3R> particlesPositions;

    float x, y, z;
    while (infile >> x >> y >> z) { particlesPositions.emplace_back(x, y, z); }
    infile.close();

    std::cout << "Number of particlesPositions: " << particlesPositions.size() << std::endl;
    std::cout << "First particle: " << particlesPositions[0] << std::endl;

    // Initialize the OpenVDB library.  This must be called at least
    // once per program and may safely be called multiple times.
    openvdb::initialize();

    // creates a SDF and a mesh from the particles (output in the results folder)
    rasterizeParticles(particlesPositions);

    // creates a point grid from the particles (output in the results folder) not compatible with Blender
    createPointGrid(particlesPositions);

    return 0;
}