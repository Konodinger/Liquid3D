#include <openvdb/openvdb.h>
#include <openvdb/Types.h>

#include <iostream>
#include <fstream>

#include "rasterize.hpp"
#include "pointGrid.hpp"

int main() {
    // Initialize the OpenVDB library.  This must be called at least
    // once per program and may safely be called multiple times.
    openvdb::initialize();

    std::vector<openvdb::Vec3R> particlesPositions;

    // read particlesPositions from file particlesPositions.txt
    // each line contains 3 float numbers separated by space
    std::ifstream infile("./particles.txt");
    if (!infile) {
        std::cout << "Error opening particlesPositions file" << std::endl;
        return 1;
    }

    float x, y, z;
    while (infile >> x >> y >> z) { particlesPositions.emplace_back(x, y, z); }

    std::cout << "Number of particlesPositions: " << particlesPositions.size() << std::endl;
    std::cout << "First particle: " << particlesPositions[0] << std::endl;

    rasterizeParticles(particlesPositions);

    createPointGrid(particlesPositions);

    return 0;
}