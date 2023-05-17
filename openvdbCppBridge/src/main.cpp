#include <openvdb/openvdb.h>
#include <openvdb/io/Stream.h>
#include <openvdb/tools/ParticlesToLevelSet.h>
#include <openvdb/Types.h>

#include <iostream>
#include <fstream>

#include "particleList.h"

using namespace openvdb;

// src: https://www.openvdb.org/documentation/doxygen/codeExamples.html#sPointsHelloWorld
// based also on: https://github.com/dneg/openvdb/blob/587c9ae84c2822bbc03d0d7eceb52898582841b9/openvdb/openvdb/unittest/TestParticlesToLevelSet.cc#L470

void writeGrid(const openvdb::GridBase::Ptr& grid, const std::string& fileName) {
    std::cout << "\nWriting \"" << fileName << "\" to file\n";
    grid->setName("TestParticlesToLevelSet");
    openvdb::GridPtrVec grids;
    grids.push_back(grid);
    openvdb::io::File file(fileName + ".vdb");
    file.write(grids);
    file.close();
}

int main() {
    // Initialize the OpenVDB library.  This must be called at least
    // once per program and may safely be called multiple times.
    initialize();

    std::vector<openvdb::Vec3R> particlesPositions;

    // read particlesPositions from file particlesPositions.txt
    // each line contains 3 float numbers separated by space
    std::ifstream infile("../particles.txt");
    if (!infile) {
        std::cout << "Error opening particlesPositions file" << std::endl;
        return 1;
    }

    float x, y, z;
    while (infile >> x >> y >> z) { particlesPositions.emplace_back(x, y, z); }

    std::cout << "Number of particlesPositions: " << particlesPositions.size() << std::endl;
    std::cout << "First particle: " << particlesPositions[0] << std::endl;

    auto particleList = new ParticleList(particlesPositions);
    std::cout << "Created OpenVDB compatible particle list" << std::endl;

    // create a level set from the particles
    const float voxelSize = 1.0f, halfWidth = 2.0f;
    openvdb::FloatGrid::Ptr ls = openvdb::createLevelSet<openvdb::FloatGrid>(voxelSize, halfWidth);
    openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid> raster(*ls);

    //raster.setGrainSize(1);//a value of zero disables threading
    raster.rasterizeSpheres(*particleList);
    raster.finalize();

    writeGrid(ls, "testRaster");

    return 0;
}