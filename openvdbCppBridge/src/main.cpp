#include <openvdb/openvdb.h>
#include <openvdb/io/Stream.h>
#include <openvdb/points/PointConversion.h>
#include <openvdb/points/PointCount.h>
#include <openvdb/tools/ParticlesToLevelSet.h>

#include <iostream>
#include <fstream>

using namespace openvdb;

// src: https://www.openvdb.org/documentation/doxygen/codeExamples.html#sPointsHelloWorld

int main() {
    // Initialize the OpenVDB library.  This must be called at least
    // once per program and may safely be called multiple times.
    initialize();

    // PART 1

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

    // The VDB Point-Partioner is used when bucketing points and requires a
    // specific interface. For convenience, we use the PointAttributeVector
    // wrapper around an stl vector wrapper here, however it is also possible to
    // write one for a custom data structure in order to match the interface
    // required.
    openvdb::points::PointAttributeVector<openvdb::Vec3R> particlePositionsWrapper(particlesPositions);

    // This method computes a voxel-size to match the number of
    // points / voxel requested. Although it won't be exact, it typically offers
    // a good balance of memory against performance.
    int pointsPerVoxel = 8;
    float voxelSize = openvdb::points::computeVoxelSize(particlePositionsWrapper, pointsPerVoxel);

    // Print the voxel-size to cout
    std::cout << "VoxelSize=" << voxelSize << std::endl;

    // Create a transform using this voxel-size.
    openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(voxelSize);

    // Create a PointDataGrid containing these four points and using the
    // transform given. This function has two template parameters, (1) the codec
    // to use for storing the position, (2) the grid we want to create
    // (ie a PointDataGrid).
    // We use no compression here for the positions.
    openvdb::points::PointDataGrid::Ptr grid = openvdb::points::createPointDataGrid<openvdb::points::NullCodec,
            openvdb::points::PointDataGrid>(particlesPositions, *transform);

    //TODO: convert the point data grid to a level set to export it to Blender
    //openvdb::FloatGrid::Ptr grid2 = openvdb::FloatGrid::create();
    //openvdb::tools::particlesToSdf(particlePositionsWrapper, *grid2);

    const char *gridName = "Particles";

    // Set the name of the grid
    grid->setName(gridName);

    // Create a VDB file object and write out the grid.
    //FIXME: the PointDataGrid is not compatible with Blender
    openvdb::io::File("particles.vdb").write({grid});

    return 0;
}