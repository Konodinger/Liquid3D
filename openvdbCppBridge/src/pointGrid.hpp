//
// Created by barth on 17/05/23.
//

#ifndef OPENVDBBRIDGE_POINTGRID_HPP
#define OPENVDBBRIDGE_POINTGRID_HPP

#include <vector>
#include <openvdb/openvdb.h>
#include <openvdb/points/PointConversion.h>

#include "utils.hpp"

/**
 * Create a point grid from a list of particles
 * @param positions the list of particle positions
 */
void createPointGrid(std::vector<openvdb::Vec3R> positions) {
    openvdb::points::PointAttributeVector<openvdb::Vec3R> particlePositionsWrapper(positions);

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
    auto grid = openvdb::points::createPointDataGrid<openvdb::points::NullCodec,
            openvdb::points::PointDataGrid>(positions, *transform);

    const char *gridName = "Particles";

    // Set the name of the grid
    grid->setName(gridName);

    // Create a VDB file object and write out the grid.
    writeGrid(grid, "pointGrid");
}

#endif //OPENVDBBRIDGE_POINTGRID_HPP
