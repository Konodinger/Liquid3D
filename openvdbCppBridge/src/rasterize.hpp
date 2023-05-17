//
// Created by barth on 17/05/23.
//

#ifndef OPENVDBBRIDGE_RASTERIZE_HPP
#define OPENVDBBRIDGE_RASTERIZE_HPP

#include <vector>
#include <openvdb/openvdb.h>
#include <openvdb/tools/ParticlesToLevelSet.h>
#include <openvdb/tools/VolumeToMesh.h>
#include "particleList.hpp"
#include "utils.hpp"

// this in inspired on https://github.com/dneg/openvdb/blob/587c9ae84c2822bbc03d0d7eceb52898582841b9/openvdb/openvdb/unittest/TestParticlesToLevelSet.cc#L470

/**
 * Rasterize a list of particles into a level set then write it to a file
 * @param positions the list of particle positions
 */
void rasterizeParticles(std::vector<openvdb::Vec3R> positions) {
    auto particleList = new ParticleList(positions);
    std::cout << "Created OpenVDB compatible particle list" << std::endl;

    // create a level set from the particles
    const float voxelSize = 1.0f, halfWidth = 2.0f;
    openvdb::FloatGrid::Ptr ls = openvdb::createLevelSet<openvdb::FloatGrid>(voxelSize, halfWidth);
    openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid> raster(*ls);

    raster.setGrainSize(1); //a value of zero disables threading
    raster.rasterizeSpheres(*particleList);
    raster.finalize();

    writeGrid(ls, "testRaster");

    // now we have to use https://github.com/AcademySoftwareFoundation/openvdb/blob/4a71881492520eb1323876becd0dad27eb0c2dcc/openvdb/openvdb/tools/VolumeToMesh.h to convert the level set to a mesh

    std::vector<openvdb::Vec3s> points;
    std::vector<openvdb::Vec4I> quads;
    std::vector<openvdb::Vec3I> triangles;
    openvdb::tools::volumeToMesh(*ls, points, triangles, quads);

    std::cout << "points: " << points.size() << std::endl;
    std::cout << "quads: " << quads.size() << std::endl;
    std::cout << "triangles: " << triangles.size() << std::endl;
}

#endif //OPENVDBBRIDGE_RASTERIZE_HPP
