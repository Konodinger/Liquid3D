//
// Created by barth on 17/05/23.
//

#ifndef OPENVDBBRIDGE_RASTERIZE_HPP
#define OPENVDBBRIDGE_RASTERIZE_HPP

#include <vector>
#include <openvdb/openvdb.h>
#include <openvdb/tools/ParticlesToLevelSet.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/points/PointConversion.h>
#include <openvdb/points/PointCount.h>
#include "particleList.hpp"
#include "utils.hpp"
#include "exportToObj.h"

// based on the slides: https://artifacts.aswf.io/io/aswf/openvdb/openvdb_toolset_2013/1.0.0/openvdb_toolset_2013-1.0.0.pdf
// code inspired on https://github.com/dneg/openvdb/blob/587c9ae84c2822bbc03d0d7eceb52898582841b9/openvdb/openvdb/unittest/TestParticlesToLevelSet.cc#L470
// see this doc: https://www.openvdb.org/documentation/doxygen/classopenvdb_1_1v10__0_1_1tools_1_1ParticlesToLevelSet.html

/**
 * Rasterize a list of particles into a level set then write it to a file
 * @param positions the list of particle positions
 */
void rasterizeParticles(std::vector<openvdb::Vec3R> positions) {
    auto particleList = new ParticleList(positions);
    std::cout << "Created OpenVDB compatible particle list" << std::endl;

    openvdb::points::PointAttributeVector<openvdb::Vec3R> particlePositionsWrapper(positions);

    // This method computes a voxel-size to match the number of
    // points / voxel requested. Although it won't be exact, it typically offers
    // a good balance of memory against performance.
    //int pointsPerVoxel = 8;
    //float voxelSize = openvdb::points::computeVoxelSize(particlePositionsWrapper, pointsPerVoxel);

    // points to levelset (https://github.com/dneg/openvdb/blob/587c9ae84c2822bbc03d0d7eceb52898582841b9/openvdb/openvdb/unittest/TestParticlesToLevelSet.cc#L471)
    // see also https://stackoverflow.com/questions/68405603/i-am-trying-to-convert-point-cloud-to-mesh-using-openvdb

    const float voxelSize = 0.5f;
    const float halfWidth = 2.0f;
    openvdb::FloatGrid::Ptr levelSet = openvdb::createLevelSet<openvdb::FloatGrid>(voxelSize, halfWidth);
    openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid> raster(*levelSet);

    raster.setGrainSize(1); //a value of zero disables threading
    raster.rasterizeSpheres(*particleList, 1);
    raster.finalize();

    std::cout << "\nParticles have been converted to a level set:" << std::endl;
    std::cout << "Memory usage: " << levelSet->tree().memUsage() << std::endl;
    std::cout << "Active voxel count: " << levelSet->activeVoxelCount() << std::endl;
    std::cout << "Bounding box: " << levelSet->evalActiveVoxelBoundingBox() << std::endl;
    levelSet->tree().print(std::cout, 4);


    // now we have to use https://github.com/AcademySoftwareFoundation/openvdb/blob/4a71881492520eb1323876becd0dad27eb0c2dcc/openvdb/openvdb/tools/VolumeToMesh.h to convert the level set to a mesh

    std::cout << "\nMeshing the level set" << std::endl;

    std::vector<openvdb::Vec3s> points;
    std::vector<openvdb::Vec4I> quads;
    std::vector<openvdb::Vec3I> triangles;
    openvdb::tools::volumeToMesh(*levelSet, points, triangles, quads, 1, 0);

    std::cout << "points: " << points.size() << std::endl;
    std::cout << "quads: " << quads.size() << std::endl;
    std::cout << "triangles: " << triangles.size() << std::endl;

    exportToObj(points, quads, triangles);

    writeGrid(levelSet, "testRaster");
}

#endif //OPENVDBBRIDGE_RASTERIZE_HPP
