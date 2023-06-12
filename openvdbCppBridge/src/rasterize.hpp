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

/**
 * Rasterize a list of particles into a level set and write it to a `.vdb` file.
 * The files are output in the `results` directory with the given name and iteration number.
 * @param positions the list of particle positions for one iteration
 * @param fileName the name of the file to output (without extension)
 * @param particleRadius the radius of the particles
 * @param voxelSize the size of the voxels
 * @param halfWidth the half width of the level set
 * @param iteration the iteration number (0 by default) to append to the file name
 * @param shouldGenerateObjFiles whether or not to generate `.obj` files (false by default)
 *
 * @see openvdb slides https://artifacts.aswf.io/io/aswf/openvdb/openvdb_toolset_2013/1.0.0/openvdb_toolset_2013-1.0.0.pdf
 * @see unit test https://github.com/dneg/openvdb/blob/587c9ae84c2822bbc03d0d7eceb52898582841b9/openvdb/openvdb/unittest/TestParticlesToLevelSet.cc#L470
 * @see openvdb doc https://www.openvdb.org/documentation/doxygen/classopenvdb_1_1v10__0_1_1tools_1_1ParticlesToLevelSet.html
 */
void rasterizeParticles(std::vector<openvdb::Vec3R> &positions, const std::string &fileName,
                        const float particleRadius, const float voxelSize, const float halfWidth,
                        const int iteration = 0, const float dt = 0.0, const bool shouldGenerateObjFiles = false) {
    auto particleList = new ParticleList(positions);
    std::cout << "Created OpenVDB compatible particle list" << std::endl;

    //openvdb::points::PointAttributeVector<openvdb::Vec3R> particlePositionsWrapper(positions);

    // This method computes a voxel-size to match the number of
    // points / voxel requested. Although it won't be exact, it typically offers
    // a good balance of memory against performance.
    //int pointsPerVoxel = 8;
    //float voxelSize = openvdb::points::computeVoxelSize(particlePositionsWrapper, pointsPerVoxel);

    // points to levelset (https://github.com/dneg/openvdb/blob/587c9ae84c2822bbc03d0d7eceb52898582841b9/openvdb/openvdb/unittest/TestParticlesToLevelSet.cc#L471)
    // see also https://stackoverflow.com/questions/68405603/i-am-trying-to-convert-point-cloud-to-mesh-using-openvdb

    openvdb::FloatGrid::Ptr levelSet = openvdb::createLevelSet<openvdb::FloatGrid>(voxelSize, halfWidth);
    levelSet->setName("levelSet" + std::to_string(iteration));
    levelSet->insertMeta("dt", openvdb::FloatMetadata(dt));

    openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid> raster(*levelSet);
    raster.setGrainSize(1); //a value of zero disables threading
    raster.rasterizeSpheres(*particleList,particleRadius);
    raster.finalize();

    std::cout << "\nParticles have been converted to a level set:" << std::endl;
    std::cout << "Memory usage: " << levelSet->tree().memUsage() << std::endl;
    std::cout << "Active voxel count: " << levelSet->activeVoxelCount() << std::endl;
    std::cout << "Bounding box: " << levelSet->evalActiveVoxelBoundingBox() << std::endl;

    writeGrid(levelSet, fileName + "SDF" + std::to_string(iteration));

    if (!shouldGenerateObjFiles) return;

    // now we have to use https://github.com/AcademySoftwareFoundation/openvdb/blob/4a71881492520eb1323876becd0dad27eb0c2dcc/openvdb/openvdb/tools/VolumeToMesh.h to convert the level set to a mesh
    std::cout << "\nMeshing the level set" << std::endl;

    std::vector<openvdb::Vec3s> points;
    std::vector<openvdb::Vec4I> quads;
    std::vector<openvdb::Vec3I> triangles;
    openvdb::tools::volumeToMesh(*levelSet, points, triangles, quads, 0, 0);

    std::cout << "points: " << points.size() << std::endl;
    std::cout << "quads: " << quads.size() << std::endl;
    std::cout << "triangles: " << triangles.size() << std::endl;

    exportToObj(points, quads, triangles, fileName + "Mesh" + std::to_string(iteration));
}

#endif //OPENVDBBRIDGE_RASTERIZE_HPP
