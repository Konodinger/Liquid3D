//
// Created by barth on 16/05/23.
//

#ifndef OPENVDBBRIDGE_PARTICLELIST_H
#define OPENVDBBRIDGE_PARTICLELIST_H

#include <openvdb/Types.h>
#include <utility>
#include <vector>
#include <cassert>

class ParticleList {
public:
    using Real = openvdb::Real;
    using PosType = openvdb::Vec3R; // required by openvdb::tools::PointPartitioner
    using AttributeType = Real; // required by openvdb::tools::ParticlesToLevelSet

    ParticleList(std::vector<openvdb::Vec3R> particlePositions) {
        this->_particlePositions = std::move(particlePositions);
    }

    // Return the total number of particles in the list.
    // Always required!
    size_t size() const {
        return this->_particlePositions.size();
    }

    // Get the world-space position of the nth particle.
    // Required by rasterizeSpheres().
    void getPos(size_t n, openvdb::Vec3R &xyz) const {
        assert(n < this->_particlePositions.size());
        xyz = this->_particlePositions[n];
    }

    // Get the world-space position and radius of the nth particle.
    // Required by rasterizeSpheres().
    void getPosRad(size_t n, openvdb::Vec3R &xyz, openvdb::Real &radius) const {
        assert(n < this->_particlePositions.size());
        xyz = this->_particlePositions[n];
        radius = 0.1; //FIXME: this will not always be the same
    }

    // Get the world-space position, radius and velocity of the nth particle.
    // Required by rasterizeTrails().
    void getPosRadVel(size_t n, openvdb::Vec3R &xyz, openvdb::Real &radius, openvdb::Vec3R &velocity) const {
        assert(n < this->_particlePositions.size());
        xyz = this->_particlePositions[n];
        radius = 0.1; //FIXME: this will not always be the same
        velocity = {0.0, 0.0, 0.0}; //FIXME: this will not always be the same
    }

    // Get the value of the nth particle's user-defined attribute (of type @c AttributeType).
    // Required only if attribute transfer is enabled in ParticlesToLevelSet.
    void getAtt(size_t n, AttributeType &att) const {
        assert(n < this->_particlePositions.size());
        att = 0.0; //FIXME: this will not always be the same
    }

private:
    std::vector<openvdb::Vec3R> _particlePositions;
};

#endif //OPENVDBBRIDGE_PARTICLELIST_H
