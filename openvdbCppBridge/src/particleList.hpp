//
// Created by barth on 16/05/23.
//

#ifndef OPENVDBBRIDGE_PARTICLELIST_HPP
#define OPENVDBBRIDGE_PARTICLELIST_HPP

#include <openvdb/Types.h>
#include <utility>
#include <vector>
#include <cassert>

/**
 * This class is required by openvdb::tools::ParticlesToLevelSet (see https://www.openvdb.org/documentation/doxygen/ParticlesToLevelSet_8h.html)
 */
class ParticleList {
protected:
    struct Particle {
        openvdb::Vec3R position, velocity;
        openvdb::Real radius;
    };

public:
    using Real = openvdb::Real;
    using PosType = openvdb::Vec3R; // required by openvdb::tools::PointPartitioner
    using AttributeType = Real; // required by openvdb::tools::ParticlesToLevelSet

    constexpr static float PARTICLE_RADIUS = 0.1;

    explicit ParticleList(std::vector<openvdb::Vec3R> &particlePositions) {
        for (auto &p: particlePositions) {
            add(p);
        }
    }

    void add(const openvdb::Vec3R &p, const openvdb::Real &r = PARTICLE_RADIUS,
             const openvdb::Vec3R &v = openvdb::Vec3R(0, 0, 0)) {
        Particle particle{};
        particle.position = p;
        particle.radius = r;
        particle.velocity = v;
        _particles.push_back(particle);
    }

    // Return the total number of particles in the list.
    // Always required!
    size_t size() const {
        return this->_particles.size();
    }

    // Get the world-space position of the nth particle.
    // Required by rasterizeSpheres().
    void getPos(size_t n, openvdb::Vec3R &xyz) const {
        assert(n < this->_particles.size());
        xyz = this->_particles[n].position;
    }

    // Get the world-space position and radius of the nth particle.
    // Required by rasterizeSpheres().
    void getPosRad(size_t n, openvdb::Vec3R &xyz, openvdb::Real &radius) const {
        assert(n < this->_particles.size());
        xyz = this->_particles[n].position;
        radius = this->_particles[n].radius;
    }

    // Get the world-space position, radius and velocity of the nth particle.
    // Required by rasterizeTrails().
    void getPosRadVel(size_t n, openvdb::Vec3R &xyz, openvdb::Real &radius, openvdb::Vec3R &velocity) const {
        assert(n < this->_particles.size());
        xyz = this->_particles[n].position;
        radius = this->_particles[n].radius;
        velocity = this->_particles[n].velocity;
    }

    // Get the value of the nth particle's user-defined attribute (of type @c AttributeType).
    // Required only if attribute transfer is enabled in ParticlesToLevelSet.
    void getAtt(size_t n, openvdb::Index32 &att) const { att = openvdb::Index32(n); }

    openvdb::CoordBBox getBBox(const openvdb::GridBase &grid) {
        openvdb::CoordBBox bbox;
        openvdb::Coord &min = bbox.min(), &max = bbox.max();
        openvdb::Vec3R pos;
        openvdb::Real rad, invDx = 1 / grid.voxelSize()[0];
        for (size_t n = 0, e = this->size(); n < e; ++n) {
            this->getPosRad(n, pos, rad);
            const openvdb::Vec3d xyz = grid.worldToIndex(pos);
            const openvdb::Real r = rad * invDx;
            for (int i = 0; i < 3; ++i) {
                min[i] = openvdb::math::Min(min[i], openvdb::math::Floor(xyz[i] - r));
                max[i] = openvdb::math::Max(max[i], openvdb::math::Ceil(xyz[i] + r));
            }
        }
        return bbox;
    }

private:
    std::vector<Particle> _particles{};
};

#endif //OPENVDBBRIDGE_PARTICLELIST_HPP
