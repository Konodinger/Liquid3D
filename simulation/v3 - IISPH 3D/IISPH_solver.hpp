#ifndef ___IISPH_SOLVER___
#define ___IISPH_SOLVER___

#include <cstdio>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cmath>

#ifdef _OPENMP

#include <omp.h>

#endif

#include "Vector.hpp"
#include "Kernel.hpp"
#include "particleInitialization.hpp"

using namespace std;

// debug index
int obsPart = 2000;

#define LCONFING_FRACTION_X 0.4f
#define LCONFING_FRACTION_Y 0.4f

class IisphSolver {
public:
    explicit IisphSolver(
            const Real dt = 0.0005, const Real nu = 0.08,
            const Real h = 0.5, const Real density = 1e3, const Real wallWeightCoef = 2.f,
            const Vec3f g = Vec3f(0, 0, -9.8), const Real initP = 0.5,
            const Real omega = 0.5, const Real pressureError = 1.) :
            _kernel(h), _rad(_kernel.supportRadius()), _dt(dt), _wallWeightCoef(wallWeightCoef),
            _nu(nu), _h(h), _d0(density), _g(g),
            _initP(initP), _omega(omega), _pressureError(pressureError) {
        _m0 = _d0 * _h * _h * _h / 8. * M_PI * 4. / 3.;
        maxVel = Vec3f();
    }

    // assume an arbitrary grid with the size of res_x*res_y; a fluid mass fill up
    // the size of f_width, f_height; each cell is sampled with 2x2 particles.
    void initScene(const Vec3f gridRes, InitType initType, bool useLConfig = false) {
        _useLConfig = useLConfig;

        _pos.clear();

        _resX = gridRes.x;
        _resY = gridRes.y;
        _resZ = gridRes.z;

        /// WALL BOUNDARY PARTICLES

        // set wall for boundary
        _left = 0.5 * _h;
        _right = static_cast<Real>(_resX) - 0.5 * _h;
        _bottom = 0.5 * _h;
        _top = static_cast<Real>(_resY) - 0.5 * _h;
        _back = 0.5 * _h;
        _front = static_cast<Real>(_resZ) - 0.5 * _h;


        const int wallDensity = 5;
        for (int x: {0, _resX - 1}) {
            for (int y = 0; y < _resY; ++y) {
                for (int z = 0; z < _resZ; ++z) {
                    for (int lx = 1; lx < wallDensity; lx++) {
                        for (int ly = 1; ly < wallDensity; ly++) {
                            for (int lz = 1; lz < wallDensity; lz++) {
                                _pos.push_back(Vec3f(x + lx / (Real) wallDensity, y + ly / (Real) wallDensity,
                                                     z + lz / (Real) wallDensity));
                            }
                        }
                    }
                    /*_pos.push_back(Vec3f(x + 0.25, y + 0.25, z + 0.25));
                    _pos.push_back(Vec3f(x + 0.75, y + 0.25, z + 0.25));
                    _pos.push_back(Vec3f(x + 0.25, y + 0.75, z + 0.25));
                    _pos.push_back(Vec3f(x + 0.75, y + 0.75, z + 0.25));
                    _pos.push_back(Vec3f(x + 0.25, y + 0.25, z + 0.75));
                    _pos.push_back(Vec3f(x + 0.75, y + 0.25, z + 0.75));
                    _pos.push_back(Vec3f(x + 0.25, y + 0.75, z + 0.75));
                    _pos.push_back(Vec3f(x + 0.75, y + 0.75, z + 0.75));*/
                }
            }
        }


        for (int x = 2; x < _resX - 2; ++x) {
            for (int y: {0, 1, _resY - 2, _resY - 1}) {
                for (int z = 0; z < _resZ; ++z) {
                    for (int lx = 1; lx < wallDensity; lx++) {
                        for (int ly = 1; ly < wallDensity; ly++) {
                            for (int lz = 1; lz < wallDensity; lz++) {
                                _pos.push_back(Vec3f(x + lx / (Real) wallDensity, y + ly / (Real) wallDensity,
                                                     z + lz / (Real) wallDensity));
                            }
                        }
                    }
                }
            }
        }


        for (int x = 1; x < _resX - 1; ++x) {
            for (int y = 1; y < _resY - 1; ++y) {
                for (int z: {0, _resZ - 1}) {
                    for (int lx = 1; lx < wallDensity; lx++) {
                        for (int ly = 1; ly < wallDensity; ly++) {
                            for (int lz = 1; lz < wallDensity; lz++) {
                                _pos.push_back(Vec3f(x + lx / (Real) wallDensity, y + ly / (Real) wallDensity,
                                                     z + lz / (Real) wallDensity));
                            }
                        }
                    }
                }
            }
        }

        // Additional obstacles

        if (_useLConfig) {
            int lLineX = int(_resX * LCONFING_FRACTION_X);
            int lLineY = int(_resY * LCONFING_FRACTION_Y);
            for (int i: {lLineX, lLineX + 1}) {
                for (int j = 2; j < lLineY + 2; ++j) {
                    for (int k = 2; k < _resZ - 2; ++k) {
                        _pos.push_back(Vec3f(i + 0.25, j + 0.25, k + 0.25));
                        _pos.push_back(Vec3f(i + 0.75, j + 0.25, k + 0.25));
                        _pos.push_back(Vec3f(i + 0.25, j + 0.75, k + 0.25));
                        _pos.push_back(Vec3f(i + 0.75, j + 0.75, k + 0.25));
                        _pos.push_back(Vec3f(i + 0.25, j + 0.25, k + 0.75));
                        _pos.push_back(Vec3f(i + 0.75, j + 0.25, k + 0.75));
                        _pos.push_back(Vec3f(i + 0.25, j + 0.75, k + 0.75));
                        _pos.push_back(Vec3f(i + 0.75, j + 0.75, k + 0.75));
                    }
                }
            }

            for (int i = 2; i < lLineX; ++i) {
                for (int j: {lLineY, lLineY + 1}) {
                    for (int k = 2; k < _resZ - 2; ++k) {
                        _pos.push_back(Vec3f(i + 0.25, j + 0.25, k + 0.25));
                        _pos.push_back(Vec3f(i + 0.75, j + 0.25, k + 0.25));
                        _pos.push_back(Vec3f(i + 0.25, j + 0.75, k + 0.25));
                        _pos.push_back(Vec3f(i + 0.75, j + 0.75, k + 0.25));
                        _pos.push_back(Vec3f(i + 0.25, j + 0.25, k + 0.75));
                        _pos.push_back(Vec3f(i + 0.75, j + 0.25, k + 0.75));
                        _pos.push_back(Vec3f(i + 0.25, j + 0.75, k + 0.75));
                        _pos.push_back(Vec3f(i + 0.75, j + 0.75, k + 0.75));
                    }
                }
            }
        }

        _nbWallParticles = _pos.size();

        obsPart += _nbWallParticles;

        /// FLUID PARTICLES

        Vec3f blockPosition = Vec3f(0.5f * gridRes.x, 0.8f * gridRes.y, 0.3f * gridRes.z);
        Vec3f blockDimensions = Vec3f(0.8f * gridRes.x, 0.15f * gridRes.y, 0.5f * gridRes.z);

        Vec3f spherePosition = Vec3f(0.5f * gridRes.x, 0.5f * gridRes.y, 0.3f * gridRes.z);
        Real sphereRadius = min(gridRes.x, min(gridRes.y, gridRes.z)) * 0.25f;

        Vec3f torusPosition = Vec3f(0.5f * gridRes.x, 0.5f * gridRes.y, 0.12f * gridRes.z);
        Real torusMajorRadius = min(gridRes.x, min(gridRes.y, gridRes.z)) * 0.25f;
        Real torusMinorRadius = torusMajorRadius / 3.0f;

        if (_useLConfig) {
            blockDimensions = Vec3f((1.0f - LCONFING_FRACTION_X) * gridRes.x * 0.3f,
                                    gridRes.y * 0.5f,
                                    0.5f * gridRes.z);
            blockPosition = Vec3f((1.0f - LCONFING_FRACTION_X) * gridRes.x / 2.0f, 0.5f * gridRes.y, 0.3f * gridRes.z);

            sphereRadius = (1.0f - LCONFING_FRACTION_X) * min(gridRes.x, min(gridRes.y, gridRes.z)) * 0.3f;
            spherePosition = Vec3f((1.0f - LCONFING_FRACTION_X) * gridRes.x / 2.0f,
                                   gridRes.y * 0.5f, 0.2f * gridRes.z);

            torusMajorRadius = 0.15f * min(gridRes.x, min(gridRes.y, gridRes.z));
            torusMinorRadius = 0.05f * min(gridRes.x, min(gridRes.y, gridRes.z));
            torusPosition = Vec3f(0.25f * gridRes.x, 0.25f * gridRes.y, 0.5f * gridRes.z);
        }

        switch (initType) {
            case InitType::BLOCK:
                initBlock(blockPosition, blockDimensions, _pos);
                break;
            case InitType::SPHERE:
                initSphere(spherePosition, sphereRadius, _pos);
                break;
            case InitType::TORUS:
                initTorus(torusPosition, torusMajorRadius, torusMinorRadius, _pos);
                break;
        }
        std::cout << "Simulating " << fluidParticleCount() << " particles of fluid" << std::endl;

        // make sure for the other particle quantities
        _vel = vector<Vec3f>(_pos.size(), Vec3f(0, 0, 0));
        _acc = vector<Vec3f>(_pos.size(), Vec3f(0, 0, 0));
        _p = vector<Real>(_pos.size(), 100000);
        for (int i = 0; i < _nbWallParticles; ++i) {
            _p[i] = 0;
        }
        _d = vector<Real>(_pos.size(), 0);
        _interVel = vector<Vec3f>(_pos.size(), Vec3f(0, 0, 0));
        _interD = vector<Real>(_pos.size(), 0);
        _a_ii = vector<Real>(_pos.size(), 0);
        _d_ii = vector<Vec3f>(_pos.size(), Vec3f(0, 0, 0));
        _c_i = vector<Vec3f>(_pos.size(), Vec3f(0, 0, 0));
        _predP = vector<Real>(_pos.size(), 0);


        for (int i = 0; i < _resX * _resY * _resZ; i++) {
            _pidxInGrid.push_back(vector<tIndex>{});
        }

        buildNeighbor();
        computeDensity();
    }

    void update(Real dt) {
        _dt = dt;
#ifdef __DEBUG1__
        cout << '.' << flush;
#endif

        buildNeighbor();
        computeDensity();

        _acc = vector<Vec3f>(_pos.size(), Vec3f(0, 0, 0));
        applyBodyForce();
        applyViscousForce();

        // IISPH part
        computeIntermediateVelocity();
        computeDiiCoeff();
        computeAiiCoeff();
        computeIntermediateDensity();

        computePressure(); // Huge loop where we compute c_i and p_i.

        // End of IISPH part
        applyPressureForce();

        updateVelocity();
        updatePosition();

        resolveCollision();

    }

    bool particleCollision(const Vec3f &part) const {
        bool collision = (part.x < _left || part.y < _bottom || part.z < _back || part.x > _right ||
                          part.y > _top || part.z > _front);

        if (_useLConfig) {
            collision = collision ||
                        (part.x > (int) ((1.0f - LCONFING_FRACTION_X) * _resX) - 0.5 * _h &&
                         part.y > (int) ((1.0f - LCONFING_FRACTION_Y) * _resY) - 0.5 * _h);
        }

        return collision;
    }

    void clampParticle(Vec3f &part) const {
        const Real newX = clamp(part.x, _left, _right);
        const Real newY = clamp(part.y, _bottom, _top);
        const Real newZ = clamp(part.z, _back, _front);

        if (newX != part.x || newY != part.y || newZ != part.z)
            std::cout << "Clamping particle" << std::endl;

        part.x = newX;
        part.y = newY;
        part.z = newZ;

        if (_useLConfig) {
            if (part.x > (int) ((1.0f - LCONFING_FRACTION_X) * _resX) - 0.5 * _h &&
                part.y > (int) ((1.0f - LCONFING_FRACTION_Y) * _resY) - 0.5 * _h) {
                if (part.x / (1.0f - LCONFING_FRACTION_X) < part.y / (1.0f - LCONFING_FRACTION_Y))
                    part.x = (1.0f - LCONFING_FRACTION_X) *
                             _resX;
                else part.y = (1.0f - LCONFING_FRACTION_Y) * _resY;

            }
        }
    }

    const CubicSpline &getKernel() const { return _kernel; }

    tIndex particleCount() const { return _pos.size(); }

    tIndex wallParticleCount() const { return _nbWallParticles; }

    tIndex fluidParticleCount() const { return _pos.size() - _nbWallParticles; }

    Real particleMass() const { return _m0; }

    const vector<tIndex> &gridParticles(const tIndex i, const tIndex j, const tIndex k) const {
        return _pidxInGrid[idx1d(i, j, k)];
    };

    const Vec3f &gravity() const { return _g; }

    const Vec3f &position(const tIndex i) const { return _pos[i]; }

    const Vec3f &velocity(const tIndex i) const { return _vel[i]; }

    const Vec3f &acceleration(const tIndex i) const { return _acc[i]; }

    Real pressure(const tIndex i) const { return _p[i]; }

    Real density(const tIndex i) const { return _d[i]; }

    const Vec3f &maxVelocity() const { return maxVel; }

    int resX() const { return _resX; }

    int resY() const { return _resY; }

    int resZ() const { return _resZ; }

    void scaleGarvity(Real d) {
        _g *= d;
    }

private:
    void buildNeighbor() {

        for (auto &pix: _pidxInGrid) {
            pix.clear();
        }

        for (tIndex i = 0; i < particleCount(); ++i) {
            const tIndex index = idx1d(floor(position(i).x), floor(position(i).y), floor(position(i).z));
            if (index > _pidxInGrid.size()) {
                std::cout << "Index out of bound: " << i << std::endl;
                std::cout << "Position: " << position(i) << std::endl;
                std::cout << "Index: " << index << std::endl;

                exit(1);
            }
            _pidxInGrid[index].push_back(i);
        }
    }

    void computeDensity() {
#pragma omp parallel for default(none)
        for (tIndex i = 0; i < particleCount(); ++i) {
            Real rho = 0.f;

            for (tIndex _kernelX = max(0, (int) floor(position(i).x - _rad));
                 _kernelX < min(resX(), (int) floor(position(i).x + _rad + 1)); ++_kernelX) {
                for (tIndex _kernelY = max(0, (int) floor(position(i).y - _rad));
                     _kernelY < min(resY(), (int) floor(position(i).y + _rad + 1)); ++_kernelY) {
                    for (tIndex _kernelZ = max(0, (int) floor(position(i).z - _rad));
                         _kernelZ < min(resZ(), (int) floor(position(i).z + _rad + 1)); ++_kernelZ) {

                        for (tIndex j: _pidxInGrid[idx1d(_kernelX, _kernelY, _kernelZ)]) {
                            rho += _m0 * _kernel.w(position(i) - position(j));
                        }
                    }
                }
            }
            _d[i] = rho;
        }
#ifdef __DEBUG3__
        cout << "Density: " << _d[obsPart] << endl;
#endif
    }

    void applyBodyForce() {
#pragma omp parallel for
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            _acc[i] += _g;
        }
    }

    void applyViscousForce() {
#pragma omp parallel for default(none)
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            Vec3f dist = Vec3f();
            Vec3f viscousAcc = Vec3f();
            for (tIndex _kernelX = max(0, (int) floor(position(i).x - _rad));
                 _kernelX < min(resX(), (int) floor(position(i).x + _rad + 1)); ++_kernelX) {
                for (tIndex _kernelY = max(0, (int) floor(position(i).y - _rad));
                     _kernelY < min(resY(), (int) floor(position(i).y + _rad + 1)); ++_kernelY) {
                    for (tIndex _kernelZ = max(0, (int) floor(position(i).z - _rad));
                         _kernelZ < min(resZ(), (int) floor(position(i).z + _rad + 1)); ++_kernelZ) {
                        for (tIndex j: _pidxInGrid[idx1d(_kernelX, _kernelY, _kernelZ)]) {
                            dist = position(i) - position(j);
                            viscousAcc += _m0 / _d[j] * (velocity(i) - velocity(j)) *
                                          (dist.dotProduct(_kernel.grad_w(dist))) / (dist.lengthSquare() + _h * _h);
                        }
                    }
                }
            }
            _acc[i] += 2 * _nu * viscousAcc;
#ifdef __DEBUG4__
            cout << 2 * _nu * viscousAcc << "    ";
#endif
        }
#ifdef __DEBUG3__
        cout << "D " << _d[obsPart] << " A " << _acc[obsPart] << " V " << velocity(obsPart) <<" P " << position(obsPart) << endl;
#endif
    }

    void computeIntermediateVelocity() {
#pragma omp parallel for
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            _interVel[i] = _vel[i] + _dt * _acc[i];
        }
    }

    void computeDiiCoeff() {
#pragma omp parallel for default(none)
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            _d_ii[i] = Vec3f(0, 0, 0);

            for (tIndex _kernelX = max(0, (int) floor(position(i).x - _rad));
                 _kernelX < min(resX(), (int) floor(position(i).x + _rad + 1)); ++_kernelX) {
                for (tIndex _kernelY = max(0, (int) floor(position(i).y - _rad));
                     _kernelY < min(resY(), (int) floor(position(i).y + _rad + 1)); ++_kernelY) {
                    for (tIndex _kernelZ = max(0, (int) floor(position(i).z - _rad));
                         _kernelZ < min(resZ(), (int) floor(position(i).z + _rad + 1)); ++_kernelZ) {
                        for (tIndex j: _pidxInGrid[idx1d(_kernelX, _kernelY, _kernelZ)]) {

                            _d_ii[i] += _m0 * _kernel.grad_w(position(i) - position(j));
                        }
                    }
                }
            }

            _d_ii[i] *= -_dt * _dt / (_d[i] * _d[i]);

        }
    }

    void computeAiiCoeff() {
#pragma omp parallel for default(none)
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            _a_ii[i] = 0;

            for (tIndex _kernelX = max(0, (int) floor(position(i).x - _rad));
                 _kernelX < min(resX(), (int) floor(position(i).x + _rad + 1)); ++_kernelX) {
                for (tIndex _kernelY = max(0, (int) floor(position(i).y - _rad));
                     _kernelY < min(resY(), (int) floor(position(i).y + _rad + 1)); ++_kernelY) {
                    for (tIndex _kernelZ = max(0, (int) floor(position(i).z - _rad));
                         _kernelZ < min(resZ(), (int) floor(position(i).z + _rad + 1)); ++_kernelZ) {
                        for (tIndex j: _pidxInGrid[idx1d(_kernelX, _kernelY, _kernelZ)]) {

                            _a_ii[i] += _m0 * (_d_ii[i] + _dt * _dt * _m0 / (_d[i] * _d[i]) *
                                                          _kernel.grad_w(position(j) - position(i))).dotProduct(
                                    _kernel.grad_w(position(i) - position(j)));
                        }
                    }
                }
            }


        }
    }

    void computeIntermediateDensity() {
#pragma omp parallel for default(none)
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            _interD[i] = _d[i];

            for (tIndex _kernelX = max(0, (int) floor(position(i).x - _rad));
                 _kernelX < min(resX(), (int) floor(position(i).x + _rad + 1)); ++_kernelX) {
                for (tIndex _kernelY = max(0, (int) floor(position(i).y - _rad));
                     _kernelY < min(resY(), (int) floor(position(i).y + _rad + 1)); ++_kernelY) {
                    for (tIndex _kernelZ = max(0, (int) floor(position(i).z - _rad));
                         _kernelZ < min(resZ(), (int) floor(position(i).z + _rad + 1)); ++_kernelZ) {
                        for (tIndex j: _pidxInGrid[idx1d(_kernelX, _kernelY, _kernelZ)]) {

                            _interD[i] += _dt * _m0 * (_interVel[i] - _interVel[j]).dotProduct(
                                    _kernel.grad_w(position(i) - position(j)));
                        }
                    }
                }
            }

        }
    }

    void computePressure() {
#pragma omp parallel for default(none)
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            _p[i] *= _initP;

        }

        int iter = 0;
        int convCriteria = true;

        while ((convCriteria) | (iter < 2)) {

            convCriteria = false;
            int nbIssue = 0;

#pragma omp parallel for default(none) shared(convCriteria)
            for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
                _c_i[i] = Vec3f(0, 0, 0);

                for (tIndex _kernelX = max(0, (int) floor(position(i).x - _rad));
                     _kernelX < min(resX(), (int) floor(position(i).x + _rad + 1)); ++_kernelX) {
                    for (tIndex _kernelY = max(0, (int) floor(position(i).y - _rad));
                         _kernelY < min(resY(), (int) floor(position(i).y + _rad + 1)); ++_kernelY) {
                        for (tIndex _kernelZ = max(0, (int) floor(position(i).z - _rad));
                             _kernelZ < min(resZ(), (int) floor(position(i).z + _rad + 1)); ++_kernelZ) {
                            for (tIndex j: _pidxInGrid[idx1d(_kernelX, _kernelY, _kernelZ)]) {

                                _c_i[i] += _p[j] / (_d[j] * _d[j]) * _kernel.grad_w(position(i) - position(j));
                            }
                        }
                    }
                }

                _c_i[i] *= -_dt * _dt * _m0;
            }

#pragma omp parallel for default(none) shared(convCriteria, nbIssue)
            for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
                Real newP = 0.f;
                if (_a_ii[i] != 0.f) {
                    _predP[i] = _d0 - _interD[i];

                    for (tIndex _kernelX = max(0, (int) floor(position(i).x - _rad));
                         _kernelX < min(resX(), (int) floor(position(i).x + _rad + 1)); ++_kernelX) {
                        for (tIndex _kernelY = max(0, (int) floor(position(i).y - _rad));
                             _kernelY < min(resY(), (int) floor(position(i).y + _rad + 1)); ++_kernelY) {
                            for (tIndex _kernelZ = max(0, (int) floor(position(i).z - _rad));
                                 _kernelZ < min(resZ(), (int) floor(position(i).z + _rad + 1)); ++_kernelZ) {
                                for (tIndex j: _pidxInGrid[idx1d(_kernelX, _kernelY, _kernelZ)]) {
                                    if (j < _nbWallParticles) { //Boundary
                                        _predP[i] -=
                                                _m0 * _wallWeightCoef *
                                                _c_i[i].dotProduct(_kernel.grad_w(position(i) - position(j)));
                                    } else {
                                        _predP[i] -= _m0 * (_c_i[i] - _d_ii[j] * _p[j] - _c_i[j] -
                                                            _dt * _dt * _m0 / (_d[i] * _d[i]) *
                                                            _kernel.grad_w(position(j) - position(i)) *
                                                            _p[i]).dotProduct(
                                                _kernel.grad_w(position(i) - position(j)));
                                    }
                                }
                            }
                        }
                    }
                    newP = max(0.f, (1 - _omega) * _p[i] + _omega / _a_ii[i] * _predP[i]);
                }

                if (abs(_p[i] - newP) > _p[i] * _pressureError) {
                    nbIssue++;
                    convCriteria = true;
                }

                _p[i] = newP;
            }

            iter++;
#ifdef __DEBUG3__
            cout << "NbIssue " << nbIssue << " InterVel " << _interVel[obsPart] << " d_ii " << _d_ii[obsPart] << " a_ii " << _a_ii[obsPart] << " InterDen " << _interD[obsPart] << " c_i " << _c_i[obsPart] << " PredP " << _predP[obsPart] <<" Pr " << _p[obsPart] << endl;
#endif
            if (iter > 500) {
                std::cout << "Pressure solver did not converge" << std::endl;
                break;
            }
        }
#ifdef __DEBUG2__
        cout << "Iterations : " << iter << endl;
#endif

    }


    void applyPressureForce() {
#pragma omp parallel for default(none)
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            Vec3f f = Vec3f(0);
            for (tIndex _kernelX = max(0, (int) floor(position(i).x - _rad));
                 _kernelX < min(resX(), (int) floor(position(i).x + _rad + 1)); ++_kernelX) {
                for (tIndex _kernelY = max(0, (int) floor(position(i).y - _rad));
                     _kernelY < min(resY(), (int) floor(position(i).y + _rad + 1)); ++_kernelY) {
                    for (tIndex _kernelZ = max(0, (int) floor(position(i).z - _rad));
                         _kernelZ < min(resZ(), (int) floor(position(i).z + _rad + 1)); ++_kernelZ) {
                        for (tIndex j: _pidxInGrid[idx1d(_kernelX, _kernelY, _kernelZ)]) {

                            f -= _m0 * (_p[i] / (_d[i] * _d[i]) + _p[j] / (_d[j] * _d[j])) *
                                 _kernel.grad_w(position(i) - position(j));
                        }
                    }
                }
            }
            _acc[i] += f;
        }

#ifdef __DEBUG2__
        cout << "Den " << _d[obsPart] << " Pr " << _p[obsPart] << " Acc " << _acc[obsPart] << " Vit " << velocity(obsPart) <<" Pos " << position(obsPart) << endl;
#endif

    }

    void updateVelocity() {
        maxVel = Vec3f(0);
#ifdef __DEBUG5__
        tIndex maxInd = 0;
#endif
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            _vel[i] += _dt * _acc[i];
            if (_vel[i].lengthSquare() > maxVel.lengthSquare()) {
                maxVel = _vel[i];
#ifdef __DEBUG5__
                maxInd = i;
#endif
            }
        }

#ifdef __DEBUG5__
        cout << "Maximum velocity: " << maxVel << " at index " << maxInd << " (norm = " << maxVel.length() << ")." << endl;
#endif
    }

    void updatePosition() {
#pragma omp parallel for default(none)
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            _pos[i] += _dt * _vel[i];
        }
    }

    // simple collision detection/resolution for each particle
    void resolveCollision() {

#pragma omp parallel for default(none) shared(_pos, _vel)
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            if (particleCollision(_pos[i])) {
                const Vec3f p0 = _pos[i];
                clampParticle(_pos[i]);
                /*if (p0 == _pos[i]) {
                    std::cout << "ALERT FLUID COLLISION IS NOT RESOLVED" << std::endl;
                }*/
                _vel[i] = (_pos[i] - p0) / _dt;
            }
        }
    }

    inline tIndex idx1d(const int i, const int j, const int k) const { return i + j * resX() + k * resX() * resY(); }

    vector<vector<tIndex>> _pidxInGrid; // will help you find neighbor particles

    const CubicSpline _kernel;
    const Real _rad;

    // simulation
    Real _dt;                     // time step

    int _resX, _resY, _resZ;             // background grid resolution

    // wall
    Real _left, _right, _bottom, _top, _back, _front;          // wall (boundary)
    bool _useLConfig = false;
    Real _wallWeightCoef;
    tIndex _nbWallParticles;          // number of particles that belong to the wall

    // SPH coefficients
    Real _nu;                     // viscosity coefficient
    Real _h;                      // particle spacing (i.e., diameter)
    Real _d0;                     // rest density
    Vec3f _g;                     // gravity
    Real _m0;                     // rest mass

    // particle data
    vector<Vec3f> _pos;      // position
    vector<Vec3f> _vel;      // velocity
    vector<Vec3f> _acc;      // acceleration
    vector<Real> _p;         // pressure
    vector<Real> _d;         // density
    Vec3f maxVel;             //maximum velocity

    //IISPH specific data
    vector<Vec3f> _interVel; // intermediate velocity (without pressure force)
    vector<Real> _interD;    // intermediate density (with _predVel)
    vector<Real> _a_ii;      // diagonal coef of the SOE
    vector<Vec3f> _d_ii;     // first level coef of the SOE
    vector<Vec3f> _c_i;      // second level coef of the SOE
    vector<Real> _predP;     // predicted pressure
    const Real _initP;                  // initial pressure for Jacobi method. Must be between 0 and 1
    const Real _omega;                  // relaxation coefficient. Must be between 0 and 1
    const Real _pressureError;           // maximum pressure variation rate serving as a limit of the Jacobi method
};

#endif /* ___IISPH_SOLVER___ */
