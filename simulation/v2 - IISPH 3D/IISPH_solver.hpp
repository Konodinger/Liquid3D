#include <stdio.h>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <math.h>

#ifdef _OPENMP

#include <omp.h>

#endif

#ifndef M_PI
#define M_PI 3.141592
#endif

#include "Vector.hpp"

using namespace std;

// debug index
int obsPart = 2000;

// for recording videos
FILE *ffmpeg;
int *movieBuffer;

// SPH Kernel function: cubic spline
class CubicSpline {
public:
    explicit CubicSpline(const Real h = 1) : _dim(3) {
        setSmoothingLen(h);
    }

    void setSmoothingLen(const Real h) {
        const Real h2 = square(h), h3 = h2 * h;
        _h = h;
        _sr = 2e0 * h;
        _c[0] = 2e0 / (3e0 * h);
        _c[1] = 10e0 / (7e0 * M_PI * h2);
        _c[2] = 1e0 / (M_PI * h3);
        _gc[0] = _c[0] / h;
        _gc[1] = _c[1] / h;
        _gc[2] = _c[2] / h;
    }

    Real smoothingLen() const { return _h; }

    Real supportRadius() const { return _sr; }

    Real f(const Real l) const {
        const Real q = l / _h;
        if (q < 1e0) return _c[_dim - 1] * (1e0 - 1.5 * square(q) + 0.75 * cube(q));
        else if (q < 2e0) return _c[_dim - 1] * (0.25 * cube(2e0 - q));
        return 0;
    }

    Real derivative_f(const Real l) const {
        const Real q = l / _h;
        if (q <= 1e0) return _gc[_dim - 1] * (-3e0 * q + 2.25 * square(q));
        else if (q < 2e0) return -_gc[_dim - 1] * 0.75 * square(2e0 - q);
        return 0;
    }

    Real w(const Vec3f &rij) const { return f(rij.length()); }

    Vec3f grad_w(const Vec3f &rij) const { return grad_w(rij, rij.length()); }

    Vec3f grad_w(const Vec3f &rij, const Real len) const {
        if (len == 0.f) {
            return Vec3f(0, 0, 0);
        }
        return derivative_f(len) * rij / len;
    }

private:
    unsigned int _dim;
    Real _h, _sr, _c[3], _gc[3];
};


class IisphSolver {
public:
    explicit IisphSolver(
            const Real dt = 0.0005, const Real nu = 0.08,
            const Real h = 0.5, const Real density = 1e3,
            const Vec3f g = Vec3f(0, -9.8, 0), const Real initP = 0.5,
            const Real omega = 0.5, const Real pressureError = 1.) :
            _kernel(h), _dt(dt), _nu(nu),
            _h(h), _d0(density), _g(g),
            _initP(initP), _omega(omega), _pressureError(pressureError) {
        _m0 = _d0 * _h * _h * _h / 8. * M_PI * 4. / 3.;
    }

    // assume an arbitrary grid with the size of res_x*res_y; a fluid mass fill up
    // the size of f_width, f_height; each cell is sampled with 2x2 particles.
    void initScene(
            const Vec3f gridRes, const Vec3f initShift, const Vec3f initBlock) {
        if ((initShift.x + initBlock.x + 1 >= gridRes.x) | (initShift.y + initBlock.y + 1 >= gridRes.y) |
            (initShift.z + initBlock.z + 1 >= gridRes.z)) {
            throw length_error("Issue with the size and position of the initial particles.");
        }
        _pos.clear();

        _resX = gridRes.x;
        _resY = gridRes.y;
        _resZ = gridRes.z;

        // set wall for boundary
        _left = 0.5 * _h;
        _right = static_cast<Real>(_resX) - 0.5 * _h;
        _bottom = 0.5 * _h;
        _top = static_cast<Real>(_resY) - 0.5 * _h;
        _back = 0.5 * _h;
        _front = static_cast<Real>(_resZ) - 0.5 * _h;


        _nbWallParticles = 0;
        for (int i: {0, _resX - 1}) {
            for (int j = 0; j < _resY; ++j) {
                for (int k = 0; k < _resZ; ++k) {
                    _pos.push_back(Vec3f(i + 0.25, j + 0.25, k + 0.25));
                    _pos.push_back(Vec3f(i + 0.75, j + 0.25, k + 0.25));
                    _pos.push_back(Vec3f(i + 0.25, j + 0.75, k + 0.25));
                    _pos.push_back(Vec3f(i + 0.75, j + 0.75, k + 0.25));
                    _pos.push_back(Vec3f(i + 0.25, j + 0.25, k + 0.75));
                    _pos.push_back(Vec3f(i + 0.75, j + 0.25, k + 0.75));
                    _pos.push_back(Vec3f(i + 0.25, j + 0.75, k + 0.75));
                    _pos.push_back(Vec3f(i + 0.75, j + 0.75, k + 0.75));
                    _nbWallParticles += 8;
                }
            }
        }


        for (int i = 1; i < _resX - 1; ++i) {
            for (int j: {0, _resY - 1}) {
                for (int k = 0; k < _resZ; ++k) {
                    _pos.push_back(Vec3f(i + 0.25, j + 0.25, k + 0.25));
                    _pos.push_back(Vec3f(i + 0.75, j + 0.25, k + 0.25));
                    _pos.push_back(Vec3f(i + 0.25, j + 0.75, k + 0.25));
                    _pos.push_back(Vec3f(i + 0.75, j + 0.75, k + 0.25));
                    _pos.push_back(Vec3f(i + 0.25, j + 0.25, k + 0.75));
                    _pos.push_back(Vec3f(i + 0.75, j + 0.25, k + 0.75));
                    _pos.push_back(Vec3f(i + 0.25, j + 0.75, k + 0.75));
                    _pos.push_back(Vec3f(i + 0.75, j + 0.75, k + 0.75));
                    _nbWallParticles += 8;
                }
            }
        }

        for (int i = 1; i < _resX - 1; ++i) {
            for (int j = 1; j < _resY - 1; ++j) {
                for (int k: {0, _resZ - 1}) {
                    _pos.push_back(Vec3f(i + 0.25, j + 0.25, k + 0.25));
                    _pos.push_back(Vec3f(i + 0.75, j + 0.25, k + 0.25));
                    _pos.push_back(Vec3f(i + 0.25, j + 0.75, k + 0.25));
                    _pos.push_back(Vec3f(i + 0.75, j + 0.75, k + 0.25));
                    _pos.push_back(Vec3f(i + 0.25, j + 0.25, k + 0.75));
                    _pos.push_back(Vec3f(i + 0.75, j + 0.25, k + 0.75));
                    _pos.push_back(Vec3f(i + 0.25, j + 0.75, k + 0.75));
                    _pos.push_back(Vec3f(i + 0.75, j + 0.75, k + 0.75));
                    _nbWallParticles += 8;
                }
            }
        }

        obsPart += _nbWallParticles;

        // sample a fluid mass
        for (int i = initShift.x; i < initBlock.x + initShift.x; ++i) {
            for (int j = initShift.y; j < initBlock.y + initShift.y; ++j) {
                for (int k = initShift.z; k < initBlock.z + initShift.z; ++k) {
                    _pos.push_back(Vec3f(i + 1.25, j + 1.25, k + 1.25));
                    _pos.push_back(Vec3f(i + 1.75, j + 1.25, k + 1.25));
                    _pos.push_back(Vec3f(i + 1.25, j + 1.75, k + 1.25));
                    _pos.push_back(Vec3f(i + 1.75, j + 1.75, k + 1.25));
                    _pos.push_back(Vec3f(i + 1.25, j + 1.25, k + 1.25));
                    _pos.push_back(Vec3f(i + 1.75, j + 1.25, k + 1.25));
                    _pos.push_back(Vec3f(i + 1.25, j + 1.75, k + 1.25));
                    _pos.push_back(Vec3f(i + 1.75, j + 1.75, k + 1.25));
                }
            }
        }



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

    void update() {
#ifdef __DEBUG1__
        cout << '.' << flush;
#endif

        buildNeighbor();
        computeDensity();

        for (auto &acc: _acc) acc.setFromFloats(0, 0, 0);

        applyBodyForce();
        applyViscousForce();

        // IISMH part
        computeIntermediateVelocity();
        computeDiiCoeff();
        computeAiiCoeff();
        computeIntermediateDensity();

        computePressure(); // Huge loop where we compute c_i and p_i

        //end of IISPH part
        applyPressureForce();

        updateVelocity();
        updatePosition();

        resolveCollision();

    }

    const CubicSpline &getKernel() const { return _kernel; }

    tIndex particleCount() const { return _pos.size(); }

    const tIndex &wallParticleCount() const { return _nbWallParticles; }

    tIndex fluidParticleCount() { return _pos.size() - _nbWallParticles; }

    const vector<tIndex> &gridParticles(const tIndex i, const tIndex j, const tIndex k) {
        return _pidxInGrid[idx1d(i, j, k)];
    };

    const Vec3f &position(const tIndex i) const { return _pos[i]; }

    const Vec3f &velocity(const tIndex i) const { return _vel[i]; }

    const Vec3f &acceleration(const tIndex i) const { return _acc[i]; }

    const Real &pressure(const tIndex i) const { return _p[i]; }

    const Real &density(const tIndex i) const { return _d[i]; }

    int resX() const { return _resX; }

    int resY() const { return _resY; }

    int resZ() const { return _resZ; }

private:
    void buildNeighbor() {

#pragma omp parallel for default(none)
        for (auto &pix: _pidxInGrid) {
            pix.clear();
        }

#pragma omp parallel for default(none)
        for (tIndex i = 0; i < particleCount(); ++i) {
            _pidxInGrid[idx1d(floor(position(i).x), floor(position(i).y), floor(position(i).z))].push_back(i);
        }
    }

    void computeDensity() {
        Real rad = _kernel.supportRadius();

#pragma omp parallel for default(none) shared(rad)
        for (tIndex i = 0; i < particleCount(); ++i) {
            Real rho = 0.f;

            tIndex minX = max(0, (int) floor(position(i).x - rad));
            tIndex maxX = min(resX(), (int) floor(position(i).x + rad + 1));

            tIndex minY = max(0, (int) floor(position(i).y - rad));
            tIndex maxY = min(resY(), (int) floor(position(i).y + rad + 1));

            tIndex minZ = max(0, (int) floor(position(i).z - rad));
            tIndex maxZ = min(resZ(), (int) floor(position(i).z + rad + 1));

            for (tIndex _kernelX = minX; _kernelX < maxX; ++_kernelX) {
                for (tIndex _kernelY = minY; _kernelY < maxY; ++_kernelY) {
                    for (tIndex _kernelZ = minZ; _kernelZ < maxZ; ++_kernelZ) {
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
#pragma omp parallel for default(none)
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            _acc[i] += _g;
        }
    }

    void applyViscousForce() {
        Real rad = _kernel.supportRadius();

#pragma omp parallel for default(none) shared(rad)
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            Vec3f dist = Vec3f();
            Vec3f viscousAcc = Vec3f();
            for (tIndex _kernelX = max(0, (int) floor(position(i).x - rad));
                 _kernelX < min(resX(), (int) floor(position(i).x + rad + 1)); ++_kernelX) {
                for (tIndex _kernelY = max(0, (int) floor(position(i).y - rad));
                     _kernelY < min(resY(), (int) floor(position(i).y + rad + 1)); ++_kernelY) {
                    for (tIndex _kernelZ = max(0, (int) floor(position(i).z - rad));
                         _kernelZ < min(resZ(), (int) floor(position(i).z + rad + 1)); ++_kernelZ) {
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
#pragma omp parallel for default(none)
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            _interVel[i] = _vel[i] + _dt * _acc[i];
        }
    }

    void computeDiiCoeff() {
        Real rad = _kernel.supportRadius();

#pragma omp parallel for default(none) shared(rad)
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            _d_ii[i] = Vec3f(0, 0, 0);

            tIndex minX = max(0, (int) floor(position(i).x - rad));
            tIndex maxX = min(resX(), (int) floor(position(i).x + rad + 1));

            tIndex minY = max(0, (int) floor(position(i).y - rad));
            tIndex maxY = min(resY(), (int) floor(position(i).y + rad + 1));

            tIndex minZ = max(0, (int) floor(position(i).z - rad));
            tIndex maxZ = min(resZ(), (int) floor(position(i).z + rad + 1));

            for (tIndex _kernelX = minX; _kernelX < maxX; ++_kernelX) {
                for (tIndex _kernelY = minY; _kernelY < maxY; ++_kernelY) {
                    for (tIndex _kernelZ = minZ; _kernelZ < maxZ; ++_kernelZ) {
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
        Real rad = _kernel.supportRadius();

#pragma omp parallel for default(none) shared(rad)
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            _a_ii[i] = 0;

            for (tIndex _kernelX = max(0, (int) floor(position(i).x - rad));
                 _kernelX < min(resX(), (int) floor(position(i).x + rad + 1)); ++_kernelX) {
                for (tIndex _kernelY = max(0, (int) floor(position(i).y - rad));
                     _kernelY < min(resY(), (int) floor(position(i).y + rad + 1)); ++_kernelY) {
                    for (tIndex _kernelZ = max(0, (int) floor(position(i).z - rad));
                         _kernelZ < min(resZ(), (int) floor(position(i).z + rad + 1)); ++_kernelZ) {
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
        Real rad = _kernel.supportRadius();

#pragma omp parallel for default(none) shared(rad)
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            _interD[i] = _d[i];

            for (tIndex _kernelX = max(0, (int) floor(position(i).x - rad));
                 _kernelX < min(resX(), (int) floor(position(i).x + rad + 1)); ++_kernelX) {
                for (tIndex _kernelY = max(0, (int) floor(position(i).y - rad));
                     _kernelY < min(resY(), (int) floor(position(i).y + rad + 1)); ++_kernelY) {
                    for (tIndex _kernelZ = max(0, (int) floor(position(i).z - rad));
                         _kernelZ < min(resZ(), (int) floor(position(i).z + rad + 1)); ++_kernelZ) {
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
        Real rad = _kernel.supportRadius();

#pragma omp parallel for default(none) shared(rad)
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            _p[i] *= _initP;

        }

        int l = 0;
        int convCriteria = true;

        while ((convCriteria) | (l < 2)) {

            convCriteria = false;
            int nbIssue = 0;

#pragma omp parallel for default(none) shared(rad, convCriteria, nbIssue)
            for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
                _c_i[i] = Vec3f(0, 0, 0);

                tIndex minX = max(0, (int) floor(position(i).x - rad));
                tIndex maxX = min(resX(), (int) floor(position(i).x + rad + 1));

                tIndex minY = max(0, (int) floor(position(i).y - rad));
                tIndex maxY = min(resY(), (int) floor(position(i).y + rad + 1));

                tIndex minZ = max(0, (int) floor(position(i).z - rad));
                tIndex maxZ = min(resZ(), (int) floor(position(i).z + rad + 1));

                for (tIndex _kernelX = minX; _kernelX < maxX; ++_kernelX) {
                    for (tIndex _kernelY = minY; _kernelY < maxY; ++_kernelY) {
                        for (tIndex _kernelZ = minZ; _kernelZ < maxZ; ++_kernelZ) {
                            for (tIndex j: _pidxInGrid[idx1d(_kernelX, _kernelY, _kernelZ)]) {
                                _c_i[i] += _p[j] / (_d[j] * _d[j]) * _kernel.grad_w(position(i) - position(j));
                            }
                        }
                    }
                }

                _c_i[i] *= -_dt * _dt * _m0;
            }

#pragma omp parallel for default(none) shared(rad, convCriteria, nbIssue)
            for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
                Real newP = 0.f;
                if (_a_ii[i] != 0.f) {
                    _predP[i] = _d0 - _interD[i];

                    tIndex minX = max(0, (int) floor(position(i).x - rad));
                    tIndex maxX = min(resX(), (int) floor(position(i).x + rad + 1));

                    tIndex minY = max(0, (int) floor(position(i).y - rad));
                    tIndex maxY = min(resY(), (int) floor(position(i).y + rad + 1));

                    tIndex minZ = max(0, (int) floor(position(i).z - rad));
                    tIndex maxZ = min(resZ(), (int) floor(position(i).z + rad + 1));

                    for (tIndex _kernelX = minX; _kernelX < maxX; ++_kernelX) {
                        for (tIndex _kernelY = minY; _kernelY < maxY; ++_kernelY) {
                            for (tIndex _kernelZ = minZ; _kernelZ < maxZ; ++_kernelZ) {
                                for (tIndex j: _pidxInGrid[idx1d(_kernelX, _kernelY, _kernelZ)]) {
                                    if (j < _nbWallParticles) { //Boundary
                                        _predP[i] -=
                                                _m0 * _c_i[i].dotProduct(_kernel.grad_w(position(i) - position(j)));
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

            l++;
#ifdef __DEBUG3__
            cout << "NbIssue " << nbIssue << " InterVel " << _interVel[obsPart] << " d_ii " << _d_ii[obsPart] << " a_ii " << _a_ii[obsPart] << " InterDen " << _interD[obsPart] << " c_i " << _c_i[obsPart] << " PredP " << _predP[obsPart] <<" Pr " << _p[obsPart] << endl;
#endif
        }
#ifdef __DEBUG2__
        cout << "Iterations : " << l << endl;
#endif

    }


    void applyPressureForce() {
        // TODO:
        Real rad = _kernel.supportRadius();

#pragma omp parallel for default(none) shared(rad)
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            Vec3f f = Vec3f(0);

            tIndex minX = max(0, (int) floor(position(i).x - rad));
            tIndex maxX = min(resX(), (int) floor(position(i).x + rad + 1));

            tIndex minY = max(0, (int) floor(position(i).y - rad));
            tIndex maxY = min(resY(), (int) floor(position(i).y + rad + 1));

            tIndex minZ = max(0, (int) floor(position(i).z - rad));
            tIndex maxZ = min(resZ(), (int) floor(position(i).z + rad + 1));

            for (tIndex _kernelX = minX; _kernelX < maxX; ++_kernelX) {
                for (tIndex _kernelY = minY; _kernelY < maxY; ++_kernelY) {
                    for (tIndex _kernelZ = minZ; _kernelZ < maxZ; ++_kernelZ) {
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
#pragma omp parallel for default(none)
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            _vel[i] += _dt * _acc[i];
        }
    }

    void updatePosition() {
#pragma omp parallel for default(none)
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            _pos[i] += _dt * _vel[i];
        }
    }

    // simple collision detection/resolution for each particle
    void resolveCollision() {
#pragma omp parallel for default(none)
        for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
            if (_pos[i].x < _left || _pos[i].y < _bottom || _pos[i].z < _back || _pos[i].x > _right ||
                _pos[i].y > _top || _pos[i].z > _front) {
                const Vec3f p0 = _pos[i];
                _pos[i].x = clamp(_pos[i].x, _left, _right);
                _pos[i].y = clamp(_pos[i].y, _bottom, _top);
                _pos[i].z = clamp(_pos[i].z, _back, _front);
                _vel[i] = (_pos[i] - p0) / _dt;
            }
        }
    }

    inline tIndex idx1d(const int i, const int j, const int k) { return i + j * resX() + k * resX() * resY(); }

    vector<vector<tIndex>> _pidxInGrid; // will help you find neighbor particles

    const CubicSpline _kernel;

    // simulation
    Real _dt;                     // time step

    int _resX, _resY, _resZ;             // background grid resolution

    // wall
    Real _left, _right, _bottom, _top, _back, _front;          // wall (boundary)
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
    vector<Real> _p;        // pressure
    vector<Real> _d;        // density

    //IISPH specific data
    vector<Vec3f> _interVel; // intermediate velocity (without pressure force)
    vector<Real> _interD;    // intermediate density (with _predVel)
    vector<Real> _a_ii;      // diagonal coef of the SOE
    vector<Vec3f> _d_ii;     // first level coef of the SOE
    vector<Vec3f> _c_i;      // second level coef of the SOE
    vector<Real> _predP;     // predicted pressure
    Real _initP;                  // initial pressure for Jacobi method. Must be between 0 and 1
    Real _omega;                  // relaxation coefficient. Must be between 0 and 1
    Real _pressureError;           // maximum pressure variation rate serving as a limit of the Jacobi method
};