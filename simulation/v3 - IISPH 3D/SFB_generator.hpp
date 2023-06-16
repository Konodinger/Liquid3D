#ifndef __SFB_GENERATOR_HPP__
#define __SFB_GENERATOR_HPP__

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <list>
#include <cmath>
#include <math.h>
#include <random>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.141592
#endif

#include "Vector.hpp"
#include "Kernel.hpp"
#include "IISPH_solver.hpp"

using namespace std;

enum sfbType {
    SPRAY,
    FOAM,
    BUBBLE
};

struct sfb {
    sfbType nature = SPRAY;
    Vec3f position = Vec3f();
    Vec3f velocity = Vec3f();
    Real lifetime = 1e-9;
    vector<tIndex> fluidNeighbor = vector<tIndex>();
};



class sfbSimulation {
public:
    explicit sfbSimulation(
        const IisphSolver * solver, const Real h = 0.5f, const Real dt = 0.0005,
                const int minBubbleNeighbor = 20, const int minFoamNeighbor = 6,
                const Real kb = 0.3f, const Real kd = 0.7f) :
            _solver(solver),
            _lKernel(h), _cKernel(h), _sr(2*h), _dt(dt),
            _minBubbleNeighbor(minBubbleNeighbor), _minFoamNeighbor(minFoamNeighbor),
            _kb(kb), _kd(kd), _partRad(h/2), _diffuseTimer(dt*5)
        {
            sfbList = list<sfb>();

        random_device _rd;
        mt19937 _gen(_rd());
        uniform_real_distribution<> _randU(0.f, 1.f);
            
        };

    void sfbStep(Real dt) {
        _dt = dt;
        //Step 1: advection and dissolution.
        computeNeighbor();
        computeVelocityPosition();

        //Step 2: formation.
        sfbFormation();
        cout << sfbList.size() << " diffuse particles." << endl;
    };

    void initScene() {
        _resX = _solver->resX();
        _resY = _solver->resY();
        _resZ = _solver->resZ();

        for (int i = 0; i < _resX * _resY * _resZ; i++) {
            _sfbIndexInGrid.push_back(vector <tIndex> {});
        }
    }

    tIndex particleCount() const { return sfbList.size(); }

    const list<sfb> *diffuseList() const { return &sfbList; }

private:
    void computeNeighbor() {
        for (sfb &diffuse : sfbList) {
            diffuse.fluidNeighbor = vector<tIndex>();

            for (tIndex _kernelX = max(0, (int) floor(diffuse.position.x - _sr));
                    _kernelX < min(_resX, (int) floor(diffuse.position.x + _sr + 1)); ++_kernelX) {
                for (tIndex _kernelY = max(0, (int) floor(diffuse.position.y - _sr));
                        _kernelY < min(_resY, (int) floor(diffuse.position.y + _sr + 1)); ++_kernelY) {
                    for (tIndex _kernelZ = max(0, (int) floor(diffuse.position.z - _sr));
                            _kernelZ < min(_resZ, (int) floor(diffuse.position.z + _sr + 1)); ++_kernelZ) {
                        for (tIndex f : (_solver->gridParticles(_kernelX, _kernelY, _kernelZ))) {

                            if (_lKernel.w(diffuse.position - _solver->position(f)) > 0.f) {
                                diffuse.fluidNeighbor.push_back(f);
                            }
                        }
                    }
                }
            }

            // if((diffuse.fluidNeighbor.size() > 0) && (diffuse.fluidNeighbor.size() != 4)) {
            //     cout << diffuse.fluidNeighbor.size() << endl;
            // }

            if (diffuse.fluidNeighbor.size() >= _minBubbleNeighbor) {
                diffuse.nature = BUBBLE;
            } else if (diffuse.fluidNeighbor.size() >= _minFoamNeighbor) {
                diffuse.nature = FOAM;
                diffuse.lifetime -= _dt;
            } else {
                diffuse.nature = SPRAY;
                diffuse.lifetime -= _dt*_sprayResistance;
            }
        }

        for (list<sfb>::iterator it = sfbList.begin(); it != sfbList.end();) {
            if ((*it).lifetime <= 0.f) {
                it = sfbList.erase(it);
            } else {
                ++it;
            }
        }

    };

    void computeVelocityPosition() {

        for (sfb diffuse : sfbList) {
            if (diffuse.nature == SPRAY) {
                diffuse.velocity += _dt*_solver->gravity();
                diffuse.position += _dt*diffuse.velocity;

            } else {
                Vec3f avLocFluVel = Vec3f();
                Real locFluDensity = 0;
                for (tIndex f : diffuse.fluidNeighbor) {
                    avLocFluVel += _solver->velocity(f) * _cKernel.w(diffuse.position - _solver->position(f));
                    locFluDensity += _cKernel.w(diffuse.position - _solver->position(f));
                }
                if (!locFluDensity) {
                    avLocFluVel /= locFluDensity;
                }

                if (diffuse.nature == BUBBLE) {
                    diffuse.velocity += _kd * (avLocFluVel - diffuse.velocity) - _dt * _kb * _solver->gravity();
                    diffuse.position += _dt*diffuse.velocity;
                } else { //If the particle is foam.
                    diffuse.position += _dt*avLocFluVel;
                }
            }
                
        }

    };

    void sfbFormation() {

        for (tIndex i = (_solver->wallParticleCount()); i < _solver->particleCount(); ++i) {
            Vec3f posI = _solver->position(i);
            Vec3f velI = _solver->velocity(i);
            Real kineticPotential = _clamp(0.5*_solver->particleMass()*velI.lengthSquare(), _tauEnergyMin, _tauEnergyMax);
            
            Real scaVelDiff = 0.f;

            Vec3f normal(0);
            Real cumul = 0.;
            //À COMPLETER
            for (tIndex _kernelX = max(0, (int) floor(posI.x - _sr));
                 _kernelX < min(_resX, (int) floor(posI.x + _sr + 1)); ++_kernelX) {
                for (tIndex _kernelY = max(0, (int) floor(posI.y - _sr));
                     _kernelY < min(_resY, (int) floor(posI.y + _sr + 1)); ++_kernelY) {
                    for (tIndex _kernelZ = max(0, (int) floor(posI.z - _sr));
                         _kernelZ < min(_resZ, (int) floor(posI.z + _sr + 1)); ++_kernelZ) {
                        for (tIndex j: _solver->gridParticles(_kernelX, _kernelY, _kernelZ)) {

                            Vec3f pos_ij = _solver->position(j) - posI;
                            if (pos_ij.length() > 1e-5) {
                                Vec3f vel_ij = _solver->velocity(j) - velI;
                                scaVelDiff += vel_ij.length() * (1 - vel_ij.normalized().dotProduct(pos_ij.normalized())) * _lKernel.w(pos_ij);

                                normal -= pos_ij;
                            }
                        }
                    }
                }
            }

            Real trappedPotential = _clamp(scaVelDiff, _tauTrappedMin, _tauTrappedMax);
            Real crestPotential = 0.f;

            int particleGen = int(_dt * kineticPotential * (trappedPotential*_generateTrapped + crestPotential*_generateCrest));
            if (i == (_solver->wallParticleCount())) {
                cout << "Energy " << kineticPotential << " Scaled velocity dif " << scaVelDiff << " Generate " << particleGen << " particles." << endl;
            }
            
            if (particleGen > 0) {
                Vec3f partXaxis(1.f, 0.f, 0.f);
                Vec3f partYaxis(0.f, 1.f, 0.f);
                Vec3f partZaxis(0.f, 0.f, 1.f);
                if (velI.length() != 0.f) {
                    Vec3f partZaxis = velI.normalized();
                    Vec3f partXaxis = partZaxis.crossProduct(Vec3f(1.f, 0.f, 0.f)).normalize();
                    if (partXaxis.length() < 1e-9) {
                        partXaxis = partZaxis.crossProduct(Vec3f(0.f, 1.f, 0.f)).normalize();
                    }
                    partYaxis = partZaxis.crossProduct(partXaxis).normalize();

                }

                for (int d = 0; d < particleGen; ++d) {
                    
                    Real randR = _randU(_gen);
                    Real randT = _randU(_gen);
                    Real randH = _randU(_gen);
                    Real r = _partRad*sqrtf(randR);
                    Real theta = randT*M_PI*2;
                    Real h = randH*_dt*velI.length();

                    sfb diffPart;
                    diffPart.position = posI + r*cosf(theta)*partXaxis + r*sinf(theta)*partYaxis + h*partZaxis;
                    diffPart.velocity = velI + r*cosf(theta)*partXaxis + r*sin(theta)*partYaxis;
                    diffPart.lifetime = _diffuseTimer;
                    sfbList.push_back(diffPart);
                }
            }

        }
    };

    void resolveCollision() {
        for (sfb diffuse : sfbList) {
            if (_solver->particleCollision(diffuse.position)) {
                const Vec3f part0 = diffuse.position;
                _solver->clampParticle(diffuse.position);
                if (part0 == diffuse.position) {
                    cout << "ALERT DIFFUSE COLLISION IS NOT RESOLVED" << endl;
                }
                diffuse.velocity = (diffuse.position - part0) / _dt;
            }
        }
    }

    int diffuseCount() {return sfbList.size();}


    inline Real _clamp(const Real value, const Real tauMin, const Real tauMax) {return (min(value, tauMax) - min(value, tauMin))/(tauMax - tauMin);}
    inline tIndex idx1d(const int i, const int j, const int k) { return i + j * _resX + k * _resX * _resY; }
    
    const IisphSolver * _solver;
    int _resX, _resY, _resZ;
    const LinearKernel _lKernel;
    const CubicSpline _cKernel;
    const Real _sr;
    Real _dt;
    list<sfb> sfbList;
    vector <vector<tIndex>> _sfbIndexInGrid; // will help you find neighbor particles.

    //Advection parameters.
    long long unsigned int _minBubbleNeighbor;
    long long unsigned int _minFoamNeighbor;
    const Real _kb;
    const Real _kd;

    //Generation parameters.
    const Real _partRad;
    const Real _tauTrappedMin = 10.;
    const Real _tauTrappedMax = 100.;
    const Real _generateTrapped = 150.;
    const Real _tauCrestMin = 0.;
    const Real _tauCrestMax = 1.;
    const Real _generateCrest = 100.;
    const Real _tauEnergyMin = 5.;
    const Real _tauEnergyMax = 50.;
    const Real _diffuseTimer;
    const Real _sprayResistance = 0.2;

    mt19937 _gen;
    uniform_real_distribution<> _randU;
};

#endif