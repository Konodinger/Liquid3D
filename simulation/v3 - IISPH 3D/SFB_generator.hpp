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

#ifdef _OPENMP
#include <omp.h>
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
            _solver(solver), _resX(_solver->resX()),
            _resY(_solver->resY()), _resZ(_solver->resZ()),
            _lKernel(h), _cKernel(h), _sr(2.f*h), _dt(dt),
            _minBubbleNeighbor(minBubbleNeighbor), _minFoamNeighbor(minFoamNeighbor),
            _kb(kb), _kd(kd)
        {
            sfbList = list<sfb>();
            for (int i = 0; i < _resX * _resY * _resZ; i++) {
                _sfbIndexInGrid.push_back(vector <tIndex> {});
            }
            
        };
    void stfStep() {
            //Step 1: advection and dissolution.
            computeNeighbor();
            cout << sfbList.size() << endl;
            computeVelocityPosition();

            //Step 2: formation.
            sfbFormation();
        };

private:
    void computeNeighbor() {
//#pragma omp parallel for
        for (sfb diffuse : sfbList) {
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

            if (diffuse.fluidNeighbor.size() >= _minBubbleNeighbor) {
                diffuse.nature = BUBBLE;
            } else if (diffuse.fluidNeighbor.size() >= _minFoamNeighbor) {
                diffuse.nature = FOAM;
                diffuse.lifetime -= _dt;
            } else {
                diffuse.nature = SPRAY;
            }
        }

        for (list<sfb>::iterator it = sfbList.begin(); it != sfbList.end();) {
            cout << (*it).nature << endl;
            if (((*it).nature == FOAM) & ((*it).lifetime <= 0.f)) {
                cout << "DEAD";
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
        sfbList.push_back(sfb());

        for (tIndex i = (_solver->wallParticleCount()); i < _solver->particleCount(); ++i) {
            Real scaVelDiff = 0.f;
            //À COMPLETER
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


    inline tIndex idx1d(const int i, const int j, const int k) { return i + j * _resX + k * _resX * _resY; }
    
    const IisphSolver * _solver;
    const int _resX, _resY, _resZ;
    const LinearKernel _lKernel;
    const CubicSpline _cKernel;
    const Real _sr;
    const Real _dt;
    list<sfb> sfbList;
    vector <vector<tIndex>> _sfbIndexInGrid; // will help you find neighbor particles.

    //Advection parameters.
    const int _minBubbleNeighbor;
    const int _minFoamNeighbor;
    const Real _kb;
    const Real _kd;
};

#endif