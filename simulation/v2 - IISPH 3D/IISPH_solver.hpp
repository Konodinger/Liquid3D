#include <stdio.h>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <math.h>
#ifdef __OPEN_MP__
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.141592
#endif

#include "Vector.hpp"

using namespace std;

// debug index
int obsPart = 82000;

// for recording videos
FILE* ffmpeg;
int* movieBuffer;

// SPH Kernel function: cubic spline
class CubicSpline {
public:
  explicit CubicSpline(const Real h=1) : _dim(3)
  {
    setSmoothingLen(h);
  }
  void setSmoothingLen(const Real h)
  {
    const Real h2 = square(h), h3 = h2*h;
    _h = h;
    _sr = 2e0*h;
    _c[0]  = 2e0/(3e0*h);
    _c[1]  = 10e0/(7e0*M_PI*h2);
    _c[2]  = 1e0/(M_PI*h3);
    _gc[0] = _c[0]/h;
    _gc[1] = _c[1]/h;
    _gc[2] = _c[2]/h;
  }
  Real smoothingLen() const { return _h; }
  Real supportRadius() const { return _sr; }

  Real f(const Real l) const
  {
    const Real q = l/_h;
    if(q<1e0) return _c[_dim-1]*(1e0 - 1.5*square(q) + 0.75*cube(q));
    else if(q<2e0) return _c[_dim-1]*(0.25*cube(2e0-q));
    return 0;
  }
  Real derivative_f(const Real l) const
  {
    const Real q = l/_h;
    if(q<=1e0) return _gc[_dim-1]*(-3e0*q+2.25*square(q));
    else if(q<2e0) return -_gc[_dim-1]*0.75*square(2e0-q);
    return 0;
  }

  Real w(const Vec3f &rij) const { return f(rij.length()); }
  Vec3f grad_w(const Vec3f &rij) const { return grad_w(rij, rij.length()); }
  Vec3f grad_w(const Vec3f &rij, const Real len) const
  {
    if (len == 0.f) {
      return Vec3f(0,0,0);
    }
    return derivative_f(len)*rij/len;
  }

private:
  unsigned int _dim;
  Real _h, _sr, _c[3], _gc[3];
};



class SphSolver {
public:
  explicit SphSolver(
    const Real nu=0.08, const Real h=0.5, const Real density=1e3,
    const Vec3f g=Vec3f(0, -9.8, 0), const Real eta=0.01, const Real gamma=7.0,
    const Real initP = 0.5, const Real omega = 0.5, const Real pressureError = 1.) :
    _kernel(h), _nu(nu), _h(h), _d0(density),
    _g(g), _eta(eta), _gamma(gamma),
    _initP(initP), _omega(omega), _pressureError(pressureError)
  {
    _dt = 0.0005;
    _m0 = _d0*_h*_h;
    _c = fabs(_g.y)/_eta;
    _k = _d0*_c*_c/_gamma;
  }

  // assume an arbitrary grid with the size of res_x*res_y; a fluid mass fill up
  // the size of f_width, f_height; each cell is sampled with 2x2 particles.
  void initScene(
    const int res_x, const int res_y, const int res_z, const int f_length, const int f_height, const int f_width)
  {
    _pos.clear();

    _resX = res_x;
    _resY = res_y;
    _resZ = res_z;

    // set wall for boundary
    _left = 0.5*_h;
    _right = static_cast<Real>(res_x) - 0.5*_h;
    _bottom = 0.5*_h;
    _top = static_cast<Real>(res_y) - 0.5*_h;
    _back = 0.5*_h;
    _front = static_cast<Real>(res_z) - 0.5*_h;


    _nbWallParticles = 0;
    for (int i : {0, res_x - 1}) {
      for (int j = 0; j<res_y; ++j) {
        for (int k = 0; k<res_z; ++k) {
          _pos.push_back(Vec3f(i+0.25, j+0.25, k+0.25));
          _pos.push_back(Vec3f(i+0.75, j+0.25, k+0.25));
          _pos.push_back(Vec3f(i+0.25, j+0.75, k+0.25));
          _pos.push_back(Vec3f(i+0.75, j+0.75, k+0.25));
          _pos.push_back(Vec3f(i+0.25, j+0.25, k+0.75));
          _pos.push_back(Vec3f(i+0.75, j+0.25, k+0.75));
          _pos.push_back(Vec3f(i+0.25, j+0.75, k+0.75));
          _pos.push_back(Vec3f(i+0.75, j+0.75, k+0.75));
          _nbWallParticles += 8;
        }
      }
    }

    
    for (int i = 1; i<res_x - 1; ++i) {
      for (int j : {0, res_y - 1}) {
        for (int k = 0; k<res_z; ++k) {
          _pos.push_back(Vec3f(i+0.25, j+0.25, k+0.25));
          _pos.push_back(Vec3f(i+0.75, j+0.25, k+0.25));
          _pos.push_back(Vec3f(i+0.25, j+0.75, k+0.25));
          _pos.push_back(Vec3f(i+0.75, j+0.75, k+0.25));
          _pos.push_back(Vec3f(i+0.25, j+0.25, k+0.75));
          _pos.push_back(Vec3f(i+0.75, j+0.25, k+0.75));
          _pos.push_back(Vec3f(i+0.25, j+0.75, k+0.75));
          _pos.push_back(Vec3f(i+0.75, j+0.75, k+0.75));
          _nbWallParticles += 8;
        }
      }
    }

    for (int i = 1; i<res_x - 1; ++i) {
      for (int j = 1; j<res_y - 1; ++j) {
        for (int k : {0, res_z - 1}) {
          _pos.push_back(Vec3f(i+0.25, j+0.25, k+0.25));
          _pos.push_back(Vec3f(i+0.75, j+0.25, k+0.25));
          _pos.push_back(Vec3f(i+0.25, j+0.75, k+0.25));
          _pos.push_back(Vec3f(i+0.75, j+0.75, k+0.25));
          _pos.push_back(Vec3f(i+0.25, j+0.25, k+0.75));
          _pos.push_back(Vec3f(i+0.75, j+0.25, k+0.75));
          _pos.push_back(Vec3f(i+0.25, j+0.75, k+0.75));
          _pos.push_back(Vec3f(i+0.75, j+0.75, k+0.75));
          _nbWallParticles += 8;
        }
      }
    }

    // sample a fluid mass
    for(int i=0; i<f_length; ++i) {
      for(int j=0; j<f_height; ++j) {
        for(int k=0; k<f_width; ++k) {
          _pos.push_back(Vec3f(i+1.25, j+1.25, k+1.25));
          _pos.push_back(Vec3f(i+1.75, j+1.25, k+1.25));
          _pos.push_back(Vec3f(i+1.25, j+1.75, k+1.25));
          _pos.push_back(Vec3f(i+1.75, j+1.75, k+1.25));
          _pos.push_back(Vec3f(i+1.25, j+1.25, k+1.25));
          _pos.push_back(Vec3f(i+1.75, j+1.25, k+1.25));
          _pos.push_back(Vec3f(i+1.25, j+1.75, k+1.25));
          _pos.push_back(Vec3f(i+1.75, j+1.75, k+1.25));
        }
      }
    }
    


    // make sure for the other particle quantities
    _vel = vector<Vec3f>(_pos.size(), Vec3f(0, 0, 0));
    _acc = vector<Vec3f>(_pos.size(), Vec3f(0, 0, 0));
    _p   = vector<Real>(_pos.size(), 100000);
    for (int i = 0; i < _nbWallParticles; ++i) {
      _p[i] = 0;
    }
    _d   = vector<Real>(_pos.size(), 0);
    _interVel = vector<Vec3f>(_pos.size(), Vec3f(0, 0, 0));
    _interD = vector<Real>(_pos.size(), 0);
    _a_ii = vector<Real>(_pos.size(), 0);
    _d_ii = vector<Vec3f>(_pos.size(), Vec3f(0, 0, 0));
    _c_i = vector<Vec3f>(_pos.size(), Vec3f(0, 0, 0));
    _predP = vector<Real>(_pos.size(), 0);


    for (int i = 0; i < res_x * res_y * res_z; i++) {
      _pidxInGrid.push_back(vector<tIndex>{});
    }
    
    buildNeighbor();
    computeDensity();
  }

  void update()
  {
    #ifdef __DEBUG1__
    cout << '.' << flush;
    #endif

    buildNeighbor();
    computeDensity();

    _acc = vector<Vec3f>(_pos.size(), Vec3f(0, 0, 0));
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

  const CubicSpline &getKernel() const {return _kernel;}
  tIndex particleCount() const { return _pos.size(); }
  const tIndex &wallParticleCount() { return _nbWallParticles; }
  tIndex fluidParticleCount() { return _pos.size() - _nbWallParticles; }
  const vector<tIndex> &gridParticles(const tIndex i, const tIndex j, const tIndex k) {return _pidxInGrid[idx1d(i, j, k)];};
  const Vec3f &position(const tIndex i) const { return _pos[i]; }
  const Vec3f &velocity(const tIndex i) const { return _vel[i]; }
  const Vec3f &acceleration(const tIndex i) const { return _acc[i]; }
  const Real &pressure(const tIndex i) const { return _p[i]; }
  const Real &density(const tIndex i) const { return _d[i]; }

  int resX() const { return _resX; }
  int resY() const { return _resY; }
  int resZ() const { return _resZ; }

private:
  void buildNeighbor()
  {
    // TODO:
    for (auto &pix : _pidxInGrid) {
        pix.clear();
      }

    for(tIndex i=0; i<particleCount(); ++i) {
        _pidxInGrid[idx1d(floor(position(i).x), floor(position(i).y), floor(position(i).z))].push_back(i);
      }
  }

  void computeDensity()
  {
    Real rad = _kernel.supportRadius();
    int nb;

#ifdef __OPEN_MP__
    #pragma omp parallel for 
#endif
    for (tIndex i = 0; i < particleCount(); ++i) {
      Real rho = 0.f;
      nb = 0;
      Real test = 0.f;

      for (tIndex _kernelX = max(0, (int) floor(position(i).x - rad)); _kernelX < min(resX(), (int) floor(position(i).x + rad + 1)); ++_kernelX) {
        for (tIndex _kernelY = max(0, (int) floor(position(i).y - rad)); _kernelY < min(resY(), (int) floor(position(i).y + rad + 1)); ++_kernelY) {
          for (tIndex _kernelZ = max(0, (int) floor(position(i).z - rad)); _kernelZ < min(resZ(), (int) floor(position(i).z + rad + 1)); ++_kernelZ) {
          
            for (tIndex j : _pidxInGrid[idx1d(_kernelX, _kernelY, _kernelZ)]) {
              nb++;
              test += _kernel.w(position(i) - position(j));
              rho += _m0 * _kernel.w(position(i) - position(j));
            }
          }
        }
      }
      #ifdef __DEBUG3__
      cout << rho << "   " << nb << "    " << test << endl;
      #endif
      _d[i] = rho;
    }
    #ifdef __DEBUG3__
    cout << " " << nb << "    " << _d[particleCount() - 1] << "   ";
    #endif
  }

  void applyBodyForce()
  {
    #ifdef __OPEN_MP__
    #pragma omp parallel for
    #endif 
    for(tIndex i=_nbWallParticles; i<particleCount(); ++i) {
      _acc[i] += _g;
    }
  }

  void applyViscousForce()
  {
    Real rad = _kernel.supportRadius();

    #ifdef __OPEN_MP__
    #pragma omp parallel for 
    #endif
    for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
      Vec3f dist = Vec3f();
      Vec3f viscousAcc = Vec3f();
      for (tIndex _kernelX = max(0, (int) floor(position(i).x - rad)); _kernelX < min(resX(), (int) floor(position(i).x + rad + 1)); ++_kernelX) {
        for (tIndex _kernelY = max(0, (int) floor(position(i).y - rad)); _kernelY < min(resY(), (int) floor(position(i).y + rad + 1)); ++_kernelY) {
          for (tIndex _kernelZ = max(0, (int) floor(position(i).z - rad)); _kernelZ < min(resZ(), (int) floor(position(i).z + rad + 1)); ++_kernelZ) {
            for (tIndex j : _pidxInGrid[idx1d(_kernelX, _kernelY, _kernelZ)]) {
              dist = position(i) - position(j);
              viscousAcc += _m0 / _d[j] * (velocity(i) - velocity(j)) * (dist.dotProduct(_kernel.grad_w(dist))) / (dist.lengthSquare() + _h*_h);
            }
          }
        }
      }
      _acc[i] += 2 * _nu * viscousAcc;
    }
    #ifdef __DEBUG3__
    cout << 2 * _nu * viscousAcc << "    ";
    cout << "D " << _d[100] << " A " << _acc[100] << " V " << velocity(100) <<" P " << position(100) << endl;
    #endif
  }

  void computeIntermediateVelocity()
  {
    #ifdef __OPEN_MP__
    #pragma omp parallel for 
    #endif
    for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
      _interVel[i] = _vel[i] + _dt * _acc[i];
    }
  }

  void computeDiiCoeff()
  {
    Real rad = _kernel.supportRadius();

    #ifdef __OPEN_MP__
    #pragma omp parallel for 
    #endif
    for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
      _d_ii[i] = Vec3f(0, 0, 0);

      for (tIndex _kernelX = max(0, (int) floor(position(i).x - rad)); _kernelX < min(resX(), (int) floor(position(i).x + rad + 1)); ++_kernelX) {
        for (tIndex _kernelY = max(0, (int) floor(position(i).y - rad)); _kernelY < min(resY(), (int) floor(position(i).y + rad + 1)); ++_kernelY) {
          for (tIndex _kernelZ = max(0, (int) floor(position(i).z - rad)); _kernelZ < min(resZ(), (int) floor(position(i).z + rad + 1)); ++_kernelZ) {
            for (tIndex j : _pidxInGrid[idx1d(_kernelX, _kernelY, _kernelZ)]) {

              _d_ii[i] += _m0 * _kernel.grad_w(position(i) - position(j));
            }
          }
        }
      }

      _d_ii[i] *= -_dt*_dt /(_d[i]*_d[i]);
      
    }
  }

  void computeAiiCoeff()
  {
    Real rad = _kernel.supportRadius();

    #ifdef __OPEN_MP__
    #pragma omp parallel for 
    #endif
    for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
      _a_ii[i] = 0;

      for (tIndex _kernelX = max(0, (int) floor(position(i).x - rad)); _kernelX < min(resX(), (int) floor(position(i).x + rad + 1)); ++_kernelX) {
        for (tIndex _kernelY = max(0, (int) floor(position(i).y - rad)); _kernelY < min(resY(), (int) floor(position(i).y + rad + 1)); ++_kernelY) {
          for (tIndex _kernelZ = max(0, (int) floor(position(i).z - rad)); _kernelZ < min(resZ(), (int) floor(position(i).z + rad + 1)); ++_kernelZ) {
            for (tIndex j : _pidxInGrid[idx1d(_kernelX, _kernelY, _kernelZ)]) {

              _a_ii[i] += _m0 * (_d_ii[i] + _dt * _dt * _m0 / (_d[i] * _d[i]) * _kernel.grad_w(position(j) - position(i))).dotProduct(_kernel.grad_w(position(i) - position(j)));
            }
          }
        }
      }

      
    }
  }

  void computeIntermediateDensity()
  {
    Real rad = _kernel.supportRadius();

    #ifdef __OPEN_MP__
    #pragma omp parallel for 
    #endif
    for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
      _interD[i] = _d[i];

      for (tIndex _kernelX = max(0, (int) floor(position(i).x - rad)); _kernelX < min(resX(), (int) floor(position(i).x + rad + 1)); ++_kernelX) {
        for (tIndex _kernelY = max(0, (int) floor(position(i).y - rad)); _kernelY < min(resY(), (int) floor(position(i).y + rad + 1)); ++_kernelY) {
          for (tIndex _kernelZ = max(0, (int) floor(position(i).z - rad)); _kernelZ < min(resZ(), (int) floor(position(i).z + rad + 1)); ++_kernelZ) {
            for (tIndex j : _pidxInGrid[idx1d(_kernelX, _kernelY, _kernelZ)]) {
              
              _interD[i] += _dt * _m0 * (_interVel[i] - _interVel[j]).dotProduct(_kernel.grad_w(position(i) - position(j)));
            }
          }
        }
      }
      
    }
  }

  void computePressure()
  {
    Real rad = _kernel.supportRadius();

    #ifdef __OPEN_MP__
    #pragma omp parallel for 
    #endif
    for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
      _p[i] *= _initP;
      
    }

    int l = 0;
    int convCriteria = true;

    while ((convCriteria) | (l < 2)) {

      convCriteria = false;

      #ifdef __OPEN_MP__
      #pragma omp parallel for
      #endif
      for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
        _c_i[i] = Vec3f(0, 0, 0);

        for (tIndex _kernelX = max(0, (int) floor(position(i).x - rad)); _kernelX < min(resX(), (int) floor(position(i).x + rad + 1)); ++_kernelX) {
          for (tIndex _kernelY = max(0, (int) floor(position(i).y - rad)); _kernelY < min(resY(), (int) floor(position(i).y + rad + 1)); ++_kernelY) {
            for (tIndex _kernelZ = max(0, (int) floor(position(i).z - rad)); _kernelZ < min(resZ(), (int) floor(position(i).z + rad + 1)); ++_kernelZ) {
              for (tIndex j : _pidxInGrid[idx1d(_kernelX, _kernelY, _kernelZ)]) {
                
                _c_i[i] += _p[j] / (_d[j] * _d[j]) * _kernel.grad_w(position(i) - position(j));
              }
            }
          }
        }

        _c_i[i] *= -_dt * _dt * _m0;
      }

      #ifdef __OPEN_MP__
      #pragma omp parallel for
      #endif
      for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
        Real newP;
        if (_a_ii[i] == 0.f) {
          newP == 0.f;
        } else {
          _predP[i] = _d0 - _interD[i];

          for (tIndex _kernelX = max(0, (int) floor(position(i).x - rad)); _kernelX < min(resX(), (int) floor(position(i).x + rad + 1)); ++_kernelX) {
            for (tIndex _kernelY = max(0, (int) floor(position(i).y - rad)); _kernelY < min(resY(), (int) floor(position(i).y + rad + 1)); ++_kernelY) {
              for (tIndex _kernelZ = max(0, (int) floor(position(i).z - rad)); _kernelZ < min(resZ(), (int) floor(position(i).z + rad + 1)); ++_kernelZ) {
                for (tIndex j : _pidxInGrid[idx1d(_kernelX, _kernelY, _kernelZ)]) {
                  if (j < _nbWallParticles) { //Boundary
                    _predP[i] -= _m0 * _c_i[i].dotProduct(_kernel.grad_w(position(i) - position(j)));
                  } else {
                    _predP[i] -= _m0 * (_c_i[i] - _d_ii[j] * _p[j] - _c_i[j] - _dt * _dt * _m0 / (_d[i] * _d[i]) * _kernel.grad_w(position(j) - position(i)) * _p[i]).dotProduct(_kernel.grad_w(position(i) - position(j)));
                  }
                }
              }
            }
          }
          newP = max(0.f, (1 - _omega) * _p[i] + _omega / _a_ii[i] * _predP[i]);
        }

        if (abs(_p[i] - newP) > _p[i] * _pressureError) {
          convCriteria = true;
        }

        _p[i] = newP;
        
      }

      l++;
      #ifdef __DEBUG3__
      cout << "InterVel " << _interVel[obsPart] << " d_ii " << _d_ii[obsPart] << " a_ii " << _a_ii[obsPart] << " InterDen " << _interD[obsPart] << " c_i " << _c_i[obsPart] << " PredP " << _predP[obsPart] <<" Pr " << _p[obsPart] << endl;
      #endif
    }
    #ifdef __DEBUG2__
    cout << "Iterations : " << l << endl;
    #endif

  }


  void applyPressureForce()
  {
    // TODO:
    Real rad = _kernel.supportRadius();

    #ifdef __OPEN_MP__
    #pragma omp parallel for
    #endif 
    for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
      Vec3f f = Vec3f(0);
      for (tIndex _kernelX = max(0, (int) floor(position(i).x - rad)); _kernelX < min(resX(), (int) floor(position(i).x + rad + 1)); ++_kernelX) {
        for (tIndex _kernelY = max(0, (int) floor(position(i).y - rad)); _kernelY < min(resY(), (int) floor(position(i).y + rad + 1)); ++_kernelY) {
          for (tIndex _kernelZ = max(0, (int) floor(position(i).z - rad)); _kernelZ < min(resZ(), (int) floor(position(i).z + rad + 1)); ++_kernelZ) {
            for (tIndex j : _pidxInGrid[idx1d(_kernelX, _kernelY, _kernelZ)]) {
              
              f -= _m0 * (_p[i]/(_d[i] * _d[i]) + _p[j]/(_d[j] * _d[j])) * _kernel.grad_w(position(i) - position(j));
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

  void updateVelocity()
  {
    #ifdef __OPEN_MP__
    #pragma omp parallel for 
    #endif
    for(tIndex i=_nbWallParticles; i<particleCount(); ++i) {
      _vel[i] += _dt * _acc[i];
    }
  }

  void updatePosition()
  {
    #ifdef __OPEN_MP__
    #pragma omp parallel for 
    #endif
    for(tIndex i=_nbWallParticles; i<particleCount(); ++i) {
      _pos[i] += _dt * _vel[i];
    }
  }

  // simple collision detection/resolution for each particle
  void resolveCollision()
  {
    vector<tIndex> need_res;
    #ifdef __OPEN_MP__
    #pragma omp parallel for
    #endif 
    for(tIndex i=_nbWallParticles; i<particleCount(); ++i) {
      if(_pos[i].x<_left || _pos[i].y<_bottom || _pos[i].z<_back|| _pos[i].x>_right || _pos[i].y>_top || _pos[i].z>_front)
        need_res.push_back(i);
    }

    #ifdef __OPEN_MP__
    #pragma omp parallel for
    #endif
    for(
      vector<tIndex>::const_iterator it=need_res.begin();
      it<need_res.end();
      ++it) {
      const Vec3f p0 = _pos[*it];
      _pos[*it].x = clamp(_pos[*it].x, _left, _right);
      _pos[*it].y = clamp(_pos[*it].y, _bottom, _top);
      _pos[*it].z = clamp(_pos[*it].z, _back, _front);
      _vel[*it] = (_pos[*it] - p0)/_dt;
    }
  }

  inline tIndex idx1d(const int i, const int j, const int k) { return i + j*resX() + k*resX()*resY(); }

  const CubicSpline _kernel;

  // particle data
  vector<Vec3f> _pos;      // position
  vector<Vec3f> _vel;      // velocity
  vector<Vec3f> _acc;      // acceleration
  vector<Real>  _p;        // pressure
  vector<Real>  _d;        // density
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



  vector< vector<tIndex> > _pidxInGrid; // will help you find neighbor particles

  // simulation
  Real _dt;                     // time step

  int _resX, _resY, _resZ;             // background grid resolution

  // wall
  Real _left, _right, _bottom, _top, _back, _front;          // wall (boundary)
  tIndex _nbWallParticles;          // number of particles that belong to the wall

  // SPH coefficients
  Real _nu;                     // viscosity coefficient
  Real _d0;                     // rest density
  Real _h;                      // particle spacing (i.e., diameter)
  Vec3f _g;                     // gravity

  Real _m0;                     // rest mass
  Real _k;                      // EOS coefficient

  Real _eta;
  Real _c;                      // speed of sound
  Real _gamma;                  // EOS power factor
};