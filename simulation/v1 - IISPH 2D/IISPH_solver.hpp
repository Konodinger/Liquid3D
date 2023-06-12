#include <stdio.h>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <math.h>
#include <omp.h>

#ifndef M_PI
#define M_PI 3.141592
#endif

#include "Vector.hpp"

using namespace std;

// debug index
int obsPart = 685;

// for recording videos
FILE* ffmpeg;
int* movieBuffer;

// SPH Kernel function: cubic spline
class CubicSpline {
public:
  explicit CubicSpline(const Real h=1) : _dim(2)
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

  Real w(const Vec2f &rij) const { return f(rij.length()); }
  Vec2f grad_w(const Vec2f &rij) const { return grad_w(rij, rij.length()); }
  Vec2f grad_w(const Vec2f &rij, const Real len) const
  {
    if (len == 0.f) {
      return Vec2f(0,0);
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
    const Vec2f g=Vec2f(0, -9.8), const Real eta=0.01, const Real gamma=7.0,
    const Real initP = 0.5, const Real omega = 0.5, const Real pressureError = 1.) :
    _kernel(h), _rad(_kernel.supportRadius()), _nu(nu), _h(h), _d0(density),
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
    const int res_x, const int res_y, const int f_width, const int f_height)
  {
    _pos.clear();

    _resX = res_x;
    _resY = res_y;

    // set wall for boundary
    _l = 0.5*_h;
    _r = static_cast<Real>(res_x) - 0.5*_h;
    _b = 0.5*_h;
    _t = static_cast<Real>(res_y) - 0.5*_h;


    _nbWallParticles = 0;
    for (int j = 0; j<res_y; ++j) {
      for (int i : {0, res_x - 1}) {
        _pos.push_back(Vec2f(i+0.25, j+0.25));
        _pos.push_back(Vec2f(i+0.75, j+0.25));
        _pos.push_back(Vec2f(i+0.25, j+0.75));
        _pos.push_back(Vec2f(i+0.75, j+0.75));
        _nbWallParticles += 4;
      }
    }

    
    for (int i = 1; i<res_x - 1; ++i) {
      for (int j : {0, res_y - 1}) {
        _pos.push_back(Vec2f(i+0.25, j+0.25));
        _pos.push_back(Vec2f(i+0.75, j+0.25));
        _pos.push_back(Vec2f(i+0.25, j+0.75));
        _pos.push_back(Vec2f(i+0.75, j+0.75));
        _nbWallParticles += 4;
      }
    }
    
    //Add a small wall. Beware to not superpose it with the fluid particles.
    /*for (int j = 1; j<3; ++j) {
      int i = 20;
      _pos.push_back(Vec2f(i+0.25, j+0.25));
      _pos.push_back(Vec2f(i+0.75, j+0.25));
      _pos.push_back(Vec2f(i+0.25, j+0.75));
      _pos.push_back(Vec2f(i+0.75, j+0.75));
      _nbWallParticles += 4;
    }*/

    //Add a pipe. Beware to not superpose it with the fluid particles.
    /*for (int i = 1; i<res_x - 3; ++i) {
      for (int j : {10, 14, 20}) {
        _pos.push_back(Vec2f(i+0.25, j+0.25));
        _pos.push_back(Vec2f(i+0.75, j+0.25));
        _pos.push_back(Vec2f(i+0.25, j+0.75));
        _pos.push_back(Vec2f(i+0.75, j+0.75));
        _nbWallParticles += 4;
      }
    }
    for (int i = 3; i<res_x-1; ++i) {
      for (int j : {12, 17}) {
        _pos.push_back(Vec2f(i+0.25, j+0.25));
        _pos.push_back(Vec2f(i+0.75, j+0.25));
        _pos.push_back(Vec2f(i+0.25, j+0.75));
        _pos.push_back(Vec2f(i+0.75, j+0.75));
        _nbWallParticles += 4;
      }
    }*/

    //Add a U pipe. Beware to not superpose it with the fluid particles.
    /*for (int j = 5; j<res_y-1; ++j) {
      for (int i : {7, 11}) {
        _pos.push_back(Vec2f(i+0.25, j+0.25));
        _pos.push_back(Vec2f(i+0.75, j+0.25));
        _pos.push_back(Vec2f(i+0.25, j+0.75));
        _pos.push_back(Vec2f(i+0.75, j+0.75));
        _nbWallParticles += 4;
      }
    }

    for (int i = 8; i<11; ++i) {
      int j = 5;
      _pos.push_back(Vec2f(i+0.25, j+0.25));
        _pos.push_back(Vec2f(i+0.75, j+0.25));
        _pos.push_back(Vec2f(i+0.25, j+0.75));
        _pos.push_back(Vec2f(i+0.75, j+0.75));
        _nbWallParticles += 4;
    }*/


    // sample a fluid mass
    for(int j=0; j<f_height; ++j) {
      for(int i=0; i<f_width; ++i) {
        _pos.push_back(Vec2f(i+1.25, j+1.25));
        _pos.push_back(Vec2f(i+1.75, j+1.25));
        _pos.push_back(Vec2f(i+1.25, j+1.75));
        _pos.push_back(Vec2f(i+1.75, j+1.75));
      }
    }
    


    // make sure for the other particle quantities
    _vel = vector<Vec2f>(_pos.size(), Vec2f(0, 0));
    _acc = vector<Vec2f>(_pos.size(), Vec2f(0, 0));
    _p   = vector<Real>(_pos.size(), 100000);
    for (int i = 0; i < _nbWallParticles; ++i) {
      _p[i] = 0;
    }
    _d   = vector<Real>(_pos.size(), 0);
    _interVel = vector<Vec2f>(_pos.size(), Vec2f(0, 0));
    _interD = vector<Real>(_pos.size(), 0);
    _a_ii = vector<Real>(_pos.size(), 0);
    _d_ii = vector<Vec2f>(_pos.size(), Vec2f(0, 0));
    _c_i = vector<Vec2f>(_pos.size(), Vec2f(0, 0));
    _predP = vector<Real>(_pos.size(), 0);

    _col = vector<float>(_pos.size()*4, 1.0); // RGBA
    _vln = vector<float>(_pos.size()*4, 0.0); // GL_LINES

    for (int i = 0; i < res_x * res_y; i++) {
      _pidxInGrid.push_back(vector<tIndex>{});
    }
    
    buildNeighbor();
    computeDensity();
    initColor();
  }

  void update()
  {
    cout << '.' << flush;

    buildNeighbor();
    computeDensity();

    _acc = vector<Vec2f>(_pos.size(), Vec2f(0, 0));
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

    updateColor();
  }

  const CubicSpline &getKernel() const {return _kernel;}
  tIndex particleCount() const { return _pos.size(); }
  tIndex wallParticleCount() { return _nbWallParticles; }
  tIndex fluidParticleCount() { return _pos.size() - _nbWallParticles; }
  const vector<tIndex> &gridParticles(const tIndex i, const tIndex j) {return _pidxInGrid[idx1d(i, j)];};
  const Vec2f &position(const tIndex i) const { return _pos[i]; }
  const Vec2f &velocity(const tIndex i) const { return _vel[i]; }
  const Vec2f &acceleration(const tIndex i) const { return _acc[i]; }
  Real pressure(const tIndex i) const { return _p[i]; }
  Real density(const tIndex i) const { return _d[i]; }
  const float &color(const tIndex i) const { return _col[i]; }
  const float &vline(const tIndex i) const { return _vln[i]; }

  int resX() const { return _resX; }
  int resY() const { return _resY; }

private:
  void buildNeighbor() {

    for (auto &pix : _pidxInGrid) {
        pix.clear();
      }

    for(tIndex i=0; i<particleCount(); ++i) {
        _pidxInGrid[idx1d(floor(position(i).x), floor(position(i).y))].push_back(i);
      }
  }

  void computeDensity()
  {
    int nb;

    //#pragma omp parallel for 
    for (tIndex i = 0; i < particleCount(); ++i) {
      Real rho = 0.f;
      nb = 0;
      Real test = 0.f;

      for (tIndex _kernelX = max(0, (int) floor(position(i).x - _rad)); _kernelX < min(resX(), (int) floor(position(i).x + _rad + 1)); ++_kernelX) {
        for (tIndex _kernelY = max(0, (int) floor(position(i).y - _rad)); _kernelY < min(resY(), (int) floor(position(i).y + _rad + 1)); ++_kernelY) {
          
          for (tIndex j : _pidxInGrid[idx1d(_kernelX, _kernelY)]) {
            nb++;
            test += _kernel.w(position(i) - position(j));
            rho += _m0 * _kernel.w(position(i) - position(j));
          }
        }
      }
      //cout << rho << "   " << nb << "    " << test << endl;
      _d[i] = rho;
    }
    //cout << " " << nb << "    " << _d[particleCount() - 1] << "   ";
  }

  void applyBodyForce()
  {
    //#pragma omp parallel for 
    for(tIndex i=_nbWallParticles; i<particleCount(); ++i) {
      _acc[i] += _g;
    }
  }

  void applyViscousForce()
  {

    //#pragma omp parallel for 
    for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
      Vec2f dist = Vec2f();
      Vec2f viscousAcc = Vec2f();
      for (tIndex _kernelX = max(0, (int) floor(position(i).x - _rad)); _kernelX < min(resX(), (int) floor(position(i).x + _rad + 1)); ++_kernelX) {
        for (tIndex _kernelY = max(0, (int) floor(position(i).y - _rad)); _kernelY < min(resY(), (int) floor(position(i).y + _rad + 1)); ++_kernelY) {
          for (tIndex j : _pidxInGrid[idx1d(_kernelX, _kernelY)]) {
            dist = position(i) - position(j);
            viscousAcc += _m0 / _d[j] * (velocity(i) - velocity(j)) * (dist.dotProduct(_kernel.grad_w(dist))) / (dist.lengthSquare() + _h*_h);
          }
        }
      }
      _acc[i] += 2 * _nu * viscousAcc;
    }
    //cout << 2 * _nu * viscousAcc << "    ";
    //cout << "D " << _d[100] << " A " << _acc[100] << " V " << velocity(100) <<" P " << position(100) << endl;
  }

  void computeIntermediateVelocity()
  {
    //#pragma omp parallel for 
    for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
      _interVel[i] = _vel[i] + _dt * _acc[i];
    }
  }

  void computeDiiCoeff()
  {
    //#pragma omp parallel for 
    for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
      _d_ii[i] = Vec2f(0, 0);

      for (tIndex _kernelX = max(0, (int) floor(position(i).x - _rad)); _kernelX < min(resX(), (int) floor(position(i).x + _rad + 1)); ++_kernelX) {
        for (tIndex _kernelY = max(0, (int) floor(position(i).y - _rad)); _kernelY < min(resY(), (int) floor(position(i).y + _rad + 1)); ++_kernelY) {
          for (tIndex j : _pidxInGrid[idx1d(_kernelX, _kernelY)]) {

            
            _d_ii[i] += _m0 * _kernel.grad_w(position(i) - position(j));
          }
        }
      }

      _d_ii[i] *= -_dt*_dt /(_d[i]*_d[i]);
      
    }
  }

  void computeAiiCoeff()
  {
    //#pragma omp parallel for 
    for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
      _a_ii[i] = 0;

      for (tIndex _kernelX = max(0, (int) floor(position(i).x - _rad)); _kernelX < min(resX(), (int) floor(position(i).x + _rad + 1)); ++_kernelX) {
        for (tIndex _kernelY = max(0, (int) floor(position(i).y - _rad)); _kernelY < min(resY(), (int) floor(position(i).y + _rad + 1)); ++_kernelY) {
          for (tIndex j : _pidxInGrid[idx1d(_kernelX, _kernelY)]) {

            
            _a_ii[i] += _m0 * (_d_ii[i] + _dt * _dt * _m0 / (_d[i] * _d[i]) * _kernel.grad_w(position(j) - position(i))).dotProduct(_kernel.grad_w(position(i) - position(j)));
          }
        }
      }

      
    }
  }

  void computeIntermediateDensity()
  {
    //#pragma omp parallel for 
    for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
      _interD[i] = _d[i];

      for (tIndex _kernelX = max(0, (int) floor(position(i).x - _rad)); _kernelX < min(resX(), (int) floor(position(i).x + _rad + 1)); ++_kernelX) {
        for (tIndex _kernelY = max(0, (int) floor(position(i).y - _rad)); _kernelY < min(resY(), (int) floor(position(i).y + _rad + 1)); ++_kernelY) {
          for (tIndex j : _pidxInGrid[idx1d(_kernelX, _kernelY)]) {
            
            _interD[i] += _dt * _m0 * (_interVel[i] - _interVel[j]).dotProduct(_kernel.grad_w(position(i) - position(j)));
          }
        }
      }
      
    }
  }

  void computePressure()
  {
    //#pragma omp parallel for 
    for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
      _p[i] *= _initP;
      
    }

    int iter = 0;
    int convCriteria = true;

    while ((convCriteria) | (iter < 2)) {

      convCriteria = false;

      //#pragma omp parallel for
      for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
        _c_i[i] = Vec2f(0, 0);

        for (tIndex _kernelX = max(0, (int) floor(position(i).x - _rad)); _kernelX < min(resX(), (int) floor(position(i).x + _rad + 1)); ++_kernelX) {
          for (tIndex _kernelY = max(0, (int) floor(position(i).y - _rad)); _kernelY < min(resY(), (int) floor(position(i).y + _rad + 1)); ++_kernelY) {
            for (tIndex j : _pidxInGrid[idx1d(_kernelX, _kernelY)]) {
              
              _c_i[i] += _p[j] / (_d[j] * _d[j]) * _kernel.grad_w(position(i) - position(j));
            }
          }
        }

        _c_i[i] *= -_dt * _dt * _m0;
      }

      //#pragma omp parallel for
      for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
        Real newP;
        if (_a_ii[i] == 0.f) {
          newP = 0.f;
        } else {
          _predP[i] = _d0 - _interD[i];

          for (tIndex _kernelX = max(0, (int) floor(position(i).x - _rad)); _kernelX < min(resX(), (int) floor(position(i).x + _rad + 1)); ++_kernelX) {
            for (tIndex _kernelY = max(0, (int) floor(position(i).y - _rad)); _kernelY < min(resY(), (int) floor(position(i).y + _rad + 1)); ++_kernelY) {
              for (tIndex j : _pidxInGrid[idx1d(_kernelX, _kernelY)]) {
                if (j < _nbWallParticles) { //Boundary
                  _predP[i] -= _m0 * _c_i[i].dotProduct(_kernel.grad_w(position(i) - position(j)));
                } else {
                  _predP[i] -= _m0 * (_c_i[i] - _d_ii[j] * _p[j] - _c_i[j] - _dt * _dt * _m0 / (_d[i] * _d[i]) * _kernel.grad_w(position(j) - position(i)) * _p[i]).dotProduct(_kernel.grad_w(position(i) - position(j)));
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

      iter++;
      cout << "InterVel " << _interVel[obsPart] << " d_ii " << _d_ii[obsPart] << " a_ii " << _a_ii[obsPart] << " InterDen " << _interD[obsPart] << " c_i " << _c_i[obsPart] << " PredP " << _predP[obsPart] <<" Pr " << _p[obsPart] << endl;
    }

    cout << "Iterations : " << iter << endl;

  }


  void applyPressureForce()
  {
    //#pragma omp parallel for 
    for (tIndex i = _nbWallParticles; i < particleCount(); ++i) {
      Vec2f f = Vec2f(0);
      for (tIndex _kernelX = max(0, (int) floor(position(i).x - _rad)); _kernelX < min(resX(), (int) floor(position(i).x + _rad + 1)); ++_kernelX) {
        for (tIndex _kernelY = max(0, (int) floor(position(i).y - _rad)); _kernelY < min(resY(), (int) floor(position(i).y + _rad + 1)); ++_kernelY) {
          for (tIndex j : _pidxInGrid[idx1d(_kernelX, _kernelY)]) {
            
            f -= _m0 * (_p[i]/(_d[i] * _d[i]) + _p[j]/(_d[j] * _d[j])) * _kernel.grad_w(position(i) - position(j));
          }
        }
      }
      _acc[i] += f;
    }
    
    cout << "Den " << _d[obsPart] << " Pr " << _p[obsPart] << " Acc " << _acc[obsPart] << " Vit " << velocity(obsPart) <<" Pos " << position(obsPart) << endl;
    
  }

  void updateVelocity()
  {
    //#pragma omp parallel for 
    for(tIndex i=_nbWallParticles; i<particleCount(); ++i) {
      _vel[i] += _dt * _acc[i];
    }
  }

  void updatePosition()
  {
    //#pragma omp parallel for 
    for(tIndex i=_nbWallParticles; i<particleCount(); ++i) {
      _pos[i] += _dt * _vel[i];
    }
  }

  // simple collision detection/resolution for each particle
  void resolveCollision()
  {
    vector<tIndex> need_res;
    //#pragma omp parallel for 
    for(tIndex i=_nbWallParticles; i<particleCount(); ++i) {
      if(_pos[i].x<_l || _pos[i].y<_b || _pos[i].x>_r || _pos[i].y>_t)
        need_res.push_back(i);
    }

    //#pragma omp parallel for
    for(
      vector<tIndex>::const_iterator it=need_res.begin();
      it<need_res.end();
      ++it) {
      const Vec2f p0 = _pos[*it];
      _pos[*it].x = clamp(_pos[*it].x, _l, _r);
      _pos[*it].y = clamp(_pos[*it].y, _b, _t);
      _vel[*it] = (_pos[*it] - p0)/_dt;
    }
  }

  void initColor()
  {
    //#pragma omp parallel for 
    for(tIndex i=0; i<_nbWallParticles; ++i) {
      _col[i*4+0] = 0.8;
      _col[i*4+1] = 0.3;
      _col[i*4+2] = _d[i]/_d0;
    }

    //#pragma omp parallel for 
    for(tIndex i=_nbWallParticles; i<particleCount(); ++i) {
      _col[i*4+0] = 0.6;
      _col[i*4+1] = 0.6;
      _col[i*4+2] = _d[i]/_d0;
    }

    // _col[obsPart*4] = 1.f;
    // _col[obsPart*4+1] = 0.f;
  }

  void updateColor()
  {
    //#pragma omp parallel for 
    for(tIndex i=0; i<_nbWallParticles; ++i) {
      _col[i*4+2] = _d[i]/_d0;
    }

    //#pragma omp parallel for 
    for(tIndex i=_nbWallParticles; i<particleCount(); ++i) {
      _col[i*4+2] = _d[i]/_d0;
    }
  }

  void updateVelLine()
  {
    //#pragma omp parallel for 
    for(tIndex i=_nbWallParticles; i<particleCount(); ++i) {
      _vln[i*4+0] = _pos[i].x;
      _vln[i*4+1] = _pos[i].y;
      _vln[i*4+2] = _pos[i].x + _vel[i].x;
      _vln[i*4+3] = _pos[i].y + _vel[i].y;
    }
  }

  inline tIndex idx1d(const int i, const int j) { return i + j*resX(); }

  const CubicSpline _kernel;
  const Real _rad;

  // particle data
  vector<Vec2f> _pos;      // position
  vector<Vec2f> _vel;      // velocity
  vector<Vec2f> _acc;      // acceleration
  vector<Real>  _p;        // pressure
  vector<Real>  _d;        // density
  //IISPH specific data
  vector<Vec2f> _interVel; // intermediate velocity (without pressure force)
  vector<Real> _interD;    // intermediate density (with _predVel)
  vector<Real> _a_ii;      // diagonal coef of the SOE
  vector<Vec2f> _d_ii;     // first level coef of the SOE
  vector<Vec2f> _c_i;      // second level coef of the SOE
  vector<Real> _predP;     // predicted pressure
  const Real _initP;                  // initial pressure for Jacobi method. Must be between 0 and 1
  const Real _omega;                  // relaxation coefficient. Must be between 0 and 1
  const Real _pressureError;           // maximum pressure variation rate serving as a limit of the Jacobi method



  vector< vector<tIndex> > _pidxInGrid; // will help you find neighbor particles

  vector<float> _col;    // particle color; just for visualization
  vector<float> _vln;    // particle velocity lines; just for visualization

  // simulation
  Real _dt;                     // time step

  int _resX, _resY;             // background grid resolution

  // wall
  Real _l, _r, _b, _t;          // wall (boundary)
  tIndex _nbWallParticles;          // number of particles that belong to the wall

  // SPH coefficients
  Real _nu;                     // viscosity coefficient
  Real _d0;                     // rest density
  Real _h;                      // particle spacing (i.e., diameter)
  Vec2f _g;                     // gravity

  Real _m0;                     // rest mass
  Real _k;                      // EOS coefficient

  Real _eta;
  Real _c;                      // speed of sound
  Real _gamma;                  // EOS power factor
};