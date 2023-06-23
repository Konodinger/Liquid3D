#ifndef __KERNEL_HPP__
#define __KERNEL_HPP__

#include <cmath>
#include <math.h>

#ifndef M_PI
#define M_PI 3.141592
#endif

#include "Vector.hpp"



//Linear kernel for diffuse particule generation.
class LinearKernel {
public:
    explicit LinearKernel(const Real h = 0.5f) : _h(2.f*h) {};
    Real supportRadius() const { return _h; }
    Real w(const Vec3f &rij) const { return max(0.f, 1 - (rij.length())/_h); }

private:
    Real _h;
};

// Cubic spline kernel for SPH.
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

#endif