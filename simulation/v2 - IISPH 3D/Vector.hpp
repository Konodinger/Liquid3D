// ----------------------------------------------------------------------------
// Vector.hpp
//
//  Created on: 22 Jan 2021
//      Author: Kiwon Um
//        Mail: kiwon.um@telecom-paris.fr
//
// Description: Vector class
//
// Copyright 2021-2023 Kiwon Um
//
// The copyright to the computer program(s) herein is the property of Kiwon Um,
// Telecom Paris, France. The program(s) may be used and/or copied only with
// the written permission of Kiwon Um or in accordance with the terms and
// conditions stipulated in the agreement/contract under which the program(s)
// have been supplied.
// ----------------------------------------------------------------------------

#ifndef _VECTOR_HPP_
#define _VECTOR_HPP_

#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

typedef float Real;
typedef long int tIndex;

inline Real square(const Real a) { return a*a; }
inline Real cube(const Real a) { return a*a*a; }
inline Real clamp(const Real v, const Real vmin, const Real vmax)
{
  if(v<vmin) return vmin;
  if(v>vmax) return vmax;
  return v;
}

template<typename T>
class Vector3 {
public:
  enum { D = 3 };

  typedef T ValueT;

  union {
    struct { T x; T y; T z; };
    struct { T i; T j; T k; };
    T v[D];
  };

  explicit Vector3(const T &value=0) : x(value), y(value), z(value) {}
  Vector3(const T &a, const T &b, const T &c) : x(a), y(b), z(c) {}

  // assignment operators
  Vector3& operator+=(const Vector3 &r) { x+=r.x; y+=r.y; z+=r.z; return *this; }
  Vector3& operator-=(const Vector3 &r) { x-=r.x; y-=r.y; z-=r.z; return *this; }
  Vector3& operator*=(const Vector3 &r) { x*=r.x; y*=r.y; z*=r.z; return *this; }
  Vector3& operator/=(const Vector3 &r) { x/=r.x; y/=r.y; z/=r.z; return *this; }

  Vector3& operator+=(const T *r) { x+=r[0]; y+=r[1]; z+=r[3]; return *this; }
  Vector3& operator-=(const T *r) { x-=r[0]; y-=r[1]; z-=r[3]; return *this; }
  Vector3& operator*=(const T *r) { x*=r[0]; y*=r[1]; z*=r[3]; return *this; }
  Vector3& operator/=(const T *r) { x/=r[0]; y/=r[1]; z/=r[3]; return *this; }

  Vector3& operator+=(const T s) { x+=s; y+=s; z+=s; return *this; }
  Vector3& operator-=(const T s) { x-=s; y-=s; z-=s; return *this; }
  Vector3& operator*=(const T s) { x*=s; y*=s; z*=s; return *this; }
  Vector3& operator/=(const T s) {
    const T d=static_cast<T>(1)/s; return operator*=(d);
  }

  // unary operators
  Vector3 operator+() const { return *this; }
  Vector3 operator-() const { return Vector3(-x, -y, -z); }

  // binary operators
  Vector3 operator+(const Vector3 &r) const { return Vector3(*this)+=r; }
  Vector3 operator-(const Vector3 &r) const { return Vector3(*this)-=r; }
  Vector3 operator*(const Vector3 &r) const { return Vector3(*this)*=r; }
  Vector3 operator/(const Vector3 &r) const { return Vector3(*this)/=r; }

  Vector3 operator+(const T *r) const { return Vector3(*this)+=r; }
  Vector3 operator-(const T *r) const { return Vector3(*this)-=r; }
  Vector3 operator*(const T *r) const { return Vector3(*this)*=r; }
  Vector3 operator/(const T *r) const { return Vector3(*this)/=r; }

  Vector3 operator+(const T s) const { return Vector3(*this)+=s; }
  Vector3 operator-(const T s) const { return Vector3(*this)-=s; }
  Vector3 operator*(const T s) const { return Vector3(*this)*=s; }
  Vector3 operator/(const T s) const { return Vector3(*this)/=s; }

  // comparison operators
  bool operator==(const Vector3 &r) const { return ((x==r.x) && (y==r.y) && (z==r.z)); }
  bool operator!=(const Vector3 &r) const { return !((*this)==r); }
  bool operator<(const Vector3 &r) const { return (x!=r.x) ? x<r.x : ((y!=r.y) ? y<r.y : z<r.z); }

  // cast operator
  template<typename T2>
  operator Vector3<T2>() const
  {
    return Vector3<T2>(static_cast<T2>(x), static_cast<T2>(y), static_cast<T2>(z));
  }

  const T& operator[](const tIndex i) const { assert(i<D); return v[i]; }
  T& operator[](const tIndex i)
  {
    return const_cast<T &>(static_cast<const Vector3 &>(*this)[i]);
  }

  // special calculative functions

  Vector3& lowerValues(const Vector3 &r)
  {
    x = min(x, r.x); y = min(y, r.y); z = min(z, r.z); return *this;
  }
  Vector3& upperValues(const Vector3 &r)
  {
    x = max(x, r.x); y = max(y, r.y); z = max(z, r.z); return *this;
  }

  template<typename T2>
  Vector3& increase(const tIndex di, const T2 d)
  {
    v[di] += static_cast<T>(d); return (*this);
  }
  template<typename T2>
  Vector3 increased(const tIndex di, const T2 d) const
  {
    return Vector3(*this).increase(di, d);
  }

  Vector3& normalize() { return (x==0 && y==0 && z==0) ? (*this) : (*this)/=length(); }
  Vector3 normalized() const { return Vector3(*this).normalize(); }

  T dotProduct(const Vector3 &r) const { return x*r.x + y*r.y + z*r.z; }

  T length() const { return sqrt(lengthSquare()); }
  T lengthSquare() const { return x*x + y*y + z*z; }
  T distanceTo(const Vector3 &t) const { return (*this-t).length(); }
  T distanceSquareTo(const Vector3 &t) const
  {
    return (*this-t).lengthSquare();
  }

  T mulAll() const { return x*y*z; }
  T sumAll() const { return x+y+z; }

  tIndex minorAxis() const { return (fabs(y)<fabs(x)) ? ((fabs(z)<fabs(y)) ? 2 : 1) : ((fabs(z)<fabs(x)) ? 2 : 0); }
  tIndex majorAxis() const { return (fabs(y)>fabs(x)) ? ((fabs(z)>fabs(y)) ? 2 : 1) : ((fabs(z)>fabs(x)) ? 2 : 0); }

  T minValue() const { return min(x, y, z); }
  T maxValue() const { return max(x, y, z); }
  T minAbsValue() const { return min(fabs(x), fabs(y), fabs(z)); }
  T maxAbsValue() const { return max(fabs(x), fabs(y), fabs(z)); }

  Vector3& rotate(const Real radian) { return *this = rotated(radian); }
  Vector3 rotated(const Real radianT, const Real radianO) const
  {
    const Real cost=cos(radianT);
    const Real sint=sin(radianT);
    const Real coso=cos(radianO);
    const Real sino=sin(radianO);
    return Vector3(cost*coso*x - sint*y + cost*sino*z,
                  sint*coso*x + coso*y - sint*sino*z,
                  -sino*x + coso*z);
  }

  Vector3& reflect(const Vector3 &n) { return *this = reflected(n); }
  Vector3 reflected(const Vector3 &n) const
  {
    return *this - 2.0*(this->dotProduct(n))*n;
  }

  Vector3& mirror(const Vector3 &n) { return *this = mirrored(n); }
  Vector3 mirrored(const Vector3 &n) const
  {
    return -(*this) + 2.0*(this->dotProduct(n))*n;
  }

  Vector3& project(const Vector3 &n) { return *this = projected(n); }
  Vector3 projected(const Vector3 &n) const
  {
    return n * (this->dotProduct(n));
  }

  Vector3& reject(const Vector3 &n) { return *this = rejected(n); }
  Vector3 rejected(const Vector3 &n) const { return *this - projected(n); }

  friend istream& operator>>(istream &in, Vector3 &vec)
  {
    return (in >> vec.x >> vec.y >> vec.z);
  }
  friend ostream& operator<<(ostream &out, const Vector3 &vec)
  {
    return (out << vec.x << " " << vec.y << " " << vec.z);
  }
};
typedef Vector3<Real> Vec3f;
inline const Vec3f operator*(const Real s, const Vec3f &r) { return r*s; }

#endif  /* _VECTOR_HPP_ */
