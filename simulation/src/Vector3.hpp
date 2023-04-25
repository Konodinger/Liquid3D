#ifndef _VECTOR3_HPP_
#define _VECTOR3_HPP_

typedef float Real;
typedef long int tIndex;

#include <iostream>

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

  Vector3& operator+=(const T *r) { x+=r[0]; y+=r[1]; z+=r[2]; return *this; }
  Vector3& operator-=(const T *r) { x-=r[0]; y-=r[1]; z-=r[2]; return *this; }
  Vector3& operator*=(const T *r) { x*=r[0]; y*=r[1]; z*=r[2]; return *this; }
  Vector3& operator/=(const T *r) { x/=r[0]; y/=r[1]; z/=r[2]; return *this; }

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

  Vector3 operator+(const T *r) {
    z = clamp(z, vmin, vmax);
    return *this;
  }

// dot product
T dot(const Vector3 &r) const
{
return xr.x + yr.y + z*r.z;
}

// cross product
Vector3 cross(const Vector3 &r) const
{
return Vector3(yr.z - zr.y, zr.x - xr.z, xr.y - yr.x);
}

// magnitude
T magnitude() const
{
return std::sqrt(square(x) + square(y) + square(z));
}

// squared magnitude
T squaredMagnitude() const
{
return square(x) + square(y) + square(z);
}

// normalize
Vector3 normalize()
{
T mag = magnitude();
if (mag != 0)
{
*this /= mag;
}
return *this;
}

// limit magnitude
Vector3 limit(const T &max)
{
T magSq = squaredMagnitude();
if (magSq > square(max))
{
*this *= max / std::sqrt(magSq);
}
return *this;
}

// clamp values
Vector3 clamp(const T &vmin, const T &vmax)
{
x = ::clamp(x, vmin, vmax);
y = ::clamp(y, vmin, vmax);
z = ::clamp(z, vmin, vmax);
return *this;
}

// binary operators
template<typename T>
Vector3<T> operator+(const Vector3<T> &l, const Vector3<T> &r)
{
return Vector3<T>(l) += r;
}

template<typename T>
Vector3<T> operator-(const Vector3<T> &l, const Vector3<T> &r)
{
return Vector3<T>(l) -= r;
}

template<typename T>
Vector3<T> operator*(const Vector3<T> &l, const Vector3<T> &r)
{
return Vector3<T>(l) *= r;
}

template<typename T>
Vector3<T> operator/(const Vector3<T> &l, const Vector3<T> &r)
{
return Vector3<T>(l) /= r;
}

template<typename T>
Vector3<T> operator+(const Vector3<T> &v, const T *r)
{
return Vector3<T>(v) += r;
}

template<typename T>
Vector3<T> operator-(const Vector3<T> &v, const T *r)
{
return Vector3<T>(v) -= r;
}

template<typename T>
Vector3<T> operator*(const Vector3<T> &v, const T *r)
{
return Vector3<T>(v) *= r;
}

template<typename T>
Vector3<T> operator/(const Vector3<T> &v, const T *r)
{
return Vector3<T>(v) /= r;
}

template<typename T>
Vector3<T> operator+(const T *l, const Vector3<T> &v)
{
return Vector3<T>(v) += l;
}

template<typename T>
Vector3<T> operator-(const T *l, const Vector3<T> &v)
{
return Vector3<T>(v) -= l;
}

template<typename T>
Vector3<T> operator*(const T *l, const Vector3<T> &v)
{
return Vector3<T>(v) *= l;
}

template<typename T>
Vector3<T> operator/(const T *l, const Vector3<T> &v)
{
return Vector3<T>(v) /=l;
}

// scalar operations
template<typename T>
Vector3<T> operator+(const Vector3<T> &v, const T &s)
{
return Vector3<T>(v) += s;
}

template<typename T>
Vector3<T> operator-(const Vector3<T> &v, const T &s)
{
return Vector3<T>(v) -= s;
}

template<typename T>
Vector3<T> operator*(const Vector3<T> &v, const T &s)
{
return Vector3<T>(v) *= s;
}

template<typename T>
Vector3<T> operator/(const Vector3<T> &v, const T &s)
{
return Vector3<T>(v) /= s;
}

template<typename T>
Vector3<T> operator+(const T &s, const Vector3<T> &v)
{
return Vector3<T>(v) += s;
}

template<typename T>
Vector3<T> operator-(const T &s, const Vector3<T> &v)
{
return Vector3<T>(s - v.x, s - v.y, s - v.z);
}

template<typename T>
Vector3<T> operator*(const T &s, const Vector3<T> &v)
{
return Vector3<T>(v) *= s;
}

template<typename T>
Vector3<T> operator/(const T &s, const Vector3<T> &v)
{
return Vector3<T>(s / v.x, s / v.y, s / v.z);
}

// comparison operators
template<typename T>
bool operator==(const Vector3<T> &l, const Vector3<T> &r)
{
return l.x == r.x && l.y == r.y && l.z == r.z;
}

template<typename T>
bool operator!=(const Vector3<T> &l, const Vector3<T> &r)
{
return !(l == r);
}

template<typename T>
bool operator<(const Vector3<T> &l, const Vector3<T> &r)
{
return l.squaredMagnitude() < r.squaredMagnitude();
}

template<typename T>
bool operator<=(const Vector3<T> &l, const Vector3<T> &r)
{
return l.squaredMagnitude() <= r.squaredMagnitude();
}

template<typename T>
bool operator>(const Vector3<T> &l, const Vector3<T> &r)
{
return l.squaredMagnitude() > r.squaredMagnitude();
}

template<typename T>
bool operator>=(const Vector3<T> &l, const Vector3<T> &r)
{
return l.squaredMagnitude() >= r.squaredMagnitude();
}

// output stream operator
template<typename T>
std::ostream& operator<<(std::ostream &os, const Vector3<T> &v) {
os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
return os;
}
};

#endif // MYLIB_MATH_VECTOR3_H
