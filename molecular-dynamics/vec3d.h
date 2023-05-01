
#ifndef __VECTORS_H__
#define __VECTORS_H__

#include <cmath>

struct vec3d {
  double x, y, z;

  vec3d() {};

  vec3d(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;
  }

  vec3d(const vec3d &other) {
    x = other.x;
    y = other.y;
    z = other.z;
  }

  vec3d operator-() const { return vec3d(-x, -y, -z); }

  vec3d& operator=(const vec3d &other) {
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
  }

  vec3d& operator+=(const vec3d &other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }

  vec3d& operator-=(const vec3d &other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
  }

  vec3d& operator*=(double a) {
    x *= a;
    y *= a;
    z *= a;
    return *this;
  }

  vec3d& operator/=(double a) {
    x /= a;
    y /= a;
    z /= a;
    return *this;
  }

  double norm() { return sqrt(x*x + y*y + z*z); }

};

vec3d operator+(const vec3d &a, const vec3d &b) {return vec3d(a) += b;}
vec3d operator-(const vec3d &a, const vec3d &b) {return vec3d(a) -= b;}
vec3d operator*(const vec3d &a, double b) {return vec3d(a) *= b;}
vec3d operator*(double a, const vec3d &b) {return vec3d(b) *= a;}
vec3d operator/(const vec3d &a, double b) {return vec3d(a) /= b;}

#endif