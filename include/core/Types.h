#pragma once

#include <array>
#include <vector>
#include <cmath>

namespace bndf {

struct Vec3 {
  double x{0}, y{0}, z{0};
};

inline Vec3 operator+(const Vec3& a, const Vec3& b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
inline Vec3 operator*(double s, const Vec3& a) { return {s*a.x, s*a.y, s*a.z}; }
inline Vec3 operator*(const Vec3& a, double s) { return s*a; }
inline double dot(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline double norm(const Vec3& a) { return std::sqrt(dot(a,a)); }
inline Vec3 normalize(const Vec3& a) { double n = norm(a); return (n>0)? (1.0/n)*a : a; }

} // namespace bndf

