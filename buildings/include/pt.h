#pragma once

#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <vector>

using namespace std;

static inline double frandom(void) {
  return (double)(random()) / INT_MAX;
}

static inline double frandom(double _min, double _max) {
  return _min + (_max - _min) * frandom();
}

class Pt {
  public:

    double x;
    double y;
    double z;

    Pt() {
      set(0,0,0);
    }

    Pt(double v) {
      set(v, v, v);
    }

    Pt(double x, double y, double z=0.) {
      set(x, y, z);
    }

    void set(double x, double y, double z) {
      this->x = x;
      this->y = y;
      this->z = z;
    }

    void load(const double coords[3]) {
      x = coords[0];
      y = coords[1];
      z = coords[2];
    }

    static Pt randoml(double min_len=0, double max_len=1)
    {
      Pt r = random();
      return r.unit(Pt(1, 0, 0)) * (min_len + (max_len - min_len) * frandom());
    }

    static Pt random(double min_c=-1, double max_c=1)
    {
      return Pt(min_c + (max_c - min_c) * frandom(),
                min_c + (max_c - min_c) * frandom(),
                min_c + (max_c - min_c) * frandom());
    }

    bool zero() const {
      return (!x) && (!y) && (!z);
    }

    double len() const {
      return sqrt(x*x + y*y + z*z);
    }

    double volume() const {
      return fabs(x * y * z);
    }

    Pt unit(const Pt &if_zero=Pt()) const {
      double l = len();
      if (l < 1.e-6)
        return if_zero;
      return (*this) / l;
    }

		double cart_len() const {
			return max(fabs(x), max(fabs(y), fabs(z)));
		}

    void rot_about(const Pt &axis, double rad) {
			*this = rotated_about(axis, rad);
		}

    Pt rotated_about(const Pt &axis, double rad) const {
#if 1
      Pt u = axis.unit();
      double cos_rad = cos(rad);
      double sin_rad = sin(rad);
			return
        u * (dot(u) * (1. - cos_rad))
        + (*this)*cos_rad
        + u.cross(*this) * sin_rad;
#else
      double xx, yy, zz;
      double l = axis.len();
      double u = axis.x / l;
      double v = axis.y / l;
      double w = axis.z / l;
      double cos_rad = cos(rad);
      double sin_rad = sin(rad);
      double f1 = (u*x + v*y + w*z)*(-cos_rad + 1);
      xx = u*f1 + x*cos_rad + (- w*y + v*z)*sin_rad;
      yy = v*f1 + y*cos_rad + (  w*x - u*z)*sin_rad;
      zz = w*f1 + z*cos_rad + (- v*x + u*y)*sin_rad;
			return Pt(xx, yy, zz);
#endif
    }

    /* Euler vector: rotate about given vector by the length of the vector in
     * radians. */
    void rot_e(const Pt &e) {
      double rad = e.len();
      if (fabs(rad) < 1e-6)
        return;
      Pt u = e / rad;
      double cos_rad = cos(rad);
      double sin_rad = sin(rad);
      *this =
        u * (dot(u) * (1. - cos_rad))
        + (*this)*cos_rad
        + u.cross(*this) * sin_rad;
    }


    void scale3(Pt &p) {
      x *= p.x;
      y *= p.y;
      z *= p.z;
    }

    Pt &set_min(double c) {
      x = min(x, c);
      y = min(y, c);
      z = min(z, c);
      return *this;
    }

    Pt &set_max(double c) {
      x = max(x, c);
      y = max(y, c);
      z = max(z, c);
      return *this;
    }

    void set_min(const Pt &p) {
      x = min(x, p.x);
      y = min(y, p.y);
      z = min(z, p.z);
    }

    void set_max(const Pt &p) {
      x = max(x, p.x);
      y = max(y, p.y);
      z = max(z, p.z);
    }

    Pt &limit(double c_min, double c_max) {
      x = min(c_max, max(c_min, x));
      y = min(c_max, max(c_min, y));
      z = min(c_max, max(c_min, z));
      return *this;
    }

    static bool wrap_coord(double &x, double c, double dist) {
      bool changed = false;
      double d = x - c;
      if (fabs(d) > dist) {
        d -= dist * 2. * (d>0? 1. : -1.) * ((int)(fabs(d)/(dist*2)) + 1);
        changed = true;
      }
      x = d;
      return changed;
    }

    int wrap_cube(const Pt &center, double dist) {
      return
        (wrap_coord(x, center.x, dist)
         || wrap_coord(y, center.y, dist)
         || wrap_coord(z, center.z, dist)) ? 1 : 0;
    }

    int wrap_sphere(const Pt &center, double dist) {
      Pt d = (center - (*this));
      if (d.len() > dist) {
        *this += d.unit() * (dist*2);
        return 1;
      }
      return 0;
    }

    double dot(const Pt &p) const
    {
      return
        x * p.x + y * p.y + z * p.z;
    }

    double udot(const Pt &p) const
    {
      Pt u = unit();
      Pt pu = p.unit();
      return
        min(1.0, max(-1.0, u.x * pu.x + u.y * pu.y + u.z * pu.z));
    }

    Pt cross(const Pt &p) const
    {
      return Pt(y * p.z - z * p.y,
                z * p.x - x * p.z,
                x * p.y - y * p.x);
    }

    double min_angle(const Pt &p) const
    {
      return acos(this->udot(p));
    }

    double angle(const Pt &p, const Pt &n) const
    {
      double ma = acos(this->udot(p));
      if (n.udot( this->cross(p) ) < 0)
        return - ma;
      return ma;
    }

    Pt uv() const
    {
      Pt p = unit();
      return Pt(
                .5 + ( atan2(p.z, p.x) / (2.*M_PI) ),
                .5 - ( asin(p.y) / M_PI ),
                0);
    }

    Pt& operator+=(const Pt &p) {
      x += p.x;
      y += p.y;
      z += p.z;
      return *this;
    }

    Pt& operator-=(const Pt &p) {
      x -= p.x;
      y -= p.y;
      z -= p.z;
      return *this;
    }

    Pt& operator+=(const double f) {
      x += f;
      y += f;
      z += f;
      return *this;
    }

    Pt& operator-=(const double f) {
      x -= f;
      y -= f;
      z -= f;
      return *this;
    }

    Pt& operator*=(const double f) {
      x *= f;
      y *= f;
      z *= f;
      return *this;
    }

    Pt& operator/=(const double f) {
      x /= f;
      y /= f;
      z /= f;
      return *this;
    }

    Pt operator+(const Pt &p) const {
      Pt x(*this);
      x += p;
      return x;
    }

    Pt operator-(const Pt &p) const {
      Pt x(*this);
      x -= p;
      return x;
    }

    Pt operator*(const double f) const {
      Pt x(*this);
      x *= f;
      return x;
    }

    Pt operator/(const double f) const {
      Pt x(*this);
      x /= f;
      return x;
    }

    Pt operator+(const double f) const {
      Pt x(*this);
      x += f;
      return x;
    }

    Pt operator-(const double f) const {
      Pt x(*this);
      x -= f;
      return x;
    }

    bool operator==(const Pt& other) const {
      return (x == other.x) && (y == other.y) && (z == other.z);
    }

    bool operator!=(const Pt& other) const {
      return (x != other.x) || (y != other.y) || (z != other.z);
    }

    void print() const
    {
      printf("Pt(%5.3f, %5.3f, %5.3f)\n", x, y, z);
    }

    Pt &operator=(double v)
    {
      x = y = z = v;
      return *this;
    }

    Pt &set(const double v[3])
    {
      x = v[0];
      y = v[1];
      z = v[2];
      return *this;
    }

    Pt scaled(const Pt &factors) const
    {
      return Pt(x * factors.x, y * factors.y, z * factors.z);
    }

    Pt unscaled(const Pt &factors) const
    {
      return Pt(fabs(factors.x) < 1e-6? 0 : x / factors.x,
                fabs(factors.y) < 1e-6? 0 : y / factors.y,
                fabs(factors.z) < 1e-6? 0 : z / factors.z);
    }

    Pt without(const Pt &axis, bool neg_means_zero=false) const
    {
      return (*this) - project(axis, neg_means_zero);
    }

    Pt project(const Pt &axis, bool neg_means_zero=true) const
    {
      double d = dot(axis);
      if (neg_means_zero && (d < 0))
        return Pt();
      return axis * d;
    }
};


