#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <cstdlib>
#include <limits.h>
#include <time.h>
#include <unistd.h>

#include <math.h>
#include <stdio.h>

#include <vector>
using namespace std;

#include "textures.h"
#include "foreach.h"
#include "pt.h"

struct Matrix33 {
    double a, b, c,
           d, e, f,
           g, h, i;

    static Matrix33 from_rot3(Pt rot3) {
      double cx = cos(rot3.x);
      double cy = cos(rot3.y);
      double cz = cos(rot3.z);
      double sx = sin(rot3.x);
      double sy = sin(rot3.y);
      double sz = sin(rot3.z);
      Matrix33 mx = {
        1, 0, 0,
        0, cx, sx,
        0, -sx, cx
      };

      Matrix33 my {
        cy, 0, -sy,
        0, 1, 0,
        sy, 0, cy
      };

      Matrix33 mz {
        cz, sz, 0,
        -sz, cz, 0,
        0, 0, 1
      };

      return my * mx * mz;
    }

    Matrix33 operator*(const Matrix33 &o) const
    {
      return Matrix33{
        a*o.a + b*o.d + c*o.g, a*o.b + b*o.e + c*o.h, a*o.c + b*o.f + c*o.i,
        d*o.a + e*o.d + f*o.g, d*o.b + e*o.e + f*o.h, d*o.c + e*o.f + f*o.i,
        g*o.a + h*o.d + i*o.g, g*o.b + h*o.e + i*o.h, g*o.c + h*o.f + i*o.i
      };
    }

    Pt operator*(const Pt &p) const
    {
      /*     a b c   x    a*x + b*y + c*z
             d e f * y =  d*x + e*y + f*z
             g h i   z    g*x + h*y + i*z */
      return Pt( a*p.x + b*p.y + c*p.z,
                 d*p.x + e*p.y + f*p.z,
                 g*p.x + h*p.y + h*p.z);
    }

    double det() const
    {
      return a*e*i + b*f*g + c*d*h - g*e*c - h*f*a - i*d*b;
    }

    Matrix33 inverse() const
    {
      double dd = det();
      if (fabs(dd) < 1e-6)
        return Matrix33();
      return Matrix33{
        e*i - f*h, c*h - b*i, b*f - c*e,
        f*g - d*i, a*i - c*g, c*d - a*f, 
        d*h - e*g, b*g - a*h, a*e - b*d
      } / dd;
    }

    Matrix33 operator/(double x) const
    {
      return Matrix33{
          a/x, b/x, c/x,
          d/x, e/x, f/x,
          g/x, h/x, i/x
      };
    }

};

class Plane {
  public:
    Pt pos;
    Pt n;

    Plane(const Pt &pos, const Pt &n)
    {
      this->pos = pos;
      this->n = n;
    }

    Plane(const Pt &p0, const Pt &p1, const Pt &p2)
    {
      pos = p0;
      n = (p1 - p0).cross(p2 - p0);
    }
};

class Line {
  public:
    Pt pos;
    Pt len;

    Line(const Pt &pos, const Pt &len)
    {
      this->pos = pos;
      this->len = len;
    }
};

bool intersect(const Plane &plane, const Line &line,
               double *at=NULL, Pt *intersection=NULL)
{
  double D = plane.n.dot(line.len);
  if (fabs(D) < 1e-6)
    return false;

  double N = plane.n.dot(plane.pos - line.pos);
  double _at = N / D;

  if (at)
    *at = _at;

  if ((_at < -1e-6) || (_at > (1. + 1e-6)))
    return false;

  if (intersection)
    *intersection = line.pos + line.len * _at;

  return true;
}

bool intersect(const Pt &p0, const Pt &p1, const Pt &p2, Line &line,
               double *at=NULL, Pt *intersection=NULL,
               bool parallelogram=false)
{
  double _at;
  Pt _intersection;

  Pt u = p1 - p0;
  Pt v = p2 - p0;

  if (! intersect(Plane(p0, u.cross(v)), line, &_at, &_intersection))
    return false;

  if (at)
    *at = _at;

  if (intersection)
    *intersection = _intersection; 

  Pt w = _intersection - p0;

  double udotv = u.dot(v);
  double wdotv = w.dot(v);
  double vdotv = v.dot(v);
  double wdotu = w.dot(u);
  double udotu = u.dot(u);

  double D = udotv * udotv - udotu * vdotv;

  double uu = ( udotv * wdotv - vdotv * wdotu ) / D;
  if ((uu < -1e-6) || (uu > (1. + 1e-6)))
    return false;

  double vv = ( udotv * wdotu - udotu * wdotv ) / D;
  if ((vv < -1e-6) || (vv > (1. + 1e-6)))
    return false;

  return parallelogram || ((vv + uu) < (1. + 1e-6));
}


class Color {
  public:
    double r;
    double g;
    double b;
    double a;

    Color() {
      set(1,1,1,1);
    }

    void set(double r, double g, double b, double a=1) {
      this->r = r;
      this->g = g;
      this->b = b;
      this->a = a;
    }

    void load(const double color[4]) {
      this->r = color[0];
      this->g = color[1];
      this->b = color[2];
      this->a = color[3];
    }

    void blend(Color &c, double fade) {
      double rfade = 1.0 - fade;
      r = rfade * r  +  fade * c.r;
      g = rfade * g  +  fade * c.g;
      b = rfade * b  +  fade * c.b;
      a = rfade * a  +  fade * c.a;
    }

    void random(double cmin=.3, double cmax=1.) {
      r = frandom();
      g = frandom();
      b = frandom();
      a = 1;
      double is_intens = max(.01, (r + g + b) / 3.);
      double intens_change = (cmin + (cmax-cmin)*frandom()) / is_intens;

      r = intens_change * r;
      g = intens_change * g;
      b = intens_change * b;
    }

    void glColor(double alpha=1., double greying=0) {
      alpha *= a;
      if (greying > .01)
        ::glColor4f(( r * (1.-greying) + greying * .33),
                    ( g * (1.-greying) + greying * .33),
                    ( b * (1.-greying) + greying * .33),
                    ( alpha * (1.-greying) + greying * .33)
                    );
      else
        ::glColor4f(r, g, b, alpha);
    }

    void glMaterial() {
      GLfloat mat_specular[] = { (float)r, (float)g, (float)b, (float)a };
      GLfloat mat_shininess[] = { .1 };
      glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
      glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    }

    void print() {
      printf("r%f g%f b%f a%f\n", r,g,b,a);
    }

    Color& operator=(const Pt &p) {
      r = p.x;
      g = p.y;
      b = p.z;
    }
};



class PtGl : public Pt {
  public:
    PtGl() : Pt() {}
    PtGl(const Pt &p) : Pt(p) {}
    PtGl(double v) : Pt(v) {}
    PtGl(double x, double y, double z=0.) : Pt(x,y,z) {}
    PtGl& operator=(const Pt &p) {
      Pt::operator=(p);
      return *this;
    }
    PtGl &operator=(double v)
    {
      Pt::operator=(v);
      return *this;
    }

    void glVertex3d() const
    {
      ::glVertex3d(x, y, z);
    }

    void glTranslated() const
    {
      ::glTranslated(x, y, z);
    }

    void glRotated() const
    {
      // last, turn left/right
      if (y) {
        ::glRotated(y * (180./M_PI), 0,1,0);
      }
      // second, turn nose up or down
      if (x) {
        ::glRotated(x * (180./M_PI), 1,0,0);
      }
      // first, roll left or right
      if (z) {
        ::glRotated(z * (180./M_PI), 0,0,1);
      }
    }

    void glNormal() const
    {
      ::glNormal3d(x, y, z);
    }


    void glScaled() const
    {
      ::glScaled(x, y, z);
    }

    void glTexCoord()
    {
      ::glTexCoord2d(x, y);
    }
};

class Point : public PtGl{
  public:

    Color c;
    PtGl t;

    Point() : PtGl() {
      c.random();
    }

    Point(double x, double y, double z) : PtGl(x, y, z) {
      c.random();
    }

    Point(const Pt &p)
    {
      *this = p;
      c.random();
    }

    void draw(double alpha=1., double greying=0) {
      c.glColor(alpha, greying);
      glVertex3d();
    }

    void random() {
      Pt::random();
      c.random();
    }

    Point &operator=(const Pt &p)
    {
      Pt::operator=(p);
      return *this;
    }

};



class FilterContext {
  public:
    int point_i;
    int color_i;
    int face_i;
    int particle_i;
    int cloud_i;
    int frame_i;

    GLenum what;

    FilterContext() {
      clear();
      frame_i = 0;
    }

    void clear() {
      point_i = 0;
      color_i = 0;
      face_i = 0;
      particle_i = 0;
      cloud_i = 0;
      what = -1;
    }
};

class Filter {
  public:
    FilterContext *fc;

    Filter() {
      fc = NULL;
    }

    virtual void ign() {};
};

typedef enum {
  level_unknown = 0,
  level_camera,
  level_cloud,
  level_particle
} level_e;


class PointFilter : public Filter {
  public:
    virtual void point(Pt &p) = 0;
};

class ColorFilter : public Filter {
  public:
    virtual void color(Color &c) = 0;
};

class TranslateFilter : public Filter {
  public:
    virtual void translate(Pt &p, level_e l) = 0;
};

class RotateFilter : public Filter {
  public:
    virtual void rotate(Pt &p, level_e l) = 0;
};

class ScaleFilter : public Filter {
  public:
    virtual void scale(Pt &p, level_e l) = 0;
};

template <typename V, typename T>
void vector_rm(V &v, T &item) {
  typename V::iterator i = v.begin();
  for (; i != v.end(); i++) {
    if (*i == item) {
      v.erase(i);
      return;
    }
  }
}

/* Laugh if you might, yet this gives a nice overview of drawing primitives. */
class Draw {
  public:
    static void point(PtGl &p) {
      p.glVertex3d();
    }

    static void color(Color &c, double opacity=1.) {
      c.glColor(opacity);
    }

    static void color_point(Point &p, double opacity=1.) {
      color(p.c, opacity);
      point(p);
    }

    static void texture_coord(double cx, double cy) {
      ::glTexCoord2d(cx, cy);
    }

    static void translate(PtGl &p) {
      p.glTranslated();
    }

    static void rotate(PtGl &p) {
      p.glRotated();
    }

    static void scale(PtGl &p) {
      p.glScaled();
    }

    static void begin(GLenum what) {
      glBegin(what);
    }

    static void end() {
      glEnd();
    }
};

class DrawBank {
  public:
    FilterContext fc;

    vector<PointFilter*> point_filters;
    vector<ColorFilter*> color_filters;
    vector<TranslateFilter*> translate_filters;
    vector<RotateFilter*> rotate_filters;
    vector<ScaleFilter*> scale_filters;

    void start() {
      fc.clear();
    }

    void add(Filter &f) {
      add(&f);
    }

    void add(Filter *f) {
      if (dynamic_cast<PointFilter*>(f)) {
        f->fc = &fc;
        point_filters.push_back((PointFilter*)f);
      }
      else if (dynamic_cast<ColorFilter*>(f)) {
        f->fc = &fc;
        color_filters.push_back((ColorFilter*)f);
      }
      else if (dynamic_cast<TranslateFilter*>(f)) {
        f->fc = &fc;
        translate_filters.push_back((TranslateFilter*)f);
      }
      else if (dynamic_cast<RotateFilter*>(f)) {
        f->fc = &fc;
        rotate_filters.push_back((RotateFilter*)f);
      }
      else if (dynamic_cast<ScaleFilter*>(f)) {
        f->fc = &fc;
        scale_filters.push_back((ScaleFilter*)f);
      }
    }

    void remove(Filter &f) {
      remove(&f);
    }

    void remove(Filter *f) {
      if (dynamic_cast<PointFilter*>(f)) {
        f->fc = &fc;
        PointFilter *ff = (PointFilter*)f;
        vector_rm<vector<PointFilter*>,PointFilter*>(point_filters, ff);
      }
      else if (dynamic_cast<ColorFilter*>(f)) {
        f->fc = &fc;
        ColorFilter *ff = (ColorFilter*)f;
        vector_rm<vector<ColorFilter*>,ColorFilter*>(color_filters, ff);
      }
      else if (dynamic_cast<TranslateFilter*>(f)) {
        f->fc = &fc;
        TranslateFilter *ff = (TranslateFilter*)f;
        vector_rm<vector<TranslateFilter*>,TranslateFilter*>(translate_filters, ff);
      }
      else if (dynamic_cast<RotateFilter*>(f)) {
        f->fc = &fc;
        RotateFilter *ff = (RotateFilter*)f;
        vector_rm<vector<RotateFilter*>,RotateFilter*>(rotate_filters, ff);
      }
      else if (dynamic_cast<ScaleFilter*>(f)) {
        f->fc = &fc;
        ScaleFilter *ff = (ScaleFilter*)f;
        vector_rm<vector<ScaleFilter*>,ScaleFilter*>(scale_filters, ff);
      }
    }

    void point(Point &p) {
      if (point_filters.size()) {
        Point q = p;
        for (int i = 0; i < point_filters.size(); i++)
          point_filters[i]->point(q);
        q.glVertex3d();
      }
      else {
        p.glVertex3d();
      }
      fc.point_i ++;
    }

    void color(Color &c) {
      if (color_filters.size()) {
        Color d = c;
        for (int i = 0; i < color_filters.size(); i++)
          color_filters[i]->color(d);
        d.glColor();
      }
      else
        c.glColor();
      fc.color_i ++;
    }

    void color_point(Point &p) {
      color(p.c);
      point(p);
    }

    void texture_coord(double cx, double cy) {
      ::glTexCoord2d(cx, cy);
    }

    void translate(Point &p, level_e l) {
      if (translate_filters.size()) {
        Point q = p;
        for (int i = 0; i < translate_filters.size(); i++)
          translate_filters[i]->translate(q, l);
        q.glTranslated();
      }
      else
        p.glTranslated();
    }

    void rotate(Point &p, level_e l) {
      if (rotate_filters.size()) {
        Point q = p;
        for (int i = 0; i < rotate_filters.size(); i++)
          rotate_filters[i]->rotate(q, l);
        q.glRotated();
      }
      else
        p.glRotated();
    }

    void scale(Point &p, level_e l) {
      if (scale_filters.size()) {
        Point q = p;
        for (int i = 0; i < scale_filters.size(); i++)
          scale_filters[i]->scale(q, l);
        q.glScaled();
      }
      else
        p.glScaled();
    }

    void begin(GLenum what) {
      fc.what = what;
      glBegin(what);
    }

    void end() {
      glEnd();
      fc.what = -1;
    }

    void end_of_face() {
      fc.face_i ++;
    }

    void end_of_particle() {
      fc.particle_i ++;
    }

    void end_of_cloud() {
      fc.face_i ++;
    }

    void end_of_frame() {
      fc.frame_i ++;
      fc.clear();
    }
};



class DrawAs {
  public:
    double alpha;
    
    DrawAs() {
      alpha = 1;
    }

    virtual void draw(vector<Point> &points, DrawBank &d) = 0;
};


class AsLines : public DrawAs {
  public:
    AsLines() : DrawAs() {}
    virtual void draw(vector<Point> &points, DrawBank &d)
    {
      int l;

      glLineWidth(2);
      glEnable( GL_LINE_SMOOTH );
      //glEnable( GL_POLYGON_SMOOTH );
      glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
      //glHint( GL_POLYGON_SMOOTH_HINT, GL_NICEST );

      d.begin(GL_LINES);

      l = points.size();
      for (int i = 1; i < l; i++) {
        d.color_point(points[i-1]);
        d.color_point(points[i]);
        d.end_of_face();
      }
      d.color_point(points.back());
      d.color_point(points.front());
      d.end_of_face();

      d.end();
    }
};

class AsPoints : public DrawAs {
  public:
    AsPoints() : DrawAs() {}
    virtual void draw(vector<Point> &points, DrawBank &d)
    {
      int l;

      d.begin(GL_POINTS);

      l = points.size();
      for (int i = 0; i < l; i++) {
        d.color_point(points[i]);
        d.end_of_face();
      }

      d.end();
    }
};

class AsTriangles : public DrawAs {
  public:
    virtual void draw(vector<Point> &points, DrawBank &d)
    {
      int l;
      d.begin(GL_TRIANGLES);

      l = points.size();
      for (int i = 2; i < l; i++) {
        d.color_point(points[i-2]);
        d.color_point(points[i-1]);
        d.color_point(points[i]);
        d.end_of_face();
      }

      d.end();
    }
};

class AsTets : public DrawAs {
  public:
    int pitch;

    AsTets() {
      pitch = 2;
    }

    virtual void draw(vector<Point> &points, DrawBank &d)
    {
      int l;
#if 1
      d.begin(GL_TRIANGLES);

      l = points.size();
      for (int i = 3; i < l; i += pitch) {
        d.color_point(points[i-3]);
        d.color_point(points[i-2]);
        d.color_point(points[i-1]);
        d.end_of_face();

        d.color_point(points[i-3]);
        d.color_point(points[i-2]);
        d.color_point(points[i]);
        d.end_of_face();

        d.color_point(points[i-3]);
        d.color_point(points[i-1]);
        d.color_point(points[i]);
        d.end_of_face();

        d.color_point(points[i-2]);
        d.color_point(points[i-1]);
        d.color_point(points[i]);
        d.end_of_face();
      }
#else

      glBegin(GL_TRIANGLE_STRIP);

      l = points.size();
      for (int i = 3; i < l; i += 2) {
        points[i-3].draw(alpha);
        points[i-2].draw(alpha);
        points[i-1].draw(alpha);

        points[i].draw(alpha);

        points[i-3].draw(alpha);

        points[i-2].draw(alpha);

      }
#endif

      glEnd();
    }
};


class AsPoly : public DrawAs {
  public:
    AsPoly() : DrawAs() {}
    virtual void draw(vector<Point> &points, DrawBank &d)
    {
      int l;
      d.begin(GL_POLYGON);

      l = points.size();
      for (int i = 0; i < l; i ++) {
        d.color_point(points[i]);
      }
      d.end_of_face();

      d.end();
    }
};

class AsQuads : public DrawAs {
  public:
    int pitch;

    AsQuads(){
      pitch = 4;
    }

    virtual void draw(vector<Point> &points, DrawBank &d)
    {
      int l;
      d.begin(GL_QUADS);

      l = points.size();
      for (int i = 3; i < l; i += pitch) {
        d.color_point(points[i-3]);
        d.color_point(points[i-2]);
        d.color_point(points[i-1]);
        d.color_point(points[i]);
        d.end_of_face();
      }

      glEnd();
    }
};

class AsTexturePlanes : public DrawAs {
  public:
    int pitch;
    Texture *texture;

    AsTexturePlanes() {
      pitch = 4;
      texture = NULL;
    }

    virtual void draw(vector<Point> &points, DrawBank &d)
    {
      int l;
      texture->begin();
      d.begin(GL_QUADS);

      l = points.size();
      for (int i = 3; i < l; i += pitch) {
        d.texture_coord(0, 0);
        d.color_point(points[i-3]);
        d.texture_coord(1, 0);
        d.color_point(points[i-2]);
        d.texture_coord(1, 1);
        d.color_point(points[i-1]);
        d.texture_coord(0, 1);
        d.color_point(points[i]);
        d.end_of_face();
      }

      glEnd();
      texture->end();
    }
};


class Placed {
  public:
    Point pos;
    Point rot3;
    Point scale;
    level_e level;

    Placed() {
      pos.set(0, 0, 0);
      rot3.set(0, 0, 0);
      scale.set(1, 1, 1);
      level = level_unknown;
    }

    void placement(DrawBank &d) {
      d.translate(pos, level);
      d.rotate(rot3, level);
      d.scale(scale, level);
    }

    void draw(DrawAs &as, DrawBank &d) { 
      glPushMatrix();
      placement(d);

      _draw(as, d);

      glPopMatrix();
    }

    virtual void _draw(DrawAs &as, DrawBank &d) = 0;

};


class Particle : public Placed {
  public:
    vector<Point> points;

    Particle() {
      level = level_particle;
    }

    virtual void _draw(DrawAs &as, DrawBank &d) {
      as.draw(points, d);
      d.end_of_particle();
    }

    Point &add_point() {
      points.resize(points.size() + 1);
      return points.back();
    }
};


class Cloud : public Placed {
  public:
    vector<Particle> particles;

    Cloud()
    {
      level = level_cloud;
    }

    Particle &add_particle() {
      particles.resize( particles.size() + 1 );
      return particles.back();
    }

    virtual void _draw(DrawAs &as, DrawBank &d) { 
      for (int i = 0; i < particles.size(); i++) {
        particles[i].draw(as, d);
      }
      d.end_of_cloud();
    }

    void step() {
    }

};

class ParticleGenesis {
  public:
    virtual void generate(Cloud &in) = 0;
};

class PointGenesis {
  public:
    virtual void generate(Particle &in) = 0;
};

class RandomPoints : public PointGenesis {
  public:
    int n;
    RandomPoints() {
      n = 4;
    }

    virtual void generate(Particle &in) {
      for (int i = 0; i < n; i++) {
        Point &p = in.add_point();
        p.random();
        p /= 2;
      }
    }
};

class Block : public PointGenesis {
  public:
    Block() {
    }

    virtual void generate(Particle &in) {
#define P(x,y,z) {\
          in.add_point().set(x, y, z); \
        }

#define A  P(-.5, -.5, -.5);
#define B  P(-.5, -.5,  .5);
#define C  P(-.5,  .5,  .5);
#define D  P(-.5,  .5, -.5);
#define E  P( .5, -.5, -.5);
#define F  P( .5, -.5,  .5);
#define G  P( .5,  .5,  .5);
#define H  P( .5,  .5, -.5);

      A B C D
      A B F E
      B C G F
      E F G H
      C D H G
      D A E H

#undef A
#undef B
#undef C
#undef D
#undef E
#undef F
#undef G
#undef H
#undef P
    }
};

class RandomParticles : public ParticleGenesis {
  public:
    int n;
    double pos_range;
    double scale_min, scale_max;
    Pt dir_range;
    PointGenesis *point_genesis;

    RandomParticles() {
      defaults();
      point_genesis = NULL;
    }

    RandomParticles(PointGenesis *point_genesis) {
      defaults();
      this->point_genesis = point_genesis;
    }

    void defaults() {
      n = 50;
      pos_range = 10.;
      scale_min = -1.;
      scale_max = 1.;
      dir_range.set(360., 360., 360.);
    }

    virtual void generate(Cloud &in) {
      for (int i = 0; i < n; i++) {
        Particle &p = in.add_particle();
        p.pos.random();
        p.pos *= pos_range;

        p.scale.random();
        p.scale *= scale_max - scale_min;
        p.scale += scale_min;

        p.rot3.random();
        p.rot3.scale3(dir_range);

        if (point_genesis) {
          point_genesis->generate(p);
        }

      }
    }
};



class Param {
  public:
    double val;
    double want_val;

    double slew;
    double change;

    bool do_limit_min;
    double val_min;

    bool do_limit_max;
    double val_max;


    Param(double val = 0, double slew=0.85){
      this->val = known_want_val = want_val = val;
      this->slew = slew;
      change = 0;
      do_limit_min = false;
      val_min = -1;
      do_limit_max = false;
      val_max = 1;
    }

    void step() {
      want_val += change;
      if (do_limit_min)
        want_val = max(val_min, want_val);
      if (do_limit_max)
        want_val = min(val_max, want_val);
      val = slew * val + (1. - slew) * want_val;
    }

    bool changed() {
      return known_want_val != want_val;
    }

    void set_unchanged() {
      known_want_val = want_val;
    }

    void limit(double min_val, double max_val) {
      val_min = min_val;
      do_limit_min = true;
      val_max = max_val;
      do_limit_max = true;
    }

    void limit_min(double min_val) {
      val_min = min_val;
      do_limit_min = true;
    }

    void limit_max(double max_val) {
      val_max = max_val;
      do_limit_max = true;
    }

    Param& operator=(double v) {
      want_val = v;
      return *this;
    }

    Param& operator=(Param &v) {
      want_val = v.want_val;
      return *this;
    }

    Param& operator+=(Param &v) {
      want_val += v.want_val;
      return *this;
    }

    Param& operator+=(double v) {
      want_val += v;
      return *this;
    }

    operator double() const {
      return val;
    }


  private:
    double known_want_val;
};

class Mass {
  public:
    double m;
    Pt v;
    Pt v_ang;

    Mass() 
      :m(1), v(0, 0, 0)
    {}

    void accelerate(const Pt &F, double dt)
    {
      v += F * (dt / m);
    }

    Pt impulse() const
    {
      return v * m;
    }
};


/**
adapted from http://www.euclideanspace.com/physics/dynamics/collision/threed/index.htm

This function calulates the velocities after a 3D collision vaf, vbf, waf and wbf from information about the colliding bodies
@param double e coefficient of restitution which depends on the nature of the two colliding materials
@param double ma total mass of body a
@param double mb total mass of body b
@param matrix Ia inertia tensor for body a in absolute coordinates (if this is known in local body coordinates it must
                 be converted before this is called).
@param matrix Ib inertia tensor for body b in absolute coordinates (if this is known in local body coordinates it must
                 be converted before this is called).
@param vector ra position of collision point relative to centre of mass of body a in absolute coordinates (if this is
                 known in local body coordinates it must be converted before this is called).
@param vector rb position of collision point relative to centre of mass of body b in absolute coordinates (if this is
                 known in local body coordinates it must be converted before this is called).
@param vector n normal to collision point, the line along which the impulse acts.
@param vector vai initial velocity of centre of mass on object a
@param vector vbi initial velocity of centre of mass on object b
@param vector wai initial angular velocity of object a
@param vector wbi initial angular velocity of object b
@param vector vaf final velocity of centre of mass on object a
@param vector vbf final velocity of centre of mass on object a
@param vector waf final angular velocity of object a
@param vector wbf final angular velocity of object b
*/
void collide(double restitution,
             double ma,
             double mb,
             const Matrix33 &Ia,
             const Matrix33 &Ib,
             const Pt &ra,
             const Pt &rb,
             const Pt &n,
             const Pt &vai,
             const Pt &vbi,
             const Pt &wai,
             const Pt &wbi,
             Pt &vaf,
             Pt &vbf,
             Pt &waf,
             Pt &wbf)
{
  Matrix33 inv_Ia = Ia.inverse();
  Matrix33 inv_Ib = Ib.inverse();

  Pt anga_a = inv_Ia * n.cross(ra);
  Pt angv_a = anga_a.cross(ra);  // calculate the linear velocity of collision point on a due to rotation of a

  Pt anga_b = inv_Ib * n.cross(rb);
  Pt angv_b = anga_b.cross(rb);  // calculate the linear velocity of collision point on b due to rotation of b

  double scalar = 1/ma + angv_a.dot(n) + 1/mb + angv_b.dot(n);

  double Jmod = (restitution+1) * (vai-vbi).len() / scalar;
  Pt J = n * Jmod;
  vaf = vai - J * (1/ma);
  vbf = vbi - J * (1/mb);
  waf = wai - anga_a;
  wbf = wbi - anga_b;
}

// vim: expandtab shiftwidth=2 tabstop=2
