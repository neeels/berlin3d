int W = 800;
int H = 500;//1200;

#include <stdlib.h>
#include <time.h>

#include <all_sections.h>

#define DRAW_SPHERES 0
#if DRAW_SPHERES
#include <GL/glut.h>
#endif

#define Pf(V) printf(#V "=%f\n", (float)V)
#define Pi(V) printf(#V "=%d\n", (int)V)

#include "draw.h"
//#include "audio.h"

static AsQuads as_quads;
static AsLines as_lines;


bool running = true;
volatile int frames_rendered = 0;
volatile int avg_frame_period = 0;
volatile double dt;
#define AVG_SHIFTING 3
float want_fps = 25;

struct Face {
  static const unsigned int N = 3;
  unsigned int point_idx[Face::N];
  PtGl n;

  void set(unsigned int i0, unsigned int i1, unsigned int i2, const Pt &normal)
  {
    point_idx[0] = i0;
    point_idx[1] = i1;
    point_idx[2] = i2;
    n = normal.unit();
  }

  void set(unsigned int i0, unsigned int i1, unsigned int i2, const vector<Point> &points)
  {
    Pt a = points[i1] - points[i0];
    Pt b = points[i2] - points[i0];
    set(i0, i1, i2, b.cross(a));
  }

  void load(const unsigned int idxs[Face::N], const Pt &normal)
  {
    memcpy(point_idx, idxs, sizeof(point_idx));
    n = normal.unit();
  }
};


struct Orientation {
  Pt nose;
  Pt top;

  Orientation()
  {
    clear();
  }

  void clear()
  {
    nose.set(0, 0, -1);
    top.set(0, 1, 0);
  }

  Pt right() {
    return nose.cross(top);
  }

  void rotate(double roll_x, double roll_y, double roll_z)
  {
      nose.rot_about(top, roll_y);
      Pt r = right();
      nose.rot_about(r, roll_x);
      top.rot_about(r, roll_x);
      top.rot_about(nose, roll_z);
  }

  void rotate(Pt r3)
  {
    rotate(r3.x, r3.y, r3.z);
  }

  void rotate_e(Pt e)
  {
    nose.rot_e(e);
    top.rot_e(e);
  }

  void glRotated() const {
    PtGl(rot3()).glRotated();
  }

  void check_normals() {
      nose = nose.unit();
      top = top.unit();

      double is_ortho = nose.dot(top);
      if (fabs(is_ortho) > 1.e-6)
        printf("ORTHO! %g\n", is_ortho);
  }


  Pt rot3() const {
    double l = nose.len();
    if (fabs(l) < 1.e-6) {
      return Pt();
    }

    /* "first roll left/right about the nose axis by angle az (around z axis),
        then turn the nose up/down by angle ax (around the x axis),
        then turn left/right by angle ay (around the vertical y axis)."

        ^ y
        |
        | az
     ax +------> x
       / ay
      v
      z

      This vector is the result of those rotations, the top vector points
      "upwards", and we're trying to find out ax, ay, and az.
     */

    Pt u = nose.unit();

    /* ax is the angle between the y=0 plane and this vector.
                                              
                  .                             
                 /|                                
               /  |                               
      y ^    /    |
        |  /      |
        |/        |
        +__ ax    |. . . . .  --> x
       /   --__   |
      v        --_|
      z

      */
    double ax = asin(u.y);


    /* ay is the angle between the z axis and the projection of this vector
       onto the y=0 plane. 
                  .                             
                 /                                 
               /                                  
      y ^    /     
        |  /       
        |/         
        +__. . . . . . . . .  --> x
       /   --__    
      /  ay    --_ 
     /-------------
    v
    z

    */
    double ay;
    bool pos = u.z > 1e-6;
    bool neg = u.z < -1e-6;
    if (!(pos || neg))
      ay = (u.x > 0? 1. : -1.) * M_PI/2;
    else
      ay = atan(u.x / fabs(u.z));

    if (neg)
      ay = M_PI - ay;

    /* Construct the top vector that would have az == 0, according to ax and
       ay: it is this vector but with an additional ax turn of pi/2 (90°).

       So it is a unit vector, starting from directly on y axis (pointing upward),
       tilted "backward" (towards negative z axis) by ax,
       and then turned around the y axis by ay.

       Turns out the zero-top's new y coordinate is the same as the length of
       this vector's projection onto the y=0 plane: ly = sqrt(x²+z²)

                    | y
             .___   |         . 
              \  ---+ly      /|                                
               \    ^      /  |                               
                \   |    /    |
                 \ax|  /      |
                  \ |/        |
          . . . . . +__ax     |. . . . .  --> x
                   /   --__   |
                  /        --_|
                 /            > ly
                v             
                z

       Also, the new distance from the z axis == this vector's y coordinate:

                    | 
          "y".___   |         .y
              \  -->+    /   /|                                
               \    |   /  /  |                               
                \   |  / /    |
             -_<-\ay+>//      |
                -.\ |/        |
          . . . . . +__       |. . . . .  --> x
                   /   --__   |
                  /<-ay--> --_V
                 /            -
                v               
                z

      So use sin and cos to get x and z coords.
     */

    Pt top_zero(-u.y * sin(ay),
                sqrt(u.x*u.x + u.z*u.z),
                -u.y * cos(ay));

    /* Find the angle between the "zero" top and the user supplied top vector.
     * Both top and top_zero must be in the plane perpendicular to this vector.
     * Remove any component from top that's pointing in this vector's dir. */
    Pt topu = top.unit();
    if (fabs(u.dot(topu)) > 1e-6) {
      printf("ORTHO! %g " __FILE__ " line %d\n", fabs(u.dot(topu)), __LINE__);
    }
    //topu -= u * topu.project(u);

    /* Angle of the "zero" top to the plane-ized and unit-ized user supplied
     * top, angle sign relative to this vector. */
    double az = top_zero.angle(topu, u);

    /* I measured the angle from the z axis rightwards. That's counter the
     * canonical angle direction OpenGL uses, so let's correct for that... */
    return Pt(M_PI - ax, ay, - az);
  }

  Matrix33 rot_matrix() const
  {
    return Matrix33::from_rot3(rot3());
  }

  bool operator!=(const Orientation &other) const {
    return (nose != other.nose) && (top != other.top);
  }

  bool operator==(const Orientation &other) const {
    return !(operator!=(other));
  }
};

Textures *gl_textures = NULL;

class Visible {
  private:
    double _radius;
    vector<Point> rotated_points;
    Orientation rotated_ori;
    bool shape_changed;
  public:
    vector<Point> points;
    vector<Face> faces;
    Texture *texture;
    PtGl scale;
    PtGl pos;
    Orientation ori;

    double opacity;

    Visible()
      :scale(1, 1, 1),
       shape_changed(false),
       _radius(0),
       opacity(1),
       texture(0)
    {
      scale.set(1, 1, 1);
    }

    void clear()
    {
      points.resize(0);
      faces.resize(0);
    }

    Point &add_point()
    {
      points.resize(points.size() + 1);
      shape_changed = true;
      return points.back();
    }

    Face &add_face()
    {
      faces.resize(faces.size() + 1);
      return faces.back();
    }

    void color_scheme(const Pt &base_rgb=Pt::random(0.1, 1), double d=.02)
    {
      foreach (p, points) {
        p->c = (base_rgb * frandom(1.-d, 1.+d)).set_min(1);
      }
    }

    void draw()
    {
      double opacity = this->opacity;
      if (opacity < 1e-6)
        return;

      opacity = max(0., min(1., opacity));

      glPushMatrix();

      pos.glTranslated();
   
      ori.glRotated();
      scale.glScaled();

#if DRAW_SPHERES
      glutSolidSphere(.5, 20, 16);
#else
      if (texture)
        texture->begin();

      int l;

#define LINE 0
#if LINE

      l = faces.size();
      for (int i = 0; i < l; i++) {
        Draw::begin(GL_LINE_STRIP);
        Face &face = faces[i];
        face.n.glNormal();
        for (int j = 0; j < (Face::N    +1     ); j++) {
          Point &p = points[face.point_idx[j%3]];
          if (texture)
            p.t.glTexCoord();
          Draw::color_point(p, opacity);
        }
        Draw::end();
      }
#else
      Draw::begin(GL_TRIANGLES);

      l = faces.size();
      for (int i = 0; i < l; i++) {
        Face &face = faces[i];
        face.n.glNormal();
        for (int j = 0; j < (Face::N); j++) {
          Point &p = points[face.point_idx[j%3]];
          if (texture)
            p.t.glTexCoord();
          Draw::color_point(p, opacity);
        }
      }
      Draw::end();
#endif


      if (texture)
        texture->end();
#endif
      glPopMatrix();
    }

    void load_points(const double points[][3], int n_points)
    {
      for (int i = 0; i < n_points; i++)
        add_point().load(points[i]);
    }

    void load_faces(const double points[][3], int n_points,
                    const unsigned int faces[][3], int n_faces)
    {
      load_points(points, n_points);
      for (int i = 0; i < n_faces; i++) {
        Pt a = this->points[faces[i][0]] - this->points[faces[i][1]];
        Pt b = this->points[faces[i][1]] - this->points[faces[i][2]];
        add_face().load(faces[i], b.cross(a).unit());
      }
    }

    void load_colors(const double colors[][4], int n_colors, bool cycle=true)
    {
      int p = 0;
      while (p < points.size()) {
        for (int i = 0; i < n_colors; i++, p++) {
          points[p].c.load(colors[i]);
        }
        if (! cycle)
          break;
      }
    }

    double radius()
    {
      if (shape_changed) {
        shape_changed = false;
        _radius = 0;
        foreach(p, points) {
          _radius = max(_radius, p->scaled(scale).len());
        }
      }
      return _radius;
    }

    bool radius_overlaps(Visible &v)
    {
      return (pos - v.pos).len() < (radius() + v.radius());
    }

    Pt face_center(unsigned int face_idx) const
    {
      return face_center(faces[face_idx]);
    }

  private:
    Pt face_center(const Face &f) const
    {
      Pt c;
      for (int i = 0; i < Face::N; i++) {
        c += points[f.point_idx[i]];
      }
      return c / Face::N;
    }

    void update_rotated_points()
    {
      if (rotated_ori != ori) {
        rotated_ori = ori;
        Pt rot3 = rotated_ori.rot3();

        rotated_points.resize(points.size());
        for (int i = 0; i < rotated_points.size(); i++) {
          rotated_points[i] = points[i];
        }
      }
    }
};

void make_rectangle(Visible &v, bool clear=true)
{
  if (clear)
    v.clear();

  static double points[][3] = {
    { .5,  .5, 0.},
    { .5, -.5, 0.},
    {-.5, -.5, 0.},
    {-.5,  .5, 0.}
  };

  static unsigned int faces[][3] = {
    {0, 1, 2},
    {2, 3, 0}
  };

  v.load_faces(points, ARRAY_SIZE(points),
               faces, ARRAY_SIZE(faces));
}

void make_block(Visible &v, bool clear=true)
{
  if (clear)
    v.clear();

  static double points[][3] = {
    {-.5, -.5,  .5},
    {-.5, -.5, -.5},
    {-.5,  .5, -.5},
    {-.5,  .5,  .5},
    { .5, -.5, -.5},
    { .5, -.5,  .5},
    { .5,  .5,  .5},
    { .5,  .5, -.5},
  };

  // side: a rectangle from two triangles
#define side(A, B, C, D) \
  {A, B, C}, \
  {C, D, A} 

  static unsigned int faces[][3] = {
    side(0, 1, 2, 3),
    side(1, 4, 7, 2),
    side(4, 5, 6, 7),
    side(5, 0, 3, 6),
    side(3, 2, 7, 6),
    side(5, 4, 1, 0)
  };
#undef side

  v.load_faces(points, ARRAY_SIZE(points),
               faces, ARRAY_SIZE(faces));
}

void make_icosahedron(Visible &v, bool clear=true)
{
  if (clear)
    v.clear();

  // icosahedron facts...
  //
  // outer radius
  // R = a * (sqrt(10. + 2.*sqrt(5)) / 4.).
  // inner radius
  // r = a * ((sqrt(3.) * (3. + sqrt(5.))) / 12.);
  //
  // icosahedron cartesian coords
  // (±1, 0 ,±φ)
  // (0, ±φ, ±1)
  // (±φ, ±1, 0)
  // with  φ = (1. + sqrt(5.))/2. = 1.618033988749895 (golden ratio)
  //
  // yields a side length a = 2
  //
  // I'd like a unit outer diameter, R = .5;
  // wanted side length A from R:
  //    .5 = A * (sqrt(10. + 2.*sqrt(5)) / 4.)
  // =>  A = .5 * 4. / sqrt(10. + 2.*sqrt(5))
  //
  // With above cartesian coords, we'd get a == 2.
  // So any wanted coordinate C from above coordinate c would be:
  //   C = c * ( A / 2.)
  //
  // thus construct from
  // (±(A/2.), 0 ,±(φ * A/2.))
  //
  // A/2 = .5 * .5 * 4. / sqrt(10. + 2.*sqrt(5))
  //     = 1. / sqrt(10. + 2.*sqrt(5))
  //     = 0.2628655560595668
  //
  // φ * A/2. = ((1. + sqrt(5.)) / 2) * A/2
  //          = ((1. + sqrt(5.)) / 2) * (1. / sqrt(10. + 2.*sqrt(5)))
  //          = (1. + sqrt(5.)) / (2. * sqrt(10. + 2.*sqrt(5)))
  //          = 0.42532540417602

#define u 0.2628655560595668
#define g 0.42532540417602

  static double points[][3] = {
    { 0.,  g,  u},
    { 0.,  g, -u},
    { 0., -g, -u},
    { 0., -g,  u},
    {  u, 0.,  g},
    { -u, 0.,  g},
    { -u, 0., -g},
    {  u, 0., -g},
    {  g,  u, 0.},
    {  g, -u, 0.},
    { -g, -u, 0.},
    { -g,  u, 0.},
  };

  static unsigned int faces[][3] = {
    {0, 1, 8},
    {8, 1, 7},
    {7, 1, 6},
    {7, 6, 2},
    {2, 6, 10},
    {2, 10, 3},
    {3, 10, 5},
    {3, 5, 4},
    {4, 5, 0},
    {4, 0, 8},

    {8, 7, 9},
    {9, 7, 2},
    {9, 2, 3},
    {9, 3, 4},
    {9, 4, 8},

    {1, 0, 11},
    {1, 11, 6},
    {6, 11, 10},
    {10, 11, 5},
    {5, 11, 0},
  };
#undef u
#undef g

  v.load_faces(points, ARRAY_SIZE(points),
               faces, ARRAY_SIZE(faces));

  if (0)
  foreach(p, v.points) {
    *p = p->unit();
  }
}

void make_star(Visible &v, double ratio = 1.5, int n = 1)
{
  make_icosahedron(v);

  foreach (p, v.points) {
    *p /= ratio;
  }

  for (unsigned int ni = 0; ni < n; ni++) {

    int faces_n = v.faces.size();
    for (int fi = 0; fi < faces_n; fi++) {
      Face &f = v.faces[fi];

      v.add_point() = v.face_center(fi).unit() / 2;
      int p = v.points.size() - 1;

      int a = f.point_idx[0];
      int b = f.point_idx[1];
      int c = f.point_idx[2];

      f.set(a, b, p, v.points);
      v.add_face().set(b, c, p, v.points);
      v.add_face().set(c, a, p, v.points);
    }
  }
}

void make_sphere(Visible &v, unsigned int n)
{
  make_icosahedron(v);

  for (unsigned int ni = 0; ni < n; ni++) {

    int faces_n = v.faces.size();
    for (int fi = 0; fi < faces_n; fi++) {
      Face &face = v.faces[fi];

      /*
                a p0                p0             a
                 /\                /\             /\
                /  \              /  \           /  \
               /    \   ===>    p5----p3        f----d
              /      \          / \  / \       / \  / \
             /________\        /___\/___\     /___\/___\
          c p2       b p1     p2    p4   p1  c    e    b
       */

      int a = face.point_idx[0];
      int b = face.point_idx[1];
      int c = face.point_idx[2];
      Pt p0 = v.points[a];
      Pt p1 = v.points[b];
      Pt p2 = v.points[c];

      Pt p3 = (p0 + p1).unit()/2;
      Pt p4 = (p1 + p2).unit()/2;
      Pt p5 = (p2 + p0).unit()/2;

      v.add_point() = p3;
      int d = v.points.size() - 1;
      v.add_point() = p4;
      int e = v.points.size() - 1;
      v.add_point() = p5;
      int f = v.points.size() - 1;

      face.set(c, f, e, v.points);
      v.add_face().set(f, a, d, v.points);
      v.add_face().set(e, f, d, v.points);
      v.add_face().set(e, d, b, v.points);
    }
  }
}

void color_grey(Visible &v, double intens)
{
  Pt g = intens;
  for (int i = 0; i < v.points.size(); i ++) {
    v.points[i].c = g * frandom(0.9, 1.1);
  }
}

class Ufo : public Visible {
  public:

    void *user_data;

    Mass mass;
    bool fixed_position;
    bool fixed_angle;

    Ufo() :
      fixed_position(false),
      fixed_angle(false)
    {}

    bool moving(double min_len=1e-6)
    {
      return (!fixed_position) && (mass.v.len() > min_len);
    }

    bool rotating(double min_len=1e-6)
    {
      return (!fixed_angle) && (mass.v_ang.len() > min_len);
    }

    bool animate(double min_len=1e-6)
    {
      return moving(min_len) || rotating(min_len);
    }

    virtual void step() {
      if (!fixed_position)
        pos += mass.v * dt;
      if (!fixed_angle)
        ori.rotate_e(mass.v_ang * dt);
    };

    bool collide(Ufo &o)
    {
      if (fixed_position) {
        if (o.fixed_position)
          return false;
        // make sure the fixed_position one is always the other one.
        return o.collide(*this);
      }

      if (!radius_overlaps(o))
        return false;

      Pt dir = o.pos - pos;
      Pt n = dir.unit();

      /*
      printf("\n");
      Pt i1 = mass.impulse() + o.mass.impulse();
      printf("imp\t");mass.impulse().print();
      printf("o.imp\t");o.mass.impulse().print();
      printf("i1\t"); i1.print();
      printf("m %f  o.m %f\n", mass.m, o.mass.m);
      printf("v\t");mass.v.print();
      printf("o.v\t"); o.mass.v.print();
      printf("n\t"); n.print();
      printf("vdiff\t"); (o.mass.v - mass.v).print();
      printf("vdiffp\t"); (o.mass.v - mass.v).project(n).print();
      printf("vdiffl\t%f\n", (o.mass.v - mass.v).project(n).len());
      */
      Pt J = n * (2.* (mass.v - o.mass.v).project(n).len() / (1./mass.m + 1./o.mass.m));

      mass.v -= J / mass.m;
      mass.v *= .59;

      if (! o.fixed_position) {
        o.mass.v += J / o.mass.m;
        o.mass.v *= .59;
      }

      /*
      printf("J\t"); J.print();
      Pt i2 = mass.impulse() + o.mass.impulse();
      printf("v\t");mass.v.print();
      printf("o.v\t"); o.mass.v.print();
      printf("imp\t");mass.impulse().print();
      printf("o.imp\t");o.mass.impulse().print();
      printf("i2\t"); i2.print();
      printf("i1 %6.2f  i2 %6.2f\n",
             i1.len(), i2.len());
             */

      Pt at = ((pos + n * radius()) + (o.pos - n * o.radius())) / 2;
      /*
      printf("pos\t");pos.print();
      printf("at\t");at.print();
      printf("o.pos\t");o.pos.print();
*/
      if (!o.fixed_position) {
        pos = at - n * (1.01 * radius());
        o.pos = at + n * (1.01 * o.radius());
      }
      else
        pos = at - n * (2.02 * radius());

      return true;
    }

    void play_bump(double vol=1) const
    {
      double v = max(.01, scale.volume());
      double dens = max(.05, min(10., sqrt(mass.m / v))) * frandom(.9, 1.1);
      double l = scale.len();

			/*
      Mix *m = new Mix();
      m->add(new Sine(dens*100./max(.01,scale.x/l), vol, frandom()*2.*M_PI));
      m->add(new Sine(dens*100./max(.01,scale.y/l), vol, frandom()*2.*M_PI));
      m->add(new Sine(dens*100./max(.01,scale.z/l), vol, frandom()*2.*M_PI));

      //Audio.play(new Envelope(max(.01,.02/dens), max(.01,.03/dens), max(.01,1./dens), m));
			*/
    }
};

void wrap_val(double &v, double ref)
{
  double d = v - ref;
  if (d > .5)
    v -= 1.;
  else
  if (d < -.5)
    v += 1.;
}

void wrap_uv(Pt &uv, const Pt &rel)
{
  wrap_val(uv.x, rel.x);
  wrap_val(uv.y, rel.y);
}

class Backdrop : public Visible {
  public:

    Backdrop() {
      pos = 0;
      scale = 1e3;

      make_sphere(*this, 2);

      foreach (p, points) {
        p->t = p->uv();
        p->c = 1;
        // reverse the faces so they go inward.
        p->x = - p->x;
      }

      foreach (f, faces) {
        // wrap texture coordinates
        Point &p0 = points[f->point_idx[0]];

        for (int pi = 1; pi < Face::N; pi++) {
          Point &p = points[f->point_idx[pi]];
          Pt t = p.t;
          wrap_uv(t, p0.t);
          if (t != p.t) {
            Point &new_point = add_point();
            new_point = p;
            new_point.t = t;
            int new_idx = points.size() - 1;
            f->point_idx[pi] = new_idx;
          }
        }

        // fix top/bottom u s (v = 0 or v = 1)
        for (int pi = 0; pi < Face::N; pi++) {
          Point &p = points[f->point_idx[pi]];

          if ((p.t.y < (1. - 1e-6)) && (p.t.y > 1.e-6))
            continue;

          //printf("\n  p.t "); p.t.print();

          double u_avg = 0.;
          for (int pj = 0; pj < Face::N; pj++) {
            if (pj == pi)
              continue;
            Point &q = points[f->point_idx[pj]];
            //q.t.print();
            u_avg += q.t.x;
          }
          u_avg /= Face::N - 1;

          p.t.x = u_avg;
          //printf("=> p.t "); p.t.print();
          //printf("\n");
        }
      }



      //texture = gl_textures->load("backdrop.jpg");
      //texture = gl_textures->load("test_backdrop.png");
    }
};

class Fly : public Ufo {
  public:
  Pt top_lag;

  Param top_angle;
  Param roll_x;
  Param roll_y;
  Param roll_z;

  Param propulsion_forward;
  Param propulsion_break;

  Param wings;

  bool do_cruise;
  Param cruise_v;

  double engines_strength;

	Fly() : Ufo() {
    ori.nose.set(0, 0, -1);
    ori.top.set(0, 1, 0);
    top_lag = ori.top;
    roll_x.slew = .8;
    roll_y.slew = .8;
    roll_z.slew = .8;
    propulsion_forward.slew = .1;
    propulsion_break.slew = .1;
    cruise_v.slew = 0;
    wings.limit(0, 1);
    wings = .1;

    do_cruise = true;
    cruise_v = .05;
    mass.v.set(0, 0, -.05);
    mass.m = 5000;
    engines_strength = mass.m * 20;

    make_block(*this);
    color_scheme(Pt(.2, .6, .2));

    int l = points.size();
    for (int i = 0; i < l; i++) {
      Point &p = points[i];
      if (p.z < 0) {
        p.x *= .1;
        p.y = max(min(p.y, .01), -.01);
      }
      else {
        p.y = max(min(p.y, .03), -.03);
      }
    }
	}

  void clear()
  {
    mass.v = 0;
    mass.v_ang = 0;
    pos = 0;
    cruise_v = 0;
    ori.clear();
    top_lag = ori.right();
  }

  virtual void step() {
    mass.v_ang *= .5;
      roll_x.step();
      roll_y.step();
      roll_z.step();
      top_angle.step();
      propulsion_forward.step();
      propulsion_break.step();
      wings.step();
      cruise_v.step();

      ori.rotate(roll_x, roll_y, roll_z);
      ori.check_normals();

      const int div = 10;
      top_lag = (top_lag * (div-1) + ori.top) / div;


      Pt wings_force;

      if (dt < 1.e-5)
        return;

      if (wings > 1.e-5) {
        Pt v_want = ori.nose * mass.v.len();
        wings_force = (v_want - mass.v) * mass.m * wings / dt;
      }

      double want_propulsion_forward = propulsion_forward;
      if (do_cruise) {
        if ((want_propulsion_forward > 1e-6) || (propulsion_break > 1e-6)) {
          /* user is changing speed. Take current speed as new desired speed. */
          cruise_v = mass.v.len();
        }
        else {
          double diff = cruise_v - mass.v.len();
          if (diff > 1e-6) {
            /* we're slowing down below cruising speed. Hit the accelerator a bit. */
            want_propulsion_forward = (6. * mass.m * diff) / engines_strength; // get percentile
            if (want_propulsion_forward < propulsion_forward)
              want_propulsion_forward = propulsion_forward;
          }
        }
      }

      // 1. is pedal to the metal
      want_propulsion_forward = min(1., want_propulsion_forward);

      Pt forward_want = ori.nose * (want_propulsion_forward * engines_strength);
      Pt break_want = mass.v * mass.m * (-propulsion_break);

      Pt propulsion = forward_want + break_want;
      if (propulsion.len() > engines_strength)
        propulsion = propulsion.unit() * engines_strength;

      mass.accelerate(wings_force + propulsion, dt);

      Ufo::step();
  }

};


struct Camera {
  Pt at;
  Pt from;
  Pt top;
  Backdrop backdrop;

  void look_at(const Pt &look_at, const Pt &from_rel, const Pt &top) {
    at = look_at;
    from = at + from_rel;
    this->top = top;
    backdrop.pos = from;
  }

  void look(const Pt &pos, const Pt &nose, const Pt &top) {
    at = pos + nose;
    from = pos;
    this->top = top;
    backdrop.pos = from;
  }

  void gluLookAt() {
    ::gluLookAt(from.x, from.y, from.z,
                at.x, at.y, at.z,
                top.x, top.y, top.z);
  }
};

class Light {
  public:
    Pt anchor;
    Pt ofs;
    GLuint id;
    bool do_wrap;

    void on() {
      glEnable(id);
    }

    void off() {
      glDisable(id);
    }

    void init(const Pt &at, bool do_wrap=false, const Pt &anchor=Pt())
    {
      this->anchor = anchor;
      this->ofs = at - anchor;
      this->do_wrap = do_wrap;
    }

    void step()
    {
      Pt pos = anchor + ofs;
      glLight(GL_POSITION, pos.x, pos.y, pos.z);
    }

    void glLight(GLuint do_id, float r, float g, float b)
    {
      GLfloat v[] = { r, g, b, 1.0 };
      glLightfv(id, do_id, v);
    }

    void ambient(double v)
    {
      ambient(v, v, v);
    }

    void diffuse(double v)
    {
      diffuse(v, v, v);
    }

    void specular(double v)
    {
      specular(v, v, v);
    }


    void ambient(float r, float g, float b)
    {
      glLight(GL_AMBIENT, r, g, b);
    }

    void diffuse(float r, float g, float b)
    {
      glLight(GL_DIFFUSE, r, g, b);
    }

    void specular(float r, float g, float b)
    {
      glLight(GL_SPECULAR, r, g, b);
    }

    void atten_const(double f)
    {
      glLightf(id, GL_CONSTANT_ATTENUATION, f);
    }

    void atten_lin(double f)
    {
      glLightf(id, GL_LINEAR_ATTENUATION, f);
    }

    void atten_quadr(double f)
    {
      glLightf(id, GL_QUADRATIC_ATTENUATION, f);
    }

    void atten(double consd=0, double lin=0, double quadr=0)
    {
      glLightf(id, GL_CONSTANT_ATTENUATION, consd);
      glLightf(id, GL_LINEAR_ATTENUATION, lin);
      glLightf(id, GL_QUADRATIC_ATTENUATION, quadr);
    }

};


class Lights {
  public:
    static vector<Light> lights;
    static const int GL_LIGHT_ID_LIST[6];

    Lights()
    {
      clear();
      add_default_light();
    }

    void add_default_light()
    {
      /* a default light. */
      Light &l = add();
      l.init(Pt(5e8, 5e8, -1e8));
      l.ambient(.5);
      l.diffuse(1);
      l.specular(1);
      l.atten(1, 0, 0);
    }

    void step()
    {
      foreach(l, lights) {
        l->step();
      }
    }

    Light &add()
    {
      int lights_n = lights.size();
      if (lights_n >= ARRAY_SIZE(GL_LIGHT_ID_LIST)) {
        printf("Too many lights: %d\n", lights_n + 1);
        exit(1);
      }
      lights.resize(lights_n + 1);
      Light &l = lights.back();
      l.id = GL_LIGHT_ID_LIST[lights_n];
      l.on();
      return l;
    }

    void clear()
    {
      foreach(l, lights) {
        l->off();
      }
      lights.resize(0);
    }

    void wrap(const Pt &center, double dist) {
      foreach(l, lights) {
        if (l->do_wrap)
          l->anchor.wrap_cube(center, dist);
        else
          l->anchor -= center;
      }
    }
};

const int Lights::GL_LIGHT_ID_LIST[6] = {
      GL_LIGHT0,
      GL_LIGHT1,
      GL_LIGHT2,
      GL_LIGHT3,
      GL_LIGHT4,
      GL_LIGHT5
    };
vector<Light> Lights::lights;

class World {
  public:
    Lights lights;
    vector<Ufo*> ufos;
    Pt wrap_ofs;

    World() {
    }

    Pt absolute_pos(const Pt &p)
    {
      return wrap_ofs + p;
    }

    Pt relative_pos(const Pt &abs_p)
    {
      return abs_p - wrap_ofs;
    }

    void clear()
    {
      ufos.resize(0);
      wrap_ofs = 0;
      lights.clear();
    }

    void add(Ufo &u) {
      ufos.resize(ufos.size() + 1);
      ufos.back() = &u;
    }

    void step() {
      foreach(u, ufos) {
        (*u)->step();
      }
    }

    void draw() {
      lights.step();
      foreach(u, ufos) {
        (*u)->draw();
      }
    }

    void wrap(const Pt &center, double dist) {
      wrap_ofs += center;
      lights.wrap(center, dist);
      foreach(u, ufos) {
        (*u)->pos.wrap_cube(center, dist);
      }
    }

};

class Game {
  public:
    int level;
    World &world;
    Camera cam;
    double r_visible;
    bool done;
    int won;
    double redraw_dist;

    Game (World &w) :
      level(0),
      world(w),
      r_visible(15),
      done(false),
      won(0),
      redraw_dist(0)
    {
    }

    virtual void on_joy_axis(int axis, double axis_val) {};
    virtual void on_joy_button(int button, bool down) {};
    virtual void on_key(int keysym, bool down) {};

    void win()
    {
      won = 1;
      done = true;
    }

    void lose()
    {
      won = -1;
      done = true;
    }

    void quit()
    {
      won = 0;
      done = true;
    }

    virtual void init_params()
    {
      level = max(0, 101 + level) % 101;
      printf("Level %d\n", level);
    }

    /* Set up local objects and settings for set difficulty level. */
    virtual void init()
    {
      init_params();
      won = 0;
    }

    /* Clear out the world, then link everything in place. */
    virtual void play()
    {
      done = false;
      init_gl();
    }

    virtual void step() = 0;
    virtual void osd_draw() {};

    void run() {
      world.step();
      step();
      draw();
    }

    virtual void on_collision(Ufo &u, Ufo &v) {
      u.play_bump(min(1., 1./((u.pos - cam.from).len()))/3);
      v.play_bump(min(1., 1./((v.pos - cam.from).len()))/3);
    };

    void init_gl()
    {
      glMatrixMode( GL_PROJECTION );
      glLoadIdentity();
      printf("r_visible %f\n", r_visible);
      gluPerspective(80, (double)W/H, .5, r_visible);

      glClearColor (0.0, 0.0, 0.0, 0.0);
      glShadeModel (GL_SMOOTH);

      glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE ) ;
      glEnable(GL_COLOR_MATERIAL);
      glEnable(GL_NORMALIZE);

      glEnable(GL_DEPTH_TEST);
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glMatrixMode( GL_MODELVIEW );
    }


    virtual void draw_scene()
    {
      world.draw();

      if (redraw_dist > .1) {
        Pt dir = cam.at - cam.from;
        PtGl p[] = {
          PtGl(-redraw_dist, -redraw_dist, -redraw_dist),
          PtGl(           0, -redraw_dist, -redraw_dist),
          PtGl( redraw_dist, -redraw_dist, -redraw_dist),
          PtGl(-redraw_dist,            0, -redraw_dist),
          PtGl(           0,            0, -redraw_dist),
          PtGl( redraw_dist,            0, -redraw_dist),
          PtGl(-redraw_dist,  redraw_dist, -redraw_dist),
          PtGl(           0,  redraw_dist, -redraw_dist),
          PtGl( redraw_dist,  redraw_dist, -redraw_dist),

          PtGl(-redraw_dist, -redraw_dist,            0),
          PtGl(           0, -redraw_dist,            0),
          PtGl( redraw_dist, -redraw_dist,            0),
          PtGl(-redraw_dist,            0,            0),

          PtGl( redraw_dist,            0,            0),
          PtGl(-redraw_dist,  redraw_dist,            0),
          PtGl(           0,  redraw_dist,            0),
          PtGl( redraw_dist,  redraw_dist,            0),

          PtGl(-redraw_dist, -redraw_dist,  redraw_dist),
          PtGl(           0, -redraw_dist,  redraw_dist),
          PtGl( redraw_dist, -redraw_dist,  redraw_dist),
          PtGl(-redraw_dist,            0,  redraw_dist),
          PtGl(           0,            0,  redraw_dist),
          PtGl( redraw_dist,            0,  redraw_dist),
          PtGl(-redraw_dist,  redraw_dist,  redraw_dist),
          PtGl(           0,  redraw_dist,  redraw_dist),
          PtGl( redraw_dist,  redraw_dist,  redraw_dist),
        };

        for (int i = 0; i < ARRAY_SIZE(p); i++) {
          if (! p[i].project(dir).zero()) {
            glPushMatrix();
            p[i].glTranslated();
            world.draw();
            glPopMatrix();
          }
        }
      }
    }

    void draw()
    {
      glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );


      glLoadIdentity( );
      
      cam.gluLookAt();

      glDisable(GL_LIGHTING);
      //glDisable(GL_DEPTH_TEST);
      cam.backdrop.draw();
      //glEnable(GL_DEPTH_TEST);
      glEnable(GL_LIGHTING);

      draw_scene();

      glDisable(GL_LIGHTING);
      osd_draw();
      glFlush();
      SDL_GL_SwapBuffers();
    }

    void collide() {
      for (int i = 0; i < world.ufos.size(); i++) {
        for (int j = i+1; j < world.ufos.size(); j++) {
          Ufo *a = world.ufos[i];
          Ufo *b = world.ufos[j];
          if (a->moving() || b->moving()) {
            if (a->collide(*b))
              on_collision(*a, *b);
          }
        }
      }
    }
};

class BlockSpace : public Game {
  public:
    double r_wrap;
    double r_action;
    double max_block_size;
    unsigned int blocks_count;
    double reverse_roll_x;

    double v_drag;
    double v_ang_drag;

    Fly fly;
    vector<Ufo> debris;

    BlockSpace(World &w)
      : Game(w),
        r_wrap(30),
        r_action(10),
        max_block_size(2),
        blocks_count(23),
        reverse_roll_x(1.)
    {
      //redraw_dist = world_r * 2.;//10.0;
        r_visible = 1000;
    }

    virtual void init_params()
    {
      Game::init_params();

      v_drag = 0;
      v_ang_drag = 0;

      r_wrap = 3. * r_action;

      blocks_count = (1./3) * (r_action * r_action * r_action) / (max_block_size * max_block_size * max_block_size);

      fly.clear();
    }

    virtual void init()
    {
      Game::init();
      make_blocks();
      cam.backdrop.color_scheme(Pt(0x7f, 0xbf, 0xff)/255);
      cam.backdrop.texture = NULL;
    }

    void make_blocks()
    {
      debris.resize(blocks_count);
      foreach(b, debris) {
        make_block(*b);
        b->color_scheme();
        b->ori.rotate(Pt::random() * 2 * M_PI);
        b->mass.v = 0;
        b->mass.v_ang = Pt::random();
        b->pos = Pt::random(-r_action, r_action);
        b->scale = Pt::random(max_block_size / 10, max_block_size);

        double dens = frandom(.5, 1.5);
        b->mass.m = dens * b->scale.volume();
      }
    }

    virtual void add_lights()
    {
      world.lights.add_default_light();
    }

    virtual void play()
    {
      Game::play();

      world.clear();
      add_lights();

      foreach(b, debris) {
        world.add(*b);
      }
      world.add(fly);

      osd_init();
    }

    virtual void step()
    {
      if ((fabs(v_drag) > 1e-6) || (fabs(v_ang_drag) > 1e-6)) {
        double v_f = 1. - v_drag;
        double v_ang_f = 1. - v_ang_drag;
        foreach(u, world.ufos) {
          (*u)->mass.v *= v_f;
          (*u)->mass.v_ang *= v_ang_f;
        }
#if 0
        fly.mass.v *= v_f;
        fly.mass.v_ang *= v_ang_f;
#endif
      }

      collide();

      {
        Pt wrap_center = fly.pos;
        if (wrap_center.len() > 1)
          world.wrap(wrap_center, r_wrap);
      }

      Pt dir = fly.mass.v.unit();
      if (! dir.zero())
        dir = (dir + fly.ori.nose.unit()) / 2;
      else
        dir = fly.ori.nose;

      cam.look_at(fly.pos + fly.top_lag,
                  dir.unit() * (-3),
                  fly.top_lag);

      cam.backdrop.pos = world.wrap_ofs * -1;
      osd_update();
    }

    /* OSD */
    Ufo speed;
    bool draw_want_speed;
    Ufo want_speed;

    void osd_init()
    {
      make_block(speed);
      speed.pos.set(0, -.68, -1);
      make_block(want_speed);
      want_speed.pos = speed.pos - Pt(0, 0, .002);
    }

    void osd_update() {
      double v = sqrt(fly.mass.v.len()) / 10;
      speed.scale.set(.001 + v/10, .02, .001);
      double c[][4] = {{1, min(.8*v, 1.), .3, .7}, {.95, min(.8*v * 1.1, 1.), .28, .7}};
      speed.load_colors(c, ARRAY_SIZE(c));

      if (fly.do_cruise) {
        v = sqrt(fly.cruise_v) / 10;
        want_speed.scale.set(.001 + v/10, .02, .001);
        double d[][4] = {{1, min(v, 1.), .0, 1}, {.95, min(v * 1.1, 1.), .0, 1}};
        want_speed.load_colors(d, ARRAY_SIZE(d));
        draw_want_speed = true;
      }
      else
        draw_want_speed = false;
    }

    virtual void osd_draw()
    {
      glLoadIdentity();

      speed.draw();
      if (draw_want_speed)
        want_speed.draw();
    }


    /* User input */

    virtual void on_joy_axis(int axis, double axis_val)
    {
      double v;

      switch(axis)
      {
      default:
        printf("%d %f\n", axis, (float)axis_val);
        break;

      case 0:
        fly.roll_y = -(axis_val*axis_val*axis_val) / 5;
        fly.roll_z = (axis_val*axis_val*axis_val) / 50;
        break;

      case 1:
      case 4:
        fly.roll_x = reverse_roll_x * (axis_val*axis_val*axis_val) / 5;
        break;

      case 3:
        fly.roll_z = (axis_val*axis_val*axis_val) / 5;
        break;

      case 5:
        // analog trigger ... -1 == not pressed, 0 = half, 1 = full
        v = (axis_val + 1) / 2;
        fly.propulsion_forward = v;
        break;
      case 2:
        // analog trigger ... -1 == not pressed, 0 = half, 1 = full
        fly.propulsion_break = (axis_val + 1) / 2;
        break;
      }
    }

    void on_joy_button(int button, bool down)
    {
      switch (button) {
      case 3:
        if (! down) {
          if (fly.wings < 1.e-5) {
            fly.wings = .1;
            fly.do_cruise = true;
          }
          else
          if (fly.wings < .9) {
            fly.wings = 1;
            fly.do_cruise = true;
          }
          else {
            fly.wings = 0;
            fly.do_cruise = false;
          }
        }
        break;

      case 1:
        fly.propulsion_forward = down? 1 : 0;
        break;

      case 0:
        fly.propulsion_break = down? 1 : 0;
        break;

      case 2:
        if (down)
          fly.play_bump(1.);
        break;

      case 6:
        if (down) {
          reverse_roll_x = (reverse_roll_x > 0)? -1. : 1.;
          printf("roll_x %f\n", reverse_roll_x);
        }
        break;

      default:
        printf("button %d\n", button);
        break;
      }
    }

    virtual void on_key(int keysym, bool down)
    {
      switch (keysym) {
      case SDLK_UP:
      case 'w':
      case 'k':
        fly.roll_x = down? -.1 : 0;
        break;

      case SDLK_LEFT:
      case 'a':
      case 'h':
        fly.roll_y = down? .1 : 0;
        break;

      case SDLK_DOWN:
      case 's':
      case 'j':
        fly.roll_x = down? .1 : 0;
        break;

      case SDLK_RIGHT:
      case 'd':
      case 'l':
        fly.roll_y = down? -.1 : 0;
        break;

      case 'q':
        fly.roll_z = down? -.1 : 0;
        break;
      case 'e':
        fly.roll_z = down? .1 : 0;
        break;

      case SDLK_PAGEUP:
      case '+':
      case '=':
      case ' ':
      case '1':
        fly.propulsion_forward = down? 1 : 0;
        break;

      case SDLK_PAGEDOWN:
      case '-':
      case 'b':
      case '2':
        fly.propulsion_break = down? 1 : 0;
        break;

      default:
        printf("key %d %c\n", keysym, (char)keysym);
        break;

      }
    }
};


class MoveAllBlocks : public BlockSpace {
  public:
    int inanimate;

    MoveAllBlocks(World &w) :BlockSpace(w)
    {}

    virtual void init_params()
    {
      r_action = 5. + .3 * level;
      BlockSpace::init_params();

      if (level == 0) {
        blocks_count = 2;
      }
      else
        blocks_count = min(blocks_count, (unsigned int)2 + 5 * level);

      v_drag = .00003 * level;
      v_ang_drag = .00003 * level;
    }

    virtual void init()
    {
      BlockSpace::init();

      for (int i = 0; i < debris.size(); i++) {
        Ufo &b = debris[i];
        b.mass.v_ang = 0;

        if (i == 0) {
          b.pos.set(0, 0, -5);
          //b.scale = 1.0;
          //b.mass.m = 1;
        }

        color_grey(b, .1);
      }
    }

    virtual void play()
    {
      BlockSpace::play();

      move_osd_init();
    }

    void count_inanimate()
    {
      inanimate = 0;
      foreach(_u, world.ufos) {
        Ufo &u = **_u;
        if (! u.animate(1e-2)) {
          if (u.animate(0)) {
            u.mass.v = 0;
            u.mass.v_ang = 0;
            color_grey(u, .1);
          }
          inanimate ++;
        }
      }
    }

    virtual void step()
    {
      count_inanimate();
      if (inanimate == 0)
        win();

      BlockSpace::step();

      move_osd_update();
    }

    virtual void on_collision(Ufo &u, Ufo &v) {
      BlockSpace::on_collision(u, v);
      u.mass.v_ang = Pt::random();
      v.mass.v_ang = Pt::random();
      u.color_scheme();
      v.color_scheme();
    };


    // OSD

    int was_inanimate;
    double show_got_one;
    Ufo got_one;
    Ufo bar_all;
    Ufo bar_inanimate;

    void move_osd_init()
    {
      make_block(got_one);
      make_block(bar_all);
      make_block(bar_inanimate);
      bar_all.pos.set(0, -.71, -1.004);
      bar_all.scale.set(.0001 + 2., .03, .001);

      bar_inanimate.pos.set(0, -.71, -1.002);

      got_one.pos.set(0, -.71, -1);
      show_got_one = 0;

      double c[][4] = {{.1, 0.7, .1, .4},};
      bar_all.load_colors(c, ARRAY_SIZE(c));
      double d[][4] = {{.7, 0.1, .1, .4},};
      bar_inanimate.load_colors(d, ARRAY_SIZE(c));

      count_inanimate();
      was_inanimate = inanimate;
    }

    void move_osd_update()
    {
      if (show_got_one > 0) {
        show_got_one -= dt;
        double c[][4] = {{.5, 1, .5, min(show_got_one, 1.)}, {.4, 1, .4, min(show_got_one, 1.)}};
        got_one.load_colors(c, ARRAY_SIZE(c));
      };
      if (was_inanimate != inanimate) {
        show_got_one = 1.1;
        was_inanimate = inanimate;
      }
      bar_inanimate.scale.set(.0001 + 2. * inanimate / world.ufos.size(), .03, .001);
      got_one.scale.set(.0001 + 2. * inanimate / world.ufos.size(), .03, .001);
    }

    virtual void osd_draw()
    {
      BlockSpace::osd_draw();

      bar_all.draw();
      bar_inanimate.draw();
      if (show_got_one > 0) {
        got_one.draw();
      }
    }
};


class TurnAllOn : public BlockSpace {
  public:
    bool fly_starts_on;
    int off_count;

    TurnAllOn(World &w) :
      BlockSpace(w),
      off_count(0),
      fly_starts_on(true)
    {}

    void off(Ufo &b)
    {
      b.color_scheme(Pt(.5, .1, .1));
      b.user_data = 0;
      off_count ++;
      b.mass.v_ang = 0;
    }

    void on(Ufo &b)
    {
      b.color_scheme(Pt(0, 1, 0));
      b.user_data = (void*)1;
      if (off_count)
        off_count --;
      b.mass.v_ang = .7;
    }

    void toggle(Ufo &b)
    {
      if (b.user_data)
        off(b);
      else
        on(b);
    }

    virtual void init_params()
    {
      r_action = 5. + .3 * level;
      BlockSpace::init_params();

      if (level == 0) {
        blocks_count = 2;
        v_drag = .01;
      }
      else {
        v_drag = .01 / (level * level);

        blocks_count = min(blocks_count, (unsigned int)3 * level);
      }
      fly_starts_on = ! (blocks_count & 1);
    }

    virtual void init()
    {
      BlockSpace::init();

      off_count = 0;
      off(fly);
      foreach(b, debris) {
        make_icosahedron(*b);
        b->mass.v = 0;
        off(*b);
      }

      if (fly_starts_on)
        on(fly);

      for (int i = 0; i < debris.size(); i++) {
        Ufo &b = debris[i];

        if (i == 0) {
          b.pos.set(0, 0, -5);
          b.scale = 1.0;
          b.mass.m = 1;
        }
        else
        if (i == 1) {
          b.pos.set(3, 3, -5);
          b.scale = 1.001;
          b.mass.m = 1.5;
        }
        else
        if (i == 2) {
          b.pos.set(-2, -2, -5);
          b.scale = 1.001;
          b.mass.m = 1.5;
        }

        b.scale = 3;
        b.mass.m = 200;
        b.ori.clear();
        //else
         // break;
      }

      cam.backdrop.color_scheme(Pt(1), 0);
      cam.backdrop.texture = gl_textures->load("backdrop_siegessaeule.jpg");
      //cam.backdrop.texture = gl_textures->load("backdrop.jpg");
    }

    virtual void play()
    {
      BlockSpace::play();
      turnon_osd_init();
    }

    virtual void step()
    {
      if (off_count == 0)
        win();

      BlockSpace::step();

      turnon_osd_update();
    }

    virtual void on_collision(Ufo &u, Ufo &v) {
      BlockSpace::on_collision(u, v);
      u.mass.v_ang = Pt::random();
      v.mass.v_ang = Pt::random();
      toggle(u);
      toggle(v);
    };

    virtual void add_lights()
    {
      {
        Light &l = world.lights.add();
        l.init(Pt(0, 200, 1000));
        l.ambient(.5);
        l.diffuse(1);
        l.specular(1);
        l.atten(1, 0, 0);
      }
    }


    // OSD

    int was_off_count;
    double show_got_one;
    Ufo got_one;
    Ufo bar_all;
    Ufo bar_off;

    void turnon_osd_init()
    {
      bar_all.pos.set(0, -.71, -1.004);
      bar_all.scale.set(.0001 + 2., .03, .001);

      bar_off.pos.set(0, -.71, -1.002);

      got_one.pos.set(0, -.71, -1);
      show_got_one = 0;

      double c[][4] = {{.1, 0.7, .1, .4},};
      bar_all.load_colors(c, ARRAY_SIZE(c));
      double d[][4] = {{.7, 0.1, .1, .4},};
      bar_off.load_colors(d, ARRAY_SIZE(c));

      was_off_count = off_count;
    }

    void turnon_osd_update()
    {
      if (show_got_one > 0) {
        show_got_one -= dt;
        double c[][4] = {{.5, 1, .5, min(show_got_one, 1.)}, {.4, 1, .4, min(show_got_one, 1.)}};
        got_one.load_colors(c, ARRAY_SIZE(c));
      };
      if (was_off_count != off_count) {
        show_got_one = 1.1;
        was_off_count = off_count;
      }
      bar_off.scale.set(.0001 + 2. * off_count / world.ufos.size(), .03, .001);
      got_one.scale.set(.0001 + 2. * off_count / world.ufos.size(), .03, .001);
    }

    virtual void osd_draw()
    {
      BlockSpace::osd_draw();

      bar_all.draw();
      bar_off.draw();
      if (show_got_one > 0) {
        got_one.draw();
      }
    }
};



class FindTheLight : public BlockSpace {
  public:
    Pt absolute_light_pos;

    FindTheLight(World &w) :BlockSpace(w)
    {
    }

    Ufo &marker()
    {
      return debris.back();
    }

    virtual void init_params()
    {
      r_action = 20.;
      BlockSpace::init_params();
      r_wrap = r_action;
      blocks_count = 152. - (150. * level)/100;
      printf("level %d blocks_count %d\n", level, blocks_count);
    }

    virtual void init()
    {
      BlockSpace::init();

      if (level == 0)
        absolute_light_pos.set(0, 0, -10);
      else
        absolute_light_pos = Pt::randoml(10 * level, 10 * level);

      foreach (b, debris) {
        b->mass.v = 0;
        b->mass.v_ang = 0;
      }

      Ufo &m = marker();

      make_star(m, 1.5);

      m.mass.v_ang = 1;
      m.scale = 1;
      m.pos = absolute_light_pos;
      m.color_scheme(Pt(1), 0);
      m.fixed_position = true;

      cam.backdrop.texture = gl_textures->load("backdrop.jpg");
      //cam.backdrop.texture = gl_textures->load("test_backdrop.png");

      if (cam.backdrop.texture)
        cam.backdrop.color_scheme(Pt(1), 0);
      else
        cam.backdrop.color_scheme(Pt(.05));
    }

    virtual void add_lights()
    {
      double intens = .3;
      if (level > 0)
        intens /= sqrt(level);

      {
        Light &l = world.lights.add();
        l.init(absolute_light_pos + Pt(1));
        l.ambient(intens/5);
        l.diffuse(intens);
        l.specular(intens);
        l.atten(0, 0, .01);
      }

      {
        Light &l = world.lights.add();
        l.init(absolute_light_pos - Pt(1));
        l.ambient(intens/5);
        l.diffuse(intens);
        l.specular(intens);
        l.atten(0, 0, .01);
      }

      {
        Light &l = world.lights.add();
        l.init(absolute_light_pos + Pt(1, -1, 1));
        l.ambient(.0);
        l.diffuse(.1);
        l.specular(0);
        l.atten_const(0.8);
      }
    }

    virtual void on_collision(Ufo &u, Ufo &v) {
      u.play_bump(min(1., 1./((u.pos - cam.from).len()))/3);
      v.play_bump(min(1., 1./((v.pos - cam.from).len()))/3);
      Ufo &m = marker();
      if ( ( (&u == &fly) || (&v == &fly) )
           && ( ( (&u == &m) || (&v == &m) )) )
        win();
    };

    virtual void step()
    {
      BlockSpace::step();

      Ufo &m = marker();
      m.pos = world.relative_pos(absolute_light_pos);

      double d = (m.pos - cam.from).len();
      m.opacity = 3. - 2. * (1. + .2 * sqrt(level)) * d / r_wrap;
    }

};

#include "city.h"
class Bezirk;

struct Building {
  const BuildingData *b;
  Color c;

  void draw(const Bezirk &bz);
};

class Bezirk {
  public:
    const DwellingData *d;
    vector<Building> buildings;
		vector<Texture*> textures;

    void setup(const DwellingData &data)
    {
      d = &data;
			printf("%s\n", d->name);

      printf("mid ");
      d->mid.print();

      buildings.resize(d->n_buildings);
      int i = 0;
      foreach(b, buildings) {
        b->b = &(d->buildings[i]);
        b->c.set(1, 1, 1);//random();
        i++;
      }

      Pt X(-541.03, -1665.42, 0.000);

      textures.resize(data.n_images);
			for (i = 0; i < data.n_images; i++) {
        textures[i] = NULL;
      }

      foreach (b, buildings) {
        if ((b->b->pos - X).cart_len() < 300) {
          int wall_idx = b->b->walls_start_idx;
          int wall_n = b->b->n_walls;
          while (wall_n--) {
            int img_idx = data.wall_indices[wall_idx++];
            
            const char *path = data.images[img_idx];
            if (path)
              textures[img_idx] = gl_textures->load(path);

            while (data.wall_indices[wall_idx++] != -2);
          }

        }
			}
#if 0
      bboxes is vector<Ufo>
      bboxes.resize(bezirk.n_buildings);
      for (int i = 0; i < bboxes.size(); i++) {
        Ufo &b = bboxes[i];

        const BuildingData *bldg = &bezirk.buildings[i++];
        b.pos = bldg->pos;
        b.scale = bldg->box_size;
        b.pos = b.pos.scaled(scale);
        b.scale = b.scale.scaled(scale);
        b.scale.set_max(0.01);
      }
#endif

    }

    void draw(const Pt &around, double dist)
    {
      foreach (b, buildings) {
        if ((b->b->pos - around).len() < dist)
          b->draw(*this);
      }
    }
};

void Building::draw(const Bezirk &bz) {
	const DwellingData &d = *bz.d;

	int wall_idx = b->walls_start_idx;
	int normal_idx = b->normals_start_idx;
	int normal_idx_end = normal_idx + b->n_walls;
	for (; normal_idx < normal_idx_end; normal_idx ++) {

		int img_idx = d.wall_indices[wall_idx++];
		int tex_points_idx = d.wall_indices[wall_idx++];

		Texture *texture = NULL;
		if ((img_idx >= 0) && (tex_points_idx >= 0))
			texture = bz.textures[img_idx];

		if (texture)
			texture->begin();

		//Draw::begin(GL_LINE_STRIP);
		//Draw::begin(GL_TRIANGLES);//GL_POLYGON);//GL_LINE_STRIP);

		Draw::begin(GL_POLYGON);

		const double *n = d.normals[normal_idx];
		::glNormal3d(n[0], n[1], n[2]);

		c.glColor();

		int point_idx;
		while ((point_idx = d.wall_indices[wall_idx++]) >= 0) {
      if (img_idx < 0)
        continue;
			if (texture) {
				const double *tp = d.tex_coords[tex_points_idx++];
				PtGl(tp[0], tp[1]).glTexCoord();
			}
			const double *point = d.points[point_idx];
			::glVertex3d(point[0], point[1], point[2]);
		}

		Draw::end();

    if (texture)
      texture->end();
	}
}

class Berlin : public Game {
  public:
    vector<Bezirk> bezirke;

    Pt pos;
    Ufo map;
    Orientation ori;
    double cam_angle;
    PtGl scale;
    Pt center;

    int points_count;

    Param move_forward;
    Param move_sideways;
    Param move_up;
    Param roll_z;
    Param cam_nod;

    Berlin(World &w) :
      Game(w),
      scale(.1),
      center(393100.658211, 5818560.009784, 94.150000)
    {
      r_visible = 100e3;
      points_count = 2e6;

      make_rectangle(map);

      Pt texture_coords[] = {
        {0, 0}, {0, 1}, {1, 1}, {1, 0}
      };

      for (int i = 0; i < 4; i++) {
        map.points[i].t = texture_coords[i];
      }

      map.texture = gl_textures->load("berlin-germany.jpg");

      populate();
    }

    void add(const DwellingData &d)
    {
      bezirke.resize(bezirke.size() + 1);
      Bezirk &b = bezirke.back();
      b.setup(d);
    }

    virtual void populate();

    virtual void add_lights()
    {
      {
        Light &l = world.lights.add();
        l.init(Pt(-1e5, -3e5, 5e5));
        l.ambient(.5);
        l.diffuse(1);
        l.specular(1);
        l.atten(1, 0, 0);
      }
    }

    virtual void init()
    {
      Game::init();

      pos.set(0, 0, 50);
      ori.nose.set(0, 1, 0);
      ori.top.set(0, 0, 1);
      cam_angle = -.3;

      map.ori.clear();
      /* xberg.png 
      map.pos.set(321.58213, 76.41885, -70);
      map.scale.set(8903.50546, 8336.67821, 0.001);
      map.ori.top.rot_about(map.ori.nose, 0.022361);
      */

      map.pos.set(-249.119087,-1717.710267,-70.000000);
      map.scale.set(73718.681667,73718.681667,73718.681667);
      map.ori.top.rot_about(map.ori.nose, -0.034907);



      //cam.backdrop.texture = gl_textures->load("backdrop_bat.jpg");
      //cam.backdrop.texture = gl_textures->load("backdrop_stone.jpg");
      //cam.backdrop.texture = gl_textures->load("backdrop_western.jpg");
      cam.backdrop.texture = gl_textures->load("backdrop.jpg");
      cam.backdrop.scale = 2e3;

      if (cam.backdrop.texture)
        cam.backdrop.color_scheme(Pt(1), 0);
      else
        cam.backdrop.color_scheme(Pt(0x7f, 0xbf, 0xff)/2550);

      if (map.texture)
        map.color_scheme(Pt(1), 0);
      else
        map.color_scheme(Pt(0x7f, 0xbf, 0xff)/255);
    }


    virtual void play()
    {
      Game::play();

      world.clear();
      add_lights();

      if (1)
      world.add(map);

    }

    virtual void step()
    {
      move_forward.step();
      move_up.step();
      move_sideways.step();
      roll_z.step();
      cam_nod.step();

      Pt z = {0, 0, 1};
      pos +=
        ori.nose.without(z) * move_forward
        + z * (move_up * max(1.,fabs(pos.z)))
        + ori.right() * move_sideways;

      ori.rotate(0, -roll_z, 0);
      ori.check_normals();

      cam_angle += cam_nod;
      cam_angle = min(M_PI/2 - .01, max(-M_PI/2 + .01, cam_angle));
      Pt cam_dir = ori.nose.rotated_about(ori.right(), cam_angle);

      cam.look(pos, cam_dir, ori.top);
      if (0) {
        printf("pos\n");
        pos.print();
        printf("%f\n", cam_angle);
        ori.nose.print();
      }

    }

    virtual void draw_scene()
    {
      glEnable(GL_LIGHTING);
      glPushMatrix();
      scale.glScaled();

      Game::draw_scene();

      Pt unscaled_pos = pos.unscaled(scale).without(Pt(0, 0, 1));
      double max_dist = 20 / scale.x;

      foreach (b, bezirke) {
        b->draw(unscaled_pos, max_dist);
      }

      glPopMatrix();
    }

    /* User input */

    virtual void on_joy_axis(int axis, double axis_val)
    {
      double v;
      double f = 10;

      switch(axis)
      {
      default:
        printf("%d %f\n", axis, (float)axis_val);
        break;

      case 0:
        move_sideways = f*(axis_val*axis_val*axis_val);
        break;

      case 1:
        move_forward = -f*(axis_val*axis_val*axis_val);
        break;

      case 4:
        move_up = (axis_val*axis_val*axis_val)/5;
        break;

      case 3:
        roll_z = (axis_val*axis_val*axis_val)/5;
        break;
      }
    }

    virtual void on_key(int keysym, bool down)
    {
      double zf = max(1., pos.z/10);
      switch (keysym) {
      case 'w':
      //case 'k':
        move_forward = down? zf : 0;
        break;

      case 's':
      //case 'j':
        move_forward = down? -zf : 0;
        break;

      case 'a':
      //case 'h':
        move_sideways = down? -zf : 0;
        break;

      case 'd':
      //case 'l':
        move_sideways = down? zf : 0;
        break;

      case SDLK_UP:
        cam_nod = down? -.1 : 0;
        break;

      case SDLK_DOWN:
        cam_nod = down? .1 : 0;
        break;

      case SDLK_LEFT:
        roll_z = down? -.15 : 0;
        break;

      case SDLK_RIGHT:
        roll_z = down? .15 : 0;
        break;

      case 'e':
        move_up = down? -.1 : 0;
        break;
      case 'q':
        move_up = down? .1 : 0;
        break;

      default:
        printf("key %d %c\n", keysym, (char)keysym);
        break;


      case SDLK_PAGEUP:
      case '+':
      case '=':
      case '1':
        if (down)
          points_count ++;
        break;

      case SDLK_PAGEDOWN:
      case '-':
      case '2':
        points_count --;
        break;

      case ' ':
        if (down) {
          Pt z(0, 0, 1);
          Pt ofs = (pos - map.pos).without(z);
          ofs.rot_about(z, - map.ori.rot3().z);
          ofs.x /= map.scale.x;
          ofs.y /= map.scale.y;

          printf("\nofs @ pos\n");
          ofs.print();
          pos.without(z).print();
          fflush(stdout);
        }
        break;

      }
   
      zf = pos.z/50;
      if (down){ 
        bool print = true;
        switch(keysym) {
        case 'i':
          map.pos.y += zf;
          break;
        case 'k':
          map.pos.y -= zf;
          break;
        case 'j':
          map.pos.x -= zf;
          break;
        case 'l':
          map.pos.x += zf;
          break;
        case 'o':
          map.scale *= 1. + .003 * zf;
          break;
        case 'u':
          map.scale /= 1. + .003 * zf;
          break;
        case 'p':
          map.scale.x *= 1. + .003 * zf;
          break;
        case 'y':
          map.scale.x /= 1. + .003 * zf;
          break;
        case 'h':
          map.ori.top.rot_about(map.ori.nose, -M_PI/360);
          break;
        case ';':
          map.ori.top.rot_about(map.ori.nose, M_PI/360);

          break;
        default:
          print = false;
          break;
        }

        if (print) {
          printf("      map.pos.set(%f,%f,%f);\n", map.pos.x, map.pos.y, map.pos.z);
          printf("      map.scale.set(%f,%f,%f);\n", map.scale.x, map.scale.y, map.scale.z);
          printf("      map.ori.top.rot_about(map.ori.nose, %f);\n", -map.ori.rot3().z);
        }
      }
    }
};


void Berlin::populate()
{

}

class Games {
  public:
    World world;

    vector<Game *> games;
    int active_idx;

    Games()
    {
      games.push_back(new Berlin(world));
      games.push_back(new FindTheLight(world));
      games.push_back(new TurnAllOn(world));
      games.push_back(new MoveAllBlocks(world));

      active_idx = 0;
    }

    Game &game()
    {
      return *games[active_idx];
    }

    void start()
    {
      foreach(i, games) {
        Game *g = *i;
        g->init();
      }

      game().play();
    }

    void next_game()
    {
      active_idx = (active_idx + 1) % games.size();
      game().play();
    }

    void prev_game()
    {
      active_idx = (active_idx + games.size() - 1) % games.size();
      game().play();
    }

    void run()
    {
      game().run();
      if (game().done) {
        game().level ++;
        game().init();
        game().play();
      }
    }

};

//char *audio_path = NULL;

SDL_Surface *screen = NULL;


typedef struct {
  int random_seed;
  bool start_blank;
} init_params_t;

init_params_t ip;

int main(int argc, char *argv[])
{
  bool usage = false;
  bool error = false;

  int c;

  //const char *audio_out_path = NULL;

  ip.random_seed = time(NULL);

  while (1) {
    c = getopt(argc, argv, "hf:g:r:"); //A:");
    if (c == -1)
      break;
   
    switch (c) {
      case 'g':
        {
          char arg[strlen(optarg) + 1];
          strcpy(arg, optarg);
          char *ch = arg;
          while ((*ch) && ((*ch) != 'x')) ch ++;
          if ((*ch) == 'x') {
            *ch = 0;
            ch ++;
            W = atoi(arg);
            H = atoi(ch);

          }
          else {
            fprintf(stderr, "Invalid -g argument: '%s'\n", optarg);
            exit(-1);
          }
        }
        break;

      case 'f':
        want_fps = atof(optarg);
        break;

      case 'r':
        ip.random_seed = atoi(optarg);
        break;

      //case 'A':
      //  audio_out_path = optarg;
      //  break;

      case '?':
        error = true;
      case 'h':
        usage = true;
        break;

    }
  }

  if (usage) {
    if (error)
      printf("\n");
    printf(
"scop3 v0.1\n"
"(c) 2014 Neels Hofmeyr <neels@hofmeyr.de>\n"
"Published under the GNU General Public License v3.\n\n"
"Scop3 produces a mesmerizing animation controlled by any game controller.\n"
"\n"
"Usage example:\n"
"  scop3 -g 320x200 -f 25\n"
"\n"
"Options:\n"
"\n"
"  -g WxH   Set window width and height in number of pixels.\n"
"           Default is '-g %dx%d'.\n"
"  -f fps   Set desired framerate to <fps> frames per second. The framerate\n"
"           may slew if your system cannot calculate fast enough.\n"
"           If zero, run as fast as possible. Default is %.1f.\n"
"  -r seed  Supply a random seed to start off with.\n"
//"  -A file  Write raw audio data to file.\n"
, W, H, want_fps
);
    if (error)
      return 1;
    return 0;
  }


	/*
  if (audio_out_path) {
#if 0
    if (access(out_stream_path, F_OK) == 0) {
      fprintf(stderr, "file exists, will not overwrite: %s\n", out_stream_path);
      exit(1);
    }
#endif
    Audio.write_to(audio_out_path);
  }
	*/

  const int maxpixels = 1e4;

#if DRAW_SPHERES
  int _argc = 1;
   glutInit(&_argc, argv);
#endif

  if ((W < 3) || (W > maxpixels) || (H < 3) || (H > maxpixels)) {
    fprintf(stderr, "width and/or height out of bounds: %dx%d\n", W, H);
    exit(1);
  }

  SDL_Event event;

  SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_JOYSTICK);

  const int n_joysticks = SDL_NumJoysticks();
  SDL_Joystick **joysticks = NULL;

  if (n_joysticks) {
    SDL_JoystickEventState(SDL_ENABLE);

    joysticks = (SDL_Joystick**)malloc(sizeof(SDL_Joystick*) * n_joysticks);
    
    int i;
    for (i = 0; i < n_joysticks; i++)
    {
      printf("%2d: '%s'\n", i, SDL_JoystickName(i));

      SDL_Joystick *j = SDL_JoystickOpen(i);
      printf("    %d buttons  %d axes  %d balls %d hats\n",
             SDL_JoystickNumButtons(j),
             SDL_JoystickNumAxes(j),
             SDL_JoystickNumBalls(j),
             SDL_JoystickNumHats(j)
             );
      joysticks[i] = j;
    }
  }

  atexit(SDL_Quit);
  SDL_WM_SetCaption("fly", NULL);
  screen = SDL_SetVideoMode(W,H, 32, SDL_OPENGL);// | SDL_FULLSCREEN);
  SDL_ShowCursor(SDL_DISABLE);

  printf("seed %d\n", (int)ip.random_seed);
  srandom(ip.random_seed);

	int maxTextureSize;
	glGetIntegerv(GL_MAX_TEXTURE_SIZE, &maxTextureSize);
	printf("max texture size: %d\n", maxTextureSize);

	Textures _textures(20000);
	gl_textures = &_textures;

  //Audio.start();
  //Audio.play(new Sine(140, 0.01));

  Games games;
  games.start();

  games.game().draw();

  Uint32 last_time = SDL_GetTicks();
  Uint32 current_time,elapsed_time;
  Uint32 start_time;

  float want_frame_period = (want_fps > .1? 1000. / want_fps : 0);
  float last_ticks = (float)SDL_GetTicks() - want_frame_period;

  while (running)
  {
    games.run();
    frames_rendered ++;

    {
      static int last_ticks2 = 0;

      int t = SDL_GetTicks();
      int elapsed = t - last_ticks2;

      last_ticks2 = t;

      avg_frame_period -= avg_frame_period >>AVG_SHIFTING;
      avg_frame_period += elapsed;

      dt = 1.e-3 * elapsed;
    }

    while (running) {
      SDL_Event event;
      while (running && SDL_PollEvent(&event)) 
      {
        switch(event.type)
        {

        case SDL_JOYHATMOTION:
          printf("%d %d\n", event.jhat.hat, event.jhat.value);

          switch (event.jhat.value) {
          case 1:
            games.game().level ++;
            games.game().init();
            games.game().play();
            break;
          case 4:
            games.game().level --;
            games.game().init();
            games.game().play();
            break;
          }
          break;

        case SDL_JOYAXISMOTION:
          games.game().on_joy_axis(event.jaxis.axis,
                                   ((double)event.jaxis.value) / 32768.);
          break;


        case SDL_JOYBUTTONDOWN:
          if (event.jbutton.button == 4) {
            games.prev_game();
            break;
          }
          if (event.jbutton.button == 5) {
            games.next_game();
          }
          if (event.jbutton.button == 7) {
            games.game().init();
            games.game().play();
          }
        case SDL_JOYBUTTONUP:
          games.game().on_joy_button(event.jbutton.button,
                                     event.type == SDL_JOYBUTTONDOWN);
          break;

        case SDL_JOYBALLMOTION:  /* Handle Joyball Motion */
          printf("%2d: ball %d += %d, %d\n",
                 event.jball.which, event.jball.ball,
                 event.jball.xrel, event.jball.yrel);
          break;

        case SDL_QUIT:
          running = false;
          break;

        case SDL_KEYDOWN:
          {
            int c = event.key.keysym.sym;

            switch(c) {
            case SDLK_ESCAPE:
              printf("Escape key. Stop.\n");
              running = false;
              break;

            case 13:
              SDL_WM_ToggleFullScreen(screen);
              break;

            case SDLK_TAB:
              games.next_game();
              break;

            default:
              games.game().on_key(c, true);
              break;
            }
          }
          break;

        case SDL_KEYUP:
          {
            int c = event.key.keysym.sym;
            games.game().on_key(c, false);
          }
          break;
        }
      } // while sdl poll event

      if (want_frame_period) {
        int elapsed = SDL_GetTicks() - last_ticks;
        if (elapsed >= want_frame_period) {
          last_ticks += want_frame_period * (int)(elapsed / want_frame_period);
          break;
        }
        SDL_Delay(min((int)want_frame_period - elapsed, 5));
      }

    } // while running, for event polling / idle waiting

  } // while running

  //Audio.stop();

  running = false;

  printf("\n");
  printf("%d frames rendered\n", frames_rendered);

  return 0;
}

