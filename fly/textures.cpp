#include "textures.h"
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <GL/glu.h>

#include <cstring>
#include <cstdlib>

#define dbg(args...) {printf(args); fflush(stdout);}

struct Preload {
  Texture *for_texture;
  SDL_Surface *surface;
};

void Texture::commit(Textures &textures)
{
  bool l = loaded_or_loading();
  if (_want_loaded && ! l) {
    try_preload(textures);
  }
  else
    if ((! _want_loaded) && l)
      unload();
}

void Texture::try_preload(Textures &textures)
{
  if (!path) {
    dbg("!path\n");
    return;
  }
  if (loaded_or_loading()) {
    return;
  }
  _loaded = textures.unused_slot();
  if (! _loaded) {
    return;
  }
  _loaded->taken = this;
  if (! preload(textures, false)) {
    _want_loaded = false;
    unload();
  }
}


SDL_Surface * flip_surface(SDL_Surface * surface);

bool Texture::preload(Textures &textures, bool use_mip_map)
{
  GLuint id = _loaded->id;
  SDL_Surface * picture_surface = NULL;
  picture_surface = IMG_Load(path);
  if (picture_surface == NULL) {
    printf("Failed to load %s\n", path);
    return false;
  }

  Uint32 rmask, gmask, bmask, amask;
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
  rmask = 0xff000000;
  gmask = 0x00ff0000;
  bmask = 0x0000ff00;
  amask = 0x000000ff;
#else
  rmask = 0x000000ff;
  gmask = 0x0000ff00;
  bmask = 0x00ff0000;
  amask = 0xff000000;
#endif

  SDL_PixelFormat format = *(picture_surface->format);
  format.BitsPerPixel = 32;
  format.BytesPerPixel = 4;
  format.Rmask = rmask;
  format.Gmask = gmask;
  format.Bmask = bmask;
  format.Amask = amask;

  SDL_Surface *gl_surface;
  gl_surface = SDL_ConvertSurface(picture_surface, &format, SDL_SWSURFACE);

  SDL_Surface *gl_flipped_surface;
  gl_flipped_surface = flip_surface(gl_surface);

  SDL_FreeSurface(gl_surface);
  SDL_FreeSurface(picture_surface);

  Preload *p = new Preload();

  p->for_texture = this;
  p->surface = gl_flipped_surface;

  SDL_SemWait(textures.load_mutex);
  textures.pending_loads.push_back(p);
  SDL_SemPost(textures.load_mutex);
  textures.bump();

  return true;
}

void Texture::load(Preload *p)
{
  SDL_Surface *gl_flipped_surface = p->surface;

  if ((!_loaded) || (_loaded->taken != this)) {
    SDL_FreeSurface(gl_flipped_surface);
    return;
  }

  this->native_w = gl_flipped_surface->w;
  this->native_h = gl_flipped_surface->h;

  glBindTexture(GL_TEXTURE_2D, _loaded->id);

#if 0
  if (use_mip_map)
  {
    gluBuild2DMipmaps(GL_TEXTURE_2D, 4, gl_flipped_surface->w,
                      gl_flipped_surface->h, GL_RGBA,GL_UNSIGNED_BYTE,
                      gl_flipped_surface->pixels);

    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,
                    GL_LINEAR_MIPMAP_LINEAR);
  }
  else
#endif
  {
    glTexImage2D(GL_TEXTURE_2D, 0, 4, gl_flipped_surface->w,
                 gl_flipped_surface->h, 0, GL_RGBA,GL_UNSIGNED_BYTE,
                 gl_flipped_surface->pixels);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  }
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);

  _loaded->loaded = true;

  SDL_FreeSurface(gl_flipped_surface);
}

void Texture::want_loaded(Textures &textures)
{
  if (! _want_loaded) {
    _want_loaded = true;
  int tt = SDL_GetTicks();
    SDL_SemWait(textures.preload_mutex);
  int t2 = SDL_GetTicks() - tt;
  if (t2 > 0)
    dbg("Texture::want_loaded waited %d\n", t2);

		textures.pending_preloads.push_back(this);
    SDL_SemPost(textures.preload_mutex);
    textures.bump();
  }
}

void Texture::want_unloaded(Textures &textures)
{
  if (_want_loaded) {
    _want_loaded = false;
    textures.bump();
  }
}

Textures::Textures(int n) :
  all_taken(0)
{
  bumper = SDL_CreateSemaphore(0);
  preload_mutex = SDL_CreateSemaphore(1);
  load_mutex = SDL_CreateSemaphore(1);

  printf("Initializing %d texture slots\n", n);
  slots.resize(n);
  GLuint id[n];

  glGenTextures(n, id);

  for (int i = 0; i < n; i++) {
    slots[i].id = id[i];
  }

  tail = 0;
}

TextureSlot *Textures::unused_slot()
{
  if (tail >= slots.size())
    tail = 0;

  if (slots[tail].taken) {
    if (! all_taken) {
      for (int i = 0; i < slots.size(); i++) {
        if (! slots[i].taken) {
          tail = i;
          break;
        }
      }
    }
    if (slots[tail].taken) {
      all_taken ++;
      printf("Out of texture slots! (%d + %d)\n", slots.size(), all_taken);
      return NULL;
    }
  }

  slots[tail].loaded = false;
  return &slots[tail ++];
}


void Textures::textures_thread()
{
  SDL_SemPost(bumper);

  running = true;
  while (running) {
    SDL_SemWait(bumper);
    while (SDL_SemTryWait(bumper) == 0);

    SDL_SemWait(preload_mutex);

    foreach (ti, pending_preloads) {
      Texture *t = *ti;
      SDL_SemPost(preload_mutex);
      if (t) {
        t->commit(*this);
      }
      SDL_SemWait(preload_mutex);
      ti = pending_preloads.erase(ti);
    }

    SDL_SemPost(preload_mutex);
  }
}

int Textures::do_pending_loads()
{
  if (! pending_loads.size())
    return 0;
  int tt = SDL_GetTicks();
  SDL_SemWait(load_mutex);
  int t2 = SDL_GetTicks() - tt;
  if (t2 > 0)
    dbg("do_pending_loads waited %d\n", t2);
	int count = 0;
  foreach (pi, pending_loads) {
    Preload *p = *pi;

    p->for_texture->load(p);
    pi = pending_loads.erase(pi);

    delete p;
		count ++;
  }
  SDL_SemPost(load_mutex);
	return count;
}


int takeScreenshot(const char * filename)
{
    GLint viewport[4];
    Uint32 rmask, gmask, bmask, amask;
    SDL_Surface * picture, * finalpicture;

    glGetIntegerv(GL_VIEWPORT, viewport);

#if SDL_BYTEORDER == SDL_BIG_ENDIAN

    rmask = 0xff000000;
    gmask = 0x00ff0000;
    bmask = 0x0000ff00;
    amask = 0x000000ff;
#else

    rmask = 0x000000ff;
    gmask = 0x0000ff00;
    bmask = 0x00ff0000;
    amask = 0xff000000;
#endif

    picture = SDL_CreateRGBSurface(SDL_SWSURFACE,viewport[2],viewport[3], 32,
                                   rmask, gmask, bmask, amask);
    SDL_LockSurface(picture);
    glReadPixels(viewport[0],viewport[1],viewport[2],viewport[3],GL_RGBA,
                 GL_UNSIGNED_BYTE,picture->pixels);
    SDL_UnlockSurface(picture);

    finalpicture = flip_surface(picture);

    if (SDL_SaveBMP(finalpicture, filename))
    {
        return -1;
    }
    SDL_FreeSurface(finalpicture);
    SDL_FreeSurface(picture);

    return 0;
}

SDL_Surface * flip_surface(SDL_Surface * surface)
{
    int current_line,pitch;
    SDL_Surface * flipped_surface = SDL_CreateRGBSurface(SDL_SWSURFACE,
                                   surface->w,surface->h,
                                   surface->format->BitsPerPixel,
                                   surface->format->Rmask,
                                   surface->format->Gmask,
                                   surface->format->Bmask,
                                   surface->format->Amask);



    SDL_LockSurface(surface);
    SDL_LockSurface(flipped_surface);

    pitch = surface->pitch;
    for (current_line = 0; current_line < surface->h; current_line ++)
    {
        memcpy(&((unsigned char* )flipped_surface->pixels)[current_line*pitch],
               &((unsigned char* )surface->pixels)[(surface->h - 1
                                                    - current_line)*pitch],
               pitch);
    }

    SDL_UnlockSurface(flipped_surface);
    SDL_UnlockSurface(surface);
    return flipped_surface;
}

void drawAxis(double scale)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glPushMatrix();
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(2);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glScaled(scale,scale,scale);
    glBegin(GL_LINES);
    glColor3ub(255,0,0);
    glVertex3i(0,0,0);
    glVertex3i(1,0,0);
    glColor3ub(0,255,0);
    glVertex3i(0,0,0);
    glVertex3i(0,1,0);
    glColor3ub(0,0,255);
    glVertex3i(0,0,0);
    glVertex3i(0,0,1);
    glEnd();
    glPopMatrix();
    glPopAttrib();
}

