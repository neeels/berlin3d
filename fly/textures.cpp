#include "textures.h"
#include <SDL/SDL.h>
#include <SDL/SDL_image.h>
#include <GL/glu.h>

#include <cstring>
#include <cstdlib>

#if 0
#ifndef GL_CLAMP_TO_EDGE
#define GL_CLAMP_TO_EDGE 0x812F
#endif
#endif

GLuint load_texture(const char * filename,bool useMipMap);
SDL_Surface * flip_surface(SDL_Surface * surface);

bool Texture::load(GLuint id, const char *path, bool use_mip_map)
{
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

  this->native_w = gl_flipped_surface->w;
  this->native_h = gl_flipped_surface->h;

  glBindTexture(GL_TEXTURE_2D, id);

  if (use_mip_map)
  {
    gluBuild2DMipmaps(GL_TEXTURE_2D, 4, gl_flipped_surface->w,
                      gl_flipped_surface->h, GL_RGBA,GL_UNSIGNED_BYTE,
                      gl_flipped_surface->pixels);

    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,
                    GL_LINEAR_MIPMAP_LINEAR);
  }
  else
  {
    glTexImage2D(GL_TEXTURE_2D, 0, 4, gl_flipped_surface->w,
                 gl_flipped_surface->h, 0, GL_RGBA,GL_UNSIGNED_BYTE,
                 gl_flipped_surface->pixels);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  }
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);

  SDL_FreeSurface(gl_flipped_surface);
  SDL_FreeSurface(gl_surface);
  SDL_FreeSurface(picture_surface);

  this->id = id;
  this->path = path;
  loaded = true;
}

Textures::Textures(int n)
{
	printf("Initializing %d textures\n", n);
  textures.resize(n);
  GLuint id[n];
  glGenTextures(n, id);

  for (int i = 0; i < n; i++) {
    textures[i].id = id[i];
		printf("texture id %d\n", (int)id[i]);
  }
}

Texture *Textures::load(const char *path)
{
	static int out_of_textures = 0;
  Texture *first_empty = NULL;

	if (!out_of_textures) { // HACK only when sure that no texture is used twice
  foreach(t, textures) {
    if (!(t->loaded)) {
      if (! first_empty)
        first_empty = &(*t);
      continue;
    }
    if ((*t) == path)
      return &(*t);
  }
	}

  if (! first_empty) {
		out_of_textures ++;
    printf("Out of textures! (%d + %d)\n", textures.size(), out_of_textures);
    return NULL;
  }

  Texture &t = *first_empty;
  t.load(t.id, path);

  if (t.loaded) {
    static int count = 0; // ick!
    count ++;
    printf("%p: Loaded texture nr %d (id=%d): %s\n", &t, count, (int)t.id, path);
  }

  return &t;
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
               &((unsigned char* )surface->pixels)[(surface->h - 1  -
                                                    current_line)*pitch],
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

