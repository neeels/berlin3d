#pragma once


#include <vector>
using namespace std;

#include <GL/gl.h>
#include <SDL/SDL.h>

#include <foreach.h>
#include <list>

class Textures;
class Texture;

struct TextureSlot {
  GLuint id;
  Texture *taken;
};

class Texture {
  public:
    const char *path;
    TextureSlot *_loaded;

    int native_w, native_h;

    Texture() :
      path(NULL),
      _want_loaded(false),
      _loaded(NULL)
    {}

    void begin()
    {
      glEnable(GL_TEXTURE_2D);
      glBindTexture(GL_TEXTURE_2D, _loaded->id);
    }

    void end()
    {
      glDisable(GL_TEXTURE_2D);
    }

    bool operator == (const char *path)
    {
      return (this->path == path)
        || ((this->path != NULL)
            && (path != NULL)
            && (strcmp(this->path, path) == 0));
    }

    void try_load(Textures &textures);

    void try_load(Textures &textures, const char *path)
    {
      this->path = path;
			_want_loaded = true;
      try_load(textures);
    }

    void unload()
    {
      if (loaded()) {
        _loaded->taken = NULL;
        _loaded = NULL;
      }
      native_w = native_h = -1;
      _want_loaded = false;
    }

    bool loaded() const
    {
      return _loaded && (_loaded->taken == this);
    }

    void want_loaded(Textures &textures);
    void want_unloaded(Textures &textures);

    bool load(Textures &textures, bool use_mip_map=false);

    void commit(Textures &textures)
    {
			if (_want_loaded) {
				printf("commit %d %s\n", (int)_want_loaded, path);
				fflush(stdout);
			}
      bool l = loaded();
      if (_want_loaded && ! l)
        try_load(textures);
      else
      if ((! _want_loaded) && l)
        unload();
    }

    bool _want_loaded;
};


class Textures {
  public:
    SDL_sem *draw_mutex;

    /* Before setting up textures, OpenGL needs to know how many there will
     * be... */
    Textures(int n, SDL_sem *draw_mutex);

    vector<TextureSlot> slots;
		list<Texture*> pending;

    TextureSlot *unused_slot();

    void some_unloaded()
    {
      all_taken = 0;
    }

    void bump() {
			printf("bumpin\n");
      SDL_SemPost(bumper);
    }
    SDL_sem *bumper;

  private:
    unsigned int tail;
    int all_taken;
    Textures() {}
};

