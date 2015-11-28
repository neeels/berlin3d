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
  bool loaded;
  Texture *taken;
};

struct Preload;

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

    void try_preload(Textures &textures);

    void try_preload(Textures &textures, const char *path)
    {
      this->path = path;
			_want_loaded = true;
      try_preload(textures);
    }

    void unload()
    {
      if (loaded_or_loading()) {
        _loaded->taken = NULL;
        _loaded = NULL;
      }
      native_w = native_h = -1;
      _want_loaded = false;
    }


    bool loaded_or_loading() const
    {
      return _loaded && (_loaded->taken == this);
    }

    bool loaded() const
    {
      return loaded_or_loading() && _loaded->loaded;
    }

    void want_loaded(Textures &textures);
    void want_unloaded(Textures &textures);

    bool preload(Textures &textures, bool use_mip_map=false);

    void commit(Textures &textures);

    bool _want_loaded;

    void load(Preload *p);
};


class Textures {
  public:
    SDL_sem *preload_mutex;
    SDL_sem *load_mutex;
    bool running;

    /* Before setting up textures, OpenGL needs to know how many there will
     * be... */
    Textures(int n);

    vector<TextureSlot> slots;
		list<Texture*> pending_preloads;
		list<Preload*> pending_loads;

    TextureSlot *unused_slot();

    void some_unloaded()
    {
      all_taken = 0;
    }

    void bump() {
      SDL_SemPost(bumper);
    }
    SDL_sem *bumper;

    void textures_thread();
    int do_pending_loads();

  private:
    unsigned int tail;
    int all_taken;
    Textures() {}
};

