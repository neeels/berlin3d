#pragma once


#include <vector>
using namespace std;

#include <GL/gl.h>
#include <SDL/SDL.h>

#include <foreach.h>

class Textures;
class Texture;

struct TextureSlot {
	GLuint id;
	Texture *taken;
};

class Texture {
  public:
    const char *path;
		bool want_loaded;

    int native_w, native_h;

    Texture() :
			path(NULL),
			want_loaded(false),
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

		void want(const char *path)
		{
			this->path = path;
			want_loaded = true;
		}

		void unload()
		{
			if (loaded()) {
				_loaded->taken = NULL;
				_loaded = NULL;
			}
			native_w = native_h = -1;
		}

		bool loaded() const
		{
			return _loaded && (_loaded->taken == this);
		}

	private:
    TextureSlot *_loaded;
    bool load(bool use_mip_map=false);

};


class Textures {
  public:

    /* Before setting up textures, OpenGL needs to know how many there will
     * be... */
    Textures(int n);

    vector<TextureSlot> slots;

		TextureSlot *unused_slot();

		void some_unloaded()
		{
			all_taken = 0;
		}

  private:
		unsigned int tail;
		int all_taken;
    Textures() {}
};

