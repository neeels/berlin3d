#pragma once


#include <vector>
using namespace std;

#include <GL/gl.h>
#include <SDL/SDL.h>

#include <foreach.h>

class Texture {
  public:
    bool loaded;
    GLuint id;
    const char *path;
    int native_w, native_h;

    Texture() :
      loaded(false),
      id(0)
    {}

    void begin() {
      glEnable(GL_TEXTURE_2D);
      glBindTexture(GL_TEXTURE_2D, id);
    }

    void end() {
      glDisable(GL_TEXTURE_2D);
    }

    bool operator == (const char *path)
    {
      return (this->path == path)
        || ((this->path != NULL)
            && (path != NULL)
            && (strcmp(this->path, path) == 0));
    }

    /* rather invoke via Textures::load(). */
    bool load(GLuint id, const char *path, bool use_mip_map=false);

    void unload()
    {
      loaded = false;
      path = NULL;
      native_w = native_h = -1;
    }
};



class Textures {
  private:
    Textures() {}

  public:

    /* Before setting up textures, OpenGL needs to know how many there will
     * be... */
    Textures(int n);
    Texture *load(const char *path);

    vector<Texture> textures;
};

