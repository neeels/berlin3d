# berlin3d

This code is the result of me playing with the Berlin 3D data published at
http://www.businesslocationcenter.de/downloadportal

It may be full of unused code, suboptimal implementations and rough edges, but
it's fun to experiment with. Tested only on GNU/Linux so far.

The Berlin data is not included here. You need to download this separately, and
possibly select the suburbs you're interested in.

Basically, this would download, parse/build everything and run the program at
current resolution:

    cd buildings
    ./get-src-citygml.sh
    make
    cd ../satbild
    ./get-src-ecw.sh
    make
    cd ../fly
    make
    # ln -s ~/some.jpg ./backdrop.jpg
    ./run_fly.sh

Current keyboard controls:

    W,A,S,D: move forward, left, backward, right
    Q,E: move up, down
    arrows: look around
    1: start line
    2: line waypoint
    3: end line
    4: clear line
    G: toggle ground
    F: toggle buildings fly
    Esc: quit

Obviously, quite a lot can be improved here...
- drop building duplications
- combined make strucure
- docs
- vertex arrays
- load textures via pixel buffer objects
- use a proper 3D lib instead of own playful code
- portability

Patches welcome!

neels@hofmeyr.de

