
run: fly
  #./run_fly.sh
	./fly -g 1920x900

F = -g -std=c++11
INC = ../include

berlin_objs = ../buildings/cpp-sections/*.o

fly: fly.cpp textures.o draw.h $(INC)/foreach.h $(INC)/pt.h $(berlin_objs)
	g++ $F -g -o fly fly.cpp textures.o $(berlin_objs) -lGL -lGLU -lSDL -lSDL_image -lm 
  
%.o: %.cpp
	g++ $F -c $< 

# vim: noexpandtab

