
run: fly
  #./run_fly.sh
	./fly -g 1920x900

INC = ../include
INC_SECTIONS = ../buildings/cpp-sections

F = -g -std=c++11 -mcmodel=medium -I$(INC) -I$(INC_SECTIONS) 

berlin_objs = ../buildings/cpp-sections/*.o

fly: fly.cpp textures.o draw.h $(INC)/foreach.h $(INC)/pt.h $(berlin_objs)
	g++ $F -o fly fly.cpp textures.o $(berlin_objs) -lGL -lGLU -lSDL -lSDL_image -lm 
  
%.o: %.cpp
	g++ $F -c $< 

# vim: noexpandtab

