# HW5 - Raytracer makefile

CC = g++
ifeq ($(shell sw_vers 2>/dev/null | grep Mac | awk '{ print $$2}'),Mac)
CFLAGS = -g -DOSX -O3 -fopenmp
INCFLAGS = -I./include/ -I./glm
LDFLAGS = -lm -lstdc++ -L./osxlib/ -lfreeimage
else
CFLAGS = -g -DGL_GLEXT_PROTOTYPES 
INCFLAGS = -I./glm-0.9.2.7 -I./include/ -I/usr/X11R6/include -I/sw/include \
		-I/usr/sww/include -I/usr/sww/pkg/Mesa/include
LDFLAGS = -L/usr/X11R6/lib -L/sw/lib -L/usr/sww/lib \
		-L/usr/sww/bin -L/usr/sww/pkg/Mesa/lib -lX11 -lfreeimage
endif

RM = /bin/rm -f 
all: raytracer
raytracer: main.o Scene.o Camera.o Screen.o Raytracer.o Triangle.o Sphere.o Primitive.o Transform.o AABB.o
	$(CC) $(CFLAGS) -o raytracer main.o Scene.o Camera.o Screen.o Raytracer.o Triangle.o Sphere.o Primitive.o Transform.o AABB.o $(INCFLAGS) $(LDFLAGS) 
main.o: main.cpp Scene.cpp Scene.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c main.cpp
Scene.o: Scene.cpp Scene.h Camera.h Screen.h Triangle.h Sphere.h Primitive.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c Scene.cpp
Camera.o: Camera.cpp Camera.h Ray.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c Camera.cpp
Screen.o: Screen.cpp Screen.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c Screen.cpp
Raytracer.o: Raytracer.cpp Raytracer.h Ray.h Scene.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c Raytracer.cpp
Primitive.o: Primitive.cpp Primitive.h Shape.h Material.h Intersection.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c Primitive.cpp
Triangle.o: Triangle.cpp Triangle.h Shape.h Intersection.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c Triangle.cpp
Sphere.o: Sphere.cpp Sphere.h Shape.h Intersection.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c Sphere.cpp
Transform.o: Transform.cpp Transform.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c Transform.cpp
AABB.o: AABB.cpp AABB.h
	$(CC) $(CFLAGS) $(INCFLAGS) -c AABB.cpp
clean: 
	$(RM) *.o transforms *.png
