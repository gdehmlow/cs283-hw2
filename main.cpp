/*
	HW5 - Raytracer
	Greg Dehmlow
	cs184-ao
*/

#include <iostream>
#include "Scene.h"

int main(int argc, char* argv[])
{
    int val;
    if (argc == 2) {
        Scene scene(argv[1]);
        val = scene.renderImage();
    } else {
        std::cerr << "Wrong number of arguments.\n";
        val = -1;
    }
    return val;
}
