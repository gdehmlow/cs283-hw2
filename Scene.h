/*
    Scene.h

    Reads a scene file, creates primitives, sets up the camera and screen, and
    renders the scene.
*/

#ifndef I_SCENE
#define I_SCENE

#include <sstream>
#include <stack>
#include <string>
#include <vector>

#include <glm/glm.hpp>

#include "Camera.h"
#include "Screen.h"
#include "Raytracer.h"
#include "Primitive.h"
#include "Light.h"
#include "AABB.h"

class Scene {
    public:
        friend class Raytracer;
        Scene(const char* filename);
        int renderImage();

    private:
        void raytrace();
        bool readvals(std::stringstream &s, const int numvals, float* values);
        void parseFile(const char* filename);
        void rightMultiply(const glm::mat4& M, std::stack<glm::mat4>& transfstack);
        void createGrid();

        Camera camera;
        Screen screen;

        glm::vec3 gridStart;
        glm::vec3 gridEnd;
        int gridSize;
        // 3d array of int arrays used to store prims at each voxel in the grid
        std::vector<std::vector<int> > grid;
        AABB* mommaBox;

        int gi;
        int gidepth;
        std::vector<Primitive> primitiveList;
        std::vector<Light> lightList;
        glm::vec3 attenuation;
        glm::vec3 ambient;
        std::string outputFilename;
        std::stack <glm::mat4> transfstack; 
};

#endif
