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
#include "Sampler.h"
#include "Raytracer.h"
#include "Primitive.h"
#include "Light.h"
#include "AABB.h"

enum TraceType { RAY, PATH };
enum SamplingType { UNIFORM, IMPORTANCE };

class Scene {
    public:
        friend class Raytracer;
        Scene(const char* filename);
        ~Scene();
        int renderImage();

    private:
        void raytrace();
        bool readvals(std::stringstream &s, const int numvals, float* values);
        void parseFile(const char* filename);
        void rightMultiply(const glm::mat4& M, std::stack<glm::mat4>& transfstack);
        void createGrid();

        Camera camera;
        Screen screen;

        std::vector<Primitive> primitiveList;
        std::vector<Light> lightList;
        std::vector<AreaLight*> areaLightList;
        glm::vec3 attenuation;
        glm::vec3 ambient;
        std::string outputFilename;
        std::stack<glm::mat4> transfstack; 

        // AABB related variables
        glm::vec3 gridStart;
        glm::vec3 gridEnd;
        int gridSize;
        std::vector<std::vector<int> > grid;
        AABB* mommaBox;

        // Rendering settings
        int maxDepth;
        int samplesPerPixel;
        bool directLighting;
        bool indirectLighting;
        SamplerType samplerType;    // How we shoot rays through pixels
        TraceType traceType;        // Raytracing or Pathtracing
        SamplingType samplingType;  // How we weight indirect lighting rays

};

#endif
