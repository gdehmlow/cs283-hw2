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
        bool readvals(std::stringstream &s, const int numvals, double* values);
        void parseFile(const char* filename);
        void rightMultiply(const glm::dmat4& M, std::stack<glm::dmat4>& transfstack);

        Camera camera;
        Screen screen;

        std::vector<Primitive> primitiveList;
        std::vector<Light> lightList;
        std::vector<AreaLight*> areaLightList;
        glm::dvec3 attenuation;
        glm::dvec3 ambient;
        std::string outputFilename;
        std::stack<glm::dmat4> transfstack; 

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
