/*
    Light.h
*/

#ifndef I_LIGHT
#define I_LIGHT

#include <glm/glm.hpp>

enum LightType { POINT, DIRECTIONAL, AREA };
enum AreaType { RECT };

// Raytracing //

typedef struct _Light {
    glm::dvec4 posdir;
    glm::dvec3 color;
    LightType type;

    AreaType areaType;
    int numSamples; // Per square side! So numSamples = 3 -> 9 samples taken. 
    glm::dvec4 upStep;
    glm::dvec4 rightStep;
} Light;


// Pathtracing //

class AreaLight {
    public:
        virtual void getSample(const glm::dvec3& position, const glm::dvec3& normal, glm::dvec3& lightIntensity, 
                                     glm::dvec3& incidentRay) = 0;
};

/*
    Quad light is defined py a position, an upVec, and a rightVec.

    ^ = = = = = = o 
    | <-upVec     |
    |             |
    *------------>o
    ^        ^
    position rightVec
*/

class QuadLight : public AreaLight {
    public:
        QuadLight(glm::dvec3 position, glm::dvec3 upVec, glm::dvec3 rightVec, glm::dvec3 color);
        ~QuadLight();
        void getSample(const glm::dvec3& position, const glm::dvec3& normal, glm::dvec3& lightIntensity, 
                             glm::dvec3& lightRay);
    private:
        glm::dvec3 position;
        glm::dvec3 upVec;
        glm::dvec3 rightVec;
        glm::dvec3 color;
        double area;
        glm::dvec3 normal;
};

/*
    Circle light defined by a position, a normal, and a radius
*/
/*
class CircleLight : public Arealight {
    public:
        CircleLight(glm::vec4 position, glm::vec3 normal, float radius, glm::vec3 color);
        ~CircleLight();
        void getSample(const glm::vec3& position, const glm::vec3& normal, 
                             glm::vec3& lightIntensity, glm::vec3& incidentRay);
    private:
        glm::vec4 position;
        glm::vec3 normal;
        float radius;
        glm::vec3 color;
        float area;
};
*/

#endif