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
    glm::vec4 posdir;
    glm::vec3 color;
    LightType type;

    AreaType areaType;
    int numSamples; // Per square side! So numSamples = 3 -> 9 samples taken. 
    glm::vec4 upStep;
    glm::vec4 rightStep;
} Light;


// Pathtracing //

class AreaLight {
    public:
        virtual void getSample(const glm::vec3& position, const glm::vec3& normal, 
                                     glm::vec3& lightIntensity, glm::vec3& incidentRay) = 0;

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
        QuadLight(glm::vec3 position, glm::vec3 upVec, glm::vec3 rightVec, glm::vec3 color);
        ~QuadLight();
        void getSample(const glm::vec3& position, const glm::vec3& normal, 
                             glm::vec3& lightIntensity, glm::vec3& incidentRay);
    private:
        glm::vec3 position;
        glm::vec3 upVec;
        glm::vec3 rightVec;
        glm::vec3 color;
        float area;
        glm::vec3 normal;
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