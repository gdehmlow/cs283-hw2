/*
    Light.cpp

    Implements area light sources:
        QuadLight
*/

#include <algorithm>
#include "Light.h"

typedef glm::vec3 vec3;
typedef glm::vec4 vec4;

QuadLight::QuadLight(vec3 position, vec3 upVec, vec3 rightVec, vec3 color)
{
    this->position  = position;
    this->upVec     = upVec;
    this->rightVec  = rightVec;
    this->color     = color;

    this->normal = glm::normalize(glm::cross(rightVec, upVec));
}

QuadLight::~QuadLight()
{
}

/*
    Takes in a point on some surface and returns (via lightIntensity) the incoming radiance at that point from a random 
    point on the QuadLight.

    Used when the DirectLighting setting is on, otherwise lighting occurs only through indirect bounces that hit the
    area light geometry.

    lightIntensity = L_o(x'_i, w'_i)*G(x,x'_i)*V(x,x'_i)*A, found in the lecture 10 slides. 

    Note: the visibility term is calculated in the DirectLighting routine
*/
void QuadLight::getSample(const vec3& position, const vec3& normal, 
                                vec3& lightIntensity, vec3& incidentRay)
{
    // Generates random point on the light
    float u = rand() / double(RAND_MAX);
    float v = rand() / double(RAND_MAX);
    vec3 positionOnLight = this->position + this->upVec*u + this->rightVec*v;

    // Creates incidentRay for point on surface
    incidentRay = positionOnLight - position;
    float incidentLength = glm::length(incidentRay);

    // Calculates G(x,x'_i) = cos(θ)cos(θ'_i)/|x-x'_i|^2
    float g = std::max(glm::dot(incidentRay,normal),0.0f) * std::max(glm::dot(-incidentRay,this->normal),0.0f) /
              std::pow(incidentLength, 4.0f);

    lightIntensity = this->color * g * area;
}
