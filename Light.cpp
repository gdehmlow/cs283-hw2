/*
    Light.cpp

    Implements area light sources:
        QuadLight
*/

#include <iostream>
#include <algorithm>
#include "Light.h"

typedef glm::dvec3 dvec3;
typedef glm::dvec4 dvec4;

QuadLight::QuadLight(dvec3 position, dvec3 upVec, dvec3 rightVec, dvec3 color)
{
    this->position  = position;
    this->upVec     = upVec;
    this->rightVec  = rightVec;
    this->color     = color;

    this->normal = glm::normalize(glm::cross(upVec, rightVec));
    this->area = glm::length(upVec) * glm::length(rightVec);
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
void QuadLight::getSample(const dvec3& position, const dvec3& normal, 
                                dvec3& lightIntensity, dvec3& incidentRay)
{
    double u = rand() / double(RAND_MAX);
    double v = rand() / double(RAND_MAX);
    dvec3 positionOnLight = this->position + this->upVec*u + this->rightVec*v;
    incidentRay = positionOnLight - position;
    double incidentLength = glm::length(incidentRay);

    // Calculates G(x,x'_i) = cos(θ)cos(θ'_i)/|x-x'_i|^2
    double g = std::max(glm::dot(incidentRay,normal) / incidentLength,0.0) * 
               std::max(glm::dot(-incidentRay,this->normal) / incidentLength,0.0) /
               std::pow(incidentLength, 2.0);

    //std::cout << std::max(glm::dot(incidentRay,normal),0.0) / incidentLength << "\n";
    //std::cout << std::max(glm::dot(-incidentRay,this->normal),0.0) / incidentLength << "\n";
    //std::cout << "<" << this->color.x << ", " << this->color.y << ", " << this->color.z << ">, " << g << ", " << area << "\n";

    lightIntensity = this->color * g * area;
}
