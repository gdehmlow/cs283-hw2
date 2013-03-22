/*
    Primitive.cpp
*/

#include <iostream>
#include "Primitive.h"

Primitive::Primitive(Shape* shape)
{
    this->shape = shape;
    this->material = new Material();
}

Primitive::~Primitive()
{
}

int Primitive::intersectionPoint(Ray* ray, Intersection* intersection)
{
    return shape->intersectionPoint(ray, intersection);
}

int Primitive::getIntersectionPoint(const Ray& ray, Intersection& intersection, const double dt)
{
    return shape->getIntersectionPoint(ray, intersection, dt);
}

bool Primitive::doesRayIntersect(const Ray& ray, const float tmax, const double dt)
{
    // This seems hacky.
    if (this->material->type == EMISSIVE) {
        return false;
    } else {
        return shape->doesRayIntersect(ray, tmax, dt);
    }
}

bool Primitive::doesIntersect(Ray* ray, float tmax)
{
    // So does this XD
    if (this->material->type == EMISSIVE) {
        return false;
    } else {
        return shape->doesIntersect(ray, tmax);
    }
}

void Primitive::setAmbient(const glm::dvec3& ambient)
{
    this->ambient = ambient;
}

void Primitive::setMaterial(const Material& material)
{
    this->material->diffuse     = glm::vec3(material.diffuse);
    this->material->specular    = glm::vec3(material.specular);
    this->material->emission    = glm::vec3(material.emission);
    this->material->shininess   = material.shininess;
    this->material->alpha       = material.alpha;
    this->material->rindex      = material.rindex;
    this->material->type        = material.type;
}

void Primitive::setTransformation(glm::dmat4& transformation)
{
    this->transformation = transformation;
    this->inverseTransformation = glm::inverse(transformation);
    glm::mat3 tempTransformation;
    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
            tempTransformation[i][k] = transformation[i][k];
            this->inverseTranspose[i][k] = this->inverseTransformation[k][i];
        }
    }
}
