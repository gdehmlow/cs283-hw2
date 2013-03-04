/*
    Primitive.cpp
*/

#include <iostream>
#include "Primitive.h"

Primitive::Primitive(Shape* shape)
{
    this->shape = shape;
    Material* tempMaterial = new Material();
    this->material = tempMaterial;
}

int Primitive::intersectionPoint(Ray* ray, Intersection* intersection)
{
    return shape->intersectionPoint(ray, intersection);
}

bool Primitive::doesIntersect(Ray* ray, float tmax)
{
    return shape->doesIntersect(ray, tmax);
}

void Primitive::setAmbient(glm::vec3& ambient)
{
    this->ambient = ambient;
}

void Primitive::setMaterial(Material* material)
{
    this->material->diffuse = glm::vec3(material->diffuse);
    this->material->specular = glm::vec3(material->specular);
    this->material->emission = glm::vec3(material->emission);
    this->material->shininess = material->shininess;
    this->material->alpha = material->alpha;
    this->material->rindex = material->rindex;
}

void Primitive::setTransformation(glm::mat4& transformation)
{
    this->transformation = transformation;
    this->inverseTransformation = glm::inverse(transformation);
    glm::mat3 tempTransformation;
    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
            tempTransformation[i][k] = transformation[i][k];
            this->inverseTranspose[i][k] = this->inverseTransformation[k][i]; // time to go kill myself ._.
        }
    }
    this->shape->createAABB(tempTransformation);
}

AABB* Primitive::getAABB() 
{
    return this->shape->getAABB();
}
