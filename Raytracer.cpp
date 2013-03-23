/*
    Raytracer.cpp
*/

#define GLM_SWIZZLE

#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <glm/gtx/component_wise.hpp>
#include "Raytracer.h"
#include "Transform.h"

typedef glm::dmat3 dmat3;
typedef glm::dmat4 dmat4;
typedef glm::dvec3 dvec3;
typedef glm::dvec4 dvec4;

Raytracer::Raytracer(Scene* scene)
{
    this->scene = scene;
    this->maxDepth = scene->maxDepth;    
}

Raytracer::~Raytracer()
{
}

///// PATHTRACER
dvec3 Raytracer::pathTraceRay(const Ray& ray, int depth, double weight, int bounce)
{
    // Grab nearest intersection
    int closestIntersectionIndex;
    Intersection closestIntersection; 
    bool haveAnIntersection = findClosestIntersection(ray, closestIntersectionIndex, closestIntersection);
    
    dvec3 color = dvec3(0.0);

    if (!haveAnIntersection) {
        return color;
    } else {

        // Grab the object that our intersection is on
        Primitive intersectedObject = scene->primitiveList[closestIntersectionIndex];
        
        // Get intersection from object to world space
        Intersection surfaceIntersection;
        surfaceIntersection.position = closestIntersection.position * intersectedObject.transformation;
        surfaceIntersection.normal   = glm::normalize(closestIntersection.normal * intersectedObject.inverseTranspose);

        color = intersectedObject.material->emission;

        ///// RUSSIAN ROULETTE
        // If the weight of the ray falls below 0.01, 50% chance of terminating ray -> double weight of survivors
        double terminatingP = 0.01;
        double vodka = 1.0;
        if (weight <= 0.01) {
            double meteoRandom = rand() / double(RAND_MAX);
            if (meteoRandom > terminatingP) {
                return color / terminatingP;
            } else {
                vodka = 1.0 / (1.0 - terminatingP);
            }
        }

        // Diffuse surfaces get processed normally
        if (intersectedObject.material->type == LAMBERTIAN) {
            if (scene->directLighting) {
                color += directLighting(ray, surfaceIntersection, intersectedObject) * weight * vodka;
            }

            if (scene->indirectLighting) {
                Ray randomRay; randomRay.t = ray.t;
                sampleCosineWeightedHemisphere(randomRay, surfaceIntersection); 
                //sampleUniformHemisphere(randomRay, surfaceIntersection); 
                dvec3 traceColor = pathTraceRay(randomRay, ++depth, weight * glm::dot(intersectedObject.material->diffuse, 
                                                                                      dvec3(0.333)), ++bounce);
                color += traceColor * intersectedObject.material->diffuse * weight * vodka;
            }
        }

        else if (intersectedObject.material->type == GLOSSY) {            
            if (scene->directLighting) {
                color += directLighting(ray, surfaceIntersection, intersectedObject) * weight * vodka;
            }

            if (scene->indirectLighting) {
                ////// IMPORTANCE SAMPLING
                // The probability of sampling diffuse interreflection is, how I chose, the sum of the diffuse terms over the
                // sum of the diffuse and the specular terms of the material reflectance.

                // dvec3 colorSpace = dvec3(0.2989, 0.5866, 0.1145);  Here are some other possibilities for determining the
                // double third = 1.0 / 3.0;                          probability for sampling the diffuse interreflection
                // dvec3 averager = dvec3(third, third, third);       versus the specular interreflection.

                double diffuseSum = glm::dot(intersectedObject.material->diffuse, dvec3(1.0));
                double specularSum = glm::dot(intersectedObject.material->specular, dvec3(1.0));

                double probDiffuse = diffuseSum / (diffuseSum + specularSum);
                double impRand = rand() / double(RAND_MAX);
                Ray randomRay; randomRay.t = ray.t;

                // Diffuse
                if (impRand < probDiffuse) {
                    sampleCosineWeightedHemisphere(randomRay, surfaceIntersection); 
                    //sampleUniformHemisphere(randomRay, surfaceIntersection); 
                    dvec3 traceColor = pathTraceRay(randomRay, ++depth, weight * glm::dot(intersectedObject.material->diffuse, 
                                                                                          dvec3(0.333)), ++bounce);
                    color += traceColor * intersectedObject.material->diffuse * weight * vodka / probDiffuse;
                }

                // Specular
                else {
                    sampleSpecularLobe(randomRay, surfaceIntersection, dvec3(glm::normalize(glm::reflect(ray.direction, surfaceIntersection.normal)).xyz),
                                       intersectedObject.material->shininess); 
                    //sampleUniformHemisphere(randomRay, surfaceIntersection); 
                    dvec3 traceColor = pathTraceRay(randomRay, ++depth, weight * glm::dot(intersectedObject.material->specular, 
                                                                                          dvec3(0.333)), ++bounce);
                    color += traceColor * intersectedObject.material->specular * weight * vodka / (1.0 - probDiffuse);
                }
            }
        }

        else if (intersectedObject.material->type == REFLECTIVE) {
            Ray reflectedRay;
            reflectedRay.direction = glm::reflect(ray.direction, surfaceIntersection.normal);
            reflectedRay.position  = surfaceIntersection.position + reflectedRay.direction * 0.001;
            reflectedRay.t = ray.t;
            color += pathTraceRay(reflectedRay, ++depth, weight * glm::dot(intersectedObject.material->specular, dvec3(0.333)), bounce) 
                     * intersectedObject.material->specular * weight * vodka;
        }

        else if (intersectedObject.material->type == EMISSIVE) {
            // Since we're sampling the direct lighting at each intersection when it's on, we don't want to factor this into
            // indirect lighting when the ray randomly hits the light source
            if (bounce > 0 && scene->directLighting) {
                return dvec3(0.0);
            }
        }

        return color;
    }
}

bool Raytracer::findClosestIntersection(const Ray& ray, int& closestIntersectionIndex, Intersection& closestIntersection)
{
    double minimumT = DBL_MAX;
    std::vector<Primitive>::iterator it;
    Ray tempRay;
    Intersection tempIntersect;
    int i = 0;
    bool didIntersect = false;

    // Used for flipping the normal for transmissive objects
    int tempFlip;
    int actualFlip;

    // Iterate over all the objects in the scene and find the nearest intersection
    for (it = scene->primitiveList.begin(); it < scene->primitiveList.end(); it++) {
        tempRay.position  = ray.position  * it->inverseTransformation;
        tempRay.direction = ray.direction * it->inverseTransformation;

        if (tempFlip = it->getIntersectionPoint(tempRay, tempIntersect, ray.t)) {
            if (tempIntersect.t < minimumT && tempIntersect.t >= 0.0) {
                actualFlip = tempFlip;
                minimumT = tempIntersect.t;
                closestIntersectionIndex = i;
                closestIntersection.position = tempIntersect.position;
                closestIntersection.normal = tempIntersect.normal;
                closestIntersection.t = tempIntersect.t;
                didIntersect = true;
            }
        }
        ++i;
    }

    // Flip transmissive norm? Nah no transmissive materials

    return didIntersect;
}

dvec3 Raytracer::directLighting(const Ray& ray, const Intersection& surfaceIntersection, 
                               const Primitive& intersectedObject)
{
    std::vector<Primitive>::iterator it;
    std::vector<AreaLight*>::iterator lit; 
    Ray shadowRay;
    Ray tempRay;
    double tmax;
    bool visible = true;

    dvec3 color = dvec3(0.0);
    dvec3 lightIntensity = dvec3(0.0);
    dvec3 rayToLight = dvec3(0.0);

    // Loop over all the area lights
    for (lit = scene->areaLightList.begin(); lit < scene->areaLightList.end(); lit++) {

        // Get sample from the area light source (puts result into rayToLight)
        (*lit)->getSample(surfaceIntersection.position.xyz, surfaceIntersection.normal.xyz, lightIntensity, rayToLight);

        tmax = glm::length(rayToLight); // So we don't count intersections behind the light
        shadowRay.direction = dvec4(glm::normalize(rayToLight),1.0);
        shadowRay.position  = surfaceIntersection.position + 0.001 * shadowRay.direction;

        // Loop through each object again to see if our light ray intersects with anything
        for (it = scene->primitiveList.begin(); it < scene->primitiveList.end(); it++) {
            // Transform shadowRay into object space
            tempRay.position = shadowRay.position * it->inverseTransformation;
            tempRay.direction = shadowRay.direction * it->inverseTransformation;

            if (it->doesRayIntersect(tempRay, tmax, ray.t)) {
                visible = false;
                break;
            }
        }

        // If the light is unobstructed from the intersection point, we're good to go with shading
        if (visible) {
            // Diffuse term
            color += intersectedObject.material->diffuse * lightIntensity / (double)M_PI;

            // Specular term (ignore if just diffuse to speed it up)
            if (intersectedObject.material->type == GLOSSY) {
                dvec3 half = glm::normalize(glm::normalize(rayToLight) + glm::normalize(ray.position.xyz - surfaceIntersection.position.xyz));
                double nDotH = std::pow(std::max(glm::dot(surfaceIntersection.normal.xyz, half), 0.0),
                                        intersectedObject.material->shininess);
                color += intersectedObject.material->specular * nDotH * lightIntensity;
            }            
            double attenuation;
            double r = tmax;
            attenuation = (1.0/(scene->attenuation[0] + scene->attenuation[1]*r + scene->attenuation[2]*r*r));
            color *= attenuation;
        }
    }

    return color;
}

void Raytracer::sampleUniformHemisphere(Ray& ray, const Intersection& surfaceIntersection) 
{
    double d = rand() / double(RAND_MAX);
    double e = rand() / double(RAND_MAX);

    double theta = 2.0 * M_PI * d;
    double phi = acos(2.0 * e - 1.0);
    dvec3 randomDirection = dvec3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));

    dvec3 normal = dvec3(surfaceIntersection.normal.xyz);
    double a = rand() / double(RAND_MAX);
    double b = rand() / double(RAND_MAX);
    double c = rand() / double(RAND_MAX);

    dvec3 ran = dvec3(a,b,c);
    dvec3 x = glm::normalize(glm::cross(ran, normal));
    dvec3 y = glm::normalize(glm::cross(normal, x));
    dmat3 rotationMatrix = dmat3(x, y, normal);
    randomDirection = rotationMatrix*randomDirection;

    ray.position = surfaceIntersection.position + dvec4(randomDirection,1.0) * 0.001;
    ray.direction = dvec4(randomDirection,1.0);

    return;
}

void Raytracer::sampleCosineWeightedHemisphere(Ray& ray, const Intersection& surfaceIntersection) 
{
    double z = rand() / double(RAND_MAX);
    double phi = 2.0 * M_PI * rand() / double(RAND_MAX);
    double theta = acos(sqrt(z));

    // This is sampled over the hemisphere around <0,0,1>
    dvec3 randomDirection = dvec3(sin(theta)*cos(phi), sin(theta)*sin(phi), z);

    // So we align it with the intersection normal
    dvec3 normal = dvec3(surfaceIntersection.normal.xyz);
    double a = rand() / double(RAND_MAX);
    double b = rand() / double(RAND_MAX);
    double c = rand() / double(RAND_MAX);

    dvec3 ran = dvec3(a,b,c);
    dvec3 x = glm::normalize(glm::cross(ran, normal));
    dvec3 y = glm::normalize(glm::cross(normal, x));
    dmat3 rotationMatrix = dmat3(x, y, normal);
    randomDirection = rotationMatrix*randomDirection;

    ray.position = surfaceIntersection.position + dvec4(randomDirection,1.0) * 0.001;
    ray.direction = dvec4(randomDirection,1.0);
    return;
}

// Practically the same thing as the cosine weighted hemisphere.. DRY takes a hit to the face
// Thx to Jensen et al. for being straight up bosses and mathing so I don't have to
void Raytracer::sampleSpecularLobe(Ray& ray, const Intersection& surfaceIntersection, 
                                   const dvec3& reflection, const double shininess) 
{
    double z = rand() / double(RAND_MAX);
    double alpha = acos(pow(z, 1.0/(1.0 + shininess)));
    double phi = 2.0 * M_PI * rand() / double(RAND_MAX);

    // This is sampled over the hemisphere around <0,0,1>
    dvec3 randomDirection = dvec3(sin(alpha)*cos(phi), sin(alpha)*sin(phi), z);

    // So we align it with the intersection normal
    double a = rand() / double(RAND_MAX);
    double b = rand() / double(RAND_MAX);
    double c = rand() / double(RAND_MAX);

    dvec3 ran = dvec3(a,b,c);
    dvec3 x = glm::normalize(glm::cross(ran, reflection));
    dvec3 y = glm::normalize(glm::cross(reflection, x));
    dmat3 rotationMatrix = dmat3(x, y, reflection);
    randomDirection = rotationMatrix*randomDirection;

    ray.position = surfaceIntersection.position + dvec4(randomDirection,1.0) * 0.001;
    ray.direction = dvec4(randomDirection,1.0);
    return;
}

    ////////////
 //////////////////
///     RIP     ///
///  IN PIECES  ///
///  400 lines  ///
///  of   code  ///
///  2012-2013  ///
///             ///
//~~~~~~~~~~~~~~~//