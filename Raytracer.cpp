/*
    Raytracer.cpp
*/

#define GLM_SWIZZLE

#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <glm/gtx/component_wise.hpp>
#include "Raytracer.h"
#include "Intersection.h"
#include "Primitive.h"

typedef glm::mat3 mat3;
typedef glm::mat4 mat4; 
typedef glm::vec3 vec3; 
typedef glm::vec4 vec4; 

int Raytracer::traceRay(Scene* scene, Ray* ray, int depth, glm::vec3& color, float rayRIndex)
{
    int intersectPrimIndex = -1;
    Ray* tempRay = new Ray;                 // Temporary ray for finding intersection
    Ray* shadowRay = new Ray;               // Ray for finding visibility term
    std::vector<Primitive>::iterator it; 
    std::vector<Light>::iterator lit;

    // Iterators are kind of stupid IMO
    int i = 0;
    int mini = 0;
    float mint = 999999999.0;
    float tempt = 0.0f;
    float tmax;
    float tempFlip; // These two are for flipping the normal for refraction
    float normFlip; //
    vec3 direction;
    vec3 position;
    vec4 L;
    vec4 V;
    bool didIntersect = false;
    Intersection* tempIntersect = new Intersection;
    Intersection* lightIntersect = new Intersection;
    Intersection* intersect = new Intersection;

    // Get the closest intersection - naive
    for (it = scene->primitiveList.begin(); it < scene->primitiveList.end(); it++) {

        // Transform ray to object space
        tempRay->position = ray->position * it->inverseTransformation;
        tempRay->direction = ray->direction * it->inverseTransformation;

        // If we intersect, get the intersection point.
        // If it's the closest one yet, save its spot in the primitive list for lighting.
        if (tempFlip = it->intersectionPoint(tempRay, tempIntersect)) {
            if (tempIntersect->t < mint && tempIntersect->t >= 0.0) {
                // Commence gross variable setting.
                normFlip = tempFlip;
                mint = tempIntersect->t;
                mini = i;
                intersect->position = tempIntersect->position;
                intersect->normal = tempIntersect->normal;
                intersect->t = tempIntersect->t; 
                didIntersect = true;
            }
        }

        i++;
    }

    // Lighting calculations
    if (didIntersect) {
        if (depth == 0) {
            intersectPrimIndex = mini;
        }

        // Move intersection back into world space
        lightIntersect->position = intersect->position * scene->primitiveList[mini].transformation;
        lightIntersect->normal = glm::normalize(intersect->normal*scene->primitiveList[mini].inverseTranspose);

        // Loop through the lights to calculate summed light contribution
        for (lit = scene->lightList.begin(); lit < scene->lightList.end(); lit++) {
            float visibility = 1.0;

            // Different ray directions depending on if it's a point or directional source.
            if (lit->point) {
                L = glm::normalize(lit->posdir - lightIntersect->position);
                shadowRay->direction = glm::normalize(lit->posdir - lightIntersect->position);
            } else {
                L = glm::normalize(lit->posdir);
                shadowRay->direction = glm::normalize(lit->posdir);
            }

            // Move the shadow ray away from the object a little bit
            shadowRay->position = lightIntersect->position + 0.001 * shadowRay->direction;
            // Don't count intersections past the light source
            tmax = glm::length(lit->posdir - lightIntersect->position);

            // Loop through each object again to see if our light ray intersects with anything
            int k = 0;
            for (it = scene->primitiveList.begin(); it < scene->primitiveList.end(); it++, k++) {
                // Transform shadowRay into object space
                tempRay->position = shadowRay->position * it->inverseTransformation;
                tempRay->direction = shadowRay->direction * it->inverseTransformation;

                if (it->doesIntersect(tempRay, tmax)) {
                    visibility = 0.0;
                    break;
                }
            }

            //// Shading time! ////
            // Diffuse term
            float nDotL = std::max(glm::dot(lightIntersect->normal, shadowRay->direction),0.0f);
            vec3 lambert = scene->primitiveList[mini].material->diffuse*nDotL;

            // Specular term
            V = glm::normalize(ray->position - lightIntersect->position);
            float nDotH = std::pow(std::max(glm::dot(lightIntersect->normal, glm::normalize(L+V)),0.0f),
                                   scene->primitiveList[mini].material->shininess); 
            vec3 phong = scene->primitiveList[mini].material->specular*nDotH;

            // Add it all together with attenuation
            float attenuation;
            float r;
            r = glm::length(lit->posdir - lightIntersect->position);
            attenuation = (1.0/(scene->attenuation[0] + scene->attenuation[1]*r + scene->attenuation[1]*r*r));
            color += (phong + lambert)*visibility*lit->color*attenuation*scene->primitiveList[mini].material->alpha;
        }

        if (depth < maxDepth) {
            // Only do reflections if material has specular component
            if (glm::length(scene->primitiveList[mini].material->specular) > 0.0) {
                Ray* newRay = new Ray;
                vec4 newDirection = ray->direction;
                newDirection[3] = 0.0;
                // New ray
                newRay->position = lightIntersect->position + lightIntersect->normal * 0.001;
                newRay->direction = newDirection-2*lightIntersect->normal*glm::dot(newDirection,lightIntersect->normal);
                vec3 newColor = vec3(0.0,0.0,0.0);
                traceRay(scene, newRay, ++depth, newColor, rayRIndex);
                color = color + newColor*scene->primitiveList[mini].material->specular;
            }

            // Only do refractions if material is not opaque
            if (scene->primitiveList[mini].material->alpha < 1.0) {
                Ray* newRay = new Ray;
                vec4 newDirection = ray->direction;
                newDirection[3] = 0.0;
                float n = rayRIndex / scene->primitiveList[mini].material->rindex;
                vec4 normal = lightIntersect->normal * normFlip;
                float cosI = -glm::dot(normal, newDirection);
                float cosT = 1.0 - n*n*(1.0 - cosI*cosI);
                if (cosT > 0.0f) { 
                    newRay->direction = n*ray->direction + (n*cosI - sqrtf(cosT))*normal;
                    newRay->position = lightIntersect->position + newRay->direction * 0.001;
                    vec3 newColor = vec3(0.0,0.0,0.0);
                    traceRay(scene, newRay, ++depth, newColor, scene->primitiveList[mini].material->rindex);
                    color = newColor + color*scene->primitiveList[mini].material->diffuse;
                }
            }

            // Basic Global Illumination
            // gi is how many random rays to fire off at a time.
            if (depth < scene->gidepth) {
                vec3 sumColor;
                for (int g = 0; g < 40; g++) {
                    float div = 1.0f / (RAND_MAX+1.0f);
                    float z = rand()*div+rand()*div*div-0.5;
                    float x = rand()*div+rand()*div*div-0.5;
                    float y = rand()*div+rand()*div*div-0.5;
                    Ray* newRay = new Ray;
                    newRay->position = lightIntersect->position + lightIntersect->normal * 0.001;
                    newRay->direction = glm::normalize(lightIntersect->normal+vec4(x,y,z,0));
                    vec3 newColor = vec3(0.0,0.0,0.0);
                    traceRay(scene, newRay, ++depth, newColor, rayRIndex);
                    sumColor += newColor;
                }
                color += (sumColor/((float) 60.0f));
            }
        }

        // Tack on ambient and emissive terms
        color = glm::clamp(color + scene->primitiveList[mini].ambient + scene->primitiveList[mini].material->emission,0.0,1.0);
    }

    return intersectPrimIndex;
}
