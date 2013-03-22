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

//**** Pathtracer
dvec3 Raytracer::pathTraceRay(const Ray& ray, int depth, double weight, int bounce)
{
    int closestIntersectionIndex;
    Intersection closestIntersection; 
    bool haveAnIntersection = findClosestIntersection(ray, closestIntersectionIndex, closestIntersection);
    
    dvec3 color = dvec3(0.0);

    if (!haveAnIntersection) {
        // This could be replaced with a skybox in the future
        return color;
    } else {

        // Grab the object that our intersection is on
        Primitive intersectedObject = scene->primitiveList[closestIntersectionIndex];
        
        // Get intersection from object to world space
        Intersection surfaceIntersection;
        surfaceIntersection.position = closestIntersection.position * intersectedObject.transformation;
        surfaceIntersection.normal   = glm::normalize(closestIntersection.normal * intersectedObject.inverseTranspose);

        color = intersectedObject.material->emission;

        //**** Russian roulette 
        // If the weight of the ray falls below 0.01, 50% chance of terminating ray -> double weight of survivors
        double terminatingP = 0.5;
        double survivorBonus = 1.0;
        if (weight <= 0.01) {
            double meteoRandom = rand() / double(RAND_MAX);
            if (meteoRandom > terminatingP) {
                return color * 2.0;
            } else {
                survivorBonus = 2.0;
            }
        }

        if (intersectedObject.material->type == LAMBERTIAN) {
            if (scene->directLighting) {
                color += directLighting(ray, surfaceIntersection, intersectedObject) * weight * survivorBonus;
            }

            if (scene->indirectLighting) {
                Ray randomRay;
                sampleCosineWeightedHemisphere(randomRay, surfaceIntersection); 
                //sampleUniformHemisphere(randomRay, surfaceIntersection); 
                dvec3 traceColor = pathTraceRay(randomRay, ++depth, weight * glm::dot(intersectedObject.material->diffuse, 
                                                                                      dvec3(0.333)), ++bounce);
                color += traceColor * intersectedObject.material->diffuse * weight * survivorBonus;
            }
        }

        else if (intersectedObject.material->type == GLOSSY) {            
            if (scene->directLighting) {
                color += directLighting(ray, surfaceIntersection, intersectedObject) * weight * survivorBonus;
            }

            if (scene->indirectLighting) {
                //**** Importance sampling
                // The probability of sampling diffuse interreflection is, how I chose, the sum of the diffuse terms over the
                // sum of the diffuse and the specular terms of the material reflectance.

                // dvec3 colorSpace = dvec3(0.2989, 0.5866, 0.1145);  Here are some other possibilities for determining the
                // double third = 1.0 / 3.0;                          probability for sampling the diffuse interreflection
                // dvec3 averager = dvec3(third, third, third);       versus the specular interreflection.

                double diffuseSum = glm::dot(intersectedObject.material->diffuse, dvec3(1.0));
                double specularSum = glm::dot(intersectedObject.material->specular, dvec3(1.0));

                double probDiffuse = diffuseSum / (diffuseSum + specularSum);
                double impRand = rand() / double(RAND_MAX);
                Ray randomRay;

                // Diffuse
                if (impRand < probDiffuse) {
                    sampleCosineWeightedHemisphere(randomRay, surfaceIntersection); 
                    //sampleUniformHemisphere(randomRay, surfaceIntersection); 
                    dvec3 traceColor = pathTraceRay(randomRay, ++depth, weight * glm::dot(intersectedObject.material->diffuse, 
                                                                                          dvec3(0.333)), ++bounce);
                    color += traceColor * intersectedObject.material->diffuse * weight * survivorBonus / probDiffuse;
                }

                // Specular
                else {
                    sampleSpecularLobe(randomRay, surfaceIntersection, dvec3(glm::normalize(glm::reflect(ray.direction, surfaceIntersection.normal)).xyz),
                                       intersectedObject.material->shininess); 
                    //sampleUniformHemisphere(randomRay, surfaceIntersection); 
                    dvec3 traceColor = pathTraceRay(randomRay, ++depth, weight * glm::dot(intersectedObject.material->specular, 
                                                                                          dvec3(0.333)), ++bounce);
                    color += traceColor * intersectedObject.material->specular * weight * survivorBonus / (1.0 - probDiffuse);
                }
            }
        }

        else if (intersectedObject.material->type == REFLECTIVE) {
            Ray reflectedRay;
            reflectedRay.direction = glm::reflect(ray.direction, surfaceIntersection.normal);
            reflectedRay.position  = surfaceIntersection.position + reflectedRay.direction * 0.001;
            color += pathTraceRay(reflectedRay, ++depth, weight * glm::dot(intersectedObject.material->specular, dvec3(0.333)), bounce) 
                     * weight;
        }

        else if (intersectedObject.material->type == EMISSIVE) {
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

        if (tempFlip = it->getIntersectionPoint(tempRay, tempIntersect)) {
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

            if (it->doesRayIntersect(tempRay, tmax)) {
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
    double a = rand() / double(RAND_MAX);
    double b = rand() / double(RAND_MAX);

    double theta = 2.0 * M_PI * a;
    double phi = acos(2.0 * b - 1);

    dvec3 randomDirection = dvec3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
    dvec3 normal = dvec3(surfaceIntersection.normal.xyz);

    // Since this is uniformly sampled over the sphere, we flip it to stay in the correct hemisphere
    if (glm::dot(randomDirection, normal) < 0.0) {
        randomDirection = -randomDirection;
    }

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

dvec3 Raytracer::indirectDiffuseLighting(const Intersection& surfaceIntersection, const Primitive& intersectedObject)
{        
    return dvec3(0.0);
}

    ////////////
 //////////////////
///     RIP     ///
///  IN PIECES  ///
///             ///
///  Old  code  ///
///  2012-2013  ///
///             ///
//~~~~~~~~~~~~~~~//
/*
                        // Get phong illumination
                        dvec4 V = glm::normalize(ray.position - surfaceIntersection.position);
                        double nDotH = std::pow(std::max(glm::dot(surfaceIntersection.normal, 
                                                                 glm::normalize(randomRay.direction+V)),0.0f),
                                               intersectedObject.material->shininess); 
                        double specProbability = pow(glm::dot(dvec3(randomRay.direction.xyz), perfectReflection), intersectedObject.material->shininess) *
                                                ((intersectedObject.material->shininess + 1.0f) / (2.0f * M_PI));
                        color += pathTraceRay(randomRay, ++depth, weight*specularAverage, ++bounce) * 
                                              intersectedObject.material->specular * nDotH / specProbability / p * 
                                              glm::dot(dvec3(randomRay.direction.xyz),dvec3(surfaceIntersection.normal.xyz));
                        Ray randomRay; 
                        dvec3 perfectReflection = glm::reflect(ray.direction, surfaceIntersection.normal).xyz;
                        sampleHemisphereSpecular(randomRay, surfaceIntersection, perfectReflection, 
                                                 intersectedObject.material->shininess);

        // Importance sampling: we'll either sample the diffuse or the specular component of indirect illumination
        double third = 1.0f/3.0f;
        dvec3 average = dvec3(third, third, third);
        double diffuseAverage = glm::dot(intersectedObject.material->diffuse, average);
        double specularAverage = glm::dot(intersectedObject.material->specular, average);
        //std::cout << weight << "\n";
        double sampleWeight = diffuseAverage / (diffuseAverage + specularAverage);


        double p;
        double ran = rand() / double(RAND_MAX);

        // sampleWeight is the portion of the reflectance of the material that is diffuse, so if our random number is
        // less than it then we are calculating the diffuse part, and the other specular.
        Ray randomRay; 
        if (ran < sampleWeight) {
            sampleHemisphereUniformly(randomRay, surfaceIntersection);
            color = pathTraceRay(randomRay, ++depth, weight*diffuseAverage, ++bounce) * 
                    intersectedObject.material->diffuse / sampleWeight;
        } else {


        }
/*
        dvec3 coefficients = dvec3(0.0f);
        p = coefficients.x*0.2989f + coefficients.y*0.5866f + coefficients.z*0.1145f;

        dvec3 color = dvec3(0.0f);

        if (scene->directLighting) {
            color += directLighting(ray, surfaceIntersection, intersectedObject);
            //std::cout << "color: " << color.x << "," << color.y << "," << color.z << "\n";
        }

        if (ran > p) {
            if (scene->directLighting && intersectedObject.material->type == EMISSIVE && bounce > 0) {
                return color;
            } else {
                return intersectedObject.material->emission / (1.0f - p);
            }
        } else {
            switch (intersectedObject.material->type) {
                case LAMBERTIAN: {
                    if (scene->indirectLighting) {
                        // color += indirectDiffuseLighting(surfaceIntersection, intersectedObject) / p;    
                        dvec3 randomColor = pathTraceRay(randomRay, depth + 1, weight, ++bounce);

                        color += intersectedObject.material->diffuse * glm::dot(surfaceIntersection.normal, randomRay.direction) * randomColor * 2.0f;
                    }

                    break;
                }

            }
        }
*/

/*
void Raytracer::traceRay(Ray* ray, int depth, dvec3& color, double weight, int bounce, double rayRIndex)
{
    int intersectPrimIndex = -1;
    Ray* tempRay = new Ray;                 // Temporary ray for finding intersection
    std::vector<Primitive>::iterator it; 

    int i = 0;
    int mini = 0;
    double mint = 999999999.0;
    double tempt = 0.0f;
    double tmax;
    double tempFlip; // These two are for flipping the normal for refraction
    double normFlip; //
    dvec3 direction;
    dvec3 position;
    bool didIntersect = false;
    Intersection* tempIntersect = new Intersection;
    Intersection* intersect = new Intersection;    
    Intersection* lightIntersect = new Intersection;

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

    Primitive intersectedObject = scene->primitiveList[mini];

    // Move intersection from object to world space
    lightIntersect->position    = intersect->position * intersectedObject.transformation;
    lightIntersect->normal      = glm::normalize(intersect->normal * intersectedObject.inverseTranspose);

    // Pathtracing
    if (scene->traceType == PATH) {        
        double p = std::max(intersectedObject.material->diffuse.x,std::max(intersectedObject.material->diffuse.y,intersectedObject.material->diffuse.z));
        double ran = rand() / double(RAND_MAX);

        if (p < ran) {
            if (scene->directLighting && intersectedObject.material->type == EMISSIVE && bounce > 0) {
                color = dvec3(0.0);
                return;
            }

            color = intersectedObject.material->emission*(1.0f/(1.0f-p));
            return;
        }

        if (scene->directLighting) {
            color += directLighting(ray, lightIntersect, &intersectedObject);
        }  

        dvec4 randomRayDir = dvec4(1.0);
        randomRayDir.z = rand() / double(RAND_MAX);
        double xyproj = sqrt(1 - randomRayDir.z * randomRayDir.z);
        double phi = 2 * M_PI * rand() / double(RAND_MAX);
        randomRayDir.x = xyproj * cos(phi);
        randomRayDir.y = xyproj * sin(phi);
        randomRayDir = glm::normalize(randomRayDir);
        if (glm::dot(randomRayDir, lightIntersect->normal) < 0.0f) {
            randomRayDir = -randomRayDir;
        }
        dvec3 randomColor = dvec3(0.0);
        Ray* randomRay = new Ray();
        randomRay->position = lightIntersect->position + randomRayDir * 0.001;
        randomRay->direction = randomRayDir;
        traceRay(randomRay, depth + 1, randomColor, weight, ++bounce, rayRIndex);

        if (intersectedObject.material->type == LAMBERTIAN) { 
            color += intersectedObject.material->diffuse * glm::dot(lightIntersect->normal, randomRay->direction) * randomColor / p;
        }
    } 

    // Raytracing
    else if (scene->traceType == RAY) {
        if (scene->directLighting) {
            if (intersectedObject.material->type == EMISSIVE) {
                color += intersectedObject.material->emission;
            }
            
            // Direct lighting for lambertian and glossy surfaces
            if (intersectedObject.material->type == LAMBERTIAN || intersectedObject.material->type == GLOSSY) {
                color += directLighting(ray, lightIntersect, &intersectedObject);
            }
        }  

        // Global illumination
        if (depth < maxDepth) {
            // Specular indirect lighting
            Ray* newRay = new Ray;
            dvec4 newDirection = ray->direction;
            newDirection[3] = 0.0;
            dvec3 specIndirectColor = dvec3(0.0);

            if (intersectedObject.material->type == TRANSMISSIVE) {
                double n = rayRIndex / intersectedObject.material->rindex;
                dvec4 normal = lightIntersect->normal * normFlip;
                double cosI = -glm::dot(normal, newDirection);
                double cosT = 1.0 - n*n*(1.0 - cosI*cosI);
                if (cosT > 0.0f) { 
                    newRay->direction = n*ray->direction + (n*cosI - sqrtf(cosT))*normal;
                    newRay->position = lightIntersect->position + newRay->direction * 0.001;
                    traceRay(newRay, ++depth, specIndirectColor, weight, ++bounce, intersectedObject.material->rindex);
                }
            } else {
                newRay->position = lightIntersect->position + lightIntersect->normal * 0.001;
                newRay->direction = newDirection-2*lightIntersect->normal*glm::dot(newDirection,lightIntersect->normal);
                traceRay(newRay, ++depth, specIndirectColor, weight, ++bounce, rayRIndex);
            }
            delete newRay;

            if (scene->traceType == RAY) {
                if (intersectedObject.material->type == GLOSSY) {
                    color = color + specIndirectColor*intersectedObject.material->specular;
                } else if (intersectedObject.material->type == TRANSMISSIVE || intersectedObject.material->type == REFLECTIVE) {
                    color = specIndirectColor;
                }
            } else if (scene->traceType == PATH) {
                // Diffuse indirect lighting
                double u = rand() / double(RAND_MAX);
                Ray* randomRay = new Ray;
                dvec4 randomDirection = dvec4(0.0,0.0,0.0,1.0);
                randomDirection.z = rand() / double(RAND_MAX);
                double theta = 2*M_PI*u;
                double xyproj = sqrt(1.0 - randomDirection.z*randomDirection.z);
                randomDirection.x = xyproj*cos(u);
                randomDirection.y = xyproj*sin(u);

                randomRay->direction = randomDirection;
                randomRay->position = lightIntersect->position + lightIntersect->normal * 0.001;
                dvec3 randomColor = dvec3(0.0);
                traceRay(randomRay, ++depth, randomColor, weight, ++bounce, rayRIndex);
                //color += randomColor;
            } else {
                std::cout << "???\n";
            }
        }
    }
    

    delete tempIntersect;
    delete lightIntersect;
    delete intersect;
    delete tempRay;
    return;
}

// TODO: refactor this whole mess. Lots of code repetition 
dvec3 Raytracer::directLighting(Ray* ray, Intersection* lightIntersect, Primitive* intersectedObject)
{
    dvec4 L;
    dvec4 V;
    std::vector<Primitive>::iterator it;
    std::vector<Light>::iterator lit; 
    Ray* shadowRay = new Ray;
    Ray* tempRay = new Ray;

    dvec3 color = dvec3(0.0);

    // Loop through the lights to calculate summed light contribution
    for (lit = scene->lightList.begin(); lit < scene->lightList.end(); lit++) {
        double visibility = 1.0;
        int numberOfShadowRays = 1; // Only really matters if area source

        // Visibility depends on what kind of light source we have
        if (lit->type == POINT || lit->type == DIRECTIONAL) {
            // Different ray directions depending on if it's a point or directional source.
            if (lit->type == POINT) {
                L = glm::normalize(lit->posdir - lightIntersect->position);
                shadowRay->direction = glm::normalize(lit->posdir - lightIntersect->position);
            } else if (lit->type == DIRECTIONAL) {
                L = glm::normalize(lit->posdir);
                shadowRay->direction = glm::normalize(lit->posdir);
            } 

            // Move the shadow ray away from the object a little bit
            shadowRay->position = lightIntersect->position + 0.001 * shadowRay->direction;
            // Don't count intersections past the light source
            double tmax = glm::length(lit->posdir - lightIntersect->position);

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
            // Diffuse term
            double nDotL = std::max(glm::dot(lightIntersect->normal, shadowRay->direction),0.0f);
            dvec3 lambert = intersectedObject->material->diffuse*nDotL;

            // Specular term
            V = glm::normalize(ray->position - lightIntersect->position);
            double nDotH = std::pow(std::max(glm::dot(lightIntersect->normal, glm::normalize(L+V)),0.0f),
                                   intersectedObject->material->shininess); 
            dvec3 phong = intersectedObject->material->specular*nDotH;

            // Add it all together with attenuation
            double attenuation;
            double r;
            r = glm::length(lit->posdir - lightIntersect->position);
            attenuation = (1.0/(scene->attenuation[0] + scene->attenuation[1]*r + scene->attenuation[1]*r*r));
            color += (phong + lambert)*visibility*lit->color*attenuation*intersectedObject->material->alpha;
        } else if (lit->type == AREA) {
            if (lit->areaType == RECT) {
                visibility = 0.0;
                double visAdd = 1.0 / ((double)(lit->numSamples*lit->numSamples));
                bool vis;
                dvec4 pos;
                for (int i = 0; i < lit->numSamples; ++i) {
                    for (int j = 0; j < lit->numSamples; ++j) {
                        vis = true;
                        //std::cout << "i: " << i << ", j: " << j << "\n";
                        double u = rand() / double(RAND_MAX);
                        double v = rand() / double(RAND_MAX);
                        // For each sample within the area light source, the initial position is the lower right coordinate
                        // of the current region of the source. Then we add a random amount within that region.
                        pos = lit->posdir + lit->rightStep*(double)j + lit->rightStep*u + lit->upStep*(double)i + lit->upStep*v;
                        pos[3] = 1.0;

                        L = glm::normalize(pos - lightIntersect->position);
                        shadowRay->direction = glm::normalize(pos - lightIntersect->position);
                        shadowRay->position = lightIntersect->position + 0.001 * shadowRay->direction
                        double tmax = glm::length(pos - lightIntersect->position);

                        for (it = scene->primitiveList.begin(); it < scene->primitiveList.end(); it++) {
                            // Transform shadowRay into object space
                            tempRay->position = shadowRay->position * it->inverseTransformation;
                            tempRay->direction = shadowRay->direction * it->inverseTransformation;

                            if (it->doesIntersect(tempRay, tmax)) {
                                vis = false;
                                //break;
                            }
                        }
                        if (vis) {
                            visibility += visAdd;
                            /*L = glm::normalize(pos - lightIntersect->position);
                            shadowRay->direction = glm::normalize(pos - lightIntersect->position);
                            shadowRay->position = lightIntersect->position + 0.001 * shadowRay->direction;

                            // Diffuse term
                            double nDotL = std::max(glm::dot(lightIntersect->normal, shadowRay->direction),0.0f);
                            dvec3 lambert = intersectedObject->material->diffuse*nDotL;

                            // Specular term
                            V = glm::normalize(ray->position - lightIntersect->position);
                            double nDotH = std::pow(std::max(glm::dot(lightIntersect->normal, glm::normalize(L+V)),0.0f),
                                                   intersectedObject->material->shininess); 
                            dvec3 phong = intersectedObject->material->specular*nDotH;

                            // Add it all together with attenuation
                            double attenuation;
                            double r;
                            r = glm::length(pos - lightIntersect->position);
                            attenuation = (1.0/(scene->attenuation[0] + scene->attenuation[1]*r + scene->attenuation[1]*r*r));
                            color += (phong + lambert)*visAdd*lit->color*attenuation*intersectedObject->material->alpha;
                        }
                    }
                }
                
                //visibility /= (double) (lit->numSamples*lit->numSamples);

                //std::cout << "Visibility: " << visibility << "\n";

                pos = lit->posdir + lit->rightStep*(double)lit->numSamples/2.0 + lit->upStep*(double)lit->numSamples/2.0;
                pos[3] = 1.0;

                L = glm::normalize(pos - lightIntersect->position);
                shadowRay->direction = glm::normalize(pos - lightIntersect->position);
                shadowRay->position = lightIntersect->position + 0.001 * shadowRay->direction;

                // Diffuse term
                double nDotL = std::max(glm::dot(lightIntersect->normal, shadowRay->direction),0.0f);
                dvec3 lambert = intersectedObject->material->diffuse*nDotL;

                // Specular term
                V = glm::normalize(ray->position - lightIntersect->position);
                double nDotH = std::pow(std::max(glm::dot(lightIntersect->normal, glm::normalize(L+V)),0.0f),
                                       intersectedObject->material->shininess); 
                dvec3 phong = intersectedObject->material->specular*nDotH;

                // Add it all together with attenuation
                double attenuation;
                double r;
                r = glm::length(pos - lightIntersect->position);
                attenuation = (1.0/(scene->attenuation[0] + scene->attenuation[1]*r + scene->attenuation[1]*r*r));
                color += (phong + lambert)*visibility*lit->color*attenuation*intersectedObject->material->alpha;
            }
        }

    }

    // deleted ambient and emision
    color = glm::clamp(color,0.0,1.0);
    delete shadowRay;
    delete tempRay;
    return color;
}
*/