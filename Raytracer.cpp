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

typedef glm::mat3 mat3;
typedef glm::mat4 mat4;
typedef glm::vec3 vec3;
typedef glm::vec4 vec4;

Raytracer::Raytracer(Scene* scene)
{
    this->scene = scene;
    this->maxDepth = scene->maxDepth;    
}

Raytracer::~Raytracer()
{
}

// Main path tracing routine
vec3 Raytracer::pathTraceRay(const Ray& ray, int depth, float weight, int bounce)
{
    int closestIntersectionIndex;
    Intersection closestIntersection;

    if (findClosestIntersection(ray, closestIntersectionIndex, closestIntersection)) {

        Primitive intersectedObject = scene->primitiveList[closestIntersectionIndex];
        Intersection surfaceIntersection;

        // Move intersection from object to world space
        surfaceIntersection.position = closestIntersection.position * intersectedObject.transformation;
        surfaceIntersection.normal   = glm::normalize(closestIntersection.normal * intersectedObject.inverseTranspose);
        
        // Russian roulette
        // We'll sum the reflectance coefficients of the material after first weighting them according to the luminance
        // of the NTSC (1953) color space. We'll continue the ray with that probability, otherwise we will terminate.
        // Note: this assumes the coefficients sum <= 1
        float p, ran;
        vec3 coefficients = vec3(0.0f);

        switch (intersectedObject.material->type) {
            case LAMBERTIAN:
                coefficients = intersectedObject.material->diffuse;
                break;
            case GLOSSY:
                coefficients = intersectedObject.material->diffuse + intersectedObject.material->specular;
                break;
            case REFLECTIVE:
                coefficients = intersectedObject.material->reflection;
                break;
            case TRANSMISSIVE:
                coefficients = intersectedObject.material->transmission;
                break;
        }
        p = coefficients.x*0.2989f + coefficients.y*0.5866f + coefficients.z*0.1145f;

        ran = rand() / double(RAND_MAX);
        vec3 color = vec3(0.0f);

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
                        vec4 randomRayDir = vec4(1.0);
                        randomRayDir.z = rand() / double(RAND_MAX);
                        float xyproj = sqrt(1 - randomRayDir.z * randomRayDir.z);
                        float phi = 2 * M_PI * rand() / double(RAND_MAX);
                        randomRayDir.x = xyproj * cos(phi);
                        randomRayDir.y = xyproj * sin(phi);
                        randomRayDir = glm::normalize(randomRayDir);
                        if (glm::dot(randomRayDir, surfaceIntersection.normal) < 0.0f) {
                            randomRayDir = -randomRayDir;
                        }
                        Ray randomRay;
                        randomRay.position = surfaceIntersection.position + randomRayDir * 0.001;
                        randomRay.direction = randomRayDir;
                        vec3 randomColor = pathTraceRay(randomRay, depth + 1, weight, ++bounce);

                        color += intersectedObject.material->diffuse * glm::dot(surfaceIntersection.normal, randomRay.direction) * randomColor * 2.0f;
                    }

                    break;
                }
/*
                case GLOSSY: {
                    if (scene->directLighting) {
                        color += directLighting(ray, lightIntersect, intersectedObject);
                    }

                    if (scene->indirectLighting) {
                        color += indirectDiffuseLighting(lightIntersect, intersectedObject) / p;
                    }

                    break;
                }

                case REFLECTIVE: {
                    if (scene->directLighting) {
                        color += directLighting(ray, lightIntersect, intersectedObject);
                    }

                    if (scene->indirectLighting) {
                        color += indirectDiffuseLighting(lightIntersect, intersectedObject) / p;
                    }

                    break;
                }

                case TRANSMISSIVE: {
                    if (scene->directLighting) {
                        color += directLighting(ray, lightIntersect, intersectedObject);
                    }

                    if (scene->indirectLighting) {
                        color += indirectDiffuseLighting(lightIntersect, intersectedObject) / p;
                    }

                    break;
                }
*/
            }
        }

        return color;
    } else {
        // This could be replaced with a skybox in the future
        return vec3(0.0);
    }
}

bool Raytracer::findClosestIntersection(const Ray& ray, int& closestIntersectionIndex, Intersection& closestIntersection)
{
    float minimumT = FLT_MAX;
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
        tempRay.position = ray.position * it->inverseTransformation;
        tempRay.direction = ray.direction * it->inverseTransformation;

        if (tempFlip = it->getIntersectionPoint(tempRay, tempIntersect)) {
            if (tempIntersect.t < minimumT && tempIntersect.t >= 0.0f) {
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

    // Flip transmissive norm?

    return didIntersect;
}

vec3 Raytracer::directLighting(const Ray& ray, const Intersection& surfaceIntersection, 
                               const Primitive& intersectedObject)
{
    std::vector<Primitive>::iterator it;
    std::vector<AreaLight*>::iterator lit; 
    Ray shadowRay;
    Ray tempRay;
    float tmax;
    bool visible = true;

    vec3 color = vec3(0.0f);
    vec3 lightIntensity = vec3(0.0f);
    vec3 rayToLight = vec3(0.0f);

    // Loop over all the area lights
    for (lit = scene->areaLightList.begin(); lit < scene->areaLightList.end(); lit++) {

        // Get sample from the area light source (puts result into rayToLight)
        (*lit)->getSample(surfaceIntersection.position.xyz, surfaceIntersection.normal.xyz, lightIntensity, rayToLight);

        tmax = glm::length(rayToLight); // So we don't count intersections behind the light
        shadowRay.direction = vec4(glm::normalize(rayToLight),1.0);
        shadowRay.position  = surfaceIntersection.position + 0.001f * shadowRay.direction;

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
            color += intersectedObject.material->diffuse * lightIntensity / (float)M_PI;
            //std::cout << "lightIntensity: " << lightIntensity.x << "," << lightIntensity.y << "," << lightIntensity.z << "\n";

            // Specular term
            if (intersectedObject.material->type == GLOSSY) {
                vec4 V = glm::normalize(ray.position - surfaceIntersection.position);
                float nDotH = std::pow(std::max(glm::dot(surfaceIntersection.normal, 
                                                         glm::normalize(shadowRay.direction+V)),
                                                0.0f),
                                       intersectedObject.material->shininess); 
                color += intersectedObject.material->specular*nDotH;
            }            
            float attenuation;
            float r = tmax;
            attenuation = (1.0f/(scene->attenuation[0] + scene->attenuation[1]*r + scene->attenuation[1]*r*r));
            color *= attenuation;
        }
    }

    color = glm::clamp(color,0.0f,1.0f);
    return color;
}

vec3 Raytracer::indirectDiffuseLighting(const Intersection& surfaceIntersection, const Primitive& intersectedObject)
{        

    return vec3(0.0);
}

/*
void Raytracer::traceRay(Ray* ray, int depth, vec3& color, float weight, int bounce, float rayRIndex)
{
    int intersectPrimIndex = -1;
    Ray* tempRay = new Ray;                 // Temporary ray for finding intersection
    std::vector<Primitive>::iterator it; 

    int i = 0;
    int mini = 0;
    float mint = 999999999.0;
    float tempt = 0.0f;
    float tmax;
    float tempFlip; // These two are for flipping the normal for refraction
    float normFlip; //
    vec3 direction;
    vec3 position;
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
        float p = std::max(intersectedObject.material->diffuse.x,std::max(intersectedObject.material->diffuse.y,intersectedObject.material->diffuse.z));
        float ran = rand() / double(RAND_MAX);

        if (p < ran) {
            if (scene->directLighting && intersectedObject.material->type == EMISSIVE && bounce > 0) {
                color = vec3(0.0);
                return;
            }

            color = intersectedObject.material->emission*(1.0f/(1.0f-p));
            return;
        }

        if (scene->directLighting) {
            color += directLighting(ray, lightIntersect, &intersectedObject);
        }  

        vec4 randomRayDir = vec4(1.0);
        randomRayDir.z = rand() / double(RAND_MAX);
        float xyproj = sqrt(1 - randomRayDir.z * randomRayDir.z);
        float phi = 2 * M_PI * rand() / double(RAND_MAX);
        randomRayDir.x = xyproj * cos(phi);
        randomRayDir.y = xyproj * sin(phi);
        randomRayDir = glm::normalize(randomRayDir);
        if (glm::dot(randomRayDir, lightIntersect->normal) < 0.0f) {
            randomRayDir = -randomRayDir;
        }
        vec3 randomColor = vec3(0.0);
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
            vec4 newDirection = ray->direction;
            newDirection[3] = 0.0;
            vec3 specIndirectColor = vec3(0.0);

            if (intersectedObject.material->type == TRANSMISSIVE) {
                float n = rayRIndex / intersectedObject.material->rindex;
                vec4 normal = lightIntersect->normal * normFlip;
                float cosI = -glm::dot(normal, newDirection);
                float cosT = 1.0 - n*n*(1.0 - cosI*cosI);
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
                float u = rand() / double(RAND_MAX);
                Ray* randomRay = new Ray;
                vec4 randomDirection = vec4(0.0,0.0,0.0,1.0);
                randomDirection.z = rand() / double(RAND_MAX);
                float theta = 2*M_PI*u;
                float xyproj = sqrt(1.0 - randomDirection.z*randomDirection.z);
                randomDirection.x = xyproj*cos(u);
                randomDirection.y = xyproj*sin(u);

                randomRay->direction = randomDirection;
                randomRay->position = lightIntersect->position + lightIntersect->normal * 0.001;
                vec3 randomColor = vec3(0.0);
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
vec3 Raytracer::directLighting(Ray* ray, Intersection* lightIntersect, Primitive* intersectedObject)
{
    vec4 L;
    vec4 V;
    std::vector<Primitive>::iterator it;
    std::vector<Light>::iterator lit; 
    Ray* shadowRay = new Ray;
    Ray* tempRay = new Ray;

    vec3 color = vec3(0.0);

    // Loop through the lights to calculate summed light contribution
    for (lit = scene->lightList.begin(); lit < scene->lightList.end(); lit++) {
        float visibility = 1.0;
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
            float tmax = glm::length(lit->posdir - lightIntersect->position);

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
            float nDotL = std::max(glm::dot(lightIntersect->normal, shadowRay->direction),0.0f);
            vec3 lambert = intersectedObject->material->diffuse*nDotL;

            // Specular term
            V = glm::normalize(ray->position - lightIntersect->position);
            float nDotH = std::pow(std::max(glm::dot(lightIntersect->normal, glm::normalize(L+V)),0.0f),
                                   intersectedObject->material->shininess); 
            vec3 phong = intersectedObject->material->specular*nDotH;

            // Add it all together with attenuation
            float attenuation;
            float r;
            r = glm::length(lit->posdir - lightIntersect->position);
            attenuation = (1.0/(scene->attenuation[0] + scene->attenuation[1]*r + scene->attenuation[1]*r*r));
            color += (phong + lambert)*visibility*lit->color*attenuation*intersectedObject->material->alpha;
        } else if (lit->type == AREA) {
            if (lit->areaType == RECT) {
                visibility = 0.0;
                float visAdd = 1.0 / ((float)(lit->numSamples*lit->numSamples));
                bool vis;
                vec4 pos;
                for (int i = 0; i < lit->numSamples; ++i) {
                    for (int j = 0; j < lit->numSamples; ++j) {
                        vis = true;
                        //std::cout << "i: " << i << ", j: " << j << "\n";
                        float u = rand() / double(RAND_MAX);
                        float v = rand() / double(RAND_MAX);
                        // For each sample within the area light source, the initial position is the lower right coordinate
                        // of the current region of the source. Then we add a random amount within that region.
                        pos = lit->posdir + lit->rightStep*(float)j + lit->rightStep*u + lit->upStep*(float)i + lit->upStep*v;
                        pos[3] = 1.0;

                        L = glm::normalize(pos - lightIntersect->position);
                        shadowRay->direction = glm::normalize(pos - lightIntersect->position);
                        shadowRay->position = lightIntersect->position + 0.001 * shadowRay->direction
                        float tmax = glm::length(pos - lightIntersect->position);

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
                            float nDotL = std::max(glm::dot(lightIntersect->normal, shadowRay->direction),0.0f);
                            vec3 lambert = intersectedObject->material->diffuse*nDotL;

                            // Specular term
                            V = glm::normalize(ray->position - lightIntersect->position);
                            float nDotH = std::pow(std::max(glm::dot(lightIntersect->normal, glm::normalize(L+V)),0.0f),
                                                   intersectedObject->material->shininess); 
                            vec3 phong = intersectedObject->material->specular*nDotH;

                            // Add it all together with attenuation
                            float attenuation;
                            float r;
                            r = glm::length(pos - lightIntersect->position);
                            attenuation = (1.0/(scene->attenuation[0] + scene->attenuation[1]*r + scene->attenuation[1]*r*r));
                            color += (phong + lambert)*visAdd*lit->color*attenuation*intersectedObject->material->alpha;
                        }
                    }
                }
                
                //visibility /= (float) (lit->numSamples*lit->numSamples);

                //std::cout << "Visibility: " << visibility << "\n";

                pos = lit->posdir + lit->rightStep*(float)lit->numSamples/2.0 + lit->upStep*(float)lit->numSamples/2.0;
                pos[3] = 1.0;

                L = glm::normalize(pos - lightIntersect->position);
                shadowRay->direction = glm::normalize(pos - lightIntersect->position);
                shadowRay->position = lightIntersect->position + 0.001 * shadowRay->direction;

                // Diffuse term
                float nDotL = std::max(glm::dot(lightIntersect->normal, shadowRay->direction),0.0f);
                vec3 lambert = intersectedObject->material->diffuse*nDotL;

                // Specular term
                V = glm::normalize(ray->position - lightIntersect->position);
                float nDotH = std::pow(std::max(glm::dot(lightIntersect->normal, glm::normalize(L+V)),0.0f),
                                       intersectedObject->material->shininess); 
                vec3 phong = intersectedObject->material->specular*nDotH;

                // Add it all together with attenuation
                float attenuation;
                float r;
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