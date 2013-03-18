/*
    Material.h
    
    Public class to hold material properties of object.
*/

#ifndef I_MATRL
#define I_MATRL

#include <glm/glm.hpp>

enum MaterialType { LAMBERTIAN, GLOSSY, TRANSMISSIVE, REFLECTIVE, EMISSIVE };

typedef struct _Material {
	MaterialType type;
    glm::vec3 diffuse;
    glm::vec3 specular;
    float shininess;
    float rindex;
    float alpha;

    // Light source
    glm::vec3 emission;
} Material;

#endif
