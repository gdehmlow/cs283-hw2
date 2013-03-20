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
    glm::vec3 reflection;
    glm::vec3 transmission;
    glm::vec3 emission;

    float shininess;
    float rindex;
    float alpha;
} Material;

#endif
