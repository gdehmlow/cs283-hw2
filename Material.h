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

    glm::dvec3 diffuse;
    glm::dvec3 specular;
    glm::dvec3 reflection;
    glm::dvec3 transmission;
    glm::dvec3 emission;

    double shininess;
    double rindex;
    double alpha;
} Material;

#endif
