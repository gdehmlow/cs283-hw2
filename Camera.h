/*
    Camera.h

    Generates rays.
*/

#ifndef I_CAMERA
#define I_CAMERA

#include <glm/glm.hpp>
#include "Ray.h"

class Camera {
    public:
        Camera();
        void init(float* values);
        void setWidthAndHeight(const int width, const int height);
        void generateRay(Ray& ray, const float x, const float y);

    private:
        glm::vec3 eye;
        glm::vec3 w;
        glm::vec3 u;
        glm::vec3 v;

        float fovy;
        float tanfovy;
        float tanfovx;
        float width;
        float w2;
        float height;
        float h2;
};

#endif
