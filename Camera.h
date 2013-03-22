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
        void init(double* values);
        void setWidthAndHeight(const int width, const int height);
        void generateRay(Ray& ray, const double x, const double y, const double t);

    private:
        glm::dvec3 eye;
        glm::dvec3 w;
        glm::dvec3 u;
        glm::dvec3 v;

        double fovy;
        double tanfovy;
        double tanfovx;
        double width;
        double w2;
        double height;
        double h2;
};

#endif
