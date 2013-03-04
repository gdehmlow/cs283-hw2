/*
    Camera.cpp
*/

#include <iostream>
#include "Camera.h"

typedef glm::mat3 mat3;
typedef glm::mat4 mat4; 
typedef glm::vec3 vec3; 
typedef glm::vec4 vec4; 

Camera::Camera()
{
    width = 640;
    height = 480;
    w2 = 640.0;
    h2 = 480.0;
}

void Camera::init(float* values)
{
    eye = vec3(values[0],values[1],values[2]);
    vec3 lookAt = vec3(values[3],values[4],values[5]);
    vec3 up = vec3(values[6],values[7],values[8]);

    w = glm::normalize(eye - lookAt);
    u = glm::normalize(glm::cross(up,w));
    v = glm::cross(w,u);

    fovy = values[9] * M_PI / 180.0;
    tanfovy = tan(fovy / 2.0);
}

void Camera::setWidthAndHeight(int width, int height)
{
    this->width = width;
    this->height = height;
    this->w2 = (float)width/2.0;
    this->h2 = (float)height/2.0;
    tanfovx = tanfovy * w2 / h2;
}

void Camera::generateRay(Ray* ray, float x, float y)
{
    x = x + 0.5;
    y = y + 0.5;
    float alpha = tanfovx * (x - w2) / w2;
    float beta = tanfovy * (y - h2) / h2;
    vec3 tempVec = glm::normalize(alpha*u + beta*v - w);
    ray->position = vec4(eye,1);
    ray->direction = vec4(tempVec,0);
}
