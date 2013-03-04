/*
    Transformation.cpp
*/

#include "Transform.h"

mat4 Transform::rotate(const float degrees, const vec3& axis) {
    float cos0 = cos(degrees*M_PI/180.0);
    float sin0 = sin(degrees*M_PI/180.0);
    vec3 axis2 = glm::normalize(axis);
    float ux = axis2[0];
    float uy = axis2[1];
    float uz = axis2[2];
    mat3 rmat = mat3(cos0 + ux*ux*(1-cos0), ux*uy*(1-cos0) - uz*sin0, ux*uz*(1-cos0) + uy*sin0,
                                uy*ux*(1-cos0) + uz*sin0, cos0 + uy*uy*(1-cos0), uy*uz*(1-cos0) - ux*sin0,
                                uz*ux*(1-cos0) - uy*sin0, uz*uy*(1-cos0) + ux*sin0, cos0 + uz*uz*(1-cos0)
                               );
    mat4 rotmat = mat4(rmat[0][0],rmat[0][1],rmat[0][2],0,
                       rmat[1][0],rmat[1][1],rmat[1][2],0,  
                       rmat[2][0],rmat[2][1],rmat[2][2],0,  
                       0,0,0,1);
    return rotmat;
}

mat4 Transform::scale(const float &sx, const float &sy, const float &sz) 
{
    return mat4(sx, 0, 0, 0,
                0, sy, 0, 0,
                0, 0, sz, 0,
                0, 0, 0,  1);
}

mat4 Transform::translate(const float &tx, const float &ty, const float &tz) 
{
    return mat4(1, 0, 0, tx,
                0, 1, 0, ty,
                0, 0, 1, tz,
                0, 0, 0, 1);
}

Transform::Transform() { }
Transform::~Transform() { }
