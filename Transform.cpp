/*
    Transformation.cpp
*/

#include "Transform.h"

glm::dmat4 Transform::rotate(const double degrees, const glm::dvec3& axis) {
    double cos0 = cos(degrees*M_PI/180.0);
    double sin0 = sin(degrees*M_PI/180.0);
    glm::dvec3 axis2 = glm::normalize(axis);
    double ux = axis2[0];
    double uy = axis2[1];
    double uz = axis2[2];
    glm::dmat3 rmat = glm::dmat3(cos0 + ux*ux*(1-cos0), ux*uy*(1-cos0) - uz*sin0, ux*uz*(1-cos0) + uy*sin0,
                                uy*ux*(1-cos0) + uz*sin0, cos0 + uy*uy*(1-cos0), uy*uz*(1-cos0) - ux*sin0,
                                uz*ux*(1-cos0) - uy*sin0, uz*uy*(1-cos0) + ux*sin0, cos0 + uz*uz*(1-cos0)
                               );
    glm::dmat4 rotmat = glm::dmat4(rmat[0][0],rmat[0][1],rmat[0][2],0,
                       rmat[1][0],rmat[1][1],rmat[1][2],0,  
                       rmat[2][0],rmat[2][1],rmat[2][2],0,  
                       0,0,0,1);
    return rotmat;
}

glm::dmat4 Transform::scale(const double &sx, const double &sy, const double &sz) 
{
    return glm::dmat4(sx, 0, 0, 0,
                0, sy, 0, 0,
                0, 0, sz, 0,
                0, 0, 0,  1);
}

glm::dmat4 Transform::translate(const double &tx, const double &ty, const double &tz) 
{
    return glm::dmat4(1, 0, 0, tx,
                0, 1, 0, ty,
                0, 0, 1, tz,
                0, 0, 0, 1);
}

Transform::Transform() { }
Transform::~Transform() { }
