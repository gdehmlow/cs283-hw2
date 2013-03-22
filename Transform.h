/*
    Transform.h

    Static class to handle matrix transformations
*/

#ifndef I_TRNSF
#define I_TRNSF

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

class Transform {
    public:
    	Transform();
    	virtual ~Transform();
        static glm::dmat4 rotate(const double degrees, const glm::dvec3& axis);
        static glm::dmat4 scale(const double &sx, const double &sy, const double &sz); 
        static glm::dmat4 translate(const double &tx, const double &ty, const double &tz);
};

#endif
