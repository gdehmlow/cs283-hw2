/*
    Transform.h

    Static class to handle matrix transformations
*/

#ifndef I_TRNSF
#define I_TRNSF

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

typedef glm::mat3 mat3;
typedef glm::mat4 mat4; 
typedef glm::vec3 vec3; 
typedef glm::vec4 vec4; 

class Transform {
    public:
    	Transform();
    	virtual ~Transform();
        static mat4 rotate(const float degrees, const vec3& axis);
        static mat4 scale(const float &sx, const float &sy, const float &sz); 
        static mat4 translate(const float &tx, const float &ty, const float &tz);
};

#endif
