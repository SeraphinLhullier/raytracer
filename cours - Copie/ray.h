#ifndef RAY_H
#define RAY_H

#include "vector.h"

class Ray
{
public:
    Ray(const Vector& C, const Vector& u) : C(C), u(u) {};
    Vector C, u;
};

#endif // RAY_H
