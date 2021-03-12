#ifndef SPHERE_H
#define SPHERE_H

#include "vector.h"
#include "ray.h"
#include "Object.h"

class Sphere : public Object
{
private:
    Vector O;
    double R;

public:
    Sphere(const Vector& O, double R, Vector albedo, bool isMirror = false, bool isTransparent = false) : O(O), R(R) { 
        this->albedo = albedo; 
        this->isMirror = isMirror; 
        this->isTransparent = isTransparent; 
    };

    bool intersect(const Ray& r, Vector& P, Vector& N, double& t, Vector& albedo);

    Vector getO() {
        return O;
    };

    double getR() {
        return R;
    };
};

#endif // SPHERE_H
