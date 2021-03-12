#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include "sphere.h"
#include "ray.h"
#include "vector.h"
#include "TriangleMesh.h"
#include <cmath>
#include <random>
#define M_PI 3.14159265358979323846


static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0, 1);


class Scene
{
public:
    std::vector<Object*> objects;    
    Vector L;
    double I;
    Scene();
    bool intersect(const Ray& r, Vector& P, Vector& N, Vector& albedo, bool& mirror, bool& transparent, double& t, int& objectid);
    Vector getColor(const Ray& r, int rebond, bool lastdefuse);
};


#endif // SCENE_H
