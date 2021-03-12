#ifndef OBJECT_H
#define OBJECT_H

#include "ray.h"
#include "vector.h"

class Object
{
public:
	Vector albedo;
	bool isMirror, isTransparent;

	Object() {};
	virtual bool intersect(const Ray& r, Vector& P, Vector& normale, double& t, Vector& albedo) = 0;
};

#endif