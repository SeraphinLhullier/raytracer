#pragma once

#include "vector.h"
#include "ray.h"
#include <algorithm>

class BoundingBox
{
public:
	Vector mini, maxi;
	bool intersect(const Ray& r, double& t);
};

class Noeud {
public:
	Noeud* fg;
	Noeud* fd;

	BoundingBox b;

	int debut;
	int fin;
};