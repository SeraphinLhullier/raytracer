#include "sphere.h"
#include <cmath>

/*
Prend en entrée un rayon, trois vecteurs et un double.
Renvoie true s'il a une intersection entre la sphere et le rayon.
Si c'est le cas, les entrées sont changées pour stocker le point d'intersection, la normale à ce point,
l'albedo de ce point ainsi que la distance e ce point par rapport au point de départ du rayon
*/


bool Sphere::intersect(const Ray& r, Vector& P, Vector& N, double& t, Vector& albedo) {
    double a = 1;
    double b = 2 * dot(r.u, r.C - O);
    double c = (r.C - O).sqrNorm() - R * R;

    double delta = b * b - 4 * a * c;
    if (delta < 0) return false;

    double sqDelta = sqrt(delta);
    double t2 = (-b + sqDelta) / (2 * a);
    if (t2 < 0) 
        return false;

    double t1 = (-b - sqDelta) / (2 * a);
    if (t1 > 0)
        t = t1;
    else
        t = t2;

    P = r.C + t * r.u;
    N = (P - O).getNormalized();

    albedo = this->albedo;
    return true;
}
