#include "BoundingBox.h"

/*
Prend en entrée un rayon et un double.
Renvoie true s'il a une intersection entre la boite et le rayon.
Si c'est le cas, t est changé en la distance entre le point d'intersection en le point de départ du rayon
*/
bool BoundingBox::intersect(const Ray& r, double &t) {

	double t1x = (mini[0] - r.C[0]) / r.u[0];
	double t2x = (maxi[0] - r.C[0]) / r.u[0];
	double txMin = std::min(t1x, t2x);
	double txMax = std::max(t1x, t2x);

	double t1y = (mini[1] - r.C[1]) / r.u[1];
	double t2y = (maxi[1] - r.C[1]) / r.u[1];
	double tyMin = std::min(t1y, t2y);
	double tyMax = std::max(t1y, t2y);

	double t1z = (mini[2] - r.C[2]) / r.u[2];
	double t2z = (maxi[2] - r.C[2]) / r.u[2];
	double tzMin = std::min(t1z, t2z);
	double tzMax = std::max(t1z, t2z);

	double tMax = std::min(txMax, std::min(tyMax, tzMax));
	double tMin = std::max(txMin, std::max(tyMin, tzMin));

	if (tMax < 0) {
		return false;
	};

	t = tMin;

	return tMax > tMin;


}