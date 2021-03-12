#ifndef VECTOR_H
#define VECTOR_H


class Vector
{
private:
    double coords[3];

public:
    explicit Vector(double x = 0, double y = 0, double z = 0);
    double operator[](int i) const;
    double& operator[](int i);
    double sqrNorm();
    Vector getNormalized();
    Vector getRotatedVector(const Vector& k, float theta);
};

Vector operator+(const Vector& a, const Vector& b);
Vector operator-(const Vector& a, const Vector& b);
Vector operator-(const Vector& a);
Vector operator*(const Vector& a, double b);
Vector operator*(double a, const Vector& b);
Vector operator*(const Vector& a, const Vector& b);
Vector operator/(const Vector& a, double b);
Vector operator/(double a, const Vector& b);
double dot(const Vector& a, const Vector& b);
double sqr(const double x);

Vector cross(const Vector& a, const Vector& b);


#endif // VECTOR_H
