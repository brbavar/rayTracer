#ifndef MATRIX_H
#define MATRIX_H

#include "color.h"

struct Matrix3D {
    double x, y, z;

    Matrix3D();
    Matrix3D(const Matrix3D&);
    Matrix3D(double, double, double);

    double distanceTo(Matrix3D);

    // Returns magnitude of vector represented by matrix
    double mag();

    // Returns dot product of this vector with another
    double dot(Matrix3D);

    // * serves as both scalar multiplication and cross product operator
    Matrix3D operator*(double);
    Matrix3D operator*(Matrix3D);
    Matrix3D operator+(Matrix3D);
    Matrix3D operator-(Matrix3D);

    Matrix3D normalize();

    Matrix3D invert();
};

struct Light : Matrix3D {
    Color clr;

    Light();
    Light(double, double, double, Color);
};

#endif