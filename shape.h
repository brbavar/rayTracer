#ifndef SHAPE_H
#define SHAPE_H

#include <string>
#include "ray.h"

struct Shape3D {
    Matrix3D center;
    double radius, height, width;
    Ray normal;
    Color clr;
    std::string type = "";

    Shape3D();
    Shape3D(Color);

    std::vector<Matrix3D> intersects(Ray);
};

struct Sphere : Shape3D {
    Sphere();
    Sphere(Matrix3D, double, Color);
};

struct Plane : Shape3D {
    Plane();
    Plane(Ray, Matrix3D, double, double);
};

#endif