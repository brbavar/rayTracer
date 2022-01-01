#ifndef SHAPE_H
#define SHAPE_H

#include <string>
#include "ray.h"

struct Sphere {
    Matrix center;
    Color clr;
    double radius;

    Sphere();
    Sphere(const Sphere&);
    Sphere(Matrix c, double r, Color color);
    ~Sphere();

    double intersects(Ray);
    
    bool operator==(Sphere);
    bool operator!=(Sphere);
    void operator=(const Sphere&);
};

#endif