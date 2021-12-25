#ifndef SHAPE_H
#define SHAPE_H

#include <string>
#include "ray.h"

struct Shape {
    std::string type = "";
    Matrix center;
    Color clr;
    Ray normal;
    double radius, height, width;

    Shape();
    Shape(const Shape&);
    Shape(Color);
    Shape(Matrix c, double r, Color color); // Sphere

    double intersects(Ray);
    
    bool operator==(Shape);
    bool operator!=(Shape);
    void operator=(const Shape&);
};

struct Sphere : Shape {
    Sphere();
    Sphere(Matrix, double, Color);
};

struct Plane : Shape {
    Plane();
    Plane(Ray, Matrix, double, double);
    Plane(Ray, Matrix, double, double, Color);
};

#endif