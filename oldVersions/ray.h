#ifndef RAY_H
#define RAY_H

#include "oldMatrix.h"

struct Ray : Matrix {
    Matrix unitDir;

    Ray();
    Ray(Matrix, Matrix);

    bool operator==(Ray);
    bool operator!=(Ray);
};

struct Camera : Ray {
    double picDist = 1;
    Matrix right;
    Matrix focus;

    Camera();
    Camera(Ray);
};

#endif