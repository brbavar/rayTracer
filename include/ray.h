#ifndef RAY_H
#define RAY_H

#include "matrix.h"

struct Ray : Matrix {
    Matrix unitDir;

    Ray();
    Ray(Matrix, Matrix);

    bool operator==(Ray);
    bool operator!=(Ray);
};

struct Camera : Ray {
    Camera();
    Camera(Ray);
};

#endif