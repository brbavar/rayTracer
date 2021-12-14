#ifndef RAY_H
#define RAY_H

#include "matrix.h"

struct Ray : Matrix3D {
    Matrix3D unitDir;
    
    Ray();
    Ray(Matrix3D, Matrix3D);
};

struct Camera : Ray {
    Camera();
    Camera(Ray);
};

#endif