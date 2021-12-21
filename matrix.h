#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
// #include "color.h"

struct Matrix {
    std::vector<std::vector<double> > entries;

    Matrix();
    Matrix(const Matrix&);
    Matrix(std::vector<std::vector<double> >);  // Constructs matrices of arbitrary height and width
    Matrix(double, double);  // Constructs 2D matrices
    Matrix(double, double, double);  // Constructs 3D matrices
    ~Matrix();

    /* Returns distance between point at coordinates in this->entries 
       and point at coords in member vector, entries, of matrix argument */
    double distanceTo(Matrix);

    // Returns magnitude of vector represented by matrix
    double mag();

    Matrix normalize();

    bool rowsColsEq(Matrix);

    // Returns dot product of this vector with another
    double dot(Matrix);

    // Returns cross product of this vector with another
    Matrix cross(Matrix);

    Matrix negate();

    // * overloaded for both scalar and matrix multiplication
    Matrix operator*(double);
    Matrix operator*(Matrix);
    Matrix operator+(Matrix);
    Matrix operator-(Matrix);
    bool operator==(Matrix);
    bool operator!=(Matrix);
    void operator=(const Matrix&);
};

/* struct Light : Matrix {
    Color clr;

    Light();
    Light(double, double, double, Color);
}; */

#endif