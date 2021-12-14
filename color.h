#ifndef COLOR_H
#define COLOR_H

struct Color {
    double blue, green, red, opacity;

    Color();
    Color(double, double, double, double);
};

#endif