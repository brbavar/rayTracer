#ifndef COLOR_H
#define COLOR_H

struct Color {
    double blue, green, red, opacity;

    Color();
    Color(const Color&);
    Color(double, double, double, double);
    ~Color();

    bool operator==(Color);
    bool operator!=(Color);
};

#endif