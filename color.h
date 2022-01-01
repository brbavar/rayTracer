#ifndef COLOR_H
#define COLOR_H

struct Color {
    double blue, green, red;

    Color();
    Color(const Color&);
    Color(double, double, double);
    ~Color();

    bool operator==(Color);
    bool operator!=(Color);
};

#endif