#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include "sphere.h"


Color::Color() {}

Color::Color(const Color& orig) {
    blue = orig.blue, green = orig.green, red = orig.red;
}

Color::Color(double b, double g, double r) {
    blue = b, green = g, red = r;
}

Color::~Color() {}

bool Color::operator==(Color other) {
    return blue == other.blue && green == other.green && red == other.red;
}

bool Color::operator!=(Color other) {
    return !(*this == other);
}

Matrix::Matrix() {}

Matrix::Matrix(const Matrix& orig) {
    std::vector<double> list;
    entries.reserve(orig.entries.size());
    list.reserve(orig.entries[0].size()); 

    for(auto row : orig.entries) {
        for(double entry : row)
            list.push_back(entry);
        entries.push_back(list);
        list.clear();
    }
}

Matrix::Matrix(std::vector<std::vector<double> > items) {
    std::vector<double> list;
    entries.reserve(items.size());
    list.reserve(items[0].size());
    for(auto row : items) {
        for(double item : row)
            list.push_back(item);
        entries.push_back(list);
        list.clear();
    }
}

Matrix::Matrix(double x, double y) {
    entries.reserve(1);
    std::vector<double> coords;
    coords.reserve(2);
    coords.push_back(x);
    coords.push_back(y);
    entries.push_back(coords);
}

Matrix::Matrix(double x, double y, double z) {
    entries.reserve(1);
    std::vector<double> coords;
    coords.reserve(3);
    coords.push_back(x);
    coords.push_back(y);
    coords.push_back(z);
    entries.push_back(coords);
}

Matrix::~Matrix() {}

/* Returns -1 if either matrix has more than one row of entries, in which case one of them 
   cannot be said to represent a point in n-dimensional space. They must both be capable of 
   representing points in space, so that a distance between the points may be calculated. 
   Of course, there is no such thing as a negative distance. */
double Matrix::distanceTo(Matrix other) {
    if(entries.size() == 1 && other.entries.size() == 1) {
        double sumSqrs = 0;
        for(int c = 0; c < entries[0].size(); c++)
            sumSqrs += pow(other.entries[0][c] - entries[0][c], 2);
        return sqrt(sumSqrs);
    }
    return -1;
}

/* Returns -1 if matrix does not have exactly one row of entries, in which case it cannot
   be said to represent a vector (i.e., list of numbers rather than a table of numbers), 
   which has a magnitude. Of course, there is no such thing as a negative magnitude. */
double Matrix::mag() {
    double sumSqrs = 0;
    if(entries.size() == 1) {
        for(double e : entries[0])
            sumSqrs += e * e;
        return sqrt(sumSqrs);
    }
    return -1;
}

/* Returns empty matrix if this matrix has "magnitude" <= 0, because in that case
   this matrix has no magnitude with which we can normalize it. A nonzero magnitude
   is necessary for vector normalization, which is the kind of normalization we 
   are interested in. */
Matrix Matrix::normalize() {
    Matrix result = Matrix();
    double mag = this->mag();
    if(mag > 0) {
        std::vector<double> list;
        list.reserve(entries[0].size());
        for(auto row: entries) {
            for(double entry : row)
                list.push_back(entry / mag);
            result.entries.push_back(list);
            list.clear();
        }
    }
    return result;
}

bool Matrix::rowsColsEq(Matrix other) {
    bool match = entries.size() == other.entries.size();
    for(int r = 0; match && r < entries.size(); r++)
        match = entries[r].size() == other.entries[r].size();
    return match;
}

/* The dot product requires that the two matrices being multiplied have the same number of entries. 
   If they do not, then in the case of this function, a garbage value of 0 will be returned by default. */
double Matrix::dot(Matrix other) {
    double sumProds = 0;
    if(rowsColsEq(other))
        for(int r = 0; r < entries.size(); r++)
            for(int c = 0; c < entries[r].size(); c++)
                sumProds += entries[r][c] * other.entries[r][c];
    return sumProds;
}

/* For our purposes, and by convention, the cross product is only defined for coordinate vectors
   with exacty three entries. If this matrix and the other do not meet those criteria, then an
   empty matrix is returned. */
Matrix Matrix::cross(Matrix other) {
    Matrix result = Matrix();
    std::vector<double> list;
    if(entries.size() == 1 && other.entries.size() == 1 && 
        entries[0].size() == 3 && other.entries[0].size() == 3) {
        list.reserve(3);
        for(int i = 0; i < 3; i++)
            list.push_back(entries[0][i == 2 ? 0 : i + 1] * other.entries[0][i == 0 ? 2 : i - 1] 
                - entries[0][i == 0 ? 2 : i - 1] * other.entries[0][i == 2 ? 0 : i + 1]);
        result.entries.push_back(list);
    }
    return result;
}

Matrix Matrix::negate() {
    Matrix result = Matrix();
    std::vector<double> list;
    list.reserve(entries[0].size());
    for(auto row : entries) {
        for(double entry : row)
            list.push_back(-entry);
        result.entries.push_back(list);
        list.clear();
    }
    return result;
}

Matrix Matrix::operator*(double scalar) {
    Matrix result = Matrix();
    std::vector<double> list;
    list.reserve(entries[0].size());
    for(auto row : entries) {
        for(double entry : row)
            list.push_back(entry * scalar);
        result.entries.push_back(list);
        list.clear();
    }
    return result;
}

/* Call this matrix A and the other matrix B. In order to be multiplied, A must have exactly as
   many columns as B does rows. If this condition is unsatisfied, an empty matrix is returned. 
   Multiplying A and B gives us AB, which has exactly as many rows as A and exactly as many 
   columns as B. */
Matrix Matrix::operator*(Matrix other) {
    Matrix result = Matrix();
    if(entries[0].size() == other.entries.size()) {
        // Make room in AB for exactly as many rows as A has.
        result.entries.reserve(entries.size());

        /* Give every row in AB exactly as many entries as there are columns in B. 
        That way, AB will have exactly as many columns as B. */
        std::vector<double> row(other.entries[0].size(), 0);

        for(int r = 0; r < entries.size(); r++) {
            result.entries.push_back(row);
            for(int c = 0; c < result.entries[r].size(); c++)
                for(int n = 0; n < entries[0].size(); n++)
                    result.entries[r][c] += entries[r][n] * other.entries[n][c];
        }
    }
    return result;
}

/* Returns an empty matrix if the matrices being added do not meet the requirement of having 
   equal numbers of rows and columns */
Matrix Matrix::operator+(Matrix other) {
    Matrix result = Matrix();
    if(rowsColsEq(other)) {
        std::vector<double> list;
        list.reserve(entries[0].size());
        for(int r = 0; r < entries.size(); r++) {
            for(int c = 0; c < entries[r].size(); c++)
                list.push_back(entries[r][c] + other.entries[r][c]);
            result.entries.push_back(list);
            list.clear();
        }
    }
    return result;
}

Matrix Matrix::operator-(Matrix other) {
    return *this + other.negate();
}

bool Matrix::operator==(Matrix other) {
    if(rowsColsEq(other)) {
        for(int r = 0; r < entries.size(); r++)
            for(int c = 0; c < entries[r].size(); c++)
                if(entries[r][c] != other.entries[r][c])
                    return false;
        return true;
    }
    else 
        return false;
}

bool Matrix::operator!=(Matrix other) {
    return !(*this == other);
}

void Matrix::operator=(const Matrix& orig) {
    entries.clear();
    entries.reserve(orig.entries.size());
    std::vector<double> list;
    list.reserve(orig.entries[0].size());
    for(auto row : orig.entries) {
        for(double entry : row)
            list.push_back(entry);
        entries.push_back(list);
        list.clear();
    }
}


Light::Light() {}

Light::Light(double x, double y, double z, double b) {
    entries.reserve(1);
    std::vector<double> coords;
    coords.reserve(3);
    coords.push_back(x);
    coords.push_back(y); 
    coords.push_back(z);
    entries.push_back(coords);

    boost = b;
}

Ray::Ray() {}

/* Subtract coordinates of initial point from those of non-initial point,
    then normalize, and you get unit vector pointing same direction that ray points. */
Ray::Ray(Matrix initPt, Matrix nonInitPt) {
    entries.reserve(initPt.entries.size());
    std::vector<double> coords;
    coords.reserve(initPt.entries[0].size());
    for(auto row : initPt.entries) {
        for(double entry : row)
            coords.push_back(entry);
        entries.push_back(coords);
        coords.clear();
    }
    unitDir = (nonInitPt - initPt).normalize();
}

bool Ray::operator==(Ray other) {
    Matrix A = Matrix((*this).entries);
    Matrix B = Matrix(other.entries);
    return A == B && unitDir == other.unitDir;
}

bool Ray::operator!=(Ray other) {
    return !(*this == other);
}

Camera::Camera() {}

Camera::Camera(Ray ray) {
    entries.reserve(ray.entries.size());
    std::vector<double> coords;
    coords.reserve(ray.entries[0].size());
    for(auto row : ray.entries) {
        for(double entry : row)
            coords.push_back(entry);
        entries.push_back(coords);
    }
    unitDir = ray.unitDir;
}

Sphere::Sphere() {}

Sphere::Sphere(const Sphere& orig) {
    center = orig.center, radius = orig.radius;
    clr = orig.clr;
}

Sphere::Sphere(Matrix c, double r, Color color) {
    center = c, radius = r;
    clr = color;
}

Sphere::~Sphere() {}

double Sphere::intersects(Ray ray) {
    // Solving for coefficients derived from equation for sphere
    double a = ray.unitDir.dot(ray.unitDir);
    Matrix initPt = ray - center;
    double b = 2 * initPt.dot(ray.unitDir);
    double c = initPt.dot(initPt) - pow(radius, 2);
    double dist;

    if(b * b - 4 * a * c > 0) {
        dist = (center - ray).dot(ray.unitDir);
        Matrix ptBtwn = ray + ray.unitDir * dist;
        double offset = sqrt( pow(radius, 2) - pow((center - ptBtwn).mag(), 2) );
        dist -= offset;
        if(dist == -1)
            return -1.1;
        return dist;
    }
    return -1;
}

bool Sphere::operator==(Sphere other) {
    bool sameCenter = true, sameDims = true, sameNormal = true;

    sameCenter = center == other.center;
    sameDims = radius == other.radius;

    return sameCenter && sameDims && sameNormal;
}

bool Sphere::operator!=(Sphere other) {
    return !(*this == other);
}

void Sphere::operator=(const Sphere& orig) {
    center = orig.center, radius = orig.radius;
    clr = orig.clr;
}


double setBrightness(Light, Matrix, Matrix, char*, Color&);
void render(Light, int, int, std::vector<Sphere*>);


double setBrightness(Light light, Matrix nearestPt, Matrix normalDir, char* pixel, Color& clr) {
    double brightness = 1, angle;

    Ray extendedLight = Ray(nearestPt, Matrix(light.entries));

    // Get direction of light beam aimed at point of intersection.
    Matrix lightDir = extendedLight.unitDir;

    /* Get the angle between the direction of the light that reaches this point and the direction of the normal 
        vector to the sphere at this point on its surface. That will determine the brightness of this point. */
    if(lightDir.mag() > 0) {
        angle = acos( normalDir.dot(lightDir) / (normalDir.mag() * lightDir.mag()) );
    }
    else
        return brightness;

    /* If the angle between the light and the normal at this point is greater than 1.5708 radians, it is greater 
        than 90 degrees. That means this point is all the way on the other side of the sphere, relative to the point 
        at which the light strikes it. So no light reaches this point. It will thus be fully shaded (black). 
        Otherwise, light reaches the point. So the color of the sphere becomes the color of this pixel; the 
        darkness or lightness of the pixel's color is determined by multiplying each RGB component of the sphere's 
        intrinsic color by the brightness factor calculated below. The smaller the angle, the brighter the sphere's 
        surface at this point. */
    brightness -= angle / 1.5708 - angle * light.boost;

    if(brightness < 0)
        brightness = 0;
    if(brightness > 1)
        brightness = 1;

    clr.blue *= brightness;
    clr.green *= brightness;
    clr.red *= brightness;

    return brightness;
}

void render(Light light, int picHeight, int picWidth, std::vector<Sphere*> spheres) {
    int area = picWidth * picHeight;
    int size = 4 * area;
    int fileSize = size + 54;

    // Dots per inch
    int dpi = 50;

    // Pixels per dot
    int ppd = 39.375;

    // px / in = (dots / in) * (px / dot)
    int ppi = dpi * ppd;

    char header[14] = {'B', 'M', 0, 0,   0, 0, 0, 0,   0, 0, 54, 0,   0, 0};

    header[2] = (char)fileSize;
    header[3] = (char)(fileSize>>8);
    header[4] = (char)(fileSize>>16);
    header[5] = (char)(fileSize>>24);

    char dibHeader[40] = {40, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   1, 0, 24, 0};

    dibHeader[4] = (char)picWidth;
    dibHeader[5] = (char)(picWidth>>8);
    dibHeader[6] = (char)(picWidth>>16);
    dibHeader[7] = (char)(picWidth>>24);

    dibHeader[8] = (char)picHeight;
    dibHeader[9] = (char)(picHeight>>8);
    dibHeader[10] = (char)(picHeight>>16);
    dibHeader[11] = (char)(picHeight>>24);

    dibHeader[21] = (char)size;
    dibHeader[22] = (char)(size>>8);
    dibHeader[23] = (char)(size>>16);
    dibHeader[24] = (char)(size>>24);

    dibHeader[25] = (char)ppi;
    dibHeader[26] = (char)(ppi>>8);
    dibHeader[27] = (char)(ppi>>16);
    dibHeader[28] = (char)(ppi>>24);

    dibHeader[29] = (char)ppi;
    dibHeader[30] = (char)(ppi>>8);
    dibHeader[31] = (char)(ppi>>16);
    dibHeader[32] = (char)(ppi>>24);  
    
    std::ofstream imgFile;
    imgFile.open("result.bmp");
    imgFile.write(header, 14);
    imgFile.write(dibHeader, 40);

    char pixel[] = {0,0,0};

    /* Iterate over every pixel in the image, setting its color based on both the intrinsic properties of the sphere
       and the lighting. */
    for(double y = 0; y < picHeight; y++)
        for(double x = 0; x < picWidth; x++) {
            Ray cam = Ray(Matrix(0, 210, -1000), Matrix(x, y, 0));  // Point camera at current pixel sample
            
            const double INF = std::numeric_limits<double>::infinity();
            double minDist = INF, dist;
            Sphere nearestSphere;

            for(Sphere* s : spheres) {
                dist = (*s).intersects(cam);
                if(dist != -1) {
                    Matrix nearestPt = cam + cam.unitDir * minDist;
                    minDist = std::min(minDist, dist);
                    if(minDist == dist)
                        nearestSphere = *s;
                }
            }

            /* If minDist is less than infinity, the color of the pixel will depend on whether light reaches the sphere at the point of 
                intersection. It will also depend on the angle between the light and normal at that intersection. On the other hand, if 
                minDist is still infinity, meaning cam ray intersects no sphere, then pixel's color remains unchanged. */
            if(minDist < INF) {
                Matrix nearestPt = cam + cam.unitDir * minDist;
                Matrix center = nearestSphere.center;
                Matrix normal = nearestPt - center;
                Matrix normalDir = normal.normalize();

                Color clr = nearestSphere.clr;
                double brightness = setBrightness(light, nearestPt, normalDir, pixel, clr);
                bool lit = brightness > 0;

                // If angle is small enough for light to reach intersection point, then give pixel same color as sphere.
                if(lit) {
                    pixel[0] = floor(clr.blue * 255);
                    pixel[1] = floor(clr.green * 255);
                    pixel[2] = floor(clr.red * 255);
                    imgFile.write(pixel, 3);
                    continue;
                }
            }
            imgFile.write(pixel, 3);
            for(int j = 0; j < 3; j++)
                pixel[j] = 0;
        }

    imgFile.close();
}

int main() {
    unsigned int seed = std::chrono::steady_clock::now().time_since_epoch().count();
    std::default_random_engine re(seed);
    std::uniform_real_distribution<double> clrDistro(0, 1);

    Color lightClr = Color(clrDistro(re), clrDistro(re), clrDistro(re));
    Light light = Light(400, 500, -300, 0);

    std::vector<Sphere*> spheres;
    spheres.reserve(2);

    Matrix center = Matrix(350, 210, 0);
    Color clr = Color(clrDistro(re), clrDistro(re), clrDistro(re));
    Sphere* ptr = new Sphere(center, 70, clr);
    spheres.push_back(ptr);

    center = Matrix(350, -5500, 0);
    clr = Color(clrDistro(re), clrDistro(re), clrDistro(re));
    ptr = new Sphere(center, 5640, clr);
    spheres.push_back(ptr);

    int picHeight = 500, picWidth = 700;

    render(light, picHeight, picWidth, spheres);

    for(Sphere* ptr : spheres) // Deallocating here may not be necessary
        delete ptr;
}