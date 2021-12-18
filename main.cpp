#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include "shape.h"


Color::Color() {}

Color::Color(double b, double g, double r, double o) {
    blue = b, green = g, red = r, opacity = o;
}

bool Color::operator==(Color other) {
    return blue == other.blue && green == other.green && red == other.red && opacity == other.opacity;
}

bool Color::operator!=(Color other) {
    return !(*this == other);
}

Matrix::Matrix() {}

Matrix::Matrix(const Matrix& orig) {
    std::vector<double> list;
    for(auto row : orig.entries) {
        for(double e : row)
            list.push_back(e);
        entries.push_back(list);
        list.clear();
    }
}

Matrix::Matrix(std::vector<std::vector<double> > items) {
    std::vector<double> list;
    for(auto row : items) {
        for(double i : row)
            list.push_back(i);
        entries.push_back(list);
        list.clear();
    }
}

Matrix::Matrix(double x, double y) {
    std::vector<double> coords;
    coords.push_back(x);
    coords.push_back(y);
    entries.push_back(coords);
}

Matrix::Matrix(double x, double y, double z) {
    std::vector<double> coords;
    coords.push_back(x);
    coords.push_back(y);
    coords.push_back(z);
    entries.push_back(coords);
}

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

/* Returns -1 if the matrix has more than one row of entries, in which case it cannot
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

/* Returns empty matrix if this matrix has "magnitude" <= 0, because in that 
   case this matrix has no magnitude with which we may normalize it. A magnitude
   is necessary for vector normalization, which is the kind of normalization we 
   are interested in. */
Matrix Matrix::normalize() {
    Matrix result = Matrix();
    double mag = this->mag();
    if(mag > 0) {
        std::vector<double> list;
        for(auto row : entries) {
            for(double e : row) {
                list.push_back(e / mag);
            }
            result.entries.push_back(list);
            list.clear();
        }
    }
    return result;
}

bool Matrix::rowsColsEq(Matrix other) {
    bool match = true;
    if(entries.size() != other.entries.size())
        match = false;
    for(int r = 0; r < entries.size(); r++) {
        if(!match || entries[r].size() != other.entries[r].size()) {
            match = false;
            break;
        }
    }
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
    for(auto row : entries) {
        for(double e : row)
            list.push_back(-e);
        result.entries.push_back(list);
        list.clear();
    }
    return result;
}

Matrix Matrix::operator*(double scalar) {
    Matrix result = Matrix();
    std::vector<double> list;
    for(int r = 0; r < entries.size(); r++) {
        for(double e : entries[r])
            list.push_back(e * scalar);
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

    /* Give every row in AB exactly as many entries as there are columns in B. That way,
       AB will have exactly as many columns as B. */
    std::vector<double> row(other.entries[0].size(), 0);

    // Give AB exactly as many rows as A has.
    for(int i = 0; i < entries.size(); i++)
        result.entries.push_back(row);

    if(entries[0].size() == other.entries.size())
        for(int r = 0; r < result.entries.size(); r++)
            for(int c = 0; c < result.entries[r].size(); c++)
                for(int n = 0; n < entries[0].size(); n++)
                    result.entries[r][c] += entries[r][n] * other.entries[n][c];
    return result;
}

/* Returns an empty matrix if the matrices being added do not meet the requirement of having 
   equal numbers of rows and columns */
Matrix Matrix::operator+(Matrix other) {
    Matrix result = Matrix();
    if(rowsColsEq(other)) {
        std::vector<double> list;
        int i = 0, j = 0;
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
    std::vector<double> list;
    for(auto row : orig.entries) {
        for(double e : row)
            list.push_back(e);
        entries.push_back(list);
        list.clear();
    }
}


Light::Light() {}

Light::Light(double x, double y, double z, Color color) {
    std::vector<double> coords;
    coords.push_back(x);
    coords.push_back(y); 
    coords.push_back(z);
    entries.push_back(coords);

    clr = color;
}

Ray::Ray() {}

/* Subtract coordinates of initial point from those of non-initial point,
    then normalize, and you get unit vector pointing same direction that ray points. */
Ray::Ray(Matrix initPt, Matrix nonInitPt) {
    std::vector<double> list;
    for(auto row : initPt.entries) {
        for(auto e : row)
            list.push_back(e);
        entries.push_back(list);
        list.clear();
    }
    unitDir = (nonInitPt - initPt).normalize();
}

bool Ray::operator==(Ray other) {
    return *this == other && unitDir == other.unitDir;
}

bool Ray::operator!=(Ray other) {
    return !(*this == other);
}

Camera::Camera() {
    focus = Matrix(0, 2, 0);
}

Camera::Camera(Ray ray) {
    std::vector<double> list;
    for(auto row : ray.entries) {
        for(auto e : row)
            list.push_back(e);
        entries.push_back(list);
        list.clear();
    }
    unitDir = ray.unitDir;
    right = Matrix(entries) + Matrix(1, 0, 0);
    focus = Matrix(0, 2, 0);
}

// Centered at origin of local coordinate system by default
Shape::Shape() {
    locCenter = Matrix(0, 0, 0);
}

Shape::Shape(const Shape& orig) {
    locCenter = orig.locCenter, globCenter = orig.globCenter;  // Some of these may have undefined vals,
    locRadius = orig.locRadius, globRadius = orig.globRadius;  // depending on the type of orig.
    locHeight = orig.locHeight, globHeight = orig.globHeight;  // MIGHT NEED TO CHANGE
    locWidth = orig.locWidth, globWidth = orig.globWidth;
    normal = orig.normal;
    type = orig.type;
    clr = orig.clr;
}

Shape::Shape(Color color) {
    clr = color;
}

double Shape::intersects(Ray ray) {
    if(type == "Sphere") {
        // Solving for coefficients derived from equation for sphere
        double a = ray.unitDir.dot(ray.unitDir);
        Matrix initPt = ray - globCenter;
        double b = 2 * initPt.dot(ray.unitDir);
        double c = initPt.dot(initPt) - pow(globRadius, 2);

        if(b * b - 4 * a * c > 0) {
            double dist = (globCenter - ray).dot(ray.unitDir);
            Matrix ptBtwn = ray + ray.unitDir * dist;
            double offset = sqrt( pow(globRadius, 2) - pow((globCenter - ptBtwn).mag(), 2) );
            dist -= offset;
            return dist;
        }
    }
    if(type == "Plane") {
        double dist = ray.entries[0][1] / (-ray.unitDir.entries[0][1]);
        return dist;
    }
    return -1;
}

bool Shape::operator==(Shape other) {
    if(type != other.type || clr != other.clr)
        return false;

    bool sameCenter = true, sameDims = true, sameNormal = true;

    if(type == "Sphere" || type == "Plane") {
        sameCenter = /* locCenter == other.locCenter && */ globCenter == other.globCenter;
    }
    if(type == "Sphere") {
        sameDims = /* locRadius == other.locRadius && */ globRadius == other.globRadius;
    }
    if(type == "Plane") {
        sameDims = /* locHeight == other.locHeight && locWidth == other.locWidth 
            && */ globHeight == other.globHeight && globWidth == other.globWidth;
        sameNormal = normal == other.normal;
    }

    return sameCenter && sameDims && sameNormal;
}

bool Shape::operator!=(Shape other) {
    return !(*this == other);
}

void Shape::operator=(const Shape& orig) {
    locCenter = orig.locCenter, globCenter = orig.globCenter;  // Some of these may have undefined vals,
    locRadius = orig.locRadius, globRadius = orig.globRadius;  // depending on the type of orig.
    locHeight = orig.locHeight, globHeight = orig.globHeight;  // MIGHT NEED TO CHANGE
    locWidth = orig.locWidth, globWidth = orig.globWidth;
    normal = orig.normal;
    type = orig.type;
    clr = orig.clr;
}

Sphere::Sphere() {
    type = "Sphere";
}

Sphere::Sphere(Matrix c, double r, Color color) {
    globCenter = c, globRadius = r; 
    clr = color, type = "Sphere";
}

Plane::Plane() {
    type = "Plane";
}

Plane::Plane(Ray r, Matrix c, double h, double w) {
    normal = r, globCenter = c, globHeight = h, globWidth = w;
    double aspectRatio = w / h;
    locHeight = (w > h ? 1 / aspectRatio : 1);
    locWidth = (w > h ? 1 : aspectRatio);
    type = "Plane";
}

Plane::Plane(Ray r, Matrix c, double h, double w, Color color) {
    normal = r, globCenter = c, globHeight = h, globWidth = w;
    double aspectRatio = w / h;
    locHeight = (w > h ? 1 / aspectRatio : 1);
    locWidth = (w > h ? 1 : aspectRatio);
    clr = color, type = "Plane";
}


void inputPicProps(int&, int&);
double findBrightness(Light, Matrix, Matrix, char*, Color&);
void render(Light, int, int, std::vector<Shape>);
void inputLight(Light&);
std::vector<Shape> inputShapes();


void inputPicProps(int& picHeight, int& picWidth) {
    std::cout << "Do you want the image to be 500 pixels tall and 700 pixels wide? (y/n) ";
    char ans = ' ';

    while(ans != 'y' && ans != 'Y' && ans != 'n' && ans != 'N') {
        std::cin >> ans;
        if(ans == 'n' || ans == 'N') {
            std::cout << "How tall should the image be, in pixels? (Min height = 20px, max = 2160px) ";
            int min = 20;
            int max = 2160;
            std::cin >> picHeight;
            while(std::cin.fail() || picHeight < min || picHeight > max) {
                std::cout << "Try again. You must enter an integer between " << min << " and " << max 
                    << " (" << min << ", " << min + 1 << ", " << min + 2 << ", . . . " << max - 1 << ", " << max << "). ";
                std::cin >> picHeight;
            }

            std::cout << "\nHow wide should the image be, in pixels? (Min width = 20px, max = 3840px) ";
            max = 3840;
            std::cin >> picWidth;
            while(std::cin.fail() || picHeight < min || picHeight > max) {
                std::cout << "Try again. You must enter an integer between " << min << " and " << max 
                    << " (" << min << ", " << min + 1 << ", " << min + 2 << ", . . . " << max - 1 << ", " << max << "). ";
                std::cin >> picWidth;
            }
        }
        else
            if(ans != 'y' && ans != 'Y') {
                std::cout << "Sorry, what was that? Your answer must be 'y' (yes) or 'n' (no). ";
                continue;
            }
    }
}

double findBrightness(Light light, Matrix nearestPt, Matrix normalDir, char* pixel, Color& clr) {
    Ray extendedLight = Ray(Matrix(light.entries), nearestPt);
    Matrix lightDir = extendedLight.unitDir;

    /* Get the angle between the direction of the light that reaches this point and the direction of the normal 
        vector to the shape at this point on its surface. That will determine the brightness of this point. */
    double angle = acos( normalDir.dot(lightDir) / (normalDir.mag() * lightDir.mag()) );

    /* If the angle between the light and the normal at this point is greater than 1.5708 radians, it is greater 
        than 90 degrees. That means this point is all the way on the other side of the shape, relative to the point 
        at which the light strikes it. So no light reaches this point. It will thus be fully shaded (black). 
        Otherwise, light reaches the point. So the color of the shape becomes the color of this pixel; the 
        darkness or lightness of the pixel's color is determined by multiplying each RGB component of the shape's 
        intrinsic color by the brightness factor calculated below. */
    if(angle > 1.5708) {
        pixel[0] = 0;
        pixel[1] = 0;
        pixel[2] = 0;
    }

    // The smaller the angle, the brighter the shape's surface at this point.
    double brightness = 1 - angle / 1.5708;

    clr.blue *= brightness;
    clr.green *= brightness;
    clr.red *= brightness;

    return brightness;
 }

void render(Light light, int picHeight, int picWidth, std::vector<Shape> shapes) {
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

    char pixel[3];

    /* Iterate over every pixel in the image, setting its color based on both the intrinsic properties of the shape
        and the lighting. */
    for(double z = - picHeight / 2; z < picHeight / 2; z++) {
        for(double x = - picWidth / 2; x < picWidth / 2; x++) {
            // Point camera at current pixel
            Camera cam = Camera(Ray(Matrix(0, 0, 0), Matrix(x, -100, z)));
            const double INF = std::numeric_limits<double>::infinity();
            double minDist = INF, dist;
            Shape nearestShape;

            for(Shape s : shapes) {
                dist = s.intersects(cam);
                if(dist != -1) {
                    minDist = std::min(minDist, dist);
                    if(minDist == dist)
                        nearestShape = s;
                }
            }
            
            if(minDist < INF) {
                Matrix nearestPt = cam + cam.unitDir * minDist;
                Matrix normalDir = Ray(nearestShape.globCenter, nearestPt).unitDir;
                Color clr = nearestShape.clr;
                double brightness = findBrightness(light, nearestPt, normalDir, pixel, clr);
                bool smallAngle = brightness > 0, unshaded = true;

                /* Cast a ray from the nearest point of intersection to the light source to see if any other shape
                   is in the way, preventing the light from reaching said point. */
                /* Ray shadow = Ray(nearestPt, light);
                for(Shape s : shapes) {
                    if(s != nearestShape) {
                        dist = s.intersects(shadow);
                        if(dist != -1)
                            unshaded = false;
                    }
                } */

                /* If angle is small enough for light to reach intersection point, and no shape casts shadow over that point, 
                   then give pixel same color as shape. */
                if(smallAngle /* && unshaded */) {
                    pixel[0] = floor(clr.blue * 255);
                    pixel[1] = floor(clr.green * 255);
                    pixel[2] = floor(clr.red * 255);
                    imgFile.write(pixel, 3);
                    continue;
                }
            }
            // If minDist is still infinity, meaning cam ray intersects no shape, then paint pixel black.
            pixel[0] = 0;
            pixel[1] = 0;
            pixel[2] = 0;
            imgFile.write(pixel, 3);
        }
    }

    imgFile.close();
}

void inputLight(Light& light) {
    std::cout << "\nDo you want to place the light source elsewhere than (" << light.entries[0][0] << ", " 
        << light.entries[0][1] << ", " << light.entries[0][2] << ")? (y/n) ";
    char ans = ' ';

    while(ans != 'y' && ans != 'Y' && ans != 'n' && ans != 'N') {
        std::cin >> ans;
        if(ans == 'y' || ans == 'Y') {
            std::cout << "\nWhere would you like to place the light source?\nx = ";
            std::cin >> light.entries[0][0];        
            std::cout << "y = ";
            std::cin >> light.entries[0][1];
            std::cout << "z = ";
            std::cin >> light.entries[0][2];

            std::cout << "Sounds good. What color should the light be? Give the RGB values. \n r = ";
            std::cin >> light.clr.red;
            std::cout << "g = ";
            std::cin >> light.clr.green;
            std::cout << "b = ";
            std::cin >> light.clr.blue;
        }
        else
            if(ans != 'n' && ans != 'N') {
                std::cout << "Sorry, what was that? Your answer must be 'y' (yes) or 'n' (no). ";
                continue;
            }
    }
}

std::vector<Shape> inputShapes() {
    std::vector<Shape> shapes;
    Matrix c = Matrix(0, 0, 0);
    double r;
    Color clr = Color();

    std::cout << "\nWhat type of shape would you like to add to the scene? (Options: sphere) ";
    std::string ans = "";

    while(ans != "sphere") {
        std::cin >> ans;
        if(ans == "sphere") {
            std::cout << "Nice choice. Where should the center of the sphere be?\nx = ";
            std::cin >> c.entries[0][0];
            std::cout << "y = ";
            std::cin >> c.entries[0][1];
            std::cout << "z = ";
            std::cin >> c.entries[0][2];

            std::cout << "All right. What should the radius of the sphere be?\nd = ";
            std::cin >> r;

            std::cout << "Splendid! What color should the sphere be painted? Give the RGB values. \nr = ";
            std::cin >> clr.red;
            std::cout << "g = ";
            std::cin >> clr.green;
            std::cout << "b = ";
            std::cin >> clr.blue;

            std::cout << "One more thing: How see-through should the sphere be, on a scale of 0 to 10 (0 being not "
                << "at all see-through, and 10 being as see-through as possible)? ";
            std::cin >> clr.opacity;
            clr.opacity = (10 - clr.opacity) / 10;

            shapes.push_back(Sphere(c, r, clr));
        }
        else {
            std::cout << "Sorry, what was that? Your answer must be \"sphere.\" ";
            continue;
        }
    }

    return shapes;
}

/* fresnelize() {

} */


int main() {
    int picHeight = 500, picWidth = 700;

    unsigned int seed = std::chrono::steady_clock::now().time_since_epoch().count();
    std::default_random_engine re(seed);

    std::normal_distribution<double> radDistro(500, 10);
    
    /* std::uniform_real_distribution<double> x1Distro(-949, 949);
    std::uniform_real_distribution<double> y1Distro(-749, 749); */
    
    std::uniform_real_distribution<double> xDistro(-picWidth / 2, picWidth);
    std::uniform_real_distribution<double> yDistro(50, 300); 
    std::uniform_real_distribution<double> zDistro(-picHeight / 2, picHeight);
    
    std::uniform_real_distribution<double> clrDistro(0, 1);
    
    std::uniform_int_distribution<int> shapeCountDistro(1, 10);
    std::uniform_int_distribution<int> floorDistro(0, 1);

    Color lightClr = Color(clrDistro(re), clrDistro(re), clrDistro(re), clrDistro(re));
    // Light light = Light(x1Distro(re), y1Distro(re), -zDistro(re), lightClr);
    Light light = Light(0, 0, 0, lightClr);

    std::vector<Shape> shapes;
    char ans = ' ';

    std::cout << "Would you like to design a custom image, or just generate a random one? (c, r) ";
    std::cin >> ans;
    while(ans != 'c' && ans != 'r') {
        std::cout << "You should answer either 'c' for custom or 'r' for random. ";
        std::cin >> ans;
    }
    if(ans == 'c') {
        inputPicProps(picHeight, picWidth);
        inputLight(light);
        shapes = inputShapes();
    }
    else {
        // int shapeCount = shapeCountDistro(re);
        int onFloor = 0 /* floorDistro(re) */;

        /* if(onFloor) {            
            Matrix center = Matrix(475, 0, 5);
            Matrix nonInitPt = Matrix(475, 0, 20);
            Ray normal = Ray(center, nonInitPt);
            double xDim = 950, zDim = 10;
            Color clr = Color(clrDistro(re), clrDistro(re), clrDistro(re), clrDistro(re));
            shapes.push_back(Plane(normal, center, zDim, xDim, clr));
            shapeCount++;
        } */

        // for(int i = 0; i < 1; i++) {
            /* double radius = -1;
            while(radius <= 0)
                radius = radDistro(re);
            Matrix center = Matrix(xDistro(re), yDistro(re), (onFloor ? radius : zDistro(re))); */
            Matrix center1 = Matrix(0, 400, 0);
            Matrix center2 = center1 + Matrix(-280, 0, 280);
            Matrix center3 = center1 - Matrix(-280, 0, 280);
            Color clr1 = Color(clrDistro(re), clrDistro(re), clrDistro(re), clrDistro(re));
            Color clr2 = Color(clrDistro(re), clrDistro(re), clrDistro(re), clrDistro(re));
            Color clr3 = Color(clrDistro(re), clrDistro(re), clrDistro(re), clrDistro(re));
            
            shapes.push_back(Sphere(center1, 280, clr1));
            shapes.push_back(Sphere(center2, 280, clr2));
            shapes.push_back(Sphere(center3, 280, clr3));

            /* std::cout << "x = " << center1.entries[0][0] << ", y = " << center1.entries[0][1] << ", z = " << center1.entries[0][2] 
                << ", r = " << 280 << std::endl;
            std::cout << "b = " << clr1.blue << ", g = " << clr1.green << ", r = " << clr1.red << ", o = " << clr1.opacity << std::endl; */
            // std::cout << "i = " << i << std::endl;
            // std::cout << std::endl;
        // }

        std::cout << (onFloor ? "Shapes are on floor." : "Shapes are free-floating.") << std::endl;
    }

    std::cout << "Light is at (" << light.entries[0][0] << ", " << light.entries[0][1] << ", " << light.entries[0][2] << ")." << std::endl;

    // testRender(light, picHeight, picWidth, shapes);
    render(light, picHeight, picWidth, shapes);
}