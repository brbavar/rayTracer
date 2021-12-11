#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

struct Color {
    double blue, green, red;

    Color() {}

    Color(double b, double g, double r) {
        blue = b, green = g, red = r;
    }
};

struct Matrix3D {
    double x, y, z;

    Matrix3D() {}

    Matrix3D(const Matrix3D& orig) {
        x = orig.x, y = orig.y, z = orig.z;
    }

    Matrix3D(double comp1, double comp2, double comp3) {
        x = comp1, y = comp2, z = comp3;
    }

    double distanceTo(Matrix3D other) {
        return sqrt( pow((other.x - x), 2) + pow((other.y - y), 2) + pow((other.z - z), 2) );
    }

    double getMagnitude() {
        return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    }

    double getDotProduct(Matrix3D other) {
        return x * other.x + y * other.y + z * other.z;
    }

    // * serves as both scalar multiplication and cross product operator
    Matrix3D operator*(double scalar) {
        return Matrix3D(x * scalar, y * scalar, z * scalar);
    }

    Matrix3D operator*(Matrix3D other) {
        return Matrix3D(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
    }

    Matrix3D operator+(Matrix3D other) {
        return Matrix3D(x + other.x, y + other.y, z + other.z);
    }

    Matrix3D operator-(Matrix3D other) {
        return Matrix3D(x - other.x, y - other.y, z - other.z);
    }

    Matrix3D normalize() {
        double mag = getMagnitude();
        return Matrix3D(x / mag, y / mag, z / mag);
    }

    Matrix3D invert() {
        return Matrix3D(-x, -y, -z);
    }
};

struct ClrMatrix : Color, Matrix3D {
    double shine;

    ClrMatrix() {}

    ClrMatrix(Color clr, double s) {
        blue = clr.blue, green = clr.green, red = clr.red, shine = s;
        x = blue, y = green, z = red;
    }
};

struct Light : Matrix3D {
    Light() {}

    Light(double comp1, double comp2, double comp3) {
        x = comp1, y = comp2, z = comp3;
    }
};

struct Shape3D {
    ClrMatrix clr;

    Shape3D() {}

    Shape3D(ClrMatrix color) {
        clr = color;
    }
};

struct Sphere : Shape3D {
    Matrix3D center;
    double radius;

    Sphere() {}

    Sphere(Matrix3D c, double r, ClrMatrix color) {
        center.x = c.x, center.y = c.y, center.z = c.z;
        radius = r;
        clr = color;
    }
};

/* 
namespace pointKnown {
    struct Ray : Matrix3D {
        Matrix3D unitDir;
        
        Ray() {}

        // Subtract coordinates of initial point from those of non-initial point,
        // then normalize, and you get unit vector pointing same direction that ray points
        Ray(Matrix3D initPt, Matrix3D nonInitPt) {
            x = initPt.x, y = initPt.y, z = initPt.z;
            unitDir = Matrix3D(nonInitPt.x - x, nonInitPt.y - y, nonInitPt.z - z).normalize();
        }
    };
}

namespace dirKnown {
    struct Ray : Matrix3D {
        Matrix3D unitDir;

        Ray() {}

        Ray(Matrix3D initPt, Matrix3D dir) {
            x = initPt.x, y = initPt.y, z = initPt.z;
            unitDir = dir;
        }

        Ray operator+(Matrix3D other) {                // May be able to delete this function, and
            return Ray(Matrix3D(x, y, z) + other, unitDir); // the namespaces along with it; no longer
        }                                               // adding rays anywhere, it seems

        std::vector<Matrix3D> intersections(std::vector<Shape3D> shapes) {  // Unit test this function
            std::vector<Matrix3D> points;
            for(Shape3D shape : shapes) {
                if(!intersects(shape, *this).empty()) {
                    for(Matrix3D m : intersects(shape, *this))
                        points.push_back(m);
                    continue;
                }
            }
            return points;
        }
    };
}

struct Camera : dirKnown::Ray { // pointKnown::Ray instead?
    Camera() {}

    Camera(Matrix3D initPt, Matrix3D dir) {
        Ray ray = Ray(initPt, dir);
        x = ray.x, y = ray.y, z = ray.z, unitDir = ray.unitDir;
    }
}; */

struct Ray : Matrix3D {
    Matrix3D unitDir;
    
    Ray() {}

    // Subtract coordinates of initial point from those of non-initial point,
    // then normalize, and you get unit vector pointing same direction that ray points
    Ray(Matrix3D initPt, Matrix3D nonInitPt) {
        x = initPt.x, y = initPt.y, z = initPt.z;
        unitDir = Matrix3D(nonInitPt.x - x, nonInitPt.y - y, nonInitPt.z - z).normalize();
    }
};

struct Camera : Ray {
    Camera() {}

    Camera(Ray ray) {
        x = ray.x, y = ray.y, z = ray.z, unitDir = ray.unitDir;
    }
};

/* struct PicPlane() {

}; */

void inputPicProps(int&, int&);
void savePic(std::vector<std::vector<Color> >, std::string, int, int, int);
bool intersects(Sphere sphere, Matrix3D point);
// bool intersects(PicPlane pic, Matrix3D point);
bool intersects(Shape3D shape, Matrix3D point);
void traceRays(Camera, Light, int, int, std::vector<Shape3D>);
void inputMatrix(Matrix3D&, std::string);
std::vector<Shape3D> inputShapes();

void inputPicProps(int& picHeight, int& picWidth) {
    std::cout << "Do you want the image to be 700 pixels tall and 500 pixels wide? (y/n) ";
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

void savePic(std::vector<std::vector<Color> > px, std::string fileName, int picHeight, int picWidth, int dpi) {
    std::ofstream imgFile;
    int area = picWidth * picHeight;
    int size = 4 * area;
    int fileSize = size + 54;

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

    imgFile.open(fileName);
    
    imgFile.write(header, 14);
    imgFile.write(dibHeader, 40);

    for(int y = 0; y < picHeight; y++) {
        for(int x = 0; x < picWidth; x++) {
            char pixel[3] = {floor(px[y][x].blue * 255), floor(px[y][x].green * 255), floor(px[y][x].red * 255)}; // Maybe don't multiply
            imgFile.write(pixel, 3);                                                                                // these (by 255)?
        }
    }

    imgFile.close();
}

void traceRays(Camera cam, Light light, int picHeight, int picWidth, std::vector<Shape3D> shapes) {
    // Declare 2D vect of pixels/colors
    std::vector<std::vector<Color> > px;
    for(int y = 0; y < picHeight; y++) {
        std::vector<Color> row;
        px.push_back(row);
        for(int x = 0; x < picWidth; x++) {
            /*
                Add code that changes cam ray direction for each pixel/iteration
            */
            std::vector<Matrix3D> points;
            double minDist;
            Matrix3D nearestPt;
            Shape3D nearestShape;
            for(Shape3D s : shapes) {
                if(!intersects(s, cam).empty()) { // Could do without this if, after some changes below
                    minDist = cam.distanceTo(intersects(s, cam)[0]);
                    for(Matrix3D m : intersects(s, cam)) {
                        minDist = std::min(minDist, cam.distanceTo(m));
                        if(minDist == cam.distanceTo(m)) {
                            nearestPt = m;
                            nearestShape = s;
                        }
                    }
                }
            }
            Ray shadow = Ray(nearestPt, light);
            bool shaded = false;
            for(Shape3D s : shapes) {
                if(!intersects(s, shadow).empty()) {
                    shaded = true;
                    break;
                }
            }
            Color clr = shaded ? Color(0, 0, 0) : nearestShape.clr;
            px[y].push_back(clr);
            points.clear();
        }
    }

    savePic(px, "result.bmp", picHeight, picWidth, 50);
}

void inputMatrix(Matrix3D& matrix, std::string type) {
    std::cout << "\nDo you want to place the " << (type == "Camera" ? "camera" : "light") 
        << " elsewhere than the origin (0, 0, 0)? (y/n) ";
    char ans = ' ';

    while(ans != 'y' && ans != 'Y' && ans != 'n' && ans != 'N') {
        std::cin >> ans;
        if(ans == 'y' || ans == 'Y') {
            std::cout << "\nWhere would you like to place the " 
                << (type == "Camera" ? "camera sensor" : "light source") << "?\nx = ";
            std::cin >> matrix.x;        
            std::cout << "y = ";
            std::cin >> matrix.y;
            std::cout << "z = ";
            std::cin >> matrix.z;
        }
        else
            if(ans != 'n' && ans != 'N') {
                std::cout << "Sorry, what was that? Your answer must be 'y' (yes) or 'n' (no). ";
                continue;
            }
    }
}

std::vector<Shape3D> inputShapes() {
    std::vector<Shape3D> shapes;

    std::cout << "\nWhat type of shape would you like to add to the scene? (Options: sphere) ";
    std::string ans = "";
    Matrix3D c;
    double r;
    ClrMatrix clr = ClrMatrix();

    while(ans != "sphere") {
        std::cin >> ans;
        if(ans == "sphere") {
            std::cout << "Where should the center of the sphere be?\nx = ";
            std::cin >> c.x;
            std::cout << "y = ";
            std::cin >> c.y;
            std::cout << "z = ";
            std::cin >> c.z;

            std::cout << "What should the radius of the sphere be?\nd = ";
            std::cin >> r;

            std::cout << "What color should the sphere be painted? Give the RGB values. \nr = ";
            std::cin >> clr.red;
            std::cout << "g = ";
            std::cin >> clr.green;
            std::cout << "b = ";
            std::cin >> clr.blue;

            shapes.push_back(Sphere(c, r, clr));
        }
        else {
            std::cout << "Sorry, what was that? Your answer must be \"sphere.\" ";
            continue;
        }
    }

    return shapes;
}

std::vector<Matrix3D> intersects(Sphere sphere, Ray ray) {     // Unit test this function
    std::vector<Matrix3D> points;

    double lilDist, bigDist, medDist = (sphere.center - ray).getDotProduct(ray.unitDir);
    Matrix3D ptBtwn = ray + ray.unitDir * medDist;
    double offset = sqrt( pow(sphere.radius, 2) - pow((sphere.center - ptBtwn).getMagnitude(), 2) );
    
    double lilDist = medDist - offset;
    double bigDist = medDist + offset;
    points.push_back(ray + ray.unitDir * lilDist);
    points.push_back(ray + ray.unitDir * bigDist);
    
    return points;
}

/* std::vector<Matrix3D> intersects(PicPlane plane, Ray ray) {    // Unit test this function
    
} */

std::vector<Matrix3D> intersects(Shape3D shape, Ray ray) {     // Unit test this function
    std::cout << "I don't recognize that shape. Couldn't find an intersection if I tried." << std::endl;
    std::vector<Matrix3D> emptyVec;
    return emptyVec;
}

int main() {
    std::vector<Shape3D> shapes;
    Sphere shape1 = Sphere(Matrix3D(60, 100, 21), 15, ClrMatrix(Color(0, 0, 0), 0));
    Sphere shape2 = Sphere(Matrix3D(106, 8, 700), 8, ClrMatrix(Color(0, 0, 0), 0));
    shapes.push_back(shape1);
    shapes.push_back(shape2);

    Matrix3D point = Matrix3D(106, 8, 692);

    std::cout << "Shape 1 does" << (intersects(shape1, point) ? "" : " NOT") << " intersect (" 
        << point.x << "," << point.y << "," << point.z << ")." << std::endl;
    std::cout << "Shape 2 does" << (intersects(shape2, point) ? "" : " NOT") << " intersect (" 
        << point.x << "," << point.y << "," << point.z << ")." << std::endl;

    /* int picHeight = 700, picWidth = 500;
    Light light = Light(0, 0, 0);
    Camera cam = Camera(Ray(Matrix3D(0, 0, 0), Matrix3D(0, 0, 0)));
    
    std::vector<Shape3D> shapes = {Sphere(), Sphere(), Sphere()};
    char ans = ' ';
    std::cout << "Would you like to design a custom image, or just generate a random one? (c, r) ";
    std::cin >> ans;
    while(ans != 'c' && ans != 'r') {
        std::cout << "You should answer either 'c' for custom or 'r' for random. ";
        std::cin >> ans;
    }
    if(ans == 'c') {
        inputPicProps(picHeight, picWidth);
        inputMatrix(cam, "Camera");
        inputMatrix(light, "Light");
        shapes = inputShapes();
    }
    traceRays(cam, light, picHeight, picWidth, shapes); */
}