#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>

struct Color {
    double blue, green, red, opacity;

    Color() {}

    Color(double b, double g, double r, double o) {
        blue = b, green = g, red = r, opacity = o;
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

    // Returns magnitude of vector represented by matrix
    double mag() {
        return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    }

    // Returns dot product of this vector with another
    double dot(Matrix3D other) {
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
        double mag = this->mag();
        return Matrix3D(x / mag, y / mag, z / mag);
    }

    Matrix3D invert() {
        return Matrix3D(-x, -y, -z);
    }
};

struct Light : Matrix3D {
    Color clr;

    Light() {}

    Light(double comp1, double comp2, double comp3, Color color) {
        x = comp1, y = comp2, z = comp3, clr = color;
    }
};

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

struct Shape3D {
    Matrix3D center;
    double radius, height, width;
    Ray normal;
    Color clr;
    std::string type = "";

    Shape3D() {}

    Shape3D(Color color) {
        clr = color;
    }

    std::vector<Matrix3D> intersects(Ray ray) {
        std::vector<Matrix3D> points;
        if(this->type == "Sphere") {
            double lilDist, bigDist, medDist = (this->center - ray).dot(ray.unitDir);
            Matrix3D ptBtwn = ray + ray.unitDir * medDist;
            double offset = sqrt( pow(this->radius, 2) - pow((this->center - ptBtwn).mag(), 2) );
            
            lilDist = medDist - offset;
            bigDist = medDist + offset;
            points.push_back(ray + ray.unitDir * lilDist);
            points.push_back(ray + ray.unitDir * bigDist);
    
            return points;
        }
        /* if(this->type == "Plane") {

        } */
        
        std::cout << "I don't recognize that shape. Couldn't find an intersection if I tried." << std::endl;
        return points;
    }
};

struct Sphere : Shape3D {
    Sphere() {
        type = "Sphere";
    }

    Sphere(Matrix3D c, double r, Color color) {
        center = c, radius = r, clr = color, type = "Sphere";
    }
};

struct Plane : Shape3D {
    Plane() {
        type = "Plane";
    }

    Plane(Ray r, Matrix3D c, double h, double w) {
        normal = r, center = c, height = h, width = w, type = "Plane";
    }
};

void inputPicProps(int&, int&);
void savePic(std::vector<std::vector<Color> >, std::string, int, int, int);
void traceRays(Light, int, int, std::vector<Shape3D>);
void inputLight(Light&);
std::vector<Shape3D> inputShapes();

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
            char pixel[3] = {floor(px[y][x].blue * 255), floor(px[y][x].green * 255), floor(px[y][x].red * 255)};
            imgFile.write(pixel, 3);
        }
    }

    imgFile.close();
}

void traceRays(Light light, int picHeight, int picWidth, std::vector<Shape3D> shapes /* std::vector<Shape3D> shapes */) {
    // Declare 2D vect of pixels/colors
    std::vector<std::vector<Color> > px;

    Matrix3D nonInitPt = Matrix3D(500, 400, 20);
    Matrix3D center = Matrix3D(500, 400, 0);
    Plane pic = Plane(Ray(center, nonInitPt), center, picHeight, picWidth);

    /* Iterate over every pixel in the image, setting its color based on both the intrinsic properties of the shape
        and the lighting. */
    for(int y = 0; y < picHeight; y++) {
        std::vector<Color> row;
        px.push_back(row);
        for(int x = 0; x < picWidth; x++) {
            // Point camera at current pixel
            Camera cam = Camera(Ray(Matrix3D(500, 400, -1500), Matrix3D(150.5 + x, 150.5 + y, 0)));

            std::vector<Matrix3D> points;
            double minDist;
            Matrix3D nearestPt = Matrix3D(-1, -1, -1);
            Shape3D nearestShape;
            Matrix3D normalDir;
            for(Shape3D s : shapes) {
                if(!s.intersects(cam).empty()) { // Could do without this if, after some changes below
                    minDist = cam.distanceTo(s.intersects(cam)[0]);
                    for(Matrix3D m : s.intersects(cam)) {
                        minDist = std::min(minDist, cam.distanceTo(m));
                        if(minDist == cam.distanceTo(m)) {
                            nearestPt = m;
                            normalDir = Ray(s.center, nearestPt).unitDir;
                            nearestShape = s;
                        }
                    }
                }
            }
            Ray extendedLight = Ray(Matrix3D(light.x, light.y, light.z), nearestPt);
            Matrix3D lightDir = extendedLight.unitDir;
            /* Get the angle between the direction of the light that reaches this point and the direction of the normal 
                vector to the shape at this point on its surface. That will determine the brightness of this point. */
            double angle = acos( normalDir.dot(lightDir) / (normalDir.mag() * lightDir.mag()) );
            // The smaller the angle, the brighter the shape's surface at this point.
            double brightness = 1 - angle / 1.5708;
            Color clr;
            /* If the angle between the light and the normal at this point is greater than 1.5708 radians, it is greater 
                than 90 degrees. That means this point is all the way on the other side of the shape, relative to the point 
                at which the light strikes it. So no light reaches this point. It will thus be fully shaded (black). 
                Otherwise, light reaches the point. So the color of the shape becomes the color of this pixel, but the 
                darkness or lightness of the pixel's color is determined by multiplying each RGB component of the shape's 
                intrinsic color by the brightness factor calculated above. */
            if(angle > 1.5708) {
                clr = Color(0, 0, 0, 1);
                px[y].push_back(clr);
                points.clear();
                continue;
            }
            else {
                nearestShape.clr.blue *= brightness;
                nearestShape.clr.green *= brightness;
                nearestShape.clr.red *= brightness;
            }

            bool shaded = false;
            bool background = false;

            /* If cam ray intersects no shape, set pixel to background color. Otherwise, cast a ray 
               from the nearest point of intersection to the light source to see if any other shape
               is in the way, preventing the light from reaching said point of intersection. */
            if(nearestPt.x == -1) 
                background = true;
            else {
                Ray shadow = Ray(nearestPt, light);
                bool shaded = false;
                for(Shape3D s : shapes) {
                    if(!s.intersects(shadow).empty()) {
                        shaded = true;
                        break;
                    }
                }
            }

            clr = shaded ? Color(0, 0, 0, 1) : ( background ? Color(0, 0, 0, 1) : nearestShape.clr );
            px[y].push_back(clr);
            points.clear();
        }
    }

    savePic(px, "result.bmp", picHeight, picWidth, 50);
}

void inputLight(Light& light) {
    std::cout << "\nDo you want to place the light source elsewhere than (" << light.x << ", " 
        << light.y << ", " << light.z << ")? (y/n) ";
    char ans = ' ';

    while(ans != 'y' && ans != 'Y' && ans != 'n' && ans != 'N') {
        std::cin >> ans;
        if(ans == 'y' || ans == 'Y') {
            std::cout << "\nWhere would you like to place the light source?\nx = ";
            std::cin >> light.x;        
            std::cout << "y = ";
            std::cin >> light.y;
            std::cout << "z = ";
            std::cin >> light.z;

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

std::vector<Shape3D> inputShapes() {
    std::vector<Shape3D> shapes;

    std::cout << "\nWhat type of shape would you like to add to the scene? (Options: sphere) ";
    std::string ans = "";
    Matrix3D c;
    double r;
    Color clr = Color();

    while(ans != "sphere") {
        std::cin >> ans;
        if(ans == "sphere") {
            std::cout << "Nice choice. Where should the center of the sphere be?\nx = ";
            std::cin >> c.x;
            std::cout << "y = ";
            std::cin >> c.y;
            std::cout << "z = ";
            std::cin >> c.z;

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
    unsigned int seed = std::chrono::steady_clock::now().time_since_epoch().count();
    std::default_random_engine re(seed);
    
    std::uniform_real_distribution<double> x1Distro(-949, 949);
    std::uniform_real_distribution<double> y1Distro(-749, 749);
    
    std::uniform_real_distribution<double> x2Distro(1, 949);
    std::uniform_real_distribution<double> y2Distro(1, 749);
    
    std::uniform_real_distribution<double> zDistro(1, 10);
    
    std::normal_distribution<double> radDistro(150, 90);
    
    std::uniform_real_distribution<double> clrDistro(0, 1);
    
    std::uniform_int_distribution<int> shapeCountDistro(1, 10);

    Color lightClr = Color(clrDistro(re), clrDistro(re), clrDistro(re), clrDistro(re));
    Light light = Light(x1Distro(re), y1Distro(re), -zDistro(re), lightClr);

    std::vector<Shape3D> shapes;
    int picHeight = 500, picWidth = 700;
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
        int shapeCount = shapeCountDistro(re);

        std::cout << shapeCount << std::endl;
        std::cout << std::endl;

        for(int i = 0; i < shapeCount; i++) {
            Matrix3D center = Matrix3D(x2Distro(re), y2Distro(re), zDistro(re));
            double radius = radDistro(re);
            Color clr = Color(clrDistro(re), clrDistro(re), clrDistro(re), clrDistro(re));

            shapes.push_back(Sphere(center, radius, clr));

            std::cout << "x = " << center.x << ", y = " << center.y << ", z = " << center.z << ", r = " << radius << std::endl;
            std::cout << "b = " << clr.blue << ", g = " << clr.green << ", r = " << clr.red << ", o = " << clr.opacity << std::endl;
            std::cout << "i = " << i << std::endl;
            std::cout << std::endl;
        }
    }

    std::cout << "Light is at (" << light.x << ", " << light.y << ", " << light.z << ")." << std::endl;

    traceRays(light, picHeight, picWidth, shapes);
}