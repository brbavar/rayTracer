#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include <string>
#include "matrix.h"


Matrix::Matrix() {
}

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


void inputPicProps(int&, int&);

double intersect(std::tuple<Matrix,Matrix>, 
    std::tuple<std::string,Matrix,Matrix,std::pair<Matrix,Matrix>,std::vector<double> >);

double setBrightness(std::pair<Matrix,Matrix>, Matrix, Matrix, char*, Matrix&);

void render(std::pair<Matrix,Matrix>, int, int, 
    std::vector<std::tuple<std::string,Matrix,Matrix,std::pair<Matrix,Matrix>,std::vector<double> >>);

void inputLight(std::pair<Matrix,Matrix>&);

auto inputShapes();


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

double intersect(std::tuple<Matrix,Matrix> ray, 
    std::tuple<std::string,Matrix,Matrix,std::pair<Matrix,Matrix>,std::vector<double> > shape) {
    Matrix rayOrig = std::get<0>(ray);  // Origin, or initial point, of ray
    Matrix rayDir = (std::get<1>(ray) - rayOrig).normalize();  // Normalized direction of ray
    std::string type = std::get<0>(shape);
    Matrix center = std::get<1>(shape);

    if(type == "Sphere") {
        double radius = std::get<4>(shape)[0];

        // Solving for coefficients derived from equation for sphere
        double a = rayDir.dot(rayDir);
        Matrix initPt = rayOrig - center;
        double b = 2 * initPt.dot(rayDir);
        double c = initPt.dot(initPt) - pow(radius, 2);

        if(b * b - 4 * a * c > 0) {
            double dist = (center - rayOrig).dot(rayDir);
            Matrix ptBtwn = rayOrig + rayDir * dist;
            double offset = sqrt( pow(radius, 2) - pow((center - ptBtwn).mag(), 2) );
            dist -= offset;
            return dist;
        }
    }
    if(type == "Plane") {
        Matrix normal = std::get<3>(shape).second - std::get<3>(shape).first;
        Matrix normalDir = normal.normalize();
        if(rayDir.dot(normalDir) > 0) {
            double dist = rayOrig.entries[0][1] / (-rayDir.entries[0][1]);
            return dist;
        }
    }
    return -1;
}

double setBrightness(std::pair<Matrix,Matrix> light, Matrix nearestPt, Matrix normalDir, char* pixel, Matrix& clr) {
    // Get direction of light beam aimed at nearest point of intersection.
    Matrix lightDir = (nearestPt - light.first).normalize();

    /* std::cout << "nearestPt.entries" << std::endl;
    for(auto row : nearestPt.entries)
        for(double entry : row)
            std::cout << entry << " ";
    std::cout << "\nlight.first.entries" << std::endl;
    for(auto row : light.first.entries)
        for(double entry : row)
            std::cout << entry << " ";
    std::cout << "\nlightDir.entries" << std::endl;
    for(auto row : lightDir.entries)
        for(double entry : row)
            std::cout << entry << " ";
    std::cout << "\nnormalDir.entries" << std::endl;
    for(auto row : normalDir.entries)
        for(double entry : row)
            std::cout << entry << " "; */

    double brightness = 1, angle;

    /* Get the angle between the direction of the light that reaches this point and the direction of the normal 
        vector to the shape at this point on its surface. That will determine the brightness of this point. */
    if(lightDir.mag() > 0)
        angle = acos( normalDir.dot(lightDir) / (normalDir.mag() * lightDir.mag()) );
    else
        return brightness;

    /* If the angle between the light and the normal at this point is greater than 1.5708 radians, it is greater 
        than 90 degrees. That means this point is all the way on the other side of the shape, relative to the point 
        at which the light strikes it. So no light reaches this point. It will thus be fully shaded (black). 
        Otherwise, light reaches the point. So the color of the shape becomes the color of this pixel; the 
        darkness or lightness of the pixel's color is determined by multiplying each RGB component of the shape's 
        intrinsic color by the brightness factor calculated below. The smaller the angle, the brighter the shape's 
        surface at this point. */
    brightness -= angle / 1.5708;

    for(int i = 0; i < 3; i++)
        clr.entries[0][i] *= brightness;

    return brightness;
}

void render(std::pair<Matrix,Matrix> light, int picHeight, int picWidth, 
    std::vector<std::tuple<std::string,Matrix,Matrix,std::pair<Matrix,Matrix>,std::vector<double> >> shapes) {
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

    unsigned int seed = std::chrono::steady_clock::now().time_since_epoch().count();
    std::default_random_engine re(seed);
    std::uniform_real_distribution<double> sampleDistro(0, 0.999999999999999);

    /* Iterate over every pixel in the image, setting its color based on both the intrinsic properties of the shape
       and the lighting. */
    for(double z = -picHeight / 2; z < picHeight / 2; z++)
        for(double x = -picWidth / 2; x < picWidth / 2; x++) {
            auto cam = std::make_tuple(Matrix(0, 0, 0), Matrix(x, -100, z));  // Point camera at current pixel sample
            const double INF = std::numeric_limits<double>::infinity();
            double minDist = INF, dist;
            std::string type = "";
            std::pair<Matrix,Matrix> nrml = {Matrix(0,0,0), Matrix(0,0,0)};
            std::vector<double> dimensions{0};
            auto nearestShape = std::make_tuple(type, Matrix(0,0,0), Matrix(0,0,0), nrml, dimensions);

            for(auto s : shapes) {
                dist = intersect(cam, s);
                if(dist != -1) {
                    minDist = std::min(minDist, dist);
                    if(minDist == dist)
                        nearestShape = s;
                }
            }

            /* If minDist is less than infinity, the color of the pixel will depend on whether light reaches the shape at the point of 
                intersection. It will also depend on the angle between the light and normal at that intersection. On the other hand, if 
                minDist is still infinity, meaning cam ray intersects no shape, then pixel's color remains unchanged. */
            if(minDist < INF) {
                Matrix camPos = std::get<0>(cam);
                Matrix camTarg = std::get<1>(cam);
                Matrix camRay = camTarg - camPos;
                Matrix camDir = camRay.normalize();
                Matrix nearestPt = camPos + camDir * minDist;

                Matrix center = std::get<1>(nearestShape);
                Matrix normal;
                if(std::get<0>(nearestShape) == "Sphere")
                    normal = nearestPt - center;
                if(std::get<0>(nearestShape) == "Plane") {
                    Matrix nrmlTarg = std::get<3>(nearestShape).second;
                    normal = nrmlTarg - center;
                }
                Matrix normalDir = normal.normalize();
                bool inside = camDir.dot(normal) > 0 ? true : false;  // 
                
                Matrix clr = std::get<2>(nearestShape);
                double brightness = setBrightness(light, nearestPt, normalDir, pixel, clr);
                bool smallAngle = brightness > 0, unshaded = true;

                /* Cast a ray from the nearest point of intersection to the light source to see if any other shape
                is in the way, preventing the light from reaching said point. */
                /* std::pair<Matrix,Matrix> shadow = std::pair<Matrix,Matrix>(nearestPt, light);
                for(std::tuple<std::string,Matrix,Matrix,std::pair<Matrix,Matrix>,std::vector<double> > s : shapes) {
                    if(s != nearestShape) {
                        dist = s.intersects(shadow);
                        if(dist != -1)
                            unshaded = false;
                    }
                } */

                /* If angle is small enough for light to reach intersection point, and no shape casts shadow over that point, 
                then give pixel same color as shape. */
                if(smallAngle /* && unshaded */) {
                    for(int i = 0; i < 3; i++) {
                        pixel[i] = floor(clr.entries[0][i] * 255);
                    }
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

void inputLight(std::pair<Matrix,Matrix>& light) {
    std::cout << "\nDo you want to place the light source elsewhere than (" << light.first.entries[0][0] << ", " 
        << light.first.entries[0][1] << ", " << light.first.entries[0][2] << ")? (y/n) ";
    char ans = ' ';

    while(ans != 'y' && ans != 'Y' && ans != 'n' && ans != 'N') {
        std::cin >> ans;
        if(ans == 'y' || ans == 'Y') {
            std::cout << "\nWhere would you like to place the light source?\nx = ";
            std::cin >> light.first.entries[0][0];        
            std::cout << "y = ";
            std::cin >> light.first.entries[0][1];
            std::cout << "z = ";
            std::cin >> light.first.entries[0][2];

            std::cout << "Sounds good. What color should the light be? Give the RGB values. \n r = ";
            std::cin >> light.second.entries[0][2];
            std::cout << "g = ";
            std::cin >> light.second.entries[0][1];
            std::cout << "b = ";
            std::cin >> light.second.entries[0][0];
        }
        else
            if(ans != 'n' && ans != 'N') {
                std::cout << "Sorry, what was that? Your answer must be 'y' (yes) or 'n' (no). ";
                continue;
            }
    }
}

auto inputShapes() {
    std::vector<std::tuple<std::string,Matrix,Matrix,std::pair<Matrix,Matrix>,std::vector<double> >> shapes;
    Matrix c = Matrix(0, 0, 0);
    double r;
    Matrix clr = Matrix();

    int shapeCount = -1;
    std::cout << "How many shapes would you like to put in the picture?" << std::endl;
    while(shapeCount < 0) {
        std::cin >> shapeCount;

    }

    shapes.reserve(shapeCount);

    std::cout << "\nWhat type of shape would you like to add to the scene? (Options: sphere) ";
    std::string ans = "";

    while(ans != "sphere" && ans != "Sphere") {
        std::cin >> ans;
        if(ans == "sphere" || ans == "Sphere") {
            std::cout << "Nice choice. Where should the center of the sphere be?\nx = ";
            std::cin >> c.entries[0][0];
            std::cout << "y = ";
            std::cin >> c.entries[0][1];
            std::cout << "z = ";
            std::cin >> c.entries[0][2];

            std::cout << "All right. What should the radius of the sphere be?\nd = ";
            std::cin >> r;

            std::cout << "Splendid! What color should the sphere be painted? Give the RGB values. \nr = ";
            std::cin >> clr.entries[0][2];
            std::cout << "g = ";
            std::cin >> clr.entries[0][1];
            std::cout << "b = ";
            std::cin >> clr.entries[0][0];

            std::cout << "One more thing: How see-through should the sphere be, on a scale of 0 to 10 (0 being not "
                << "at all see-through, and 10 being as see-through as possible)? ";
            std::cin >> clr.entries[0][3];
            clr.entries[0][3] = (10 - clr.entries[0][3]) / 10;

            std::pair<Matrix,Matrix> normal = {Matrix(0,0,0), Matrix(0,0,0)};
            std::vector<double> dimensions{r};
            auto shape = std::make_tuple("Sphere", c, clr, normal, dimensions);
            shapes.push_back(shape);
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
    
    std::uniform_real_distribution<double> xDistro(-picWidth / 2, picWidth);
    std::uniform_real_distribution<double> yDistro(50, 300); 
    std::uniform_real_distribution<double> zDistro(-picHeight / 2, picHeight);
    
    std::uniform_real_distribution<double> clrDistro(0, 1);
    
    std::uniform_int_distribution<int> shapeCountDistro(1, 10);
    std::uniform_int_distribution<int> floorDistro(0, 1);

    std::vector<std::vector<double> > entries;
    entries.reserve(1);
    std::vector<double> list;
    list.reserve(4);
    for(int i = 0; i < 4; i++)
        list.push_back(clrDistro(re));
    entries.push_back(list);
    Matrix lightClr = Matrix(entries);
    std::pair<Matrix,Matrix> light = {Matrix(0,0,0), lightClr};

    std::vector<std::tuple<std::string,Matrix,Matrix,std::pair<Matrix,Matrix>,std::vector<double> >> shapes;
    char ans = ' ';

    std::cout << "Would you like to design a custom image, or just generate one automatically? (c, r) ";
    std::cin >> ans;
    while(ans != 'c' && ans != 'r') {
        std::cout << "You should answer either 'c' for custom or 'r' for random. ";
        std::cin >> ans;
    }
    if(ans == 'c') {
        inputPicProps(picHeight, picWidth);
        inputLight(light);
        int size = inputShapes().size();
        shapes.reserve(size);
        for(int i = 0; i < size; i++)
            shapes[i] = inputShapes()[i];
    }
    else {
        int shapeCount = 1 /* shapeCountDistro(re) */;
        shapes.reserve(shapeCount);
        int onFloor = 1 /* floorDistro(re) */;

        if(onFloor) {            
            Matrix center = Matrix(0,400,280);
            std::pair<Matrix,Matrix> normal = {center, Matrix(0,400,0)};
            std::vector<double> dimensions;
            entries.clear();
            list.clear();
            for(int i = 0; i < 4; i++)
                list.push_back(clrDistro(re));
            entries.push_back(list);
            Matrix clr = Matrix(entries);
            auto plane = std::make_tuple("Plane", center, clr, normal, dimensions);
            shapes.push_back(plane);
            shapeCount++;
        }

        for(int i = 0; i < shapeCount; i++) {
            Matrix center = Matrix(0,400,0);

            entries.clear();
            list.clear();
            for(int i = 0; i < 4; i++)
                list.push_back(clrDistro(re));
            entries.push_back(list);
            Matrix clr = Matrix(entries);

            std::pair<Matrix,Matrix> normal = {Matrix(0,0,0), Matrix(0,0,0)};

            std::vector<double> dimensions{280};
            
            auto shape = std::make_tuple("Sphere", center, clr, normal, dimensions);
            shapes.push_back(shape);
        }

        std::cout << (onFloor ? "Shapes are on floor." : "Shapes are free-floating.") << std::endl;
    }

    std::cout << "Light is at (" << light.first.entries[0][0] << ", " << light.first.entries[0][1] << ", " << light.first.entries[0][2] << ")." << std::endl;

    render(light, picHeight, picWidth, shapes);
}