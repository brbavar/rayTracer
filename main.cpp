#include <iostream>
#include <string>

void inputPicProps(int&, int&);
void traceRays();
void inputLightLoc(int[3], char);
void inputShape();

struct Color {
    double red;
    double green;
    double blue;

    Color(double r, double g, double b) {
        this->red = r;
        this->green = g;
        this->blue = b;
    }
};

struct Camera {
    int sensor[3];

    Camera(int x, int y, int z) {
        this->sensor[0] = x;
        this->sensor[1] = y;
        this->sensor[2] = z;
    }


};

struct Sphere {
    int center[3];
    int diameter;

    Sphere(int c[3], int d) {
        for(int i = 0; i < 3; i++)
            this->center[i] = c[i];
        this->diameter = d;
    }
};

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
            while(std::cin.fail() || picHeight < min || picHeight > max)
                std::cout << "Try again. You must enter an integer between " << min << " and " << max 
                    << " (" << min << ", " << min + 1 << ", " << min + 2 << ", . . . " << max - 1 << ", " << max << "). ";

            std::cout << "\nHow wide should the image be, in pixels? (Min width = 20px, max = 3840px) ";
            max = 3840;
            std::cin >> picWidth;
            while(std::cin.fail() || picHeight < min || picHeight > max)
                std::cout << "Try again. You must enter an integer between " << min << " and " << max 
                    << " (" << min << ", " << min + 1 << ", " << min + 2 << ", . . . " << max - 1 << ", " << max << "). ";
        }
        else
            if(ans != 'y' && ans != 'Y') {
                std::cout << "Sorry, what was that? Your answer must be 'y' (yes) or 'n' (no). ";
                continue;
            }
    }
}

void traceRays(int picHeight, int picWidth) {
     for(int x = 0; x < picHeight; x++)
         for(int y = 0; y < picWidth; y++) {

         }
}

void inputLightLoc(int light[3]) {
    std::cout << "\nDo you want to place the light elsewhere than the origin (0, 0, 0)? (y/n) ";
    char ans = ' ';

    while(ans != 'y' && ans != 'Y' && ans != 'n' && ans != 'N') {
        std::cin >> ans;
        if(ans == 'y' || ans == 'Y') {
            std::cout << "\nWhere would you like to place the light source?\nx = ";
            std::cin >> light[0];        
            std::cout << "y = ";
            std::cin >> light[1];
            std::cout << "z = ";
            std::cin >> light[2];
        }
        else
            if(ans != 'n' && ans != 'N') {
                std::cout << "Sorry, what was that? Your answer must be 'y' (yes) or 'n' (no). ";
                continue;
            }
    }
}

void inputShape() {
    std::cout << "\nWhat type of shape would you like to add to the scene? (Options: sphere) ";
    std::string ans = "";
    int c[3];
    int d;

    while(ans != "sphere") {
        std::cin >> ans;
        if(ans == "sphere") {
            std::cout << "Where should the center of the sphere be?\nx = ";
            std::cin >> c[0];
            std::cout << "y = ";
            std::cin >> c[1];
            std::cout << "z = ";
            std::cin >> c[2];
            std::cout << "What should the diameter of the sphere be, in pixels?\nd = ";
            std::cin >> d;
            Sphere ball = Sphere(c, d);
        }
        else {
            std::cout << "Sorry, what was that? Your answer must be \"sphere.\" ";
            continue;
        }
    }
}

int main() {
    int picHeight = 700, picWidth = 500;
    inputPicProps(picHeight, picWidth);
    int light[] = {0, 0, 0};
    inputLightLoc(light);
    inputShape();
    traceRays(picHeight, picWidth);
}