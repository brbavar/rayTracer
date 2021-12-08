#include <iostream>
#include <string>

void inputLightLoc(int[3], char);
void inputShape();

struct Camera {
    int sensor[3];

    Camera(int x, int y, int z) {
        this->sensor[0] = x;
        this->sensor[1] = y;
        this->sensor[2] = z;
    }


};

struct Sphere {


    Sphere() {

    }
};

void inputLightLoc(int light[3]) {
    std::cout << "Do you want to place the light elsewhere than the origin (0, 0, 0)? (y/n) ";
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

        std::cout << "\nOkay, the light will be placed at " 
            << (light[0] == 0 && light[1] == 0 && light[2] == 0 ? "the origin " : "") 
            << "(" << light[0] << ", " << light[1] << ", " << light[2] << ")." << std::endl;
    }
}

void inputShape() {
    std::cout << "\nWhat type of shape would you like to add to the scene? (Options: sphere) ";
    std::string ans = "";
    

    while(ans != "sphere") {
        std::cin >> ans;
        if(ans == "sphere") {
            std::cout << ""
            std::cin
        }
        else {
            std::cout << "Sorry, what was that? Your answer must be \"sphere.\" ";
            continue;
        }

    }
}

int main() {
    int light[] = {0, 0, 0};
    inputLightLoc(light);
    std::cout << light[0] << light[1] << light[2] << std::endl;
    inputShape();
}