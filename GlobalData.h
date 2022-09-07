
#pragma once
#include <string>
#include <random>


#define EPSILON 0.00001


class RayTracer;


class GlobalData {

public:
	GlobalData(int argc, char** argv);
    ~GlobalData();

    double rand() {
#if 1
        // random seed
        static std::random_device rd;
        static std::mt19937 engine(rd());
#else
        // deterministic randomness
        static std::mt19937 engine(37);
#endif
        static std::uniform_real_distribution<double> dist(0.0, 1.0);
        return dist(engine);
    }

    float angle;

    const int width = 1920, height = 1080;
    const double invW = 1.0 / double(width), invH = 1.0 / double(height);
    const double aspectRatio = width / double(height);
    const int FPS = 60;


    int numOpalLayers = 3;
    // Range of sphere size is roughly 200-350nm
    double opalSphereSize = 220;
    // Ratio of water inside the opal
    const double waterRatio = 0.1;

    const int numBounces = 15;
    //const int numPhotons = 1000;
    const double energyThresh = pow(10, -9);

    const int fov = 60;
    const double zNear = 0.1, zFar = 100.0;

    RayTracer *raytracer;
};


extern GlobalData* GLOBAL_data;