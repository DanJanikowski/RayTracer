#pragma once

#include <cmath>
#include <iostream>
#include <iomanip>
#include <GL/glut.h>


const double PI = 3.141592653589793;
#define TO_RADS PI/180.0
#define TO_DEGS 180.0/PI


class Vec3 {
public:
    double x, y, z;
    // Vector constructors.
    Vec3() : x(0.0), y(0.0), z(0.0) {}
    Vec3(double xx) : x(xx), y(xx), z(xx) {}
    Vec3(double xx, double yy, double zz) : x(xx), y(yy), z(zz) {}
    Vec3(const Vec3& v) { x = v.x; y = v.y; z = v.z; }
    Vec3& operator = (const Vec3& v) { x = v.x; y = v.y; z = v.z; return *this; }

    double sqnorm() const { return x * x + y * y + z * z; }
    double norm() const { return sqrt(x * x + y * y + z * z); }
    Vec3& normalize() {
        double nor = sqnorm();
        if (nor > 0) {
            double invNor = 1 / sqrt(nor);
            x *= invNor, y *= invNor, z *= invNor;
        }
        return *this;
    }

    Vec3 operator - () const { return Vec3(-x, -y, -z); }
    friend Vec3 operator * (const double& c, const Vec3& v) { return Vec3(v.x * c, v.y * c, v.z * c); }
    friend Vec3 operator * (const Vec3& v, const double& c) { return Vec3(v.x * c, v.y * c, v.z * c); }
    friend Vec3 operator / (const double& c, const Vec3& v) { return Vec3(c / v.x, c / v.y, c / v.z); }
    friend Vec3 operator / (const Vec3& v, const double& c) { return Vec3(v.x / c, v.y / c, v.z / c); }
    Vec3 operator + (const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator - (const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator * (const Vec3& v) const { return Vec3(x * v.x, y * v.y, z * v.z); }
    Vec3& operator += (const Vec3& v) { x += v.x, y += v.y, z += v.z; return *this; }
    Vec3& operator -= (const Vec3& v) { x -= v.x, y -= v.y, z -= v.z; return *this; }
    Vec3& operator *= (const Vec3& v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
    Vec3& operator *= (const double& c) { x *= c, y *= c, z *= c; return *this; }

    double dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }
    Vec3 cross(const Vec3& v) const { return Vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); }
    Vec3 rotateAboutV(const Vec3& v, double theta) const {
        const Vec3 temp = *this;
        Vec3 crossed = -this->cross(v);
        return (temp * cos(theta) + crossed * sin(theta) + (v * (v.dot(temp))) * (1.0 - cos(theta)));
    }

    friend std::ostream& operator<<(std::ostream& os, const Vec3& v) {
        return os << std::setprecision(4) << "<" << v.x << ", " << v.y << ", " << v.z << ">" << std::endl;
    }
};


class Ray {
public:
    Vec3 origin, dir;

    // Wavelength of the ray (in a vacuum in NANOMETERS)
    double wavelength, energy, phase;

    Ray() : wavelength(0.7), energy(1.0), phase(0) {}

    Ray(Vec3 orig_, Vec3 dir_) :
        origin(orig_), dir(dir_), wavelength(0.7), energy(1.0), phase(0)
    {}

    Ray(Vec3 orig_, Vec3 dir_, double lambda, double energy_, double phase_) :
        origin(orig_), dir(dir_), wavelength(lambda), energy(energy_), phase(phase_)
    {}

    Ray(Ray* r) : origin(r->origin), dir(r->dir), wavelength(r->wavelength), energy(r->energy), phase(r->phase)
    {}

    void AddPhase(double phase_) { phase = fmod(phase + phase_, 1.0); }

    Vec3 PointAtParameter(double t) const {
        return origin + dir * t;
    }
};

// Like a ray but has 2 endpoints so we can draw it
class Line {
public:
    Line(const Vec3& startP, const Vec3& endP) :
        p(startP), q(endP), c(Vec3(1.0, 0, 0))
    {}

    Line(const Vec3& startP, const Vec3& dir, const double& t) :
        p(startP), q(startP + dir * t), c(Vec3(1.0, 0, 0))
    {}

    Line(const Vec3& startP, const Vec3& endP, const Vec3& color) :
        p(startP), q(endP), c(color)
    {}

    Vec3 GetColor() { return c; }

    void Draw() const {
        glColor4d(c.x, c.y, c.z, 1.0);
        glBegin(GL_LINES);
        glVertex3f(p.x, p.y, p.z);
        glVertex3f(q.x, q.y, q.z);
        glEnd();
        glClearColor(0.0, 0.0, 0.0, 0.0);
    }

private:
    Vec3 p, q, c;
};


class Photon {
public:
    Vec3 pos;

    Photon(const Vec3& pos_) : pos(pos_) {}
};