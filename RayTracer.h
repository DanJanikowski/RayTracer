
#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <GL/glut.h>

#include "VectorStuff.h"
#include "KDTree.h"


class Material;
class Hit;

class Primitive;
class Plane; 
class Sphere;

enum MatType { Normal, Opal, Water, Reflector, Collector };
class Material {
public:
    Material() {}

    Material(const Vec3& sc, const Vec3& ec) :
        surfaceColor(sc), emissionColor(ec)
    {}

    Material(const Material& m) :
        surfaceColor(m.surfaceColor), emissionColor(m.emissionColor)
    {}

    Material& operator = (const Material& m) { 
        surfaceColor = m.surfaceColor;
        emissionColor = m.emissionColor;
        return *this; 
    }

    Vec3 GetSC() const { return surfaceColor; }
    Vec3 GetEC() const { return emissionColor; }

    void SetMatType(MatType mt) { matType = mt; }
    MatType GetMatType() { return matType; }

    friend std::ostream& operator<<(std::ostream& os, const Material& m) {
        return os << "Surface Color: " << m.surfaceColor;
    }

private:
    Vec3 surfaceColor, emissionColor;           /// surface color and emission (light)
    MatType matType = Normal;
};


class Hit {
public:
    Hit() {
        t = INT_MAX;
        normal = Vec3(0.0, 0.0, 0.0);
        material = NULL;
        frontface = false;
    }
    Hit(const Hit& h) {
        t = h.t;
        normal = h.normal;
        material = h.material;
        frontface = h.frontface;
    }

    float GetT() const { return t; }
    Vec3 GetNormal() const { return normal; }
    Material* GetMaterial() const { return material; }
    bool GetFrontFace() const { return frontface; }

    // MODIFIER
    void Set(double t_, const Vec3& n, Material* m, bool frontface_) {
        t = t_;  normal = n; material = m; frontface = frontface_;
    }

    void SetNormal(const Vec3& norm) { normal = norm; }

private:
    double t;
    Vec3 normal;
    Material* material;
    bool frontface;
};




class Primitive {
public:
    Material* GetMaterial() { return &material; }
    void SetMatType(MatType mt) { material.SetMatType(mt); }
    MatType GetMatType() { return material.GetMatType(); }
    virtual bool Intersect(Ray * ray, Hit& hit) = 0;
    virtual void Draw() = 0;
protected:
    Material material;
};

class Plane : public Primitive {
public:

    Vec3 p1, p2, p3, p4;                 // corners of quad (in counter-clockwise order)
    Vec3 normal;

    Plane() : p1(Vec3(0)), p2(Vec3(0)), p3(Vec3(0)), p4(Vec3(0)) {
        material = Material(Vec3(1.0), Vec3(0));
        normal = Vec3();
    }

    Plane(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d,
        const Vec3& sc, const Vec3& ec = 0) :
        p1(a), p2(b), p3(c), p4(d)
    {
        material = Material(sc, ec);
        normal = Vec3();
    }

    // compute a ray-plane intersection using the geometric solution
    bool Intersect(Ray * ray, Hit& hit) {
        Vec3 E01, E03, P, Te, Q, E23, E21;
        double det, alpha, beta, ep = 0.0000001, t;
        E01 = p2 - p1;
        E03 = p4 - p1;
        P = ray->dir.cross(E03);
        det = E01.dot(P);
        if (fabs(det) < ep) return false;
        Te = ray->origin - p1;
        alpha = (Te.dot(P)) / det;
        if (alpha < 0 || alpha > 1) return false;
        Q = Te.cross(E01);
        beta = (ray->dir.dot(Q)) / det;
        if (beta < 0 || beta > 1) return false;

        t = (E03.dot(Q)) / det;
        if ((alpha + beta) > 1) {
            E23 = p4 - p3;
            E21 = p2 - p3;
            P = ray->dir.cross(E21);
            det = E23.dot(P);
            if (fabs(det) < ep) return false;
            Te = ray->origin - p3;
            alpha = (Te.dot(P)) / det;
            if (alpha < 0) return false;
            Q = Te.cross(E23);
            beta = (ray->dir.dot(Q)) / det;
            if (beta < 0) return false;
        }

        if (t < 0) return false;
        if (t < hit.GetT()) {
            Vec3 normal = E01.cross(E03);
            normal.normalize();
            bool frontface = normal.dot(ray->dir) < 0;
            hit.Set(t, normal * ((frontface) ? 1.0 : -1.0), &material, frontface);
        }
        return true;
    }

    Vec3 getCenter() {
        return (p1 + p2 + p3 + p4) / 4.0;
    }

    Vec3 getNormal() {
        if (normal.norm() == 0) {
            Vec3 E01 = p2 - p1;
            Vec3 E03 = p4 - p1;
            normal = E01.cross(E03);
        }
        return normal;
    }

    void Draw() {
        if (material.GetSC().sqnorm() > 0)
            glColor3f(material.GetSC().x, material.GetSC().y, material.GetSC().z);
        else
            glColor3f(material.GetEC().x, material.GetEC().y, material.GetEC().z);
        glBegin(GL_QUADS);
        glVertex3f(p1.x, p1.y, p1.z);
        glVertex3f(p2.x, p2.y, p2.z);
        glVertex3f(p3.x, p3.y, p3.z);
        glVertex3f(p4.x, p4.y, p4.z);
        glEnd();
    }
};

class Light : public Primitive {
public:

    Plane p, proj;                // corners of quad (in counter-clockwise order)

    Light(const Plane& plane, const Plane& projPlane) {
        p = plane;
        proj = projPlane;
    }

    // compute a ray-plane intersection using the geometric solution
    bool Intersect(Ray * ray, Hit& hit) {
        return p.Intersect(ray, hit);
    }

    Vec3 getCenter() {
        return p.getCenter();
    }

    Vec3 getNormal() {
        return p.getNormal();
    }

    void Draw() {
        p.Draw();
    }
};


class Sphere : public Primitive {
public:
    // Sphere variables.
    Vec3 center;                         /// position of the sphere
    double radius, radius2;              /// sphere radius and radius^2
    //Material material;

    // Sphere constructor.
    // position(c), radius(r), surface color(sc), reflectivity(refl), transparency(transp), emission color(ec)
    Sphere(const Vec3& c, const double& r,
        const Vec3& sc, const Vec3& ec = 0) :
        center(c), radius(r), radius2(r* r)
    {
        material = Material(sc, ec);
    }

    // compute a ray-sphere intersection using the geometric solution
    bool Intersect(Ray* ray, Hit& hit)  {
        Vec3 m = ray->origin - center;
        double b = m.dot(ray->dir);
        double c = m.dot(m) - radius2;
        if (c > 0 && b > 0) return false;   // False if both intersections are behind
        double discr = b * b - c;
        if (discr < 0) return false;

        double t0 = -b - sqrt(discr);
        double t1 = -b + sqrt(discr);
        bool frontface = (t0 > 0);
        double t = (frontface) ? t0 : t1;
        if (t < hit.GetT()) {
            Vec3 normal = ray->PointAtParameter(t) - center;
            normal.normalize();
            hit.Set(t, normal * ((frontface) ? 1.0 : -1.0), &material, frontface);
        }
        return true;
    }

    void Draw() {
        glPushMatrix();
        glTranslatef(center.x, center.y, center.z);
        if (material.GetSC().sqnorm() > 0)
            glColor3f(material.GetSC().x, material.GetSC().y, material.GetSC().z);
        else
            glColor3f(material.GetEC().x, material.GetEC().y, material.GetEC().z);
        glutSolidSphere(radius, 20, 20);
        glPopMatrix();

        //GLfloat mat_specular[4];
        //GLfloat mat_shininess[] = { 50.0 };
        //if (material.GetSC().sqnorm() > 0) {
        //    mat_specular[0] = material.GetSC().x;
        //    mat_specular[1] = material.GetSC().y;
        //    mat_specular[2] = material.GetSC().z;
        //}
        //else {
        //    mat_specular[0] = material.GetEC().x;
        //    mat_specular[1] = material.GetEC().y;
        //    mat_specular[2] = material.GetEC().z;
        //}
        //glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_specular);
        //glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
        //glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

        //glPushMatrix();
        //glTranslatef(center.x, center.y, center.z);
        //glutSolidSphere(radius, 20, 20);
        //glPopMatrix();
    }
};




class RayTracer {
public:
	RayTracer();
    ~RayTracer();

    void createScene();
    void GenerateTriclinicPack(const Vec3& offset, const int width, const int height, const int layers);

	void Display();
    void SetCamMovement(char key, bool down);

    void RotateCam(int dyaw, int dpitch);
    
    Ray ComputeReflectance(const Vec3& dir, const Vec3& normal, const double& n1, const double& n2, double& reflectance, bool method);
    double SilicaDispersionFunction(double lambda);
    double WaterDispersionFunction(double lambda);
    void OpalTraceRay(std::vector<Ray>& contactPhotons, const Vec3& orig, Ray* ray, Hit& hit, int bounces, int layer, bool visible = false);
    void TraceRay(std::vector<Ray>& contactPhotons, Ray* ray, Hit& hit, int bounces, bool visible = false);
    void TraceViewRay();
    void TraceLightRays();
    void DrawRayTree();
private:
    // Private Functions (Helpers)
    void RayTrace();

    double dx, dy, dz;
    double moveScale;
    Vec3 camPos;
    Vec3 camDir;
    Vec3 camUp;
    double totYaw = 0, totPitch = 0;

    double A1 = 0.6961663, A2 = pow(0.06840432, 2),
        B1 = 0.4079426, B2 = pow(0.1162414, 2),
        C1 = 0.8974794, C2 = pow(9.896161, 2);

    double X1 = 0.75831, X2 = 0.01007,
        Y1 = 0.08495, Y2 = 8.91377;


    GLuint texture;

    std::vector<unsigned char> buffer;
    std::vector<Primitive*> shapes;
    std::vector<Light*> lights;
    std::vector<Line> rayTree;
    Plane* collector;

    //KDTree* kdtree;
    std::vector<Photon> photons;

    GLuint dlist = 0;
};