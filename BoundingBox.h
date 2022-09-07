#pragma once

#include <cassert>
#include <algorithm>
#include <vector>
#include <GL/glut.h>

#include "VectorStuff.h"


// ====================================================================
// ====================================================================

class BoundingBox {

public:
    Vec3 minimum;
    Vec3 maximum;

    // ========================
    // CONSTRUCTOR & DESTRUCTOR
    BoundingBox() {
        Set(Vec3(0, 0, 0), Vec3(0, 0, 0));
    }
    BoundingBox(const Vec3& pos) {
        Set(pos, pos);
    }
    BoundingBox(const Vec3& _minimum, const Vec3& _maximum) {
        Set(_minimum, _maximum);
    }

    void getCenter(Vec3& c) const {
        c = maximum;
        c -= minimum;
        c *= 0.5f;
        c += minimum;
    }
    double maxDim() const {
        double x = maximum.x - minimum.x;
        double y = maximum.y - minimum.y;
        double z = maximum.z - minimum.z;
        return fmax(x, fmax(y, z));
    }

    // =========
    // MODIFIERS
    void Set(const BoundingBox& bb) {
        minimum = bb.minimum;
        maximum = bb.maximum;
    }
    void Set(const Vec3& _minimum, const Vec3& _maximum) {
        assert(minimum.x <= maximum.x &&
            minimum.y <= maximum.y &&
            minimum.z <= maximum.z);
        minimum = _minimum;
        maximum = _maximum;
    }
    void Extend(const Vec3 v) {
        minimum = Vec3(fmin(minimum.x, v.x),
            fmin(minimum.y, v.y),
            fmin(minimum.z, v.z));
        maximum = Vec3(fmax(maximum.x, v.x),
            fmax(maximum.y, v.y),
            fmax(maximum.z, v.z));
    }
    void Extend(const BoundingBox& bb) {
        Extend(bb.minimum);
        Extend(bb.maximum);
    }

    void Draw() {
        Vec3 p1, p2, p3, p4, p5, p6, p7, p8;

        p1 = minimum;
        p8 = maximum;
        p2 = Vec3(p8.x, p1.y, p1.z);
        p3 = Vec3(p8.x, p1.y, p8.z);
        p4 = Vec3(p1.x, p1.y, p8.z);
        p5 = Vec3(p8.x, p8.y, p1.z);
        p6 = Vec3(p1.x, p8.y, p1.z);
        p7 = Vec3(p1.x, p8.y, p8.z);

        glColor3f(0.8, 0.1, 0.4);

        glBegin(GL_QUADS);
        glVertex3f(p1.x, p1.y, p1.z);
        glVertex3f(p2.x, p2.y, p2.z);
        glVertex3f(p3.x, p3.y, p3.z);
        glVertex3f(p4.x, p4.y, p4.z);
        glEnd();
        glBegin(GL_QUADS);
        glVertex3f(p5.x, p5.y, p5.z);
        glVertex3f(p6.x, p6.y, p6.z);
        glVertex3f(p7.x, p7.y, p7.z);
        glVertex3f(p8.x, p8.y, p8.z);
        glEnd();
        glBegin(GL_QUADS);
        glVertex3f(p1.x, p1.y, p1.z);
        glVertex3f(p4.x, p4.y, p4.z);
        glVertex3f(p7.x, p7.y, p7.z);
        glVertex3f(p6.x, p6.y, p6.z);
        glEnd();
        glBegin(GL_QUADS);
        glVertex3f(p3.x, p3.y, p3.z);
        glVertex3f(p8.x, p8.y, p8.z);
        glVertex3f(p7.x, p7.y, p7.z);
        glVertex3f(p4.x, p4.y, p4.z);
        glEnd();
        glBegin(GL_QUADS);
        glVertex3f(p2.x, p2.y, p2.z);
        glVertex3f(p5.x, p5.y, p5.z);
        glVertex3f(p8.x, p8.y, p8.z);
        glVertex3f(p3.x, p3.y, p3.z);
        glEnd();
        glBegin(GL_QUADS);
        glVertex3f(p1.x, p1.y, p1.z);
        glVertex3f(p6.x, p6.y, p6.z);
        glVertex3f(p5.x, p5.y, p5.z);
        glVertex3f(p2.x, p2.y, p2.z);
        glEnd();
    }
};