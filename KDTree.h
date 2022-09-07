#pragma once

#include <cstdlib>
#include <vector>

#include "VectorStuff.h"
#include "BoundingBox.h"


// ==================================================================
// A hierarchical spatial data structure to store photons.  This data
// struture allows for fast nearby neighbor queries for use in photon
// mapping.

class KDTree {
public:

    // ========================
    // CONSTRUCTOR & DESTRUCTOR
    KDTree(const BoundingBox& _bbox, int _depth = 0) {
        bbox = _bbox;
        depth = _depth;
        child1 = NULL;
        child2 = NULL;
    }
    ~KDTree();

    // =========
    // ACCESSORS
    // boundingbox
    const Vec3& getMin() const { return bbox.minimum; }
    const Vec3& getMax() const { return bbox.maximum; }
    bool overlaps(const BoundingBox& bb) const;
    // hierarchy
    int getDepth() const { return depth; }
    bool isLeaf() const {
        if (child1 == NULL && child2 == NULL) return true;
        assert(child1 != NULL && child2 != NULL);
        return false;
    }
    const KDTree* getChild1() const { assert(!isLeaf()); assert(child1 != NULL); return child1; }
    const KDTree* getChild2() const { assert(!isLeaf()); assert(child2 != NULL); return child2; }
    // photons
    const std::vector<Photon>& getPhotons() const { return photons; }
    void CollectPhotonsInBox(const BoundingBox& bb, std::vector<Photon>& photons) const;

    // =========
    // MODIFIERS
    void AddPhoton(const Photon& p);
    bool PhotonInCell(const Photon& p);

    int numPhotons();
    int numBoxes();

    void DrawHelp() {
        if (child1 != NULL) child1->DrawHelp();
        if (child2 != NULL) child2->DrawHelp();

        bbox.Draw();
    }

    void Draw() {
        glPolygonMode(GL_FRONT, GL_LINE);
        glPolygonMode(GL_BACK, GL_LINE);
        DrawHelp();
        glPolygonMode(GL_FRONT, GL_FILL);
        glPolygonMode(GL_BACK, GL_FILL);
    }

private:

    // HELPER FUNCTION
    void SplitCell();

    // REPRESENTATION
    BoundingBox bbox;
    KDTree* child1;
    KDTree* child2;
    int split_axis;
    float split_value;
    std::vector<Photon> photons;
    int depth;
};