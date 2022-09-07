#include "kdtree.h"
#include "GlobalData.h"

#define MAX_PHOTONS_BEFORE_SPLIT 100
#define MAX_DEPTH 1

// ==================================================================
// DESTRUCTOR
// ==================================================================
KDTree::~KDTree() {
    if (!isLeaf()) {
        delete child1;
        delete child2;
    }
    else {
        // delete all the photons (this is done automatically since they
        // are stored directly in an STL vector, not using pointers)
    }
}


// ==================================================================
// HELPER FUNCTIONS

bool KDTree::PhotonInCell(const Photon& p) {
    const Vec3& p1 = bbox.minimum;
    const Vec3& p8 = bbox.maximum;
    const Vec3& pos = p.pos;
    std::cout << p1 << p8 << pos << std::endl;

    Vec3 p2, p4, p6;
    p2 = Vec3(p8.x, p1.y, p1.z);
    p4 = Vec3(p1.x, p1.y, p8.z);
    p6 = Vec3(p1.x, p8.y, p1.z);
    Vec3 u, v, w;
    u = p1 - p2;
    v = p1 - p6;
    w = p1 - p4;

    if (u.dot(pos) < u.dot(p1) + EPSILON &&
        u.dot(pos) > u.dot(p2) - EPSILON &&
        v.dot(pos) < v.dot(p1) + EPSILON &&
        v.dot(pos) > v.dot(p6) - EPSILON &&
        w.dot(pos) < w.dot(p1) + EPSILON &&
        w.dot(pos) > w.dot(p4) - EPSILON)
        return true;
    return false;
}

int KDTree::numPhotons() {
    int count = photons.size();
    if (child1 != NULL) count += child1->numPhotons();
    if (child2 != NULL) count += child2->numPhotons();
    return count;
}

int KDTree::numBoxes() {
    int count = isLeaf();
    if (child1 != NULL) count += child1->numBoxes();
    if (child2 != NULL) count += child2->numBoxes();
    return count;
}

bool KDTree::overlaps(const BoundingBox& bb) const {
    const Vec3& bb_min = bb.minimum;
    const Vec3& bb_max = bb.maximum;
    const Vec3& tmp_min = bbox.minimum;
    const Vec3& tmp_max = bbox.maximum;
    if (bb_min.x > tmp_max.x) return false;
    if (tmp_min.x > bb_max.x) return false;
    if (bb_min.y > tmp_max.y) return false;
    if (tmp_min.y > bb_max.y) return false;
    if (bb_min.z > tmp_max.z) return false;
    if (tmp_min.z > bb_max.z) return false;
    return true;
}


// ==================================================================
void KDTree::AddPhoton(const Photon& p) {
    const Vec3& position = p.pos;
    //assert(PhotonInCell(p));
    if (isLeaf()) {
        // this cell is a leaf node
        photons.push_back(p);
        if (photons.size() > MAX_PHOTONS_BEFORE_SPLIT && depth < MAX_DEPTH) {
            std::cout << "SPLITTING" << std::endl;
            SplitCell();
        }
    }
    else {
        // this cell is not a leaf node
        // decide which subnode to recurse into
        if (split_axis == 0) {
            if (position.x < split_value)
                child1->AddPhoton(p);
            else
                child2->AddPhoton(p);
        }
        else if (split_axis == 1) {
            if (position.y < split_value)
                child1->AddPhoton(p);
            else
                child2->AddPhoton(p);
        }
        else {
            assert(split_axis == 2);
            if (position.z < split_value)
                child1->AddPhoton(p);
            else
                child2->AddPhoton(p);
        }
    }
}


// ==================================================================
void KDTree::CollectPhotonsInBox(const BoundingBox& bb, std::vector<Photon>& photons) const {
    // explicitly store the queue of cells that must be checked (rather
    // than write a recursive function)
    std::vector<const KDTree*> todo;
    todo.push_back(this);
    while (!todo.empty()) {
        const KDTree* node = todo.back();
        todo.pop_back();
        if (!node->overlaps(bb)) continue;
        if (node->isLeaf()) {
            // if this cell overlaps & is a leaf, add all of the photons into the master list
            // NOTE: these photons may not be inside of the query bounding box
            const std::vector<Photon>& photons2 = node->getPhotons();
            int num_photons = photons2.size();
            for (int i = 0; i < num_photons; i++) {
                photons.push_back(photons2[i]);
            }
        }
        else {
            // if this cell is not a leaf, explore both children
            todo.push_back(node->getChild1());
            todo.push_back(node->getChild2());
        }
    }
}


// ==================================================================
void KDTree::SplitCell() {
    const Vec3& min = bbox.minimum;
    const Vec3& max = bbox.maximum;
    float dx = fabs(max.x - min.x);
    float dy = fabs(max.y - min.y);
    float dz = fabs(max.z - min.z);
    std::cout << dx << " " << dy << " " << dz << std::endl;
    // split this cell in the middle of the longest axis
    Vec3 min1, min2, max1, max2;
    if (dx >= dy && dx >= dz) { // X axis
        std::sort(photons.begin(), photons.end(), [](const Photon& p1, const Photon& p2) {
            return p1.pos.x < p2.pos.x;
            });
        split_axis = 0;
        split_value = photons.at(photons.size() / 2).pos.x + EPSILON;
        min1 = Vec3(min.x, min.y, min.z);
        max1 = Vec3(split_value, max.y, max.z);
        min2 = Vec3(split_value, min.y, min.z);
        max2 = Vec3(max.x, max.y, max.z);
    }
    else if (dy >= dx && dy >= dz) { // Y axis
        std::sort(photons.begin(), photons.end(), [](const Photon& p1, const Photon& p2) {
            return p1.pos.y < p2.pos.y;
            });
        split_axis = 1;
        split_value = photons.at(photons.size() / 2).pos.y + EPSILON;
        min1 = Vec3(min.x, min.y, min.z);
        max1 = Vec3(max.x, split_value, max.z);
        min2 = Vec3(min.x, split_value, min.z);
        max2 = Vec3(max.x, max.y, max.z);
    }
    else { // Z axis
        assert(dz >= dx && dz >= dy);
        std::sort(photons.begin(), photons.end(), [](const Photon& p1, const Photon& p2) {
            return p1.pos.z < p2.pos.z;
            });
        split_axis = 2;
        split_value = photons.at(photons.size() / 2).pos.z + EPSILON;
        min1 = Vec3(min.x, min.y, min.z);
        max1 = Vec3(max.x, max.y, split_value);
        min2 = Vec3(min.x, min.y, split_value);
        max2 = Vec3(max.x, max.y, max.z);
    }
    std::cout << split_axis << " " << min1 << min2 << std::endl;
    //if (dx >= dy && dx >= dz) {
    //    split_axis = 0;
    //    split_value = min.x + dx / 2.0;
    //    min1 = Vec3(min.x, min.y, min.z);
    //    max1 = Vec3(split_value, max.y, max.z);
    //    min2 = Vec3(split_value, min.y, min.z);
    //    max2 = Vec3(max.x, max.y, max.z);
    //}
    //else if (dy >= dx && dy >= dz) {
    //    split_axis = 1;
    //    split_value = min.y + dy / 2.0;
    //    min1 = Vec3(min.x, min.y, min.z);
    //    max1 = Vec3(max.x, split_value, max.z);
    //    min2 = Vec3(min.x, split_value, min.z);
    //    max2 = Vec3(max.x, max.y, max.z);
    //}
    //else {
    //    assert(dz >= dx && dz >= dy);
    //    split_axis = 2;
    //    split_value = min.z + dz / 2.0;
    //    min1 = Vec3(min.x, min.y, min.z);
    //    max1 = Vec3(max.x, max.y, split_value);
    //    min2 = Vec3(min.x, min.y, split_value);
    //    max2 = Vec3(max.x, max.y, max.z);
    //}
    // create two new children
    child1 = new KDTree(BoundingBox(min1, max1), depth + 1);
    child2 = new KDTree(BoundingBox(min2, max2), depth + 1);
    int num_photons = photons.size();
    std::vector<Photon> tmp = photons;
    photons.clear();
    // add all the photons to one of those children
    for (int i = 0; i < num_photons; i++) {
        const Photon& p = tmp[i];
        this->AddPhoton(p);
    }
}

// ==================================================================
