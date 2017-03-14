#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_SHAPES_HEIGHTFIELD2_H
#define PBRT_SHAPES_HEIGHTFIELD2_H

// shapes/heightfield.h*
#include "shape.h"

#include "geometry.h" // normal

// Heightfield2 Declarations
class Heightfield2 : public Shape {
public:
    // Heightfield2 Public Methods
    Heightfield2(const Transform *o2, const Transform *w2o, bool ro, int nu, int nv, const float *zs);
    ~Heightfield2();
    bool CanIntersect() const;
    void Refine(vector<Reference<Shape> > &refined) const;
    BBox ObjectBound() const;
    bool Intersect(const Ray &r, float *tHit, float *rayEpsilon,
                    DifferentialGeometry *dg) const;
    bool IntersectP(const Ray &r) const;
private:
    // Heightfield2 Private Data
    float *z;
    int nx, ny;
    Normal *n;

    // Private Functions
    bool triIntersect(const Ray &r, float *tHit, float *rayEpsilon,
                        DifferentialGeometry *dg, Point triangle[]) const;

    bool triIntersectWithShadingGeometry(const Ray &r, float *tHit, float *rayEpsilon,
                        DifferentialGeometry *dg, Point triangle[], int nxnyCoordinate[3][2]) const;
    
    bool triIntersectP(const Ray &r, Point triangle[]) const;

    
};


Heightfield2 *CreateHeightfield2Shape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params);

#endif // PBRT_SHAPES_HEIGHTFIELD2_H