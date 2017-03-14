// shapes/heightfield2.cpp*
#include "stdafx.h"
//#include "vdb.h"
#include "shapes/heightfield2.h"
#include "shapes/trianglemesh.h"
#include "paramset.h"

// perhaps, used in 2D Grids. perhaps.
#include "core/primitive.h"
#include "accelerators/grid.h"

#include "geometry.h"

#include <cmath> // floor(), ceiling() etc
#include <cstdio>



// Heightfield2 Method Definitions
Heightfield2::Heightfield2(const Transform *o2w, const Transform *w2o,
        bool ro, int x, int y, const float *zs)
    : Shape(o2w, w2o, ro) {
    nx = x;
    ny = y;
    z = new float[nx*ny];
    memcpy(z, zs, nx*ny*sizeof(float));

    n = new Normal[nx*ny];



    Vector **normalUpper = new Vector*[(nx-1)*(ny-1)];
    Vector **normalLower = new Vector*[(nx-1)*(ny-1)];



    for (int y_coord = 0; y_coord < (ny-1); y_coord++) {
        for (int x_coord = 0; x_coord < (nx-1); x_coord++) {
            Point p[4];

            p[0].x = (float)x_coord / (float)(nx-1), p[0].y = (float)y_coord / (float)(ny-1), p[0].z = z[y_coord * nx + x_coord];
            p[1].x = (float)(x_coord + 1)/ (float)(nx-1), p[1].y = (float)y_coord/ (float)(ny-1), p[1].z = z[y_coord * nx + x_coord + 1];
            p[2].x = (float)(x_coord)/ (float)(nx-1), p[2].y = (float)(y_coord + 1)/ (float)(ny-1), p[2].z = z[(y_coord+1) * nx + x_coord];
            p[3].x = (float)(x_coord + 1)/ (float)(nx-1), p[3].y = (float)(y_coord + 1)/ (float)(ny-1), p[3].z = z[(y_coord+1) * nx + x_coord + 1];

            normalUpper[y_coord * nx + x_coord] = new Vector(Cross(p[3]-p[0], p[2]-p[0]));
            normalLower[y_coord * nx + x_coord] = new Vector(Cross(p[1]-p[0], p[3]-p[0]));


            /*vdb_line(p[0].x, p[0].y, p[0].z,
                p[2].x, p[2].y, p[2].z);
            vdb_line(p[3].x, p[3].y, p[3].z,
                p[2].x, p[2].y, p[2].z);
            vdb_line(p[0].x, p[0].y, p[0].z,
                p[3].x, p[3].y, p[3].z);

            vdb_line(p[0].x, p[0].y, p[0].z,
                p[1].x, p[1].y, p[1].z);
            vdb_line(p[3].x, p[3].y, p[3].z,
                p[1].x, p[1].y, p[1].z);*/

        }
    }


    n[0].x = normalUpper[0]->x + normalLower[0]->x;
    n[0].y = normalUpper[0]->y + normalLower[0]->y;
    n[0].x = normalUpper[0]->z + normalLower[0]->z;
    n[nx-1].x = normalLower[nx-2]->x;
    n[nx-1].y = normalLower[nx-2]->y;
    n[nx-1].z = normalLower[nx-2]->z;
    n[(ny-1)*nx+0].x = normalUpper[(ny-2)*(nx-1)+0]->x;
    n[(ny-1)*nx+0].y = normalUpper[(ny-2)*(nx-1)+0]->y;
    n[(ny-1)*nx+0].z = normalUpper[(ny-2)*(nx-1)+0]->z;
    n[(ny-1)*nx+nx-1].x = normalUpper[(ny-2)*(nx-1)+nx-2]->x + normalLower[(ny-2)*(nx-1)+nx-2]->x;
    n[(ny-1)*nx+nx-1].y = normalUpper[(ny-2)*(nx-1)+nx-2]->y + normalLower[(ny-2)*(nx-1)+nx-2]->y;
    n[(ny-1)*nx+nx-1].z = normalUpper[(ny-2)*(nx-1)+nx-2]->z + normalLower[(ny-2)*(nx-1)+nx-2]->z;


    for (int x_coord = 1; x_coord < (nx-1); x_coord++) {
        n[x_coord].x = normalLower[x_coord - 1]->x + normalUpper[x_coord]->x + normalLower[x_coord]->x;
        n[x_coord].y = normalLower[x_coord - 1]->y + normalUpper[x_coord]->y + normalLower[x_coord]->y;
        n[x_coord].z = normalLower[x_coord - 1]->z + normalUpper[x_coord]->z + normalLower[x_coord]->z;
        n[(ny-1)*nx + x_coord].x = normalUpper[(ny-2)*nx + x_coord - 1]->x + normalLower[(ny-2)*nx + x_coord]->x + normalUpper[(ny-2)*nx + x_coord]->x;
        n[(ny-1)*nx + x_coord].y = normalUpper[(ny-2)*nx + x_coord - 1]->y + normalLower[(ny-2)*nx + x_coord]->y + normalUpper[(ny-2)*nx + x_coord]->y;
        n[(ny-1)*nx + x_coord].z = normalUpper[(ny-2)*nx + x_coord - 1]->z + normalLower[(ny-2)*nx + x_coord]->z + normalUpper[(ny-2)*nx + x_coord]->z;
    }
   for (int y_coord = 1; y_coord < (ny-1); y_coord++) {
        n[y_coord * nx].x = normalLower[(y_coord - 1) * nx]->x + normalUpper[(y_coord) * nx]->x + normalLower[(y_coord) * nx]->x;
        n[y_coord * nx].y = normalLower[(y_coord - 1) * nx]->y + normalUpper[(y_coord) * nx]->y + normalLower[(y_coord) * nx]->y;
        n[y_coord * nx].z = normalLower[(y_coord - 1) * nx]->z + normalUpper[(y_coord) * nx]->z + normalLower[(y_coord) * nx]->z;
        n[y_coord * nx + (nx-1)].x = normalLower[(y_coord-1) * nx + (nx-2)]->x + normalUpper[(y_coord-1) * nx + (nx-2)]->x + normalLower[(y_coord) * nx + (nx-2)]->x;
        n[y_coord * nx + (nx-1)].y = normalLower[(y_coord-1) * nx + (nx-2)]->y + normalUpper[(y_coord-1) * nx + (nx-2)]->y + normalLower[(y_coord) * nx + (nx-2)]->y;
        n[y_coord * nx + (nx-1)].z = normalLower[(y_coord-1) * nx + (nx-2)]->z + normalUpper[(y_coord-1) * nx + (nx-2)]->z + normalLower[(y_coord) * nx + (nx-2)]->z;
    }

    for (int y_coord = 1; y_coord < ny-1; y_coord++) {
        for (int x_coord = 1; x_coord < nx-1; x_coord++) {
            n[y_coord * nx + x_coord].x = normalLower[(y_coord-1) * nx + (x_coord-1)]->x + normalUpper[(y_coord-1) * nx + (x_coord-1)]->x + normalLower[(y_coord) * nx + (x_coord-1)]->x + 
                                                normalLower[(y_coord) * nx + (x_coord)]->x + normalUpper[(y_coord) * nx + (x_coord)]->x + normalUpper[(y_coord-1) * nx + (x_coord)]->x;
            n[y_coord * nx + x_coord].y = normalLower[(y_coord-1) * nx + (x_coord-1)]->y + normalUpper[(y_coord-1) * nx + (x_coord-1)]->y + normalLower[(y_coord) * nx + (x_coord-1)]->y + 
                                                normalLower[(y_coord) * nx + (x_coord)]->y + normalUpper[(y_coord) * nx + (x_coord)]->y + normalUpper[(y_coord-1) * nx + (x_coord)]->y;
            n[y_coord * nx + x_coord].z = normalLower[(y_coord-1) * nx + (x_coord-1)]->z + normalUpper[(y_coord-1) * nx + (x_coord-1)]->z + normalLower[(y_coord) * nx + (x_coord-1)]->z + 
                                                normalLower[(y_coord) * nx + (x_coord)]->z + normalUpper[(y_coord) * nx + (x_coord)]->z + normalUpper[(y_coord-1) * nx + (x_coord)]->z;
        }
    }
    /*
    vdb_color(0, 1.0, 0);
    for (int y_coord = 0; y_coord < ny; y_coord++) {
        for (int x_coord = 0; x_coord < nx; x_coord++) {
            vdb_line(x_coord / (float)(nx-1), y_coord / (float)(ny-1), z[y_coord*nx+x_coord], 
               n[y_coord * nx + x_coord].x/(float)(nx-1) + (x_coord / (float)(nx-1)),
                n[y_coord * nx + x_coord].y/(float)(ny-1) + (y_coord / (float)(ny-1)),
                 n[y_coord * nx + x_coord].z + (z[y_coord*nx+x_coord]));

        }

    }*/


 



}


Heightfield2::~Heightfield2() {
    delete[] z;
}


BBox Heightfield2::ObjectBound() const {
    float minz = z[0], maxz = z[0];
    for (int i = 1; i < nx*ny; ++i) {
        if (z[i] < minz) minz = z[i];
        if (z[i] > maxz) maxz = z[i];
    }

    // draw bounding box
    /*vdb_color(0.0, 0.0, 1.0);
    vdb_line(0, 0, minz,
            1, 0, minz);
    vdb_line(0, 0, minz,
            0, 1, minz);
    vdb_line(1, 1, minz,
            0, 1, minz);
    vdb_line(1, 1, minz,
            1, 0, minz);

    vdb_line(0, 0, maxz,
            1, 0, maxz);
    vdb_line(0, 0, maxz,
            0, 1, maxz);
    vdb_line(1, 1, maxz,
            0, 1, maxz);
    vdb_line(1, 1, maxz,
            1, 0, maxz);

    vdb_line(1, 0, maxz,
            1, 0, minz);
    vdb_line(0, 0, maxz,
            0, 0, minz);
    vdb_line(1, 1, maxz,
            1, 1, minz);
    vdb_line(0, 1, maxz,
            0, 1, minz);*/


    return BBox(Point(0,0,minz), Point(1,1,maxz));
}


bool Heightfield2::CanIntersect() const {
    return true;    // PJ1 Step 2
}
/*
    About Coordinates in the following function(s):

    nx-ny coordinate              object coordinate (0-1 coordinate, since is bounded between 0 and 1)     world coordinate
    
    maxdelta                             ray                                                                      r
    m                                    heightfield2
    deltax                               p[4]
    current                              triangle1, triangle2
    previous
    floor_p

                  
            divided by (nx-1), (ny-1), 1                                                      (*ObjectToWorld)(...)
            ------------------------->                                                        ----------------->

            multiplied by (nx-1), (ny-1), 1                                                   (*WorldToObject)(...)  
            <-------------------------                                                        <-----------------
*/

// PJ1 Step 2:
// TODO:
// 0. Ray transform
// 1. Apply DDA
// 2. Check Intersection on two triangles created by four points within a x-y pixel
bool Heightfield2::Intersect(const Ray &r, float *tHit, float *rayEpsilon,
                        DifferentialGeometry *dg) const {

    // Ray transform
    Ray ray;
    (*WorldToObject) (r, &ray);

    /*if (rayCount % 500 == 0) {
    vdb_color(1.0, 0, 0);
    vdb_point(ray.o.x, ray.o.y, ray.o.z);
    vdb_line(ray.o.x, ray.o.y, ray.o.z,
            ray.o.x + ray.d.x, ray.o.y + ray.d.y, ray.o.z + ray.d.z);  
    }*/


    float maxdelta = 0.0625; // max{|deltax|, |deltay|} = max{|deltax|, |m * deltax|} = const. , in nx * ny coordinate

    float deltax = maxdelta * ray.d.x / sqrt(ray.d.x * ray.d.x + ray.d.y * ray.d.y + ray.d.z * ray.d.z);
    float deltay = maxdelta * ray.d.y / sqrt(ray.d.x * ray.d.x + ray.d.y * ray.d.y + ray.d.z * ray.d.z);
    float deltaz = maxdelta * ray.d.z / sqrt(ray.d.x * ray.d.x + ray.d.y * ray.d.y + ray.d.z * ray.d.z);


  /*  float m = (ray.d.y * (float)(ny-1)) / (ray.d.x * (float)(nx-1));
    float deltax;
    if (m > 1.0) {deltax = maxdelta / m;}
    else {deltax = maxdelta;}
    if (ray.d.x < 0.0) {deltax = -deltax;} // determine traverse direction
    */
    

    /*if (rayCount % 450 == 0) 
        fprintf(stderr, "\nm = %f, deltax = %f\n", m, deltax);*/

    // start DDA (2D)

    Point previous; // to record previously-examined point, which records "the floor of x and y"
    previous.x = -1.0; previous.y = -1.0, previous.z = -1.0; // initialize previous

    BBox boundBox = ObjectBound();

    float rayT;
    if (boundBox.Inside(ray(ray.mint)))
        rayT = ray.mint;
    else if (!boundBox.IntersectP(ray, &rayT))
        return false;

    Point current;
    current.x = (ray.o.x + rayT * ray.d.x) * (float)(nx-1), current.y = (ray.o.y + rayT * ray.d.y) * (float)(ny-1), current.z = ray.o.z + rayT * ray.d.z;

    bool doIntersect = false; // The return flag

    for (; current.x >= boundBox.pMin.x && current.x <= (float)(nx-1) * boundBox.pMax.x && 
        current.y >= boundBox.pMin.y && current.y <= (float)(ny-1) * boundBox.pMax.y /*&&
        current.z >= boundBox.pMin.z && current.z <= boundBox.pMax.z */;current.x += deltax, current.y += deltay, current.z += deltaz) {
    
        if (floor(current.x) == previous.x && floor(current.y) == previous.y )
            continue;
/*
p2(a, b+1)    p3(a+1, b+1)
.___________.
|          /|
| tri-1  /  |
|     /     |
|  /  tri-2 |
./__________.
p0(a,b)       p1(a+1, b)


*/
        Point floor_p;
        floor_p.x = floor(current.x), floor_p.y = floor(current.y);
        Point p[4];

        p[0].x = floor_p.x / (float)(nx-1), p[0].y = floor_p.y / (float)(ny-1), p[0].z = z[(int)floor_p.y * nx + (int)floor_p.x];
        p[1].x = (floor_p.x + 1.0) / (float)(nx-1), p[1].y = floor_p.y / (float)(ny-1), p[1].z = z[(int)floor_p.y * nx + (int)floor_p.x + 1];
        p[2].x = floor_p.x / (float)(nx-1), p[2].y = (floor_p.y + 1.0) / (float)(ny-1), p[2].z = z[((int)floor_p.y + 1) * nx + (int)floor_p.x];
        p[3].x = (floor_p.x + 1.0) / (float)(nx-1), p[3].y = (floor_p.y + 1.0) / (float)(ny-1), p[3].z = z[((int)floor_p.y + 1) * nx + (int)floor_p.x + 1];


        // 3 points form a triangle, in total 2 triangles
        Point triangle1[3] = {p[0], p[2], p[3]};
        Point triangle2[3] = {p[0], p[3], p[1]};


        int nxnyCoordinate1[3][2] = { {(int)floor_p.x, (int)floor_p.y}, {(int)floor_p.x, (int)floor_p.y + 1}, {(int)floor_p.x + 1, (int)floor_p.y + 1} };
        int nxnyCoordinate2[3][2] = { {(int)floor_p.x, (int)floor_p.y}, {(int)floor_p.x + 1, (int)floor_p.y + 1}, {(int)floor_p.x + 1, (int)floor_p.y} };

        // Call a function specialized in checking triangle-ray intersection, tri-1 and tri-2 respectively
        if (triIntersectWithShadingGeometry(r, tHit, rayEpsilon, dg, triangle1, nxnyCoordinate1)) {
            doIntersect = true;
            break;
        }
        else if (triIntersectWithShadingGeometry(r, tHit, rayEpsilon, dg, triangle2, nxnyCoordinate2)) {
            doIntersect = true;
            break;
        }

/*
        if (triIntersect(r, tHit, rayEpsilon, dg, triangle1)) {
            doIntersect = true;
            break;
        }
        else if (triIntersect(r, tHit, rayEpsilon, dg, triangle2)) {
            doIntersect = true;
            break;
        }

*/


        // update previous point
        previous.x = floor_p.x;
        previous.y = floor_p.y;
    }


    return doIntersect; 
}
// Derived from trianglemesh.cpp: Triangle::Intersect() and Triangle::GetShadingGeometry()
bool Heightfield2::triIntersectWithShadingGeometry(const Ray &ray, float *tHit, float *rayEpsilon,
                        DifferentialGeometry *dg, Point triangle[], int nxnyCoordinate[3][2]) const {

    const Point p1 = (*ObjectToWorld)(triangle[0]);
    const Point p2 = (*ObjectToWorld)(triangle[1]);
    const Point p3 = (*ObjectToWorld)(triangle[2]);

    Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector s1 = Cross(ray.d, e2);
    float divisor = Dot(s1, e1);

    if (divisor == 0.)
        return false;
    float invDivisor = 1.f / divisor;

    // Compute first barycentric coordinate
    Vector s = ray.o - p1;
    float b1 = Dot(s, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;

    // Compute second barycentric coordinate
    Vector s2 = Cross(s, e1);
    float b2 = Dot(ray.d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;

    // Compute _t_ to intersection point
    float t = Dot(e2, s2) * invDivisor;
    if (t < ray.mint || t > ray.maxt)
        return false;

    // Now, the ray does intersects with the triangle.

    // Compute triangle partial derivatives
    Vector dpdu, dpdv;
    

    // Compute deltas for triangle partial derivatives
    float du1 = triangle[0].x - triangle[2].x;
    float du2 = triangle[1].x - triangle[2].x;
    float dv1 = triangle[0].y - triangle[2].y;
    float dv2 = triangle[1].y - triangle[2].y;
    Vector dp1 = p1 - p3, dp2 = p2 - p3;
    float determinant = du1 * dv2 - dv1 * du2;
    if (determinant == 0.f) {
        // Handle zero determinant for triangle partial derivative matrix
        CoordinateSystem(Normalize(Cross(e2, e1)), &dpdu, &dpdv);
    }
    else {
        float invdet = 1.f / determinant;
        dpdu = ( dv2 * dp1 - dv1 * dp2) * invdet;
        dpdv = (-du2 * dp1 + du1 * dp2) * invdet;
    }

    // Interpolate $(u,v)$ triangle parametric coordinates
    float b0 = 1 - b1 - b2;
    float tu = b0*triangle[0].x + b1*triangle[1].x + b2*triangle[2].x;
    float tv = b0*triangle[0].y + b1*triangle[1].y + b2*triangle[2].y;

    // derived from trianglemesh.cpp: Triangle::GetShadingGeometry()

    Normal ns;
    Vector ss, ts;

    ns = Normalize(b0*n[nxnyCoordinate[0][1] * nx + nxnyCoordinate[0][0]] + b1*n[nxnyCoordinate[1][1] * nx + nxnyCoordinate[1][0]] + b2*n[nxnyCoordinate[2][1] * nx + nxnyCoordinate[2][0]]);
    ss = Normalize(dpdu);

    if (ts.LengthSquared() > 0.f) {
        ts = Normalize(ts);
        ss = Cross(ts, ns);
    }
    else
        CoordinateSystem((Vector)ns, &ss, &ts);

    Normal dndu, dndv;
    // Compute $\dndu$ and $\dndv$ for triangle shading geometry

    Normal dn1 = n[nxnyCoordinate[0][1] * nx + nxnyCoordinate[0][0]] - n[nxnyCoordinate[2][1] * nx + nxnyCoordinate[2][0]];
    Normal dn2 = n[nxnyCoordinate[1][1] * nx + nxnyCoordinate[1][0]] - n[nxnyCoordinate[2][1] * nx + nxnyCoordinate[2][0]];
    if (determinant == 0.f)
        dndu = dndv = Normal(0,0,0);
    else {
        float invdet = 1.f / determinant;
        dndu = ( dv2 * dn1 - dv1 * dn2) * invdet;
        dndv = (-du2 * dn1 + du1 * dn2) * invdet;
    }
    *dg = DifferentialGeometry(ray(t), ss, ts,
                               (*ObjectToWorld)(dndu), (*ObjectToWorld)(dndv),
                               tu, tv, this);

    //end of sections added from trianglemesh.cpp: Triangle::GetShadingGeometry()


    *tHit = t;
    *rayEpsilon = 1e-3f * *tHit;

    return true;


}

// Derived from trianglemesh.cpp: Triangle::Intersect()
bool Heightfield2::triIntersect(const Ray &ray, float *tHit, float *rayEpsilon,
                        DifferentialGeometry *dg, Point triangle[]) const {

    const Point p1 = (*ObjectToWorld)(triangle[0]);
    const Point p2 = (*ObjectToWorld)(triangle[1]);
    const Point p3 = (*ObjectToWorld)(triangle[2]);

    Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector s1 = Cross(ray.d, e2);
    float divisor = Dot(s1, e1);

    if (divisor == 0.)
        return false;
    float invDivisor = 1.f / divisor;

    // Compute first barycentric coordinate
    Vector s = ray.o - p1;
    float b1 = Dot(s, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;

    // Compute second barycentric coordinate
    Vector s2 = Cross(s, e1);
    float b2 = Dot(ray.d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;

    // Compute _t_ to intersection point
    float t = Dot(e2, s2) * invDivisor;
    if (t < ray.mint || t > ray.maxt)
        return false;

    // Now, the ray does intersects with the triangle.

    // Compute triangle partial derivatives
    Vector dpdu, dpdv;
    

    // Compute deltas for triangle partial derivatives
    float du1 = triangle[0].x - triangle[2].x;
    float du2 = triangle[1].x - triangle[2].x;
    float dv1 = triangle[0].y - triangle[2].y;
    float dv2 = triangle[1].y - triangle[2].y;
    Vector dp1 = p1 - p3, dp2 = p2 - p3;
    float determinant = du1 * dv2 - dv1 * du2;
    if (determinant == 0.f) {
        // Handle zero determinant for triangle partial derivative matrix
        CoordinateSystem(Normalize(Cross(e2, e1)), &dpdu, &dpdv);
    }
    else {
        float invdet = 1.f / determinant;
        dpdu = ( dv2 * dp1 - dv1 * dp2) * invdet;
        dpdv = (-du2 * dp1 + du1 * dp2) * invdet;
    }

    // Interpolate $(u,v)$ triangle parametric coordinates
    float b0 = 1 - b1 - b2;
    float tu = b0*triangle[0].x + b1*triangle[1].x + b2*triangle[2].x;
    float tv = b0*triangle[0].y + b1*triangle[1].y + b2*triangle[2].y;



    // Fill in _DifferentialGeometry_ from triangle hit
    *dg = DifferentialGeometry(ray(t), dpdu, dpdv,
                               Normal(0,0,0), Normal(0,0,0),
                               tu, tv, this);
    *tHit = t;
    *rayEpsilon = 1e-3f * *tHit;

    return true;


}


bool Heightfield2::IntersectP(const Ray &r) const {
    Ray ray;
    (*WorldToObject) (r, &ray);

    float maxdelta = 0.0625;
    float deltax = maxdelta * ray.d.x / sqrt(ray.d.x * ray.d.x + ray.d.y * ray.d.y + ray.d.z * ray.d.z);
    float deltay = maxdelta * ray.d.y / sqrt(ray.d.x * ray.d.x + ray.d.y * ray.d.y + ray.d.z * ray.d.z);
    float deltaz = maxdelta * ray.d.z / sqrt(ray.d.x * ray.d.x + ray.d.y * ray.d.y + ray.d.z * ray.d.z);

    Point previous; 
    previous.x = -1.0; previous.y = -1.0, previous.z = -1.0;

    BBox boundBox = ObjectBound();
    float rayT;
    if (boundBox.Inside(ray(ray.mint)))
        rayT = ray.mint;
    else if (!boundBox.IntersectP(ray, &rayT))
        return false;

    Point current;
    current.x = (ray.o.x + rayT * ray.d.x) * (float)(nx-1), current.y = (ray.o.y + rayT * ray.d.y) * (float)(ny-1), current.z = ray.o.z + rayT * ray.d.z;

    bool doIntersect = false; 

    for (; current.x >= boundBox.pMin.x && current.x <= (float)(nx-1) * boundBox.pMax.x && 
        current.y >= boundBox.pMin.y && current.y <= (float)(ny-1) * boundBox.pMax.y ;current.x += deltax, current.y += deltay, current.z += deltaz) {
    
        if (floor(current.x) == previous.x && floor(current.y) == previous.y )
            continue;

        Point floor_p;
        floor_p.x = floor(current.x), floor_p.y = floor(current.y);
        Point p[4];

        p[0].x = floor_p.x / (float)(nx-1), p[0].y = floor_p.y / (float)(ny-1), p[0].z = z[(int)floor_p.y * nx + (int)floor_p.x];
        p[1].x = (floor_p.x + 1.0) / (float)(nx-1), p[1].y = floor_p.y / (float)(ny-1), p[1].z = z[(int)floor_p.y * nx + (int)floor_p.x + 1];
        p[2].x = floor_p.x / (float)(nx-1), p[2].y = (floor_p.y + 1.0) / (float)(ny-1), p[2].z = z[((int)floor_p.y + 1) * nx + (int)floor_p.x];
        p[3].x = (floor_p.x + 1.0) / (float)(nx-1), p[3].y = (floor_p.y + 1.0) / (float)(ny-1), p[3].z = z[((int)floor_p.y + 1) * nx + (int)floor_p.x + 1];

        Point triangle1[3] = {p[0], p[2], p[3]};
        Point triangle2[3] = {p[0], p[3], p[1]};

        if (triIntersectP(r, triangle1)) {
            doIntersect = true;
            break;
        }
        else if (triIntersectP(r, triangle2)) {
            doIntersect = true;
            break;
        }

        previous.x = floor_p.x;
        previous.y = floor_p.y;
    }

    return doIntersect; 
}

bool Heightfield2::triIntersectP(const Ray &ray, Point triangle[]) const {
    const Point p1 = (*ObjectToWorld)(triangle[0]);
    const Point p2 = (*ObjectToWorld)(triangle[1]);
    const Point p3 = (*ObjectToWorld)(triangle[2]);

    Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector s1 = Cross(ray.d, e2);
    float divisor = Dot(s1, e1);

    if (divisor == 0.)
        return false;
    float invDivisor = 1.f / divisor;

    // Compute first barycentric coordinate
    Vector s = ray.o - p1;
    float b1 = Dot(s, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;

    // Compute second barycentric coordinate
    Vector s2 = Cross(s, e1);
    float b2 = Dot(ray.d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;

    // Compute _t_ to intersection point
    float t = Dot(e2, s2) * invDivisor;
    if (t < ray.mint || t > ray.maxt)
        return false;
    return true;
}

void Heightfield2::Refine(vector<Reference<Shape> > &refined) const {
    int ntris = 2*(nx-1)*(ny-1);
    refined.reserve(ntris);
    int *verts = new int[3*ntris];
    Point *P = new Point[nx*ny];
    float *uvs = new float[2*nx*ny];
    int nverts = nx*ny;
    int x, y;
    // Compute heightfield vertex positions
    int pos = 0;
    for (y = 0; y < ny; ++y) {
        for (x = 0; x < nx; ++x) {
            P[pos].x = uvs[2*pos]   = (float)x / (float)(nx-1);
            P[pos].y = uvs[2*pos+1] = (float)y / (float)(ny-1);
            P[pos].z = z[pos];
            ++pos;
        }
    }

    // Fill in heightfield vertex offset array
    int *vp = verts;
    for (y = 0; y < ny-1; ++y) {
        for (x = 0; x < nx-1; ++x) {
#define VERT(x,y) ((x)+(y)*nx)
            *vp++ = VERT(x, y);
            *vp++ = VERT(x+1, y);
            *vp++ = VERT(x+1, y+1);
    
            *vp++ = VERT(x, y);
            *vp++ = VERT(x+1, y+1);
            *vp++ = VERT(x, y+1);
        }
#undef VERT
    }
    ParamSet paramSet;
    paramSet.AddInt("indices", verts, 3*ntris);
    paramSet.AddFloat("uv", uvs, 2 * nverts);
    paramSet.AddPoint("P", P, nverts);
    refined.push_back(CreateTriangleMeshShape(ObjectToWorld, WorldToObject, ReverseOrientation, paramSet));
    delete[] P;
    delete[] uvs;
    delete[] verts;
}


Heightfield2 *CreateHeightfield2Shape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params) {
    int nu = params.FindOneInt("nu", -1);
    int nv = params.FindOneInt("nv", -1);
    int nitems;
    const float *Pz = params.FindFloat("Pz", &nitems);
    Assert(nitems == nu*nv);
    Assert(nu != -1 && nv != -1 && Pz != NULL);
    return new Heightfield2(o2w, w2o, reverseOrientation, nu, nv, Pz);
}