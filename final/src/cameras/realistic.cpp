
#include "stdafx.h"
#include "cameras/realistic.h"

#include <cstdio>
#include <cstring>


#include <cmath>
#include "sampler.h"
#include "core/montecarlo.h"

//#include "vdb.h"

#define _USE_MATH_DEFINES 

#define RAYCOUNT_DEBUG 1000

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
				 float hither, float yon, 
				 float sopen, float sclose, 
				 float filmdistance, float aperture_diameter, string specfile, 
				 float filmdiag, Film *f)
	: Camera(cam2world, sopen, sclose, f) // pbrt-v2 doesnot specify hither and yon
{
   // YOUR CODE HERE -- build and store datastructures representing the given lens
   // and film placement.

// Fill in some variables

	this->aperture_diameter = aperture_diameter;
	this->filmdistance = filmdistance;
	this->filmdiag = filmdiag;
	this->hither = hither;
	this->yon = yon;

// Calculate some transformations, i.e. rastertocamera

	RasterToCamera = Translate(Vector(filmdiag / hypot(f->xResolution, f->yResolution) * f->xResolution * 0.5f, 
			- (filmdiag / hypot(f->xResolution, f->yResolution) * f->yResolution * 0.5f), 0.f)) // NDC-To-Camera
		* Scale( filmdiag / hypot(f->xResolution, f->yResolution) ,
                           filmdiag / hypot(f->xResolution, f->yResolution) , 1.f) // Raster-To-NDC
		* Scale(-1.f, 1.f, 1.f);

// Read data from lens specfile
	FILE *lensData;

	float lensdistance = 0.0;

	char inputLine[1024];
	lensTable.clear();
	lensData = fopen(specfile.c_str(), "r");
	if (lensData == NULL) perror("Error opening file");
	else {
		while (fgets (inputLine, 1024, lensData) != NULL) {
			for (int i = 0; i < (int) strlen(inputLine); i++) {
				if (inputLine[i] == '#') { 
					inputLine[i] = '\0'; 
					bool flag = true;
					for (int j = 0; j < i; j++) {
						if (inputLine[j] != ' ') {flag = false; break;}
					}
					if (flag) inputLine[0] = '\0'; 
					break; }
			}
			if (strlen(inputLine) > 0) { // assume correct syntax
				RealisticCameraLens len;
				sscanf(inputLine, "%f %f %f %f", &len.radius, &len.thick, &len.n_d, &len.aperture);
				if (len.n_d == 0.0) {
					len.n_d = 1;
					len.isAperture = true;
				}
				else {len.isAperture = false;}
				lensTable.push_back(len);
				lensdistance += len.thick;
			}
		}
		float temp = lensdistance;
		for (int i = 0; i < (int) lensTable.size(); i++) {

			lensTable[i].center_z = filmdistance + temp;
			temp -= lensTable[i].thick;
		}

		// for debug
		/*for (int i = 0; i < (int) lensTable.size(); i++) {
			printf ("r = %f, thick = %f, n_d = %f, ap = %f, z = %f\n", lensTable[i].radius, lensTable[i].thick, lensTable[i].n_d, lensTable[i].aperture, lensTable[i].center_z);
		}*/
	}
	fclose(lensData);
	lensTotalWidth = lensdistance;



// debug vdb
/*
	count = -1;
	vdb_color(0, 0, 1.0);
	vdb_line(0, 0, -10, 0, 0, lensdistance + filmdistance);
	vdb_color(1.0, 0, 0);
	vdb_line(-20, 0, 0, 20, 0, 0);
	vdb_point(20, 0, 0);
	vdb_color(0, 1.0, 0);
	vdb_line(0, -20, 0, 0, 20, 0);
	vdb_point(0, 20, 0);
	vdb_color(1, 1, 0);
*/	
}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
  // YOUR CODE HERE -- make that ray!
  // use sample->imageX and sample->imageY to get raster-space coordinates
  // of the sample point on the film.
  // use sample->lensU and sample->lensV to get a sample position on the lens

// create initial ray
	//count++; //debug vdb
	Point Pras(sample.imageX, sample.imageY, 0);	// raster space coord (imageX, imageY, 0)
	Point Pcamera;
    RasterToCamera(Pras, &Pcamera);

    float lensU, lensV;
    float lastLenRadius = lensTable[lensTable.size() - 1].aperture * 0.5;
    ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);

    lensU *= lastLenRadius;
    lensV *= lastLenRadius;
    float lensW = filmdistance;
    

    ray->o = Pcamera;
    ray->d = Normalize(Point(lensU, lensV, lensW) - Pcamera);

    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    ray->mint = 0.f;
    ray->maxt = INFINITY;

// save cos theta
    float cos_theta = ray->d.z ; // assume normalized
// iteration: trace ray len by len
    for (int i = lensTable.size() - 1; i >= 0; i--) {
    	if (lensTable[i].isAperture) {
    		if (!RunThroughAperture(ray, lensTable[i]))
    			return 0.f;
    	}
    	else {
    		Point p; // the intersect point
    		if (!Intersect(ray, i, p))
    			return 0.f;

    		float n_next = (i == 0)? 1 : lensTable[i - 1].n_d;
    		if (!Refract(ray, i, p, n_next))
    			return 0.f;
    	}
    }

// compute weight: A * cos^4 theta / filmdistance^2
    float A = 0.25 * lensTable[lensTable.size() - 1].aperture * lensTable[lensTable.size() - 1].aperture * M_PI;
    float weight = A * cos_theta * cos_theta / (filmdistance);
    weight *=  (cos_theta * cos_theta / filmdistance);

    ray->maxt = (yon - hither) / ray->d.z;
    ray->mint = 0.f;
    ray->o.z  -= (lensTotalWidth + filmdistance);
	CameraToWorld(*ray, ray);
	ray->d = Normalize(ray->d);

	/*float a = ray->o.x * ray->o.x + ray->o.y + ray->o.y;
	float r = 0.5 * lensTable[lensTable.size() - 1].aperture;
	float Z = filmdistance;

	float f = -(a * a + Z * Z - r * r);
	f /= sqrt( (a * a + Z * Z + r * r) * (a * a + Z * Z + r * r) - 4 * r * r * a * a );
	f += 1.0;
	f *= 0.5;*/


    return weight;


   return 0;
}

bool RealisticCamera::RunThroughAperture(Ray *ray, const RealisticCameraLens len) const {
// calculate parameter t
	float t = (len.center_z - ray->o.z) / ray->d.z;
	float i_x = ray->o.x + t * ray->d.x;
	float i_y = ray->o.y + t * ray->d.y;

// check if the point on the aperture plane is out of bound
	if (i_x * i_x + i_y * i_y > 0.25 * len.aperture * len.aperture )
		return false;

// debug vdb
/*
	if (count % RAYCOUNT_DEBUG == 0) {
    	vdb_color(0.5, 0.5, 0.5);
    	vdb_point(i_x, i_y, len.center_z);
    	
	}
*/
	return true;
}

//referenced from /shapes/sphere.cpp
bool RealisticCamera::Intersect(Ray *ray, const int i, Point &p) const { 
	Point o(0.f, 0.f, lensTable[i].center_z - lensTable[i].radius); //reversed

	float A = ray->d.x*ray->d.x + ray->d.y*ray->d.y + ray->d.z*ray->d.z;
	float B = 2 * ((ray->o.x - o.x) * ray->d.x + (ray->o.y - o.y) * ray->d.y + (ray->o.z - o.z) * ray->d.z );
	float C = ( (ray->o.x - o.x) * (ray->o.x - o.x) + (ray->o.y - o.y) * (ray->o.y - o.y) + (ray->o.z - o.z) * (ray->o.z - o.z) )
		- lensTable[i].radius * lensTable[i].radius;

	float t0, t1; 
	// Solve quadratic equation for _t_ values
	if (!Quadratic(A, B, C, &t0, &t1))
		return false;
	// Compute intersection distance along ray
	/*if (t0 > ray->maxt || t1 < ray->mint)
		return false;
	float thit = t0;
    if (t0 < ray->mint) {
	        thit = t1;
        if (thit > ray->maxt) return false;
    }*/

    // make sure at least one root is >= 0
    float thit = t0;
    if (t0 < 0.0) {
    	thit = t1;
    	if (thit < 0.0) return false;
    }

    // now, get the point p

    p.x = ray->o.x + thit * ray->d.x;
    p.y = ray->o.y + thit * ray->d.y;
    p.z = ray->o.z + thit * ray->d.z;

    // examine if p is on the lens, not on the sphere whihout being on the lens
    // i.e. check if (o-p)^2 - (o.z-p.z)^2 <= (aperture * 0.5)^2
    if ( (o.x - p.x) * (o.x - p.x) + (o.y - p.y) * (o.y - p.y) > lensTable[i].aperture * lensTable[i].aperture * 0.25 )
    	return false;

// debug vdb
/*
    if (count % RAYCOUNT_DEBUG == 0) {
    	vdb_color((float)i / (float) lensTable.size(), 1.0 - (float)i /(float)lensTable.size(), 0);
    	vdb_line(ray->o.x, ray->o.y, ray->o.z, 
    	p.x, p.y, p.z);
    	vdb_color(0, (float)i /(float)lensTable.size(), 1.0 - (float)i /(float)lensTable.size());
    	vdb_point(ray->o.x, ray->o.y, ray->o.z);	
	}
*/

    return true;
}

bool RealisticCamera::Refract(Ray *ray, const int i, const Point p, const float n_next) const {
	// Apply Heckbert’s Method
	float eta = lensTable[i].n_d / n_next;
	Vector I = Normalize(ray->d);
	Point o(0.f, 0.f, lensTable[i].center_z - lensTable[i].radius);
	Vector N;
	// Handle Normal direction
	if (lensTable[i].radius >= 0)
		N = Normalize(o - p);
	else
		N = Normalize(p - o);
	float c_1 = -I.x * N.x + -I.y * N.y + -I.z * N.z;
	float c_2 = 1.0 - eta * eta * (1.0 - c_1 * c_1);
	if (c_2 < 0.f) // total reflection
		return false;
	c_2 = sqrt(c_2);

	// update ray
	ray->d = eta * I + (eta * c_1 - c_2) * N;
	ray->o = p;
	return true;
}


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
	// Extract common camera parameters from \use{ParamSet}
	float hither = params.FindOneFloat("hither", -1);
	float yon = params.FindOneFloat("yon", -1);
	float shutteropen = params.FindOneFloat("shutteropen", -1);
	float shutterclose = params.FindOneFloat("shutterclose", -1);

	// Realistic camera-specific parameters
	string specfile = params.FindOneFilename("specfile", "");  // original FindOneString 
	float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
 	float fstop = params.FindOneFloat("aperture_diameter", 1.0);	
	float filmdiag = params.FindOneFloat("filmdiag", 35.0);

	Assert(hither != -1 && yon != -1 && shutteropen != -1 &&
		shutterclose != -1 && filmdistance!= -1);
	if (specfile == "") {
	    Severe( "No lens spec file supplied!\n" );
	}
	return new RealisticCamera(cam2world, hither, yon,
				   shutteropen, shutterclose, filmdistance, fstop, 
				   specfile, filmdiag, film);
}
