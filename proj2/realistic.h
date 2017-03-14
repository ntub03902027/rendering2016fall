
#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "camera.h"
#include "paramset.h"
#include "film.h"

#include <vector>

class RealisticCameraLens {
public:
	float radius;
	float thick;
	float n_d;
	float aperture;
	float center_z;
	bool isAperture;
};

// RealisticCamera Declarations
class RealisticCamera : public Camera {
public:
	// RealisticCamera Public Methods
	RealisticCamera(const AnimatedTransform &cam2world,
						float hither, float yon, float sopen,
						float sclose, float filmdistance, float aperture_diameter, string specfile,
						float filmdiag, Film *film);
	float GenerateRay(const CameraSample &sample, Ray *) const;
  
private:
	// RealisticCamera Private Methods
	std::vector< RealisticCameraLens > lensTable;

	Transform RasterToCamera;
    float aperture_diameter, filmdistance, filmdiag;
    float hither, yon;
    float lensTotalWidth;

    bool RunThroughAperture(Ray *ray, const RealisticCameraLens len) const;
    bool Intersect(Ray *ray, const int i, Point &p) const;
    bool Refract(Ray *ray, const int i, const Point p, const float n_next) const;

    //mutable int count;

};


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);


#endif	// PBRT_CAMERAS_REALISTIC_H