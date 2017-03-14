
#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_LIGHTFIELD_H
#define PBRT_CAMERAS_LIGHTFIELD_H

#include "camera.h"
#include "paramset.h"
#include "film.h"

#include <vector>

// for EXRread/write
#include <ImfInputFile.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>
#include <ImfRgbaFile.h>
#include <half.h>
#include <algorithm>

using namespace Imf;
using namespace Imath;

//#define MONTAGE

class LightfieldCameraLens {
public:
	float radius;
	float thick;
	float n_d;
	float aperture;
	float center_z;
	bool isAperture;
};

// LightfieldCamera Declarations
class LightfieldCamera : public Camera {
public:
	// LightfieldCamera Public Methods
	LightfieldCamera(const AnimatedTransform &cam2world,
						float hither, float yon, float sopen,
						float sclose, float filmdistance, float aperture_diameter, string specfile,
						float filmdiag, int microlensdiameter, std::vector<float> &alphatable, Film *f);
	float GenerateRay(const CameraSample &sample, Ray *) const;

	void ReconstructFilms(const string filename);

private:
	// LightfieldCamera Private Methods
	std::vector< LightfieldCameraLens > lensTable;

	Transform RasterToCamera;
    float aperture_diameter, filmdistance, filmdiag;
    float hither, yon;
    float lensTotalWidth;
    int microlensdiameter;
    float filmResRatio;
    std::vector<float> alpha; 


#ifdef MONTAGE
    int resX, resY;
#endif

    bool RunThroughAperture(Ray *ray, const LightfieldCameraLens len) const;
    bool Intersect(Ray *ray, const int i, Point &p) const;
    bool Refract(Ray *ray, const int i, const Point p, const float n_next) const;

    void WriteEXR(const char *name, float *pixels, int xRes, int yRes);
    bool ReadEXR(const char *name, float *&rgba, int &xRes, int &yRes, bool &hasAlpha);

    //mutable int count;

};


LightfieldCamera *CreateLightfieldCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);


#endif	// PBRT_CAMERAS_LIGHTFIELD_H