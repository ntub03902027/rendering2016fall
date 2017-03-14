
#include "stdafx.h"
#include "cameras/lightfield.h"

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <sstream>

#include <cmath>
#include "sampler.h"
#include "core/montecarlo.h"

//#define DEBUG_VDB

#ifdef DEBUG_VDB
#include "vdb.h"
#endif

#define _USE_MATH_DEFINES 

#define RAYCOUNT_DEBUG 1000

LightfieldCamera::LightfieldCamera(const AnimatedTransform &cam2world,
				 float hither, float yon, 
				 float sopen, float sclose, 
				 float filmdistance, float aperture_diameter, string specfile, 
				 float filmdiag, int microlensdiameter, std::vector<float> &alphatable, Film *f)
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
	this->microlensdiameter = microlensdiameter;

	this->alpha.clear();
	for (int i = 0; i < (int)alphatable.size(); i++)
		this->alpha.push_back(alphatable[i]);


#ifdef MONTAGE
	this->resX = f->xResolution;
	this->resY = f->yResolution;
#endif

// Calculate some transformations, i.e. rastertocamera

	this->filmResRatio = filmdiag / hypot(f->xResolution, f->yResolution);

	RasterToCamera = Translate(Vector(this->filmResRatio * f->xResolution * 0.5f, - (this->filmResRatio * f->yResolution * 0.5f), 0.f)) // NDC-To-Camera
		* Scale( this->filmResRatio, this->filmResRatio , 1.f) // Raster-To-NDC
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
				LightfieldCameraLens len;
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
#ifdef DEBUG_VDB
	//count = -1;
	vdb_color(0, 0, 1.0);
	vdb_line(0, 0, -10, 0, 0, lensdistance + filmdistance);
	vdb_color(1.0, 0, 0);
	vdb_line(-20, 0, 0, 20, 0, 0);
	vdb_point(20, 0, 0);
	vdb_color(0, 1.0, 0);
	vdb_line(0, -20, 0, 0, 20, 0);
	vdb_point(0, 20, 0);
	vdb_color(1, 1, 0);
#endif
}

float LightfieldCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
  // YOUR CODE HERE -- make that ray!
  // use sample->imageX and sample->imageY to get raster-space coordinates
  // of the sample point on the film.
  // use sample->lensU and sample->lensV to get a sample position on the lens

// create initial ray
	//count++; //debug vdb

	// lightfield camera
#ifdef MONTAGE
	int montage_resX = this->resX / this->microlensdiameter;
	int montage_resY = this->resY / this->microlensdiameter;
	int montageX = (Floor2Int(sample.imageX) / montage_resX);
	int montageY = (Floor2Int(sample.imageY) / montage_resY);
	int mlenS = Floor2Int(sample.imageX - (float)(montageX * montage_resX)) % montage_resX;
	int mlenT = Floor2Int(sample.imageY - (float)(montageY * montage_resY)) % montage_resY;
	float rearrangeX = (sample.imageX - floor(sample.imageX)) + (float) (mlenS * this->microlensdiameter) + (float) montageX;
	float rearrangeY = (sample.imageY - floor(sample.imageY)) + (float) (mlenT * this->microlensdiameter) + (float) montageY;
	Point Pinholeras( Floor2Int(rearrangeX / microlensdiameter) * microlensdiameter + microlensdiameter / 2,
						Floor2Int(rearrangeY / microlensdiameter) * microlensdiameter + microlensdiameter / 2, 0);
	Point Pras(rearrangeX, rearrangeY, 0);
	Point Pinholecamera;
	Point Pcamera;
	RasterToCamera(Pinholeras, &Pinholecamera);
	RasterToCamera(Pras, &Pcamera);
#else

	Point Pinholeras( Floor2Int(sample.imageX / microlensdiameter) * microlensdiameter + microlensdiameter / 2,
						Floor2Int(sample.imageY / microlensdiameter) * microlensdiameter + microlensdiameter / 2, 0);
	Point Pinholecamera;
	RasterToCamera(Pinholeras, &Pinholecamera);


	Point Pras(sample.imageX, sample.imageY, 0);	// raster space coord (imageX, imageY, 0)

	Point Pcamera;
    RasterToCamera(Pras, &Pcamera);
#endif

    float lensU, lensV;
    float lastLenAperture = lensTable[lensTable.size() - 1].aperture;
    //ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);



    lensU = (Pcamera.x - Pinholecamera.x) * (lastLenAperture / (microlensdiameter * filmResRatio));
    lensV = (Pcamera.y - Pinholecamera.y) * (lastLenAperture / (microlensdiameter * filmResRatio));
    float lensW = filmdistance;
    
    //ray->o = Pinholecamera;
    ray->o = Pcamera; // fix at pinhole
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

bool LightfieldCamera::RunThroughAperture(Ray *ray, const LightfieldCameraLens len) const {
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
bool LightfieldCamera::Intersect(Ray *ray, const int i, Point &p) const { 
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

bool LightfieldCamera::Refract(Ray *ray, const int i, const Point p, const float n_next) const {
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

string toString(const int i) {
	if (i == 0) return "0";
	int tens = 1;
	while (i / tens > 0 && tens < 1000000000)
		tens *= 10;
	tens /= 10;
	string result;
	result.clear();
	int remains = i;
	while (tens > 0) {
		int quotient = remains / tens;
		result.push_back('0' + quotient);
		remains = remains % tens;
		tens /= 10;
	}
	return result;
}

void LightfieldCamera::ReconstructFilms(const string filename) {
	//FILE *dump;
	//dump = fopen("dump.out", "w");


	// read lightfielddata file
	float *lightfielddata;
	int res[2];
	bool hasAlpha;
	if (!ReadEXR(filename.c_str(), lightfielddata, res[0], res[1], hasAlpha)) {
		fprintf(stderr, "lightfield.cpp: error: Couldn't read image: %s\n", filename.c_str());
		return;
	}
	/*else {
		for (int i = 0; i < res[0] * res[1]; i++) {
			fprintf(stderr, "(%d, %d) %f %f %f %f\n", i / res[0], i % res[0], lightfielddata[i*4], lightfielddata[i*4+1], lightfielddata[i*4+2], lightfielddata[i*4+3]);
		}
	}*/
	//declare outfile array, res, and other neccessary data
	float lenRadius = lensTable[lensTable.size() - 1].aperture / 2;
	float lenPointDist = lenRadius / (microlensdiameter / 2);
	//float syntheticFilmPointDist = (microlensdiameter / filmResRatio);
	float lenPointMin = - lenRadius + 0.5 * lenPointDist;

	
	int outres[2] = {res[0] / microlensdiameter, res[1] / microlensdiameter};
	float *synthesizeddata = new float[outres[0] * outres[1] * 4];

	Transform CameraToRaster = Inverse(RasterToCamera);
	float filmSynRatio = filmdiag / hypot(outres[0], outres[1]);

	Transform SynRasterToCamera;

	// calculate synthetic film planes
	for (int i = 0; i < (int)this->alpha.size(); i++) {
		// decide outfilename
		string outfilename = filename.substr(0, filename.length() - 4);
		outfilename += ("_" + toString(i) + ".exr");
		fprintf(stderr, "Generating: %s, alpha = %f\n", outfilename.c_str(), alpha[i]); 

		SynRasterToCamera = Translate(Vector(alpha[i] *filmSynRatio * outres[0] * 0.5f, - (alpha[i] *filmSynRatio * outres[1] * 0.5f), (1-alpha[i]) * filmdistance)) // NDC-To-Camera
		* Scale( alpha[i] * filmSynRatio, alpha[i] * filmSynRatio , 1.f) // Raster-To-NDC
		* Scale(-1.f, 1.f, 1.f);

		// start reconstruction

		for (int synT = 0; synT < outres[1]; synT++) {
			for (int synS = 0; synS < outres[0]; synS++) {
				// R, G, B, A: initialize
				int indexBase = (synT * outres[0] + synS)*4;
				synthesizeddata[indexBase  ] = 0.f;
				synthesizeddata[indexBase+1] = 0.f;
				synthesizeddata[indexBase+2] = 0.f;
				synthesizeddata[indexBase+3] = 1.f;
				for (int lenV = 0; lenV < microlensdiameter; lenV++) {
					for (int lenU = 0; lenU < microlensdiameter; lenU++) {
						Point synthRas((float)synS+0.5, (float)synT+0.5, 0.f);
						Point synthCamera;
						SynRasterToCamera(synthRas, &synthCamera);
						Point lenCamera(lenPointMin + lenU * lenPointDist, lenPointMin + lenV * lenPointDist, filmdistance);

						Point Pcamera(lenCamera.x + (synthCamera.x - lenCamera.x)/alpha[i], lenCamera.y + (synthCamera.y - lenCamera.y)/alpha[i], 0.f);
						Point Pras;
						CameraToRaster(Pcamera, &Pras);
						
						if (Pras.x >= 0.0 && Pras.x <= (float)res[0] && Pras.y >= 0.0 && Pras.y <= (float)res[1]) {
							float mlens = Pras.x / microlensdiameter; // e.g. 0-63
							float mlent = Pras.y / microlensdiameter; // e.g. 0-63
							float interS = mlens - floor(mlens);
							float interT = mlent - floor(mlent);

							int lightfieldIndexBase = ((Floor2Int(mlent) * microlensdiameter + lenV) * res[0] + (Floor2Int(mlens) * microlensdiameter + lenU)) * 4;
							int lightfieldIndexBase_1_0, lightfieldIndexBase_0_1, lightfieldIndexBase_1_1;
							if (Floor2Int(mlens) >= res[0] / microlensdiameter - 1) {lightfieldIndexBase_1_0 = lightfieldIndexBase;}
							else {lightfieldIndexBase_1_0 = ((Floor2Int(mlent) * microlensdiameter + lenV) * res[0] + ((Floor2Int(mlens) + 1) * microlensdiameter + lenU)) * 4;}
							if (Floor2Int(mlent) >= res[1] / microlensdiameter - 1) {lightfieldIndexBase_0_1 = lightfieldIndexBase;}
							else {lightfieldIndexBase_0_1 = (((Floor2Int(mlent) + 1) * microlensdiameter + lenV) * res[0] + (Floor2Int(mlens) * microlensdiameter + lenU)) * 4;}
							if (Floor2Int(mlens) >= res[0] / microlensdiameter - 1) {
								if (Floor2Int(mlent) >= res[1] / microlensdiameter - 1) {lightfieldIndexBase_1_1 = lightfieldIndexBase;}
								else {lightfieldIndexBase_1_1 = lightfieldIndexBase_0_1;}
							} else {
								if (Floor2Int(mlent) >= res[1] / microlensdiameter - 1) {lightfieldIndexBase_1_1 = lightfieldIndexBase_1_0;}
								else {lightfieldIndexBase_1_1 = (((Floor2Int(mlent) + 1) * microlensdiameter + lenV) * res[0] + ((Floor2Int(mlens) + 1) * microlensdiameter + lenU)) * 4;}
							}

							//fprintf(dump, "(%f, %f), (%d, %d), (%f, %f), (%f, %f), (%f, %f, %f), (%d, %d)\n", 
							//lenCamera.x, lenCamera.y, synS, synT, synthCamera.x, synthCamera.y, Pras.x, Pras.y
							//, lightfielddata[lightfieldIndexBase  ], lightfielddata[lightfieldIndexBase+1], lightfielddata[lightfieldIndexBase+2],
							//(Floor2Int(mlens) * microlensdiameter + lenU), (Floor2Int(mlent) * microlensdiameter + lenV));

							// interpolation


							synthesizeddata[indexBase  ] += ( (1.0-interS)*(1.0-interT)*lightfielddata[lightfieldIndexBase  ] + (interS)*(1.0-interT)*lightfielddata[lightfieldIndexBase_1_0  ]
								+ (1.0-interS)*(interT)*lightfielddata[lightfieldIndexBase_0_1  ] + (interS)*(interT)*lightfielddata[lightfieldIndexBase_1_1  ]) / (microlensdiameter * microlensdiameter);
							synthesizeddata[indexBase+1] += ( (1.0-interS)*(1.0-interT)*lightfielddata[lightfieldIndexBase+1] + (interS)*(1.0-interT)*lightfielddata[lightfieldIndexBase_1_0+1]
								+ (1.0-interS)*(interT)*lightfielddata[lightfieldIndexBase_0_1+1] + (interS)*(interT)*lightfielddata[lightfieldIndexBase_1_1+1]) / (microlensdiameter * microlensdiameter);
							synthesizeddata[indexBase+2] += ( (1.0-interS)*(1.0-interT)*lightfielddata[lightfieldIndexBase+2] + (interS)*(1.0-interT)*lightfielddata[lightfieldIndexBase_1_0+2]
								+ (1.0-interS)*(interT)*lightfielddata[lightfieldIndexBase_0_1+2] + (interS)*(interT)*lightfielddata[lightfieldIndexBase_1_1+2]) / (microlensdiameter * microlensdiameter);
						}
					}
				}
			}
		}
		// writefile
		WriteEXR(outfilename.c_str(), synthesizeddata, outres[0], outres[1]);
	}
	/*
	for (int i = 0; i < (int)this->alpha.size(); i++) {
		// decide outfilename
		string outfilename = filename.substr(0, filename.length() - 4);
		outfilename += ("_" + toString(i) + ".exr");
		fprintf(stderr, "%s\n", outfilename.c_str()); 

		SynRasterToCamera = Translate(Vector(alpha[i] *filmSynRatio * outres[0] * 0.5f, - (alpha[i] *filmSynRatio * outres[1] * 0.5f), (1-alpha[i]) * filmdistance)) // NDC-To-Camera
		* Scale( alpha[i] * filmSynRatio, alpha[i] * filmSynRatio , 1.f) // Raster-To-NDC
		* Scale(-1.f, 1.f, 1.f);

		int divider = (res[1] / microlensdiameter) * (res[0] / microlensdiameter) / microlensdiameter / microlensdiameter;

		for (int synT = 0; synT < outres[1]; synT++) {
			for (int synS = 0; synS < outres[0]; synS++) {
				// R, G, B, A: initialize
				int indexBase = (synT * outres[0] + synS)*4;
				synthesizeddata[indexBase  ] = 0.f;
				synthesizeddata[indexBase+1] = 0.f;
				synthesizeddata[indexBase+2] = 0.f;
				synthesizeddata[indexBase+3] = 1.f;

				for (int microT = 0; microT < res[1] / microlensdiameter; microT++) {
					for (int microS = 0; microS < res[0] / microlensdiameter; microS++) {
						Point synthRas((float)synS+0.5, (float)synT+0.5, 0.f);
						Point synthCamera;
						SynRasterToCamera(synthRas, &synthCamera);
						Point microRas( ((float)microS+0.5)*microlensdiameter, ((float)microT+0.5)*microlensdiameter, 0.f);
						Point microCamera;
						RasterToCamera(microRas, &microCamera);
						Point Plens(microCamera.x + (synthCamera.x - microCamera.x)/(1- alpha[i]), microCamera.y + (synthCamera.y - microCamera.y)/(1- alpha[i]), filmdistance);

						Point Pcamera(microCamera.x + Plens.x * microlensdiameter * filmResRatio / lenRadius * 2, microCamera.y + Plens.y * microlensdiameter * filmResRatio / lenRadius * 2, 0.f);
						Point Pras;
						CameraToRaster(Pcamera, &Pras);

						if (Pras.x >= (float)microS * microlensdiameter && Pras.x < (float)(microS+1) * microlensdiameter &&
							Pras.y >= (float)microT * microlensdiameter && Pras.y < (float)(microT+1) * microlensdiameter) {
							int lightfieldIndexBase = (Floor2Int(Pras.y) * res[0] + Floor2Int(Pras.x)) * 4;
							synthesizeddata[indexBase  ] += lightfielddata[lightfieldIndexBase  ] / (divider * (1-alpha[i]));
							synthesizeddata[indexBase+1] += lightfielddata[lightfieldIndexBase+1] / (divider * (1-alpha[i]));
							synthesizeddata[indexBase+2] += lightfielddata[lightfieldIndexBase+2] / (divider * (1-alpha[i]));
						}
					}
				}
			}
		}
		WriteEXR(outfilename.c_str(), synthesizeddata, outres[0], outres[1]);
	}*/

	delete [] synthesizeddata;
	return;

}

bool LightfieldCamera::ReadEXR(const char *name, float *&rgba, int &xRes, int &yRes, bool &hasAlpha)
{
    try {
    InputFile file(name);
    Box2i dw = file.header().dataWindow();
    xRes = dw.max.x - dw.min.x + 1;
    yRes = dw.max.y - dw.min.y + 1;

    half *hrgba = new half[4 * xRes * yRes];

    // for now...
    hasAlpha = true;
    int nChannels = 4;

    hrgba -= 4 * (dw.min.x + dw.min.y * xRes);
    FrameBuffer frameBuffer;
    frameBuffer.insert("R", Slice(HALF, (char *)hrgba,
				  4*sizeof(half), xRes * 4 * sizeof(half), 1, 1, 0.0));
    frameBuffer.insert("G", Slice(HALF, (char *)hrgba+sizeof(half),
				  4*sizeof(half), xRes * 4 * sizeof(half), 1, 1, 0.0));
    frameBuffer.insert("B", Slice(HALF, (char *)hrgba+2*sizeof(half),
				  4*sizeof(half), xRes * 4 * sizeof(half), 1, 1, 0.0));
    frameBuffer.insert("A", Slice(HALF, (char *)hrgba+3*sizeof(half),
				  4*sizeof(half), xRes * 4 * sizeof(half), 1, 1, 1.0));

    file.setFrameBuffer(frameBuffer);
    file.readPixels(dw.min.y, dw.max.y);

    hrgba += 4 * (dw.min.x + dw.min.y * xRes);
    rgba = new float[nChannels * xRes * yRes];
    for (int i = 0; i < nChannels * xRes * yRes; ++i)
	rgba[i] = hrgba[i];
    delete[] hrgba;
    } catch (const std::exception &e) {
        fprintf(stderr, "Unable to read image file \"%s\": %s", name, e.what());
        return NULL;
    }

    return rgba;
}

void LightfieldCamera::WriteEXR(const char *name, float *pixels, int xRes, int yRes) {
    Rgba *hrgba = new Rgba[xRes * yRes];
    for (int i = 0; i < xRes * yRes; ++i)
        hrgba[i] = Rgba(pixels[4*i], pixels[4*i+1], pixels[4*i+2], 1.);

    Box2i displayWindow(V2i(0,0), V2i(xRes-1, yRes-1));
    Box2i dataWindow = displayWindow;

    RgbaOutputFile file(name, displayWindow, dataWindow, WRITE_RGBA);
    file.setFrameBuffer(hrgba, 1, xRes);
    try {
        file.writePixels(yRes);
    }
    catch (const std::exception &e) {
        fprintf(stderr, "Unable to write image file \"%s\": %s", name,
                e.what());
    }

    delete[] hrgba;
}

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}


LightfieldCamera *CreateLightfieldCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
	// Extract common camera parameters from \use{ParamSet}
	float hither = params.FindOneFloat("hither", -1);
	float yon = params.FindOneFloat("yon", -1);
	float shutteropen = params.FindOneFloat("shutteropen", -1);
	float shutterclose = params.FindOneFloat("shutterclose", -1);

	// Lightfield camera-specific parameters
	string specfile = params.FindOneFilename("specfile", "");  // original FindOneString 
	float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
 	float fstop = params.FindOneFloat("aperture_diameter", 1.0);	
	float filmdiag = params.FindOneFloat("filmdiag", 35.0);
	int microlensdiameter = (int)params.FindOneFloat("microlens_diameter", 10.0);
	string rawalpha = params.FindOneString("alpha", "");
	std::vector<std::string> alphastring = split(rawalpha, ' ');
	std::vector<float> alpha;
	alpha.clear();
	for (int i = 0; i < (int) alphastring.size(); i++) {
		alpha.push_back((float)atof(alphastring[i].c_str()));
	}

	Assert(hither != -1 && yon != -1 && shutteropen != -1 &&
		shutterclose != -1 && filmdistance!= -1);
	if (specfile == "") {
	    Severe( "No lens spec file supplied!\n" );
	}
	return new LightfieldCamera(cam2world, hither, yon,
				   shutteropen, shutterclose, filmdistance, fstop, 
				   specfile, filmdiag, microlensdiameter, alpha, film);
}
