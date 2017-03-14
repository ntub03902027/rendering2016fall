// lights/medianCutEnvironmentLight.cpp*
#include "stdafx.h"
#include "lights/medianCutEnvironmentLight.h"
#include "sh.h"
#include "montecarlo.h"
#include "paramset.h"
#include "imageio.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

//#define DEBUG_PROJ3

// WARNING
// THE FOLLOWING TWO CANNOT BE BOTH DEFINED
#define SQUARE_WEIGHT_PROJ3
//#define CENTER_PROJ3
// END OF WARNING

// vdb debug
//#include "vdb.h"

// unmodified
// MedianCutEnvironmentLight Utility Classes
struct InfiniteAreaCube {
    // InfiniteAreaCube Public Methods
    InfiniteAreaCube(const MedianCutEnvironmentLight *l, const Scene *s,
                     float t, bool cv, float pe)
        : light(l), scene(s), time(t), pEpsilon(pe), computeVis(cv) { }
    Spectrum operator()(int, int, const Point &p, const Vector &w) {
        Ray ray(p, w, pEpsilon, INFINITY, time);
        if (!computeVis || !scene->IntersectP(ray))
            return light->Le(RayDifferential(ray));
        return 0.f;
    }
    const MedianCutEnvironmentLight *light;
    const Scene *scene;
    float time, pEpsilon;
    bool computeVis;
};


// unmodified
// MedianCutEnvironmentLight Method Definitions
MedianCutEnvironmentLight::~MedianCutEnvironmentLight() {
    delete distribution;
    delete radianceMap;
}

// modified in hw3
MedianCutEnvironmentLight::MedianCutEnvironmentLight(const Transform &light2world,
        const Spectrum &L, int ns, const string &texmap)
    : Light(light2world, ns) {
    int width = 0, height = 0;
    RGBSpectrum *texels = NULL;
    // Read texel data from _texmap_ into _texels_
    if (texmap != "") {
        texels = ReadImage(texmap, &width, &height);
        if (texels)
            for (int i = 0; i < width * height; ++i)
                texels[i] *= L.ToRGBSpectrum();
    }
    if (!texels) {
        width = height = 1;
        texels = new RGBSpectrum[1];
        texels[0] = L.ToRGBSpectrum();
    }

// =========================================== MODIFICATION ==================================================


    std::vector<std::vector<float> > sumTable;
// initial size of sumTable
    sumTable.resize(height);
    for (int i = 0; i < height; i++)
        sumTable[i].resize(width);
    // sumTable[0][0]
    sumTable[0][0] = texels[0].y(); // y = 0.212671f * R + 0.715160f * G + 0.072169f * B
    // fill sumTable
    for (int i = 1; i < width; i++)
        sumTable[0][i] = sumTable[0][i - 1] + texels[i].y();
    for (int i = 1; i < height; i++)
        sumTable[i][0] = sumTable[i - 1][0] + texels[width * i].y();
    for (int i = 1; i < height; i++)
        for (int j = 1; j < width; j++)
            sumTable[i][j] = sumTable[i - 1][j] + sumTable[i][j - 1] - sumTable[i - 1][j - 1] + texels[width * i + j].y();
/*
    FILE *dump;
    dump = fopen("dump.out", "w");
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            fprintf(dump, "sumTable[%4d][%4d] = %6.6f ", i, j, sumTable[i][j]);
        }
        fprintf(dump, "\n");
    }
    fclose(dump);
*/
    // fill pointLightTable

    // vdb debug
    //vdb_color(1.f, 1.f, 1.f);
    //vdb_point(0.0, 0.0, 0.0);
    //vdb_line(0.0, 0.0, 0.0, width, 0.0, 0.0);
    //vdb_line(0.0, 0.0, 0.0, 0.0, height, 0.0);
    //vdb_line(0.0, height, 0.0, width, height, 0.0);
    //vdb_line(width, 0.0, 0.0, width, height, 0.0);

    pointLightTable.clear();
    FillPointLightTable(Ceil2Int(log2f(ns)), 0, height - 1, 0, width - 1, sumTable[height - 1][width - 1], sumTable, texels, height, width);

#ifdef DEBUG_PROJ3
    fprintf(stderr, "%ld\n", pointLightTable.size());
    for (int i = 0; i < (int) pointLightTable.size(); i++ ) {
        fprintf(stderr, "%f, %f, %f\n", pointLightTable[i].lightPos.x, pointLightTable[i].lightPos.y, pointLightTable[i].lightPos.z);
    }
#endif

    intensityDivisor = (float) height * (float) width / 40.0;
    uniformPDF = 1.f / (float) pointLightTable.size();

    sumTable.clear();


// =========================================== END OF MODIFICATION ==================================================
// BELOW: UNMODIFIED FOR Le() to use

    // process Mipmap: for Le() to use
    radianceMap = new MIPMap<RGBSpectrum>(width, height, texels);
    delete[] texels;
    // Initialize sampling PDFs for infinite area light

    // Compute scalar-valued image _img_ from environment map
    float filter = 1.f / max(width, height);
    float *img = new float[width*height];
    for (int v = 0; v < height; ++v) {
        float vp = (float)v / (float)height;
        float sinTheta = sinf(M_PI * float(v+.5f)/float(height));
        for (int u = 0; u < width; ++u) {
            float up = (float)u / (float)width;
            img[u+v*width] = radianceMap->Lookup(up, vp, filter).y();
            img[u+v*width] *= sinTheta;
        }
    }

    // Compute sampling distributions for rows and columns of image
    distribution = new Distribution2D(img, width, height);
    delete[] img;
}

// Get sample point lights: core mediancut algorithm
void MedianCutEnvironmentLight::FillPointLightTable( const int depth, const int hmin, const int hmax, const int wmin, const int wmax, const float sum , const std::vector<std::vector<float> > sumTable, RGBSpectrum* texels, const int height, const int width) {
    // add point light
    if (depth == 0) {
#ifdef SQUARE_WEIGHT_PROJ3
        float gamma_square = 0.f;
        float gamma_square_width = 0.f;
        float gamma_square_height = 0.f;
        RGBSpectrum totalEnergy = RGBSpectrum(0.f);
        for (int i = hmin; i <= hmax; i++) {
            for (int j = wmin; j <= wmax; j++) {
                gamma_square += texels[width * i + j].y() * texels[width * i + j].y();
                gamma_square_width += texels[width * i + j].y() * texels[width * i + j].y() * (float) j;
                gamma_square_height += texels[width * i + j].y() * texels[width * i + j].y() * (float) i;
                totalEnergy += texels[width * i + j];
            }
        }
        float h = gamma_square_height / gamma_square;
        float w = gamma_square_width / gamma_square;
#endif

#ifdef CENTER_PROJ3
        float h = 0.5 * (float)(hmin + hmax);
        float w = 0.5 * (float)(wmin + wmax);
        RGBSpectrum totalEnergy = RGBSpectrum(0.f);
        for (int i = hmin; i <= hmax; i++) 
            for (int j = wmin; j <= wmax; j++) 
                totalEnergy += texels[width * i + j];
#endif

#if !defined(CENTER_PROJ3) && !defined(SQUARE_WEIGHT_PROJ3)
        float gamma = 0.f;
        float gamma_width = 0.f;
        float gamma_height = 0.f;
        RGBSpectrum totalEnergy = RGBSpectrum(0.f);
        for (int i = hmin; i <= hmax; i++) {
            for (int j = wmin; j <= wmax; j++) {
                gamma += texels[width * i + j].y();
                gamma_width += texels[width * i + j].y() * (float) j;
                gamma_height += texels[width * i + j].y() * (float) i;
                totalEnergy += texels[width * i + j];
            }
        }
        float h = gamma_height / gamma;
        float w = gamma_width / gamma;
#endif
        // vdb debug
        //vdb_color(1.0, 1.0, 0.0);
        //vdb_point(w, h, 0.0);
#ifdef DEBUG_PROJ3
        fprintf(stderr, "final: %f\n", totalEnergy.y());
#endif
        pointLight newLight;
        float theta = h / height * M_PI;
        float phi = w / width * 2.f * M_PI;
        float costheta = cosf(theta), sintheta = sinf(theta);
        float cosphi = cosf(phi), sinphi = sinf(phi);
        newLight.lightPos = Point(sintheta * cosphi, sintheta * sinphi, costheta);
        newLight.intensity = totalEnergy;
        pointLightTable.push_back(newLight);

    }
    // cut in median, recursively
    else {
        // cut vertically
        if (wmax - wmin >= hmax - hmin) {
            //cannot cut anymore: make (depth + 1) same light samples
            if (wmax == wmin) {
                float h = hmin;
                float w = wmin;
                RGBSpectrum totalEnergy = texels[width * hmin + wmin];
#ifdef DEBUG_PROJ3
                fprintf(stderr, "final (single point): %f, number of dulpicated lights: %d\n", totalEnergy.y(), depth + 1);
#endif          
                pointLight newLight;
                float theta = h / height * M_PI;
                float phi = w / width * 2.f * M_PI;
                float costheta = cosf(theta), sintheta = sinf(theta);
                float cosphi = cosf(phi), sinphi = sinf(phi);
                newLight.lightPos = Point(sintheta * cosphi, sintheta * sinphi, costheta);
                newLight.intensity = totalEnergy;
                for (int i = 0; i < depth + 1; i++)
                    pointLightTable.push_back(newLight);
            }
            else {
                int cutPoint = 0;
                for (int i = wmin; i <= wmax; i++) {
                    float partSum = sumTable[hmax][i] - ( (wmin > 0)? sumTable[hmax][wmin - 1] : 0.f )
                     - ( (hmin > 0)? sumTable[hmin - 1][i] : 0.f ) + ( (wmin > 0 && hmin > 0)? sumTable[hmin - 1][wmin - 1] : 0.f );
                    if (partSum >= 0.5 * sum) {
                        if (i == wmin)
                            cutPoint = wmin;
                        else if (i == wmax)
                            cutPoint = wmax - 1;
                        // compare which is nearer to the median
                        else {
                            float previousPartSum = sumTable[hmax][i - 1] - ( (wmin > 0)? sumTable[hmax][wmin - 1] : 0.f )
                            - ( (hmin > 0)? sumTable[hmin - 1][i - 1] : 0.f ) + ( (wmin > 0 && hmin > 0)? sumTable[hmin - 1][wmin - 1] : 0.f );
                            if (fabs(partSum - 0.5 * sum) <= fabs(previousPartSum - 0.5 * sum))
                                cutPoint = i;
                            else
                                cutPoint = i - 1;
                        }
                        break;
                    }
                }
                // vdb debug
                //vdb_color(0.0, 0.0, 1.0);
                //vdb_line(cutPoint + 1, hmin, 0.0, cutPoint + 1, hmax + 1, 0.0);
                // caluclate sums
                float newSum1 = sumTable[hmax][cutPoint] - ( (wmin > 0)? sumTable[hmax][wmin - 1] : 0.f )
                     - ( (hmin > 0)? sumTable[hmin - 1][cutPoint] : 0.f ) + ( (wmin > 0 && hmin > 0)? sumTable[hmin - 1][wmin - 1] : 0.f );
                float newSum2 = sumTable[hmax][wmax] - sumTable[hmax][cutPoint]
                     - ( (hmin > 0)? sumTable[hmin - 1][wmax] : 0.f ) + ( (hmin > 0)? sumTable[hmin - 1][cutPoint] : 0.f );
#ifdef DEBUG_PROJ3
                fprintf(stderr, "vert: %d, %f, %f\n", depth, newSum1, newSum2);
#endif
                FillPointLightTable(depth - 1, hmin, hmax, wmin, cutPoint, newSum1, sumTable, texels, height, width);
                FillPointLightTable(depth - 1, hmin, hmax, cutPoint + 1, wmax, newSum2, sumTable, texels, height, width);
            }
        }
        //cut horizontally
        else {
            //cannot cut anymore
            if (hmax == hmin) {
                float h = hmin;
                float w = wmin;
                RGBSpectrum totalEnergy = texels[width * hmin + wmin];
#ifdef DEBUG_PROJ3
                fprintf(stderr, "final (single point): %f, number of dulpicated lights: %d\n", totalEnergy.y(), depth + 1);
#endif          
                pointLight newLight;
                float theta = h / height * M_PI;
                float phi = w / width * 2.f * M_PI;
                float costheta = cosf(theta), sintheta = sinf(theta);
                float cosphi = cosf(phi), sinphi = sinf(phi);
                newLight.lightPos = Point(sintheta * cosphi, sintheta * sinphi, costheta);
                newLight.intensity = totalEnergy;
                for (int i = 0; i < depth + 1; i++)
                    pointLightTable.push_back(newLight);
            }
            else {
                int cutPoint = 0;
                for (int i = hmin; i <= hmax; i++) {
                    float partSum = sumTable[i][wmax] - ( (hmin > 0)? sumTable[hmin - 1][wmax] : 0.f )
                     - ( (wmin > 0)? sumTable[i][wmin - 1] : 0.f ) + ( (wmin > 0 && hmin > 0)? sumTable[hmin - 1][wmin - 1] : 0.f );
                    if (partSum >= 0.5 * sum) {
                        if (i == hmin)
                            cutPoint = hmin;
                        else if (i == hmax)
                            cutPoint = hmax - 1;
                        // compare which is nearer to the median
                        else {
                            float previousPartSum = sumTable[i - 1][wmax] - ( (hmin > 0)? sumTable[hmin - 1][wmax] : 0.f )
                                - ( (wmin > 0)? sumTable[i - 1][wmin - 1] : 0.f ) + ( (wmin > 0 && hmin > 0)? sumTable[hmin - 1][wmin - 1] : 0.f );
                            if (fabs(partSum - 0.5 * sum) <= fabs(previousPartSum - 0.5 * sum))
                                cutPoint = i;
                            else
                                cutPoint = i - 1;
                        }
                        break;
                    }
                }
                // vdb debug
                //vdb_color(0.0, 0.0, 1.0);
                //vdb_line(wmin, cutPoint + 1, 0.0, wmax + 1, cutPoint + 1, 0.0);
                // caluclate sums
                float newSum1 = sumTable[cutPoint][wmax] - ( (hmin > 0)? sumTable[hmin - 1][wmax] : 0.f )
                     - ( (wmin > 0)? sumTable[cutPoint][wmin - 1] : 0.f ) + ( (wmin > 0 && hmin > 0)? sumTable[hmin - 1][wmin - 1] : 0.f );
                float newSum2 = sumTable[hmax][wmax] - sumTable[cutPoint][wmax]
                     - ( (wmin > 0)? sumTable[hmax][wmin - 1] : 0.f ) + ( (wmin > 0)? sumTable[cutPoint][wmin - 1] : 0.f );
#ifdef DEBUG_PROJ3
                fprintf(stderr, "hori: %d, %f, %f\n", depth, newSum1, newSum2);
#endif
                FillPointLightTable(depth - 1, hmin, cutPoint, wmin, wmax, newSum1, sumTable, texels, height, width);
                FillPointLightTable(depth - 1, cutPoint + 1, hmax, wmin, wmax, newSum2, sumTable, texels, height, width);
            }
        }
    }
    return;
}

// unmodified
Spectrum MedianCutEnvironmentLight::Power(const Scene *scene) const {
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    return M_PI * worldRadius * worldRadius *
        Spectrum(radianceMap->Lookup(.5f, .5f, .5f), SPECTRUM_ILLUMINANT);
}

// unmodified
Spectrum MedianCutEnvironmentLight::Le(const RayDifferential &r) const {
    Vector wh = Normalize(WorldToLight(r.d));
    float s = SphericalPhi(wh) * INV_TWOPI;
    float t = SphericalTheta(wh) * INV_PI;
    return Spectrum(radianceMap->Lookup(s, t), SPECTRUM_ILLUMINANT);
}

// unmodified
void MedianCutEnvironmentLight::SHProject(const Point &p, float pEpsilon,
        int lmax, const Scene *scene, bool computeLightVis,
        float time, RNG &rng, Spectrum *coeffs) const {
    // Project MedianCutEnvironmentLight to SH using Monte Carlo if visibility needed
    if (computeLightVis) {
        Light::SHProject(p, pEpsilon, lmax, scene, computeLightVis,
                         time, rng, coeffs);
        return;
    }
    for (int i = 0; i < SHTerms(lmax); ++i)
        coeffs[i] = 0.f;
    int ntheta = radianceMap->Height(), nphi = radianceMap->Width();
    if (min(ntheta, nphi) > 50) {
        // Project MedianCutEnvironmentLight to SH from lat-long representation

        // Precompute $\theta$ and $\phi$ values for lat-long map projection
        float *buf = new float[2*ntheta + 2*nphi];
        float *bufp = buf;
        float *sintheta = bufp;  bufp += ntheta;
        float *costheta = bufp;  bufp += ntheta;
        float *sinphi = bufp;    bufp += nphi;
        float *cosphi = bufp;
        for (int theta = 0; theta < ntheta; ++theta) {
            sintheta[theta] = sinf((theta + .5f)/ntheta * M_PI);
            costheta[theta] = cosf((theta + .5f)/ntheta * M_PI);
        }
        for (int phi = 0; phi < nphi; ++phi) {
            sinphi[phi] = sinf((phi + .5f)/nphi * 2.f * M_PI);
            cosphi[phi] = cosf((phi + .5f)/nphi * 2.f * M_PI);
        }
        float *Ylm = ALLOCA(float, SHTerms(lmax));
        for (int theta = 0; theta < ntheta; ++theta) {
            for (int phi = 0; phi < nphi; ++phi) {
                // Add MedianCutEnvironmentLight texel's contribution to SH coefficients
                Vector w = Vector(sintheta[theta] * cosphi[phi],
                                  sintheta[theta] * sinphi[phi],
                                  costheta[theta]);
                w = Normalize(LightToWorld(w));
                Spectrum Le = Spectrum(radianceMap->Texel(0, phi, theta),
                                       SPECTRUM_ILLUMINANT);
                SHEvaluate(w, lmax, Ylm);
                for (int i = 0; i < SHTerms(lmax); ++i)
                    coeffs[i] += Le * Ylm[i] * sintheta[theta] *
                        (M_PI / ntheta) * (2.f * M_PI / nphi);
            }
        }

        // Free memory used for lat-long theta and phi values
        delete[] buf;
    }
    else {
        // Project MedianCutEnvironmentLight to SH from cube map sampling
        SHProjectCube(InfiniteAreaCube(this, scene, time, computeLightVis,
                                       pEpsilon),
                      p, 200, lmax, coeffs);
    }
}

// unmodified
MedianCutEnvironmentLight *CreateMedianCutEnvironmentLight(const Transform &light2world,
        const ParamSet &paramSet) {
    Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    string texmap = paramSet.FindOneFilename("mapname", "");
    int nSamples = paramSet.FindOneInt("nsamples", 1);
    if (PbrtOptions.quickRender) nSamples = max(1, nSamples / 4);
    return new MedianCutEnvironmentLight(light2world, L * sc, nSamples, texmap);
}

// modified, will be called by pbrt
Spectrum MedianCutEnvironmentLight::Sample_L(const Point &p, float pEpsilon,
        const LightSample &ls, float time, Vector *wi, float *pdf,
        VisibilityTester *visibility) const {



    PBRT_INFINITE_LIGHT_STARTED_SAMPLE();
    // Fill wi
    int lightNo = Floor2Int( ls.uComponent * pointLightTable.size() );
    *wi = LightToWorld(Vector(pointLightTable[lightNo].lightPos.x, pointLightTable[lightNo].lightPos.y, pointLightTable[lightNo].lightPos.z));

    // Fill PDF
    *pdf = this->uniformPDF;

    visibility->SetRay(p, pEpsilon, *wi, time);

    Spectrum Ls = Spectrum(pointLightTable[lightNo].intensity / intensityDivisor, SPECTRUM_ILLUMINANT);

    PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();

    return Ls;
}

// unmodified
float MedianCutEnvironmentLight::Pdf(const Point &, const Vector &w) const {
    PBRT_INFINITE_LIGHT_STARTED_PDF();
    Vector wi = WorldToLight(w);
    float theta = SphericalTheta(wi), phi = SphericalPhi(wi);
    float sintheta = sinf(theta);
    if (sintheta == 0.f) return 0.f;
    float p = distribution->Pdf(phi * INV_TWOPI, theta * INV_PI) /
           (2.f * M_PI * M_PI * sintheta);
    PBRT_INFINITE_LIGHT_FINISHED_PDF();
    return p;
}

// unmodified, will not be used despite the function name
Spectrum MedianCutEnvironmentLight::Sample_L(const Scene *scene,
        const LightSample &ls, float u1, float u2, float time,
        Ray *ray, Normal *Ns, float *pdf) const {
    PBRT_INFINITE_LIGHT_STARTED_SAMPLE();
    // Compute direction for infinite light sample ray

    // Find $(u,v)$ sample coordinates in infinite light texture
    float uv[2], mapPdf;
    distribution->SampleContinuous(ls.uPos[0], ls.uPos[1], uv, &mapPdf);
    if (mapPdf == 0.f) return Spectrum(0.f);

    float theta = uv[1] * M_PI, phi = uv[0] * 2.f * M_PI;
    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    Vector d = -LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                                    costheta));
    *Ns = (Normal)d;

    // Compute origin for infinite light sample ray
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    Vector v1, v2;
    CoordinateSystem(-d, &v1, &v2);
    float d1, d2;
    ConcentricSampleDisk(u1, u2, &d1, &d2);
    Point Pdisk = worldCenter + worldRadius * (d1 * v1 + d2 * v2);
    *ray = Ray(Pdisk + worldRadius * -d, d, 0., INFINITY, time);

    // Compute MedianCutEnvironmentLight ray PDF
    float directionPdf = mapPdf / (2.f * M_PI * M_PI * sintheta);
    float areaPdf = 1.f / (M_PI * worldRadius * worldRadius);
    *pdf = directionPdf * areaPdf;
    if (sintheta == 0.f) *pdf = 0.f;
    Spectrum Ls = (radianceMap->Lookup(uv[0], uv[1]), SPECTRUM_ILLUMINANT);
    PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();
    return Ls;
}


