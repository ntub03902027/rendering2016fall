#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_LIGHTS_MEDIAN_CUT_ENVIRONMENT_H
#define PBRT_LIGHTS_MEDIAN_CUT_ENVIRONMENT_H

// lights/infinite.h*
#include "pbrt.h"
#include "light.h"
#include "texture.h"
#include "shape.h"
#include "scene.h"
#include "mipmap.h"

#include <vector>

class pointLight {
public:
    Point lightPos;
    Spectrum intensity;
};

// MedianCutEnvironmentLight Declarations
class MedianCutEnvironmentLight : public Light {
public:
    // MedianCutEnvironmentLight Public Methods
    MedianCutEnvironmentLight(const Transform &light2world, const Spectrum &power, int ns,
        const string &texmap);
    ~MedianCutEnvironmentLight();
    Spectrum Power(const Scene *) const;
    bool IsDeltaLight() const { return true; }
    Spectrum Le(const RayDifferential &r) const;
    Spectrum Sample_L(const Point &p, float pEpsilon, const LightSample &ls,
        float time, Vector *wi, float *pdf, VisibilityTester *visibility) const;
    Spectrum Sample_L(const Scene *scene, const LightSample &ls, float u1, float u2,
        float time, Ray *ray, Normal *Ns, float *pdf) const;
    float Pdf(const Point &, const Vector &) const;
    void SHProject(const Point &p, float pEpsilon, int lmax, const Scene *scene,
        bool computeLightVis, float time, RNG &rng, Spectrum *coeffs) const;
private:
    // MedianCutEnvironmentLight Private Data

    //std::vector<std::vector<float> > sumTable;

    std::vector<pointLight> pointLightTable;

    float intensityDivisor;
    float uniformPDF;

    void FillPointLightTable(const int depth, const int hmin, const int hmax, const int wmin, const int wmax, const float sum, const std::vector<std::vector<float> > sumTable, RGBSpectrum* texels, const int height, const int width);

    MIPMap<RGBSpectrum> *radianceMap;
    Distribution2D *distribution;
};



MedianCutEnvironmentLight *CreateMedianCutEnvironmentLight(const Transform &light2world,
        const ParamSet &paramSet);

#endif // PBRT_LIGHTS_MEDIAN_CUT_ENVIRONMENT_H
