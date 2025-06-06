
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_SHAPES_SPHERE_H
#define PBRT_SHAPES_SPHERE_H

// shapes/sphere.h*
#include "shape.h"

#include "pbrt.h"
#include "transform.h"
#include "texture.h"

namespace pbrt {

class HeightFieldGeneration : public Shape {
public:
    // HeightField Public Methods
    HeightFieldGeneration(const Transform *ObjectToWorld, const Transform *WorldToObject,
                         bool reverseOrientation, Float width, Float height, int nOctaves, Float freq, Float amplitude)
        : Shape(ObjectToWorld, WorldToObject, reverseOrientation),
          width(width), height(height), nOctaves(nOctaves), freq(freq), amplitude(amplitude) {}

    // Declared Functions
    Bounds3f ObjectBound() const;
    bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
                   bool testAlphaTexture) const;
    bool IntersectP(const Ray &ray, bool testAlphaTexture) const;
    Float Area() const;
    Interaction Sample(const Point2f &u, Float *pdf) const;
    Float Pdf(const Interaction &ref, const Vector3f &wi) const;

    // Noise
    float OctaveNoise(Point3f p, int nOctaves, Float freq, Float amplitude) const;

private:
    // Private data members
    const Float width, height;
    const int nOctaves;
    const Float freq, amplitude;
};

// Function declaration for creating HeightFieldGeneration shape
std::shared_ptr<Shape> CreateHeightFieldGenerationShape(const Transform *o2w,
                                     const Transform *w2o,
                                     bool reverseOrientation,
                                     const ParamSet &params);


 float OctaveNoise(const Point3f &p, int nOctaves, float freq, float amplitude);

/*
 
 
 ----------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------
 
 

 
 ----------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------
 
 
 */



class Sphere : public Shape {
  public:
    // Sphere Public Methods
    Sphere(const Transform *ObjectToWorld, const Transform *WorldToObject,
           bool reverseOrientation, Float radius, Float zMin, Float zMax,
           Float phiMax)
        : Shape(ObjectToWorld, WorldToObject, reverseOrientation),
          radius(radius),
          zMin(Clamp(std::min(zMin, zMax), -radius, radius)),
          zMax(Clamp(std::max(zMin, zMax), -radius, radius)),
          thetaMin(std::acos(Clamp(std::min(zMin, zMax) / radius, -1, 1))),
          thetaMax(std::acos(Clamp(std::max(zMin, zMax) / radius, -1, 1))),
          phiMax(Radians(Clamp(phiMax, 0, 360))) {}
    Bounds3f ObjectBound() const;
    bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
                   bool testAlphaTexture) const;
    bool IntersectP(const Ray &ray, bool testAlphaTexture) const;
    Float Area() const;
    Interaction Sample(const Point2f &u, Float *pdf) const;
    Interaction Sample(const Interaction &ref, const Point2f &u,
                       Float *pdf) const;
    Float Pdf(const Interaction &ref, const Vector3f &wi) const;
    Float SolidAngle(const Point3f &p, int nSamples) const;

  private:
    // Sphere Private Data
    const Float radius;
    const Float zMin, zMax;
    const Float thetaMin, thetaMax, phiMax;
};



std::shared_ptr<Shape> CreateSphereShape(const Transform *o2w,
                                         const Transform *w2o,
                                         bool reverseOrientation,
                                         const ParamSet &params);


}  // namespace pbrt



#endif  // PBRT_SHAPES_SPHERE_H
