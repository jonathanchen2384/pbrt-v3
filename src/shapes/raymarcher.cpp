
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

// Modification to Sphere Shape to Implement RayMarching
// by Kevin M. Smith 3-2-2019
//


// shapes/sphere.cpp*
#include "shapes/raymarcher.h"
#include "sampling.h"
#include "paramset.h"
#include "efloat.h"
#include "stats.h"

namespace pbrt {

// Sphere Method Definitions
Bounds3f RayMarcher::ObjectBound() const {
    return Bounds3f(Point3f(-radius, -radius, zMin),
                    Point3f(radius, radius, zMax));
}


//  Note: create PBRT parameters for these
//

#define MAX_RAY_STEPS 1000
#define DIST_THRESHOLD .01
#define MAX_DISTANCE 100
#define NORMAL_EPS .01


//  Template Method
//
bool RayMarcher::Intersect(const Ray &r, Float *tHit,
                                  SurfaceInteraction *isect,
                                  bool testAlphaTexture) const {
   
    
    // Part 4:
    const int counter;
    printf(counter);
    counter++;
    
    /*
     Q: 
        The numbers are likely not in order.
        Why and how to we remedy this for debugging purposes?
     
     A:
        The numbers are likely not in order
        due to the multi threading process causing a race condition
        among the three lines of code.
        
        One method of resolving this issue without hindering
        the process of the threads and conditions would be
        to utilize locks to ensure the variable is incrementing correctly
        without the interference. Through the use of locks,
        deadlocks must also be taken into consideration to prevent
        the threads from getting locked altogether.
     
        For debugging purposes:
        you may also set the amount of threads used by pbrt to
        the value "1" by doing "--nthreads" in the commandline,
        so only run thread exists and the increment will
        not be affected by race conditions.
     */
    
    // Part 6
    
	bool hit = false;
    Vector3f dir = Normalize(r.d);
    
    for ( int i = 0; i < MAX_RAY_STEPS; i ++)
    {
        float dist = sdf(dir);
        if (dist < DIST_THRESHOLD)
        {
            hit = true;
            break;
        }
        
        else if (dist > DIST_THRESHOLD)
        {
            break;
        }
        
        else {
            dir = dir + r.d*dist;
        }
    }
    

    ProfilePhase p(Prof::ShapeIntersect);

    Vector3f oErr, dErr;
    Ray ray = (*WorldToObject)(r, &oErr, &dErr);

    EFloat ox(ray.o.x, oErr.x), oy(ray.o.y, oErr.y), oz(ray.o.z, oErr.z);
    EFloat dx(ray.d.x, dErr.x), dy(ray.d.y, dErr.y), dz(ray.d.z, dErr.z);
    EFloat a = dx * dx + dy * dy + dz * dz;
    EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
    EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);
    EFloat t0, t1;
    
    if (!Quadratic(a, b, c, &t0, &t1)) {
        return false;
    }
    if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= 0) {
        return false;
    }
    
    EFloat tShapeHit = t0;
    
    if (tShapeHit.LowerBound() <= 0) {
        tShapeHit = t1;
        if (tShapeHit.UpperBound() > ray.tMax) {
            return false;
        }
    }

    Point3f pHit = ray((Float)tShapeHit);

    pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
    if (pHit.x == 0 && pHit.y == 0) 
    {
        pHit.x = 1e-5f * radius;
    }

    Float phi = std::atan2(pHit.y, pHit.x);
    if (phi < 0) {
        phi += 2 * Pi;
    }
    Float u = phi / phiMax;
    Float theta = std::acos(Clamp(pHit.z / radius, -1, 1));
    Float v = (theta - thetaMin) / (thetaMax - thetaMin);

    Float zRadius = std::sqrt(pHit.x * pHit.x + pHit.y * pHit.y);
    Float invZRadius = 1 / zRadius;
    Float cosPhi = pHit.x * invZRadius;
    Float sinPhi = pHit.y * invZRadius;
    Vector3f dpdu(-phiMax * pHit.y, phiMax * pHit.x, 0);
    Vector3f dpdv = (thetaMax - thetaMin) * Vector3f(pHit.z * cosPhi, pHit.z * sinPhi, -radius * std::sin(theta));

    Vector3f d2Pduu = -phiMax * phiMax * Vector3f(pHit.x, pHit.y, 0);
    Vector3f d2Pduv = (thetaMax - thetaMin) * pHit.z * phiMax * Vector3f(-sinPhi, cosPhi, 0.);
    Vector3f d2Pdvv = -(thetaMax - thetaMin) * (thetaMax - thetaMin) * Vector3f(pHit.x, pHit.y, pHit.z);
    Float E = Dot(dpdu, dpdu);
    Float F = Dot(dpdu, dpdv);
    Float G = Dot(dpdv, dpdv);
    Vector3f N = Normalize(Cross(dpdu, dpdv));
    Float e = Dot(N, d2Pduu);
    Float f = Dot(N, d2Pduv);
    Float g = Dot(N, d2Pdvv);
    Float invEGF2 = 1 / (E * G - F * F);
    Normal3f dndu = Normal3f((f * F - e * G) * invEGF2 * dpdu + (e * F - f * E) * invEGF2 * dpdv);
    Normal3f dndv = Normal3f((g * F - f * G) * invEGF2 * dpdu + (f * F - g * E) * invEGF2 * dpdv);

    Vector3f pError = gamma(5) * Abs((Vector3f)pHit);
    *isect = (*ObjectToWorld)(SurfaceInteraction(pHit, pError, Point2f(u, v), -ray.d, dpdu, dpdv, dndu, dndv, ray.time, this));
    *tHit = (Float)tShapeHit;

	if (hit && tHit != nullptr && isect != nullptr) {
		// This where you return your SurfaceInteraction structure and your tHit
		// Important Note: You must check for null pointer as Intersect is called 
		// by IntersectP() with null values for these parameters.
        
        for ( int i = 0; i < MAX_RAY_STEPS; i ++)
        {
            float dist = sdf(dir);
            if (dist < DIST_THRESHOLD)
            {
                hit = true;
                break;
            }
            
            else if (dist > DIST_THRESHOLD)
            {
                break;
            }
            
            else {
                dir = dir + r.d*dist;
            }
        }
	}
    return hit;
}

bool RayMarcher::IntersectP(const Ray &r, bool testAlphaTexture) const {
    Ray ray = (*WorldToObject)(r);

    Float tMin, tMax;
    if (!ObjectBound().IntersectP(ray, &tMin, &tMax))
        return false;

    Float t = tMin;
    Vector3f dir = Normalize(ray.d);

    for (int i = 0; i < MAX_RAY_STEPS; i++) {
        Float dist = sdf(ray(t));
        if (dist < DIST_THRESHOLD)
            return true;
        t += dist;
        if (t >= tMax || t >= MAX_DISTANCE)
            break;
    }

    return false;
}



//  Template Method
//
Float RayMarcher::sdf(const Point3f &pos) const {
    // Part 5:
    Float distanceToCenter = ((*WorldToObject)(Vector3f(pos.x, pos.y, pos.z))).Length();
    Float signedDistance = distanceToCenter - radius;
    return signedDistance;
}


// Get Normal using Gradient (Finite Distances Methods )  - See class slides.
//  Note if the normal you calculate has zero length, return the defaultNormal
//
Vector3f RayMarcher::GetNormalRM( const Point3f &p, float eps, const Vector3f &defaultNormal) const {
    float sdfCenter = sdf(pos);
    if (sdfCenter == 0)
        return defaultNormal;

    Vector3f grad(
           sdf(Point3f(pos.x + eps, pos.y, pos.z)) - sdf(Point3f(pos.x - eps, pos.y, pos.z)),
           sdf(Point3f(pos.x, pos.y + eps, pos.z)) - sdf(Point3f(pos.x, pos.y - eps, pos.z)),
           sdf(Point3f(pos.x, pos.y, pos.z + eps)) - sdf(Point3f(pos.x, pos.y, pos.z - eps))
    );

    grad = Normalize(grad);

    return grad;
}


Float RayMarcher::Area() const { return phiMax * radius * (zMax - zMin); }

// These functions are stubbed
//
Interaction RayMarcher::Sample(const Point2f &u, Float *pdf) const {
    LOG(FATAL) << "RayMarcher::Sample not implemented.";
    return Interaction();
}

Interaction RayMarcher::Sample(const Interaction &ref, const Point2f &u,
                           Float *pdf) const {
    LOG(FATAL) << "RayMarcher::Sample not implemented.";
    return Interaction();
}

std::shared_ptr<Shape> CreateRayMarcherShape(const Transform *o2w,
                                         const Transform *w2o,
                                         bool reverseOrientation,
                                         const ParamSet &params) {
    Float radius = params.FindOneFloat("radius", 1.f); 
    Float zmin = params.FindOneFloat("zmin", -radius);
    Float zmax = params.FindOneFloat("zmax", radius);
    Float phimax = params.FindOneFloat("phimax", 360.f);
    return std::make_shared<RayMarcher>(o2w, w2o, reverseOrientation, radius, zmin,
                                    zmax, phimax);
}

}  // namespace pbrt
