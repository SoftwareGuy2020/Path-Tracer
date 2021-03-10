#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>


using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO:
  // Implement MirrorBSDF
    *pdf = 1.0;
    BSDF::reflect(wo, wi);

    return reflectance / abs_cos_theta(*wi);

}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D& wo, const Vector3D& wi) {
  return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D& h) {
  // TODO: proj3-2, part 3
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.
    Vector3D n(0, 0, 1);
    
    double thetaH = acos(cos(h.z / (h.norm() * n).norm()));
    return exp(-pow(tan(thetaH), 2.0) / pow(alpha, 2.0)) / (PI * pow(alpha, 2.0) * pow(cos(thetaH), 4.0));
}

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
  // TODO: proj3-2, part 3
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Spectrum.
    Spectrum rs = ((eta * eta + k * k) - 2 * eta * cos_theta(wi) + pow(cos_theta(wi), 2.0))
                  / ((eta * eta + k * k) + 2 * eta * cos_theta(wi) + pow(cos_theta(wi), 2.0));

    Spectrum rp = ((eta * eta + k * k) * pow(cos_theta(wi), 2.0) - 2 * eta * cos_theta(wi) + 1)
                  / ((eta * eta + k * k) * pow(cos_theta(wi), 2.0) + 2 * eta * cos_theta(wi) + 1);

    return (rs + rp) / 2.0;
    
    
}

Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  // TODO: proj3-2, part 3
  // Implement microfacet model here.
    Vector3D h = (wo + wi) / (wo + wi).norm();
    return F(wi) * G(wo, wi) * D(h) / (4 * wo.z * wi.z);
}

Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO: proj3-2, part 3
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.
  
  
  //return MicrofacetBSDF::f(wo, *wi);
    Vector2D r = sampler.get_sample();
    double r1 = r.x, r2 = r.y;
    double thetaH = atan(sqrt(-(alpha * alpha) * log(1 - r1))), phiH = 2 * PI * r2;    

    double pTheta = ((2 * sin(thetaH)) / ((alpha * alpha) * pow(cos(thetaH), 2.0))) * exp(-pow(tan(thetaH / (alpha * alpha)), 2));
    double pPhi = 0.5 / PI;
    Vector3D h(sin(thetaH) * cos(phiH), sin(thetaH) * sin(phiH), cos(thetaH));
    //Vector3D h(sin(phiH) * cos(thetaH), sin(phiH) * sin(thetaH), cos(phiH));    
    *wi = 2.0 * dot(h, wo) * h - wo;

    if (dot(*wi, h) > 0.0) {
        if (dot(Vector3D(0.0, 0.0, 1.0), *wi) > 0.0 && dot(Vector3D(0.0, 0.0, 1.0), wo) > 0.0) {
            
            //*wi = cosineHemisphereSampler.get_sample(pdf); //placeholder
            *pdf = (pTheta * pPhi / sin(thetaH)) / (4 * dot(*wi, h));
            /*std::cout << "pdf: " << *pdf << " h: " << std::endl << h << " theta_h" << thetaH << " phi_h: " << phiH << std::endl
                << " P_theta: " << pTheta << " P_phi: " << pPhi << std::endl << std::endl << " wi_x: " << wi->x << std::endl
                << " wi_y: " << wi->y << std::endl << " wi_z: " << wi->z << std::endl;*/
            return f(wo, *wi);
        }
    }
    *pdf = .000001;
    return Spectrum();
}

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO:
  // Implement RefractionBSDF
  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO:
  // Compute Fresnel coefficient and use it as the probability of reflection
  // - Fundamentals of Computer Graphics page 305
    if (!refract(wo, wi, this->ior)) {
        reflect(wo, wi);
        *pdf = 1;
        return reflectance / abs_cos_theta(*wi);
    }

    double ro = pow((this->ior - 1) / (this->ior + 1), 2);
    double r = ro + (1 - ro) * pow(1 - abs_cos_theta(wo), 5);

    if (coin_flip(r)) {
        reflect(wo, wi);
        *pdf = r;
        return r * reflectance / abs_cos_theta(*wi);
    }

    else {
        *pdf = 1 - r;
        double eta = this->ior;
        if (wo.z > 0.0) {
            eta = 1.0 / eta;
        }
        return *pdf * transmittance / abs_cos_theta(*wi) * pow(eta, 2);
    }

}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {
  // TODO:
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
    
    *wi = 2 * Vector3D(0, 0, wo.z) - wo;

}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {
  // TODO:
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
    int sign = 1;

    if (wo.z > 0.0) {
        ior = 1.0 / ior;
        sign = -1;
    }
    double internalReflection = 1.0 - pow(ior, 2) * (1.0 - pow(wo.z, 2));

    if (internalReflection < 0.0)
        return false;

    wi->x = -ior * wo.x;
    wi->y = -ior * wo.y;
    wi->z = sign * sqrt(internalReflection);

    return true;
}

} // namespace CGL
