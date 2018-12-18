#pragma once

#define PI 3.14159265358979323846

#include <Eigen/Core>
#include <cmath>

typedef Eigen::Vector3d Vec3D;

namespace SmoothKernel {
    // Refer to https://blog.csdn.net/Knight_Lyh/article/details/51597901
    double Poly6(const double &r, const double &h) {
        if (r < 0 || r >= h)
            return 0.0;
        
        double q = h * h - r * r;
        double alpha = 315.0 / (64.0 * PI * pow(h, 9));
        return alpha * pow(q, 3);
    }

    Vec3D Poly6Gradient(const double &r, const Vec3D &rij, const double &h) {
        if (r < 0 || r >= h)
            return Vec3D(0.0);

        double q = h * h - r * r;
        double alpha = (-945.0) / (32.0 * PI * pow(h, 9));
        return alpha * q * q * rij;
    }

    double Poly6Laplace(const double &r, const double &h) {
        if (r < 0 || r >= h)
            return 0.0;

        double q = h * h - r * r;
        double alpha = (-945.0) / (32.0 * PI * pow(h, 9));
        return alpha * (3 * q * q - 4 * r * r * q);
    }

    double Spiky(const double &r, const double &h) {
        if (r < 0 || r >= h)
            return 0.0;
        
        double q = h - r;
        double alpha = 15.0 / (PI * pow(h, 6));
        return alpha * pow(q, 3);
    }

    Vec3D SpikyGradient(const double &r, const Vec3D &rij, const double &h) {
        if (r < 0 || r >= h)
            return Vec3D(0.0);
        
        double q = h - r;
        double alpha = (-45.0) / (PI * pow(h, 6));
        return alpha * q * q * rij / r;
    }

    double SpikyLaplace(const double &r, const double &h) {
        if (r < 0 || r >= h)
            return 0.0;
        
        double q = h - r;
        double alpha = (-90.0) / (PI * pow(h, 6));
        return alpha * (q * q / r - q);
    }
}