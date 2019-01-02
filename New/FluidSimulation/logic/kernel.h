#include "constants.h"
#include "utils.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <QVector3D>
#include <cmath>

namespace SmoothKernel
{
    double Poly6(const QVector3D& r)
    {
        // W_Poly6(r, h) = 315 / (64pi * h^9) * (h^2 - |r|^2)^3

        double r_norm = r.length();
        double h = KERNEL_H;
        
        if (r_norm >= h)
        {
            return 0.0;
        }

        return 315.0 / (64.0 * MATH_PI * pow(h, 9)) * pow((h * h - r_norm * r_norm), 3);
    }

    double Poly6Anisotropic(const QVector3D& _r, const Eigen::Matrix3d& G) {
        // For Poly6 kernel
        // \sigma = 315 / (64 * pi)
        // P(x) = (1 - x^2)^3

        double r_norm = _r.length();
        double h = KERNEL_H;

        if (r_norm >= h)
        {
            return 0.0;
        }

        Eigen::Vector3d r = QtEigen(_r);
        double G_det = G.determinant();
        double Gr_norm = (G * r).norm();

        return 315.0 / (64.0 * MATH_PI) * G_det * pow((1 - Gr_norm * Gr_norm), 3);
    }

    QVector3D SpikyGradient(const QVector3D& r)
    {
        // grad W_Spiky(r, h) = -45 / (pi * h^6) * (h - r)^2 * hat(r)

        double r_norm = r.length();
        double h = KERNEL_H;

        if (r_norm <= EPS || r_norm >= h) {
            return QVector3D(0.0, 0.0, 0.0);
        }

        return -45.0 / (MATH_PI * pow(h, 6)) * pow((h - r_norm), 2) * r.normalized();
    }

    /*
    double Poly6(const double &r, const double &h) {
        if (r < 0 || r >= h)
            return 0.0;
        
        double q = h * h - r * r;
        double alpha = 315.0 / (64.0 * PI * pow(h, 9));
        return alpha * pow(q, 3);
    }

    QVector3D Poly6Gradient(const double &r, const QVector3D &rij, const double &h) {
        if (r < 0 || r >= h)
            return QVector3D(0.0, 0.0, 0.0);

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

    QVector3D SpikyGradient(const double &r, const QVector3D &rij, const double &h) {
        if (r < 0 || r >= h)
            return QVector3D(0.0, 0.0, 0.0);
        
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
    */
}
