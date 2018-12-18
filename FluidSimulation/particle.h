#pragma once

#include <Eigen/Core>
#include <vector>

typedef Eigen::Vector3d Vec3D;
typedef Eigen::Vector4d Vec4D;

namespace FluidSimulation {
    // Refer to http://jamesjia.com/cs184-fluidsim
    // Here we assume each particle the has same mass as 1
    class Particle {
    public:
        Vec3D position;
        Vec3D pred_position;
        Vec3D delta_position;
        Vec3D velocity;
        Vec3D vorticity;
        Vec3D extForce;
        Vec4D color;
        std::vector<Particle*> neighbors;
        double lambda;
        int id;
        static const int size = 3;
    };

}