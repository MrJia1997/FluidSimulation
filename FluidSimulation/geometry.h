#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

typedef Eigen::Vector3d Vec3D;
typedef Eigen::Vector4d Vec4D;

namespace FluidSimulation {
    // Refer to http://jamesjia.com/cs184-fluidsim
    // Here we assume each particle the has same mass 1
    class Particle {
    public:
        Vec3D position;
        Vec3D predPosition;
        Vec3D deltaPosition;
        Vec3D velocity;
        Vec3D vorticity;
        Vec3D extForce;
        //Vec4D color;
        std::vector<Particle*> neighbors;
        double lambda;
        int id;
        //static const int size = 3;

    };

    class Container {
    public:
        virtual void handleCollision(Particle & particle) = 0;
        virtual bool isInside(Vec3D position) = 0;
    };

    class CubeContainer : public Container {
    public:
        double xmin, xmax, ymin, ymax, zmin, zmax;
    public:
        CubeContainer();
        CubeContainer(Vec3D center, double size);

        void handleCollision(Particle & particle);
        bool isInside(Vec3D position);
    };
}