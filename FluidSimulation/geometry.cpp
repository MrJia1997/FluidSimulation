#include "geometry.h"

using namespace FluidSimulation;

CubeContainer::CubeContainer() :
    xmin(-1.0), xmax(1.0),
    ymin(-1.0), ymax(1.0),
    zmin(-1.0), zmax(1.0) {}

CubeContainer::CubeContainer(Vec3D center, double size) :
    xmin(center.x() - size / 2), xmax(center.x() + size / 2),
    ymin(center.y() - size / 2), ymax(center.y() + size / 2),
    zmin(center.z() - size / 2), zmax(center.z() + size / 2) {}


// Refer to https://github.com/kvmohyi/fluid-simulation-project
void CubeContainer::handleCollision(Particle & particle) {
    double lossFactor = 0.6;

    Vec3D afterDelta = particle.predPosition + particle.deltaPosition;

    if (afterDelta.x() < xmin) {
        particle.deltaPosition.x() = 2 * xmin - 2 * particle.predPosition.x() - particle.deltaPosition.x();
        //particle.velocity.x() = (-1.0) * lossFactor * particle.velocity.x();
    }
    else if (afterDelta.x() > xmax) {
        particle.deltaPosition.x() = 2 * xmax - 2 * particle.predPosition.x() - particle.deltaPosition.x();
        //particle.velocity.x() = (-1.0) * lossFactor * particle.velocity.x();
    }

    if (afterDelta.y() < ymin) {
        particle.deltaPosition.y() = 2 * ymin - 2 * particle.predPosition.y() - particle.deltaPosition.y();
        //particle.velocity.y() = (-1.0) * lossFactor * particle.velocity.y();
    }
    else if (afterDelta.y() > ymax) {
        particle.deltaPosition.y() = 2 * ymax - 2 * particle.predPosition.y() - particle.deltaPosition.y();
        //particle.velocity.y() = (-1.0) * lossFactor * particle.velocity.y();
    }

    if (afterDelta.z() < zmin) {
        particle.deltaPosition.z() = 2 * zmin - 2 * particle.predPosition.z() - particle.deltaPosition.z();
        //particle.velocity.z() = (-1.0) * lossFactor * particle.velocity.z();
    }
    else if (afterDelta.z() > zmax) {
        particle.deltaPosition.z() = 2 * zmax - 2 * particle.predPosition.z() - particle.deltaPosition.z();
        //particle.velocity.z() = (-1.0) * lossFactor * particle.velocity.z();
    }
}

bool CubeContainer::isInside(Vec3D position) {

    if (position.x() > xmin && position.x() < xmax &&
        position.y() > ymin && position.y() < ymax &&
        position.z() > zmin && position.z() < zmax)
        return true;
    else
        return false;
}

