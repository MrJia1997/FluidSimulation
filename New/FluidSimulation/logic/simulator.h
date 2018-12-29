#ifndef LOGIC_SIMULATOR_H
#define LOGIC_SIMULATOR_H

#include "geometry.h"


class Simulator
{
public:
    std::vector<Particle> particles;
    CubeContainer container;
    
    // double containerSize;
    // double defaultDensity;
    // double neighborDistance;
    // double numGrids;
    int numParticles;

    QVector3D gravity;

    double timeStep;
    double solverIteration;

public:
    Simulator();
    ~Simulator();

    void buildCubeFluid(const QVector3D& center, int numParticleX, int numParticleY, int numParticleZ);
    void clearFluid();
    double getDensity(const Particle& p);
    void simulate();
    void calcIsosurface();
};


#endif
