#pragma once

#include "geometry.h"

using std::vector;

namespace FluidSimulation {

    class Simulator {
    public:
        vector<vector<Particle>> grids;
        CubeContainer *container = nullptr;
        
        double containerSize;
        double defaultDensity;
        double neighborDistance;
        double numGrids;
        int numParticles;

        Vec3D gravity;

        double timeStep;
        double solverIteration;

    public:
        Simulator();
        ~Simulator();

        int positionToIndex(Vec3D position);
        int gridToIndex(int gridx, int gridy, int gridz);
        void buildCubeFluid(Vec3D center, double size, double particleNum);
        void clearFluid();
        // Calculate density using Poly6 kernel
        double density(Particle & particle);

        void simulate();
    };
}