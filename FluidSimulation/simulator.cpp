#include "simulator.h"
#include "kernel.h"

#include <fstream>
#include <cmath>

using namespace FluidSimulation;
using namespace SmoothKernel;
using std::vector;

Simulator::Simulator() {
    
    container = new CubeContainer();
    containerSize = 2.0;
    defaultDensity = 1000;
    neighborDistance = 0.15;
    numGrids = std::ceil(containerSize / neighborDistance);
    numParticles = 0;
    grids.resize(numGrids * numGrids * numGrids);

    gravity = Vec3D(0.0, -9.8, 0.0);

    timeStep = 0.01;
    solverIteration = 10;

    buildCubeFluid(Vec3D(0.0, 0.0, 0.0), 1, 10);
}

Simulator::~Simulator() {
    if (container != nullptr)
        delete container;
    container = nullptr;
}

int Simulator::positionToIndex(Vec3D position) {
    if (container == nullptr || !container->isInside(position))
        return -1;

    int gridx = (position.x() - container->xmin) / neighborDistance,
        gridy = (position.y() - container->ymin) / neighborDistance,
        gridz = (position.z() - container->zmin) / neighborDistance;

    return gridToIndex(gridx, gridy, gridz);
}

int Simulator::gridToIndex(int gridx, int gridy, int gridz)
{
    if (gridx < 0 || gridx >= numGrids ||
        gridy < 0 || gridy >= numGrids ||
        gridz < 0 || gridz >= numGrids)
        return -1;
    
    return gridz + numGrids * (gridy + numGrids * gridx);
}

// Here size and numParticle are for edges
// That means we totally have (numParticle ^ 3) particles
void Simulator::buildCubeFluid(Vec3D center, double size, double numParticle) {
    double xmin = center.x() - size / 2,
        ymin = center.y() - size / 2,
        zmin = center.z() - size / 2;

    double delta = size / numParticle;
    
    int iterId = numParticles;
    numParticles += numParticle * numParticle * numParticle;

    for (int x = 0; x < numParticle; x++) {
        for (int y = 0; y < numParticle; y++) {
            for (int z = 0; z < numParticle; z++) {
                Vec3D position(xmin + x * delta, ymin + y * delta, zmin + z * delta);
                Particle p;
                p.id = iterId;
                iterId++;

                p.position = position;
                grids[positionToIndex(position)].push_back(p);
            }
        }
    }
}

void Simulator::clearFluid() {
    for (auto grid : grids) {
        grid.clear();
    }
    numParticles = 0;
}

double Simulator::density(Particle & particle) {
    double rho = 0.0;

    for (auto neighbor : particle.neighbors) {
        Vec3D rij = particle.predPosition - neighbor->predPosition;
        double distance = rij.norm();
        rho += Poly6(distance, neighborDistance);
    }

    return rho;
}

void Simulator::simulate() {
    // Step 1: apply external forces
    for (auto grid : grids) {
        for (auto particle : grid) {
            particle.velocity += particle.extForce * timeStep;
            particle.predPosition = particle.position + particle.velocity * timeStep;
        }
    }
    
    // Step 2: find neighbor particles
    vector<vector<Particle>> newGrids(numGrids * numGrids * numGrids);
    for (auto grid : grids) {
        for (auto particle : grid) {
            newGrids[positionToIndex(particle.predPosition)].push_back(particle);
        }
    }
    for (auto nGrid : newGrids) {
        for (auto particle : nGrid) {
            for (int xoff = -1; xoff <= 1; xoff++) {
                for (int yoff = -1; yoff <= 1; yoff++) {
                    for (int zoff = -1; zoff <= 1; zoff++) {
                        particle.neighbors.clear();
                        Vec3D predPos = particle.predPosition;

                        int gridx = (predPos.x() - container->xmin) / neighborDistance,
                            gridy = (predPos.y() - container->ymin) / neighborDistance,
                            gridz = (predPos.z() - container->zmin) / neighborDistance;
                        
                        int index = gridToIndex(gridx + xoff, gridy + yoff, gridz + zoff);
                        if (index < 0 || index >= newGrids.size())
                            continue;

                        for (Particle other : newGrids[index]) {
                            Vec3D rij = particle.predPosition - other.predPosition;
                            if (rij.norm() < neighborDistance)
                                particle.neighbors.push_back(&other);
                        }
                    }
                }
            }
        }
    }
    // Step 3: calculate partical constraints
    int iter = 0;
    double epsilon = 1;
    while (iter < solverIteration) {
        // Calculate lambda using Spiky kernel
        for (auto nGrid : newGrids) {
            for (auto particle : nGrid) {
                double currentDensity = density(particle);
                double constraint = currentDensity / defaultDensity - 1;
                
                Vec3D gradpiCi(0.0, 0.0, 0.0);
                double denominator = 0.0;

                for (auto neighbor : particle.neighbors) {
                    if (neighbor->id == particle.id)
                        continue;
                    
                    Vec3D rij = particle.predPosition - neighbor->predPosition;
                    Vec3D grad = SpikyGradient(rij.norm(), rij, neighborDistance) / defaultDensity;
                    gradpiCi += grad;
                    denominator += grad.squaredNorm();
                }
                denominator += gradpiCi.squaredNorm();

                particle.lambda = (-constraint) / (denominator + epsilon);
            }
        }
        // Calculate deltaPosition
        double k = 0.1, n = 4, qNorm = 0.2 * neighborDistance;
        for (auto nGrid : newGrids) {
            for (auto particle : nGrid) {
                Vec3D deltaPosition(0.0, 0.0, 0.0);
                
                for (auto neighbor : particle.neighbors) {
                    Vec3D rij = particle.predPosition - neighbor->predPosition;
                    double scorr = (-k) * pow(
                        Spiky(rij.norm(), neighborDistance) / Spiky(qNorm, neighborDistance),
                        4);

                    deltaPosition += (particle.lambda + neighbor->lambda + scorr)
                        * SpikyGradient(rij.norm(), rij, neighborDistance) / defaultDensity;
                }
                
                particle.deltaPosition = deltaPosition;
                container->handleCollision(particle);
            }
        }
        // Update predict position
        for (auto nGrid : newGrids) {
            for (auto particle : nGrid) {
                particle.predPosition += particle.deltaPosition;
            }
        }

        iter++;
    }
    for (auto nGrid : newGrids) {
        for (auto particle : nGrid) {
            particle.velocity = (particle.predPosition - particle.position) / timeStep;
        }
    }
    double c = 0.01;
    for (auto nGrid : newGrids) {
        for (auto particle : nGrid) {
            Vec3D vorticity(0.0, 0.0, 0.0);
            for (auto neighbor : particle.neighbors) {
                Vec3D vij = neighbor->velocity - particle.velocity;
                Vec3D rij = particle.predPosition - neighbor->predPosition;
                vorticity += vij.cross(SpikyGradient(rij.norm(), rij, neighborDistance));
            }

            particle.vorticity = vorticity;

            // TODO: add fvorticity
            particle.extForce = gravity;

            Vec3D viscosity(0.0, 0.0, 0.0);
            for (auto neighbor : particle.neighbors) {
                Vec3D vij = neighbor->velocity - particle.velocity;
                Vec3D rij = particle.predPosition - neighbor->predPosition;
                viscosity += Spiky(rij.norm(), neighborDistance) * vij;
            }

            particle.velocity += c * viscosity;
        }
    }
    
    for (auto nGrid : newGrids) {
        for (auto particle : nGrid) {
            particle.position = particle.predPosition;
        }
    }

    grids = newGrids;
}

