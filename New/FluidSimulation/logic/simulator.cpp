#include "constants.h"
#include "kernel.h"
#include "simulator.h"

#include <omp.h>
#include <qDebug>
#include <fstream>
#include <cmath>

using namespace SmoothKernel;

Simulator::Simulator(): gravity(0.0, -9.8, 0.0)
{
    numParticles = 0;
    buildCubeFluid(QVector3D(0.0, 0.0, 0.5), 10, 10, 10);
}

Simulator::~Simulator()
{
    // do nothing
}

void Simulator::buildCubeFluid(const QVector3D& center, int numParticleX, int numParticleY, int numParticleZ)
{
    double sizeX = (numParticleX - 1) * ADJ_DISTANCE;
    double sizeY = (numParticleY - 1) * ADJ_DISTANCE;
    double sizeZ = (numParticleZ - 1) * ADJ_DISTANCE;

    double xmin = center.x() - sizeX / 2;
    double ymin = center.y() - sizeY / 2;
    double zmin = center.z() - sizeZ / 2;

    for (int x = 0; x < numParticleX; x++)
    {
        for (int y = 0; y < numParticleY; y++)
        {
            for (int z = 0; z < numParticleZ; z++)
            {
                Particle p;
                p.position = QVector3D(xmin + x * ADJ_DISTANCE, ymin + y * ADJ_DISTANCE, zmin + z * ADJ_DISTANCE);
                p.predPosition = p.position;
                p.neighbors.clear();
                p.id = numParticles;
                particles.push_back(p);
                ++numParticles;
            }
        }
    }
}

void Simulator::clearFluid()
{
    particles.clear();
    container.clear();
    numParticles = 0;
}

double Simulator::getDensity(const Particle& p)
{
    double density = 0.0;
    for (int neighbor: p.neighbors)
    {
        QVector3D r = p.predPosition - particles[neighbor].predPosition;
        density += Poly6(r);
    }
    return density;
}

void Simulator::simulate()
{
    int iter = 0, size = numParticles;

    // Step 1: apply external forces
#pragma omp parallel for
    for (iter = 0; iter < size; iter++)
    {
        particles[iter].velocity += gravity * DELTA_T;
        particles[iter].predPosition = particles[iter].position + particles[iter].velocity * DELTA_T;
        // p.position = p.predPosition;
    }

    // Step 2: find neighbor particles
    container.clear();
    int _id = -1;
    for (const Particle& p: particles)
    {
        container.add(++_id, p);
    }

#pragma omp parallel for
    for (iter = 0; iter < size; iter++)
    {
        container.findNeighbors(particles[iter], particles);
    }

    // Step 3: calculate partical constraints
    for (int _ = 0; _ < SOLVER_ITERATIONS; ++_)
    {
        // calculate lambda
        double maxdensity = -1.0;
        double mindensity = 99999.0;
        int maxneigh = 0.0;

#pragma omp parallel for
        for (iter = 0; iter < size; iter++) {
            double density = getDensity(particles[iter]);
            if (density > maxdensity) maxdensity = density;
            if (density < mindensity) mindensity = density;
            if (particles[iter].neighbors.size() > maxneigh) maxneigh = particles[iter].neighbors.size();
            double constraint = density / REST_DENSITY - 1;

            //
            // if (constraint <= 0) constraint = 0;
            //

            double denominator = 0.0;
            QVector3D gradsum;

            for (int neighbor: particles[iter].neighbors)
            {
                if (neighbor != particles[iter].id)
                {
                    QVector3D dummy = SpikyGradient(particles[iter].predPosition - particles[neighbor].predPosition) / REST_DENSITY;

                    denominator += dummy.lengthSquared();
                    gradsum += dummy;
                }
            }

            denominator += gradsum.lengthSquared();
            particles[iter].lambda = -constraint / (denominator + CFM_EPSILON);
        }

        // qDebug() << "max density = " << maxdensity << mindensity << _ << maxneigh;
        // calculate delta p & perform collision
        int parti_cnt = -1;
#pragma omp parallel for
        for (iter = 0; iter < size; iter++) {
            QVector3D deltapsum;

            for (int neighbor: particles[iter].neighbors)
            {
                if (neighbor != particles[iter].id)
                {
                    double scorr = -TENSILE_K * pow(Poly6(particles[iter].predPosition - particles[neighbor].predPosition) / Poly6(QVector3D(TENSILE_DELTA_Q, 0.0, 0.0)), TENSILE_N);
                    /*if (flag_sorr)
                    {
                        qDebug() << "sorr = " << scorr;
                        flag_sorr = false;
                    }*/
                    // scorr = 0.0;
                    deltapsum += (particles[iter].lambda + particles[neighbor].lambda + scorr) * SpikyGradient(particles[iter].predPosition - particles[neighbor].predPosition);
                }
            }

            particles[iter].deltaPosition = deltapsum / REST_DENSITY;
            /*if (p.deltaPosition.y() > 0.5)
            {
                qDebug() << "wtf";
                qDebug() << p.lambda << getDensity(p) << p.deltaPosition << p.velocity << p.id;
                while(1);
            }*/
            ++parti_cnt;
            // qDebug() << parti_cnt << p.deltaPosition << _;
            // p.deltaPosition =
        }

        // update position
#pragma omp parallel for
        for (iter = 0; iter < size; iter++) {
            // if (!flag) qDebug() << p.deltaPosition;
            particles[iter].predPosition += particles[iter].deltaPosition;
            container.handleCollision(particles[iter]);
            /*if (p.predPosition.z() < -0.2 || p.predPosition.z() > 1.2)
            {
                qDebug() << "question";
            }*/
        }
    }

    // Step 4: update velocity
#pragma omp parallel for
    for (iter = 0; iter < size; iter++) {
        particles[iter].velocity = 1.0 / DELTA_T * (particles[iter].predPosition - particles[iter].position);

        // vorticity
        QVector3D v_vorticity;
        QVector3D omega, omegax0, omegax1, omegay0, omegay1, omegaz0, omegaz1;

        for (int neighbor: particles[iter].neighbors)
        {
            QVector3D vij = particles[neighbor].velocity - particles[iter].velocity;
            QVector3D pij = particles[iter].predPosition - particles[neighbor].predPosition;

            QVector3D gradx0 = SpikyGradient(pij - QVector3D(VORTICITY_DELTA, 0.0, 0.0));
            QVector3D gradx1 = SpikyGradient(pij + QVector3D(VORTICITY_DELTA, 0.0, 0.0));

            QVector3D grady0 = SpikyGradient(pij - QVector3D(0.0, VORTICITY_DELTA, 0.0));
            QVector3D grady1 = SpikyGradient(pij + QVector3D(0.0, VORTICITY_DELTA, 0.0));

            QVector3D gradz0 = SpikyGradient(pij - QVector3D(0.0, 0.0, VORTICITY_DELTA));
            QVector3D gradz1 = SpikyGradient(pij + QVector3D(0.0, 0.0, VORTICITY_DELTA));

            omega += QVector3D::crossProduct(vij, SpikyGradient(pij));

            omegax0 += QVector3D::crossProduct(vij, gradx0);
            omegax1 += QVector3D::crossProduct(vij, gradx1);
            omegay0 += QVector3D::crossProduct(vij, grady0);
            omegay1 += QVector3D::crossProduct(vij, grady1);
            omegaz0 += QVector3D::crossProduct(vij, gradz0);
            omegaz1 += QVector3D::crossProduct(vij, gradz1);
        }

        QVector3D eta((omegax1.length() - omegax0.length()) / (VORTICITY_DELTA * 2.0),
                  (omegay1.length() - omegay0.length()) / (VORTICITY_DELTA * 2.0),
                  (omegaz1.length() - omegaz0.length()) / (VORTICITY_DELTA * 2.0));

        QVector3D f_vorticity = VORTICITY_EPSILON * QVector3D::crossProduct(eta.normalized(), omega);
        v_vorticity = f_vorticity * REST_DENSITY * DELTA_T * 0.0;

        // viscosity
        QVector3D v_viscosity;

        for (int neighbor: particles[iter].neighbors)
        {
            v_viscosity += (particles[neighbor].velocity - particles[iter].velocity) * Poly6(particles[iter].predPosition - particles[neighbor].predPosition);
        }

        v_viscosity *= VISCOSITY_C * 0.001;

        // update
        particles[iter].velocity += v_vorticity + v_viscosity;

        particles[iter].position = particles[iter].predPosition;

    }

    // set 719
    // qDebug() << "that 719" << particles[719].position << particles[719].velocity;
}

void Simulator::calcIsosurface() {
    // Implemented according to
    // Reconstructing Surfaces of Particle-Based Fluids Using Anisotropic Kernels
    // Yu J, Turk G. Siggraph 2010


}
