//#define DEBUG_MODE

#include "constants.h"
#include "kernel.h"
#include "simulator.h"
#include "marchingcube.h"

#include <omp.h>
#include <qDebug>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/SVD>

using namespace SmoothKernel;

Simulator::Simulator(): gravity(0.0, -9.8, 0.0)
{
    numParticles = 0;
    buildCubeFluid(QVector3D(0.0, 0.0, 0.5), 10, 10, 10);
    printScalar = true;
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

void Simulator::calcIsosurface(std::vector<Triangle> &mesh) {
    // Implemented according to
    // Reconstructing Surfaces of Particle-Based Fluids Using Anisotropic Kernels
    // Yu J, Turk G. Siggraph 2010

    typedef Eigen::Vector3d Vec3d;
    typedef Eigen::Matrix3d Mat3d;
    using namespace Eigen;

    auto weight = [](const Vec3d& r) {
        double r_norm = r.norm();
        if (r_norm < KERNEL_H)
            return 1 - pow(r_norm / KERNEL_H, 3);
        else
            return 0.0;
    };

    auto vvT = [](const Vec3d& v) {
        return v * v.transpose();
    };

    int i = 0, size = particles.size();

    // Calculate Gi
    std::vector<Vec3d> x_weighted(size,Vec3d::Zero());
    std::vector<Vec3d> x_bar(size, Vec3d::Zero());
    std::vector<double> density(size, 0.0);
    std::vector<Mat3d> C(size, Mat3d::Zero());
    std::vector<Mat3d> G(size, Mat3d::Zero());

    double lambda = 0.9;
#pragma omp parallel for
    for (i = 0; i < size; i++) {
        if (particles[i].neighbors.empty())
            continue;

        double weight_sum = 0.0;
        for (int neighbor : particles[i].neighbors) {
            Vec3d r = QtEigen(particles[i].position - particles[neighbor].position);
            double w = weight(r);
            weight_sum += w;
            x_weighted[i] += w * QtEigen(particles[neighbor].position);
        }
        x_weighted[i] /= weight_sum;
        x_bar[i] = (1 - lambda) * QtEigen(particles[i].position) + lambda * x_weighted[i];
    }

#pragma omp parallel for
    for (i = 0; i < size; i++) {
        if (particles[i].neighbors.empty())
            continue;

        double weight_sum = 0.0;
        for (int neighbor : particles[i].neighbors) {
            Vec3d r = QtEigen(particles[i].position - particles[neighbor].position);
            double w = weight(r);
            weight_sum += w;
            C[i] += w * vvT(QtEigen(particles[neighbor].position) - x_weighted[i]);
        }
        C[i] /= weight_sum;
    }


    double kr = 4, ks = 1400, kn = 0.5, Neps = 25;

#pragma omp parallel for
    for (i = 0; i < size; i++) {
        // Here C[i] is symmetric for sure
        // Using SVD and we will get C[i] = U(\Sigma)V', here U = V
        JacobiSVD<Mat3d> svd(C[i], ComputeFullU | ComputeFullV);
        Mat3d u = svd.matrixU(),
                eigens = svd.singularValues().asDiagonal(),
                v = svd.matrixV();

        Mat3d eigens_tilde = Mat3d::Identity();

        if (particles[i].neighbors.size() > Neps) {
            eigens_tilde(0, 0) = eigens(0, 0);
            eigens_tilde(1, 1) = std::max(eigens(1, 1), eigens(0, 0) / kr);
            eigens_tilde(2, 2) = std::max(eigens(2, 2), eigens(0, 0) / kr);
            eigens_tilde *= ks;
        }
        else {
            eigens_tilde *= kn;
        }

        G[i] = u * eigens_tilde.inverse() * v.transpose();
        G[i] /= ADJ_DISTANCE;
    }


    // Calculate scalar field
#pragma omp parallel for
    for (i = 0; i < size; i++) {
        for (auto &j : particles[i].neighbors) {
            density[i] += Poly6(particles[j].position - particles[i].position);
        }
    }

    auto ScalarFunction = [&](QVector3D pos) {
        // Find neighbors
        double scalar = 0.0;
        IntTriple index = container.indexOfPosition(pos.x(), pos.y(), pos.z());
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                for (int dz = -1; dz <= 1; dz++) {
                    IntTriple neighborIndex = index.addOffset(dx, dy, dz);
                    if (container.particleMap.count(neighborIndex) == 0) continue;
                    for (auto& j: container.particleMap[neighborIndex]) {
                        scalar += Poly6Anisotropic(pos - EigenQt(x_bar[j]), G[j]) / density[j];
                    }
                }
            }
        }
        return scalar;
    };

#ifdef DEBUG_MODE
    std::ofstream scalar_output("scalar_output.txt", std::ios::out);
    std::ofstream mesh_output("mesh_output.txt", std::ios::out);
#endif
    // Marching Cube Algorithm

    double x_max = X_MAX + MARCHINGCUBE_DISTANCE,
            x_min = X_MIN - MARCHINGCUBE_DISTANCE,
            y_max = Y_MAX + MARCHINGCUBE_DISTANCE,
            y_min = Y_MIN - MARCHINGCUBE_DISTANCE,
            z_max = Z_MAX + MARCHINGCUBE_DISTANCE,
            z_min = Z_MIN - MARCHINGCUBE_DISTANCE;

    int cube_number_x = std::floor((x_max - x_min) / MARCHINGCUBE_DISTANCE);
    int cube_number_y = std::floor((y_max - y_min) / MARCHINGCUBE_DISTANCE);
    int cube_number_z = std::floor((z_max - z_min) / MARCHINGCUBE_DISTANCE);
    int vertex_number = (cube_number_x + 1) * (cube_number_y + 1) * (cube_number_z + 1);
    int cube_number = cube_number_x * cube_number_y * cube_number_z;
    std::vector<int> scalar_storage(vertex_number, 0);

#pragma omp parallel for
    for (int i = 0; i < vertex_number; i++) {
        int ix = i / ((cube_number_y + 1) * (cube_number_z + 1)),
                iy = (i % ((cube_number_y + 1) * (cube_number_z + 1))) / (cube_number_z + 1),
                iz = i % (cube_number_z + 1);
        double x = x_min + (double)ix * MARCHINGCUBE_DISTANCE,
                y = y_min + (double)iy * MARCHINGCUBE_DISTANCE,
                z = z_min + (double)iz * MARCHINGCUBE_DISTANCE;

        QVector3D pos(x, y, z);
        scalar_storage[i] = ScalarFunction(pos);
    }

#ifdef DEBUG_MODE
    if (printScalar) {
        for (int i = 0; i < vertex_number; i++) {
            int ix = i / ((cube_number_y + 1) * (cube_number_z + 1)),
                    iy = (i % ((cube_number_y + 1) * (cube_number_z + 1))) / (cube_number_z + 1),
                    iz = i % (cube_number_z + 1);
            double x = x_min + (double)ix * MARCHINGCUBE_DISTANCE,
                    y = y_min + (double)iy * MARCHINGCUBE_DISTANCE,
                    z = z_min + (double)iz * MARCHINGCUBE_DISTANCE;

            scalar_output << x << " " << y << " " << z << " " << scalar_storage[i] << std::endl;
        }
        qDebug() << "Scalar output finished.";
    }
#endif

    double isolevel = -200.0;
    std::vector<std::vector<Triangle>> surface_meshes(cube_number, std::vector<Triangle>());

    IntTriple delta[] = {
        IntTriple(0, 0, 0),
        IntTriple(0, 1, 0),
        IntTriple(1, 1, 0),
        IntTriple(1, 0, 0),
        IntTriple(0, 0, 1),
        IntTriple(0, 1, 1),
        IntTriple(1, 1, 1),
        IntTriple(1, 0, 1)
    };

#pragma omp parallel for
    for (int i = 0; i < cube_number; i++) {
        int ix = i / (cube_number_y * cube_number_z),
                iy = (i % (cube_number_y * cube_number_z)) / cube_number_z,
                iz = i % cube_number_z;

        Cube cube;

        for (int j = 0; j < 8; j++) {
            int v_ix = ix + delta[j].x;
            int v_iy = iy + delta[j].y;
            int v_iz = iz + delta[j].z;

            int v_i = v_iz + (cube_number_z + 1) * (v_iy + (cube_number_y + 1) * v_ix);
            double x = x_min + (double)v_ix * MARCHINGCUBE_DISTANCE,
                    y = y_min + (double)v_iy * MARCHINGCUBE_DISTANCE,
                    z = z_min + (double)v_iz * MARCHINGCUBE_DISTANCE;
            cube.p[j] = QVector3D(x, y, z);
            cube.val[j] = scalar_storage[v_i];
        }

        Polygonise(cube, isolevel, surface_meshes[i]);
    }

    mesh.clear();
    for (auto surface : surface_meshes) {
        mesh.insert(mesh.end(), surface.begin(), surface.end());
    }

#ifdef DEBUG_MODE
    if (printScalar) {
        for (auto t : mesh) {
            for (int i = 0; i < 3; i++) {
                mesh_output << t.p[i].x() << " " << t.p[i].y() << " " << t.p[i].z() << " ";
            }
            mesh_output << std::endl;
        }
        qDebug() << "Mesh output finished.";
        printScalar = false;
    }

    scalar_output.close();
    mesh_output.close();
#endif


}
