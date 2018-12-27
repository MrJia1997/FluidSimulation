#ifndef LOGIC_PARTICLE_H
#define LOGIC_PARTICLE_H

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <QVector3D>
#include <vector>

class Particle
{
public:
    QVector3D position;       // x_i
    QVector3D predPosition;   // x_i_star
    QVector3D deltaPosition;  // x_i_delta
    QVector3D velocity;       // v_i
    QVector3D vorticity;
    //QVector3D extForce;
    //Vec4D color;
    std::vector<int> neighbors;
    double lambda;        // lambda_i
    int id;
    //static const int size = 3;
};

#endif
