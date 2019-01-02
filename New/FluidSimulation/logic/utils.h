#ifndef UTILS_H
#define UTILS_H

#include <Eigen/Core>
#include <QVector3D>

using Eigen::Vector3d;

Vector3d QtEigen (const QVector3D &_v);
QVector3D EigenQt (const Vector3d &_v);

#endif // UTILS_H
