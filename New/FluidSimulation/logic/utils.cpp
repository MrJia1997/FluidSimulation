#include "utils.h"

Vector3d QtEigen (const QVector3D &_v) {
    Vector3d v;
    v(0) = _v.x();
    v(1) = _v.y();
    v(2) = _v.z();
    return v;
}

QVector3D EigenQt (const Vector3d &_v) {
    QVector3D v;
    v.setX(_v(0));
    v.setY(_v(1));
    v.setZ(_v(2));
    return v;
}
