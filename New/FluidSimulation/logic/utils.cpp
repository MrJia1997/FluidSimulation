#include "utils.h"

Vector3d QtEigen (const QVector3D &_v) {
    Vector3d v;
    v(0) = _v.x();
    v(1) = _v.y();
    v(2) = _v.z();
    return v;
}
