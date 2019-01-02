#ifndef MARCHINGCUBE_H
#define MARCHINGCUBE_H

#include "geometry.h"

#include <QVector3D>
#include <vector>

// Refer to http://paulbourke.net/geometry/polygonise/


QVector3D VertexInterpolation(double isolevel, QVector3D p1, QVector3D p2, double val1, double val2);
void Polygonise(Cube cube, double isolevel, std::vector<Triangle> &triangles);


#endif // MARCHINGCUBE_H
