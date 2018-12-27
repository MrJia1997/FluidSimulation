#ifndef GLWALLS_H
#define GLWALLS_H

#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions>
#include <QVector3D>
#include <QOpenGLShaderProgram>

class GLWalls: protected QOpenGLFunctions
{
public:
    GLWalls(double _xmin, double _xmax, double _ymin, double _ymax, double _zmin, double _zmax);
    ~GLWalls();

public:
    float x[10], y[10], z[10];
    int c = 0;

private:
    float v[200];

    QOpenGLVertexArrayObject vao;
    unsigned int vbo_v;

public:
    void compile();
    void render();
};

#endif
