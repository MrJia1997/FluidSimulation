#ifndef GLSURFACE_H
#define GLSURFACE_H

#include "logic/geometry.h"

#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions>
#include <QVector3D>
#include <QOpenGLShaderProgram>

class GLSurface: protected QOpenGLFunctions
{
public:
    GLSurface();
    ~GLSurface();

public:
    std::vector<Triangle> isosurface;
    float *vs;

    int size;
private:
    QOpenGLVertexArrayObject vao;
    unsigned int vbo_v;

public:
    void compile();
    void render();
    void myupdate();
    void hack();
};

#endif // GLSURFACE_H
