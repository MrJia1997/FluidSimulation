#ifndef GLPARTICLES_H
#define GLPARTICLES_H

#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QVector3D>

class GLParticles: protected QOpenGLFunctions
{
public:
    GLParticles();
    ~GLParticles();

private:
    float* v;
    float* vt;
    float* vn;
    unsigned int* f;
    int v_cnt, vt_cnt, vn_cnt, f_cnt;

    float* v_all;
    float* vt_all;
    float* vn_all;
    unsigned int* f_all;

    QOpenGLVertexArrayObject vao;
    unsigned int vbo_v, vbo_vt, vbo_vn, ebo_f;

    QVector<QVector3D> positions;

public:
    void scale(float dd);
    void translate(float dx, float dy, float dz);
    void add(const QVector3D& pos);
    void compile();
    void render();
    void myupdate(int _id, const QVector3D& pos);
    void hack();
};

#endif
