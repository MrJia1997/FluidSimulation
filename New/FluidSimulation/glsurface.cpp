#include "glsurface.h"

GLSurface::GLSurface() : vs(nullptr), size(0) {
    initializeOpenGLFunctions();

}

GLSurface::~GLSurface() {
    if (vs != nullptr)
    delete[] vs;

    isosurface.clear();
}

void GLSurface::compile() {
    vao.create();
    vao.bind();

    glGenBuffers(1, &vbo_v);

    glBindBuffer(GL_ARRAY_BUFFER, vbo_v);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 9 * size, vs, GL_STREAM_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(0);

    vao.release();
}

void GLSurface::render() {
    vao.bind();
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glDrawArrays(GL_TRIANGLES, 0, size * 3);
    vao.release();
}

void GLSurface::myupdate() {
    if (vs != nullptr)
        delete[] vs;

    size = isosurface.size();
    vs = new float[size * 9];

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < 3; j++) {
            vs[9 * i + 3 * j] = isosurface[i].p[j].x();
            vs[9 * i + 3 * j + 1] = isosurface[i].p[j].y();
            vs[9 * i + 3 * j + 2] = isosurface[i].p[j].z();
        }
    }
}

void GLSurface::hack() {
    vao.bind();
    glBindBuffer(GL_ARRAY_BUFFER, vbo_v);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 9 * size, vs, GL_STREAM_DRAW);
    vao.release();
}
