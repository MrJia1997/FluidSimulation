#include "glwalls.h"

#include <QDebug>

GLWalls::GLWalls(double _xmin, double _xmax, double _ymin, double _ymax, double _zmin, double _zmax)
{
    initializeOpenGLFunctions();

    x[0] = (float)_xmin;
    x[1] = (float)_xmax;
    y[0] = (float)_ymin;
    y[1] = (float)_ymax;
    z[0] = (float)_zmin;
    z[1] = (float)_zmax;

    qDebug() << "z" << z[0] << z[1];

    int cnt = -1;

    for (int ot = 0; ot < 2; ++ot)
    {
        if (ot == 0)
        {
            for (int i = 0; i < 2; ++i)
            {
                v[++cnt] = x[ot];
                v[++cnt] = y[0];
                v[++cnt] = z[1];
                v[++cnt] = x[ot];
                v[++cnt] = y[1];
                v[++cnt] = z[0];
                v[++cnt] = x[ot];
                v[++cnt] = y[i];
                v[++cnt] = z[i];
            }
        }
        else
        {
            for (int i = 0; i < 2; ++i)
            {
                v[++cnt] = x[ot];
                v[++cnt] = y[0];
                v[++cnt] = z[0];
                v[++cnt] = x[ot];
                v[++cnt] = y[1];
                v[++cnt] = z[1];
                v[++cnt] = x[ot];
                v[++cnt] = y[i];
                v[++cnt] = z[1 - i];
            }
        }
    }

    for (int ot = 0; ot < 2; ++ot)
    {
        if (ot == 0)
        {
            for (int i = 0; i < 2; ++i)
            {
                v[++cnt] = x[0];
                v[++cnt] = y[ot];
                v[++cnt] = z[1];
                v[++cnt] = x[1];
                v[++cnt] = y[ot];
                v[++cnt] = z[0];
                v[++cnt] = x[i];
                v[++cnt] = y[ot];
                v[++cnt] = z[i];
            }
        }
        else
        {
            for (int i = 0; i < 2; ++i)
            {
                v[++cnt] = x[0];
                v[++cnt] = y[ot];
                v[++cnt] = z[0];
                v[++cnt] = x[1];
                v[++cnt] = y[ot];
                v[++cnt] = z[1];
                v[++cnt] = x[i];
                v[++cnt] = y[ot];
                v[++cnt] = z[1 - i];
            }
        }
    }

    for (int ot = 0; ot < 2; ++ot)
    {
        if (ot == 0)
        {
            for (int i = 0; i < 2; ++i)
            {
                v[++cnt] = x[0];
                v[++cnt] = y[1];
                v[++cnt] = z[ot];
                v[++cnt] = x[1];
                v[++cnt] = y[0];
                v[++cnt] = z[ot];
                v[++cnt] = x[i];
                v[++cnt] = y[i];
                v[++cnt] = z[ot];
            }
        }
        else
        {
            for (int i = 0; i < 2; ++i)
            {
                v[++cnt] = x[0];
                v[++cnt] = y[0];
                v[++cnt] = z[ot];
                v[++cnt] = x[1];
                v[++cnt] = y[1];
                v[++cnt] = z[ot];
                v[++cnt] = x[i];
                v[++cnt] = y[1 - i];
                v[++cnt] = z[ot];
            }
        }
    }

    for (int i = 72; i < 108; i += 3)
    {
        qDebug() << v[i] << v[i + 1] << v[i + 2];
    }
}

void GLWalls::compile()
{
    vao.create();
    vao.bind();

    glGenBuffers(1, &vbo_v);

    glBindBuffer(GL_ARRAY_BUFFER, vbo_v);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 108, v, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(0);

    vao.release();
}

void GLWalls::render()
{
    vao.bind();
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glDrawArrays(GL_TRIANGLES, 0, 36);
    vao.release();
}
