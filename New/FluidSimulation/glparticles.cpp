#include "GLParticles.h"

#include <QDebug>
#include <QFile>
#include <QTextStream>
#include <QDebug>
#include <QVector>
#include <algorithm>

GLParticles::GLParticles()
{
    initializeOpenGLFunctions();

    QFile file("./scene/sphere.obj");
    QTextStream in(&file);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        return;
    }

    QVector<float> vec_v, vec_vt, vec_vn;
    QVector<unsigned int> vec_f;

    float x_min = 1e10, x_max = -1e10;
    float y_min = 1e10, y_max = -1e10;
    float z_min = 1e10, z_max = -1e10;

    while (!in.atEnd())
    {
        QString line = in.readLine();
        QStringList opt = line.split(" ");
        if (opt[0] == 'v')
        {
            float x = opt[1].toFloat();
            float y = opt[2].toFloat();
            float z = opt[3].toFloat();

            x_min = std::min(x_min, x);
            x_max = std::max(x_max, x);
            y_min = std::min(y_min, y);
            y_max = std::max(y_max, y);
            z_min = std::min(z_min, z);
            z_max = std::max(z_max, z);

            vec_v.push_back(x);
            vec_v.push_back(y);
            vec_v.push_back(z);
        }
        else if (opt[0] == "vt")
        {
            vec_vt.push_back(opt[1].toFloat());
            vec_vt.push_back(opt[2].toFloat());
        }
        else if (opt[0] == "vn")
        {
            vec_vn.push_back(opt[1].toFloat());
            vec_vn.push_back(opt[2].toFloat());
            vec_vn.push_back(opt[3].toFloat());
        }
        else if (opt[0] == "f")
        {
            for (int i = 1; i <= 3; ++i)
            {
                QStringList elements = opt[i].split("/");
                for (int j = 0; j < 3; ++j)
                {
                    vec_f.push_back(elements[j].toUInt() - 1);
                }
            }
        }
    }

    v_cnt = vec_v.size() / 3;
    vt_cnt = vec_vt.size() / 2;
    vn_cnt = vec_vn.size() / 3;
    f_cnt = vec_f.size() / 9;

    float diff = std::max(x_max - x_min, std::max(y_max - y_min, z_max - z_min));

    for (int i = 0; i < v_cnt * 3; i += 3)
    {
        vec_v[i] = (vec_v[i] * 2.0 - x_min - x_max) / diff;
        vec_v[i + 1] = (vec_v[i + 1] * 2.0 - y_min - y_max) / diff;
        vec_v[i + 2] = (vec_v[i + 2] * 2.0 - z_min - z_max) / diff;
    }

    v = new float[f_cnt * 9];
    vt = new float[f_cnt * 6];
    vn = new float[f_cnt * 9];
    f = new unsigned int[f_cnt * 9];

    for (int i = 0, j = 0; i < f_cnt * 9; i += 3, ++j)
    {
        for (int k = 0; k < 3; ++k) v[j * 3 + k] = vec_v[vec_f[i] * 3 + k];
        for (int k = 0; k < 2; ++k) vt[j * 2 + k] = vec_vt[vec_f[i + 1] * 2 + k];
        for (int k = 0; k < 3; ++k) vn[j * 3 + k] = vec_vn[vec_f[i + 2] * 3 + k];
    }

    for (int i = 0; i < f_cnt * 9; ++i) f[i] = i;

    file.close();
}

GLParticles::~GLParticles()
{
    delete[] v;
    delete[] vt;
    delete[] vn;
    delete[] f;

    delete[] v_all;
    delete[] vt_all;
    delete[] vn_all;
    delete[] f_all;
}

void GLParticles::scale(float dd)
{
    for (int i = 0; i < f_cnt * 9; ++i)
    {
        v[i] *= dd;
    }
}

void GLParticles::translate(float dx, float dy, float dz)
{
    for (int i = 0; i < f_cnt * 9; i += 3)
    {
        v[i] += dx;
        v[i + 1] += dy;
        v[i + 2] += dz;
    }
}

void GLParticles::add(const QVector3D& pos)
{
    positions.push_back(pos);
}

void GLParticles::compile()
{
    scale(0.04);
    int cnt = positions.size();

    v_all = new float[f_cnt * 9 * cnt];
    vt_all = new float[f_cnt * 6 * cnt];
    vn_all = new float[f_cnt * 9 * cnt];
    f_all = new unsigned int[f_cnt * 9 * cnt];

    qDebug() << "good 1" << cnt << f_cnt * 9 * cnt;

    for (int ccc = 0; ccc < cnt; ++ccc)
    {
        int base_9 = f_cnt * 9 * ccc;
        int base_6 = f_cnt * 6 * ccc;

        for (int i = 0; i < f_cnt * 9; i += 3)
        {
            v_all[base_9 + i] = v[i] + positions[ccc].x();
            v_all[base_9 + i + 1] = v[i + 1] + positions[ccc].y();
            v_all[base_9 + i + 2] = v[i + 2] + positions[ccc].z();

            vn_all[base_9 + i] = vn[i];
            vn_all[base_9 + i + 1] = vn[i + 1];
            vn_all[base_9 + i + 2] = vn[i + 2];
        }

        for (int i = 0; i < f_cnt * 6; ++i)
        {
            vt_all[base_6 + i] = vt[i];
        }
    }

    for (int i = 0; i < f_cnt * 9 * cnt; ++i) f_all[i] = i;

    qDebug() << "start";

    vao.create();
    vao.bind();

    qDebug() << "create succ";

    glGenBuffers(1, &vbo_v);
    glGenBuffers(1, &vbo_vt);
    glGenBuffers(1, &vbo_vn);
    glGenBuffers(1, &ebo_f);

    glBindBuffer(GL_ARRAY_BUFFER, vbo_v);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * f_cnt * 9 * cnt, v_all, GL_STREAM_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, vbo_vt);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * f_cnt * 6 * cnt, vt_all, GL_STATIC_DRAW);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, vbo_vn);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * f_cnt * 9 * cnt, vn_all, GL_STATIC_DRAW);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(2);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_f);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * f_cnt * 9 * cnt, f_all, GL_STATIC_DRAW);

    vao.release();

    // qDebug() << "release succ";
}

void GLParticles::render()
{
    vao.bind();
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDrawElements(GL_TRIANGLES, f_cnt * 3 * positions.size(), GL_UNSIGNED_INT, 0);
    vao.release();
}

void GLParticles::myupdate(int _id, const QVector3D& pos)
{
    int base_9 = f_cnt * 9 * _id;
    for (int i = 0; i < f_cnt * 9; i += 3)
    {
        v_all[base_9 + i] = v[i] + pos.x();
        v_all[base_9 + i + 1] = v[i + 1] + pos.y();
        v_all[base_9 + i + 2] = v[i + 2] + pos.z();
    }
}

void GLParticles::hack()
{
    vao.bind();
    glBindBuffer(GL_ARRAY_BUFFER, vbo_v);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * f_cnt * 9 * positions.size(), v_all, GL_STREAM_DRAW);
    vao.release();
}
