#ifndef OPENGLWIDGET_H
#define OPENGLWIDGET_H

#include "glwalls.h"
#include "glparticles.h"
#include "glsurface.h"
#include "logic/simulator.h"

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLShader>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLTexture>
#include <QMouseEvent>
#include <QKeyEvent>

class OpenGLWidget: public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT
public:
    OpenGLWidget(QWidget* parent = 0);

public:
    QOpenGLShaderProgram* wallProgram;
    QOpenGLShaderProgram* particleProgram;
    QOpenGLShaderProgram* surfaceProgram;

    GLWalls* walls;
    GLParticles* glparticles;
    GLSurface* glsurface;
    Simulator* simulator;


    QMatrix4x4 model;
    QMatrix4x4 view;
    QMatrix4x4 projection;

    QVector2D mousePos;

    float X_STEP;
    float Y_STEP;
    float Z_STEP;

    QVector3D cameraPos;
    QVector3D cameraFront;
    QVector3D cameraUp;

    float yaw;
    float pitch;

    int cur_key;
    bool isChanging;
    int record = 0;


protected:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();

    void mousePressEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
    void mouseMoveEvent(QMouseEvent* event);

    void keyPressEvent(QKeyEvent* event);
    void keyReleaseEvent(QKeyEvent* event);

    void updateView();

public slots:
    void moveCamera();
    void getSimulate();
};

#endif
