#include "openglwidget.h"

#include "logic/constants.h"

#include <QOpenGLBuffer>
#include <QTimer>
#include <QMatrix4x4>
#include <QDebug>
#include <QtMath>
#include <algorithm>
#include <cmath>

OpenGLWidget::OpenGLWidget(QWidget *parent): QOpenGLWidget(parent)
{
    X_STEP = Y_STEP = Z_STEP = 0.025f;
    cameraPos = QVector3D(0.0f, 0.0f, 4.0f);
    cameraFront = QVector3D(0.0f, 0.0f, -1.0f);
    cameraUp = QVector3D(0.0f, 1.0f, 0.0f);

    yaw = -90.0f;
    pitch = 0.0f;
    cur_key = -1;
    isChanging = false;

    setFocusPolicy(Qt::StrongFocus);
}

void OpenGLWidget::initializeGL()
{
    initializeOpenGLFunctions();

    glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
    glEnable(GL_DEPTH_TEST);

    glShadeModel(GL_SMOOTH);

    glEnable(GL_DOUBLEBUFFER);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glClearDepth(1.0);

    QOpenGLShader *vshader = new QOpenGLShader(QOpenGLShader::Vertex, this);
    vshader->compileSourceFile("./scene/wall.vert");

    QOpenGLShader *fshader = new QOpenGLShader(QOpenGLShader::Fragment, this);
    fshader->compileSourceFile("./scene/wall.frag");

    wallProgram = new QOpenGLShaderProgram();
    wallProgram->addShader(vshader);
    wallProgram->addShader(fshader);
    wallProgram->link();

    QOpenGLShader *vshader2 = new QOpenGLShader(QOpenGLShader::Vertex, this);
    vshader2->compileSourceFile("./scene/particle.vert");

    QOpenGLShader *fshader2 = new QOpenGLShader(QOpenGLShader::Fragment, this);
    fshader2->compileSourceFile("./scene/particle.frag");

    particleProgram = new QOpenGLShaderProgram();
    particleProgram->addShader(vshader2);
    particleProgram->addShader(fshader2);
    particleProgram->link();

    qDebug() << "good1";
    simulator = new Simulator();
    /*for (int i = 0; i < 100; ++i)
    {
        qDebug() << i << ss->particles[0].position << ss->particles[0].velocity;
        ss->simulate();
    }
    qDebug() << ss->particles[0].position << ss->particles[0].velocity;
    qDebug() << "good3";*/

    walls = new GLWalls(X_MIN, X_MAX, Y_MIN, Y_MAX, Z_MIN, Z_MAX);
    walls->compile();

    glparticles = new GLParticles();
    for (Particle& p: simulator->particles)
    {
        glparticles->add(p.predPosition);
    }
    glparticles->compile();

    model.setToIdentity();

    updateView();

    projection.setToIdentity();

    QTimer* timer = new QTimer(this);
    QObject::connect(timer, SIGNAL(timeout()), this, SLOT(moveCamera()));
    timer->start(10);

    QTimer* timerSimulator = new QTimer(this);
    QObject::connect(timerSimulator, SIGNAL(timeout()), this, SLOT(getSimulate()));
    timerSimulator->start(10);
}

void OpenGLWidget::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h);

    projection.setToIdentity();
    projection.perspective(60.0f, (GLfloat)w / (GLfloat)h, 0.001f, 10.0f);
}

void OpenGLWidget::paintGL()
{
    glShadeModel(GL_SMOOTH);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    model.setToIdentity();

    updateView();

    QMatrix4x4 MVP = projection * view * model;

    wallProgram->setUniformValue("MVP", MVP);
    wallProgram->setUniformValue("Eye", cameraPos);

    particleProgram->setUniformValue("MVP", MVP);
    particleProgram->setUniformValue("Kd", QVector3D(0.0, 0.4, 0.8));
    particleProgram->setUniformValue("Ld", QVector3D(1.0, 1.0, 1.0));
    particleProgram->setUniformValue("LightPosition", QVector3D(0.0, 0.0, 100.0));

    particleProgram->bind();
    glparticles->render();
    // particleProgram->release();

    // wallProgram->bind();
    walls->render();
    // wallProgram->release();

    // wallProgram->release();


}

void OpenGLWidget::mousePressEvent(QMouseEvent *event)
{
    if (event->buttons() == Qt::LeftButton)
    {
        mousePos = QVector2D(event->pos());
        isChanging = true;
    }
}

void OpenGLWidget::mouseReleaseEvent(QMouseEvent* event)
{
    if (event->buttons() == Qt::LeftButton)
    {
        isChanging = false;
    }
}

void OpenGLWidget::mouseMoveEvent(QMouseEvent *event)
{
    if (!isChanging) return;

    QVector2D mousePosDiff = QVector2D(event->pos()) - mousePos;
    yaw += mousePosDiff.x() / 8.0;
    pitch -= mousePosDiff.y() / 8.0;

    if (pitch > 89.0f)
    {
        pitch = 89.0f;
    }
    else if (pitch < -89.0f)
    {
        pitch = -89.0f;
    }

    cameraFront = QVector3D(qCos(qDegreesToRadians(yaw)) * qCos(qDegreesToRadians(pitch)),
                            qSin(qDegreesToRadians(pitch)),
                            qSin(qDegreesToRadians(yaw)) * qCos(qDegreesToRadians(pitch))).normalized();

    mousePos = QVector2D(event->pos());
    this->update();
}

void OpenGLWidget::keyPressEvent(QKeyEvent* event)
{
    int k = event->key();
    cur_key = k;
}

void OpenGLWidget::keyReleaseEvent(QKeyEvent *event)
{
    int k = event->key();
    if (cur_key == k)
    {
        cur_key = -1;
    }
}

void OpenGLWidget::updateView()
{
    view.setToIdentity();
    view.lookAt(cameraPos, cameraFront + cameraPos, cameraUp);
    this->update();
}

void OpenGLWidget::moveCamera()
{
    if (cur_key == -1)
    {
        return;
    }

    QVector3D cameraPosNew;
    if (cur_key == Qt::Key_W)
    {
        cameraPosNew = cameraPos + cameraFront * X_STEP;
    }
    else if (cur_key == Qt::Key_S)
    {
        cameraPosNew = cameraPos - cameraFront * X_STEP;
    }
    else if (cur_key == Qt::Key_A)
    {
        cameraPosNew = cameraPos - QVector3D::crossProduct(cameraFront, cameraUp).normalized() * X_STEP;
    }
    else if (cur_key == Qt::Key_D)
    {
        cameraPosNew = cameraPos + QVector3D::crossProduct(cameraFront, cameraUp).normalized() * X_STEP;
    }
    else if (cur_key == Qt::Key_Q)
    {
        cameraPosNew = cameraPos + cameraUp * X_STEP;
    }
    else if (cur_key == Qt::Key_E)
    {
        cameraPosNew = cameraPos - cameraUp * X_STEP;
    }
    else
    {
        return;
    }
    cameraPos = cameraPosNew;
    this->update();
}

void OpenGLWidget::getSimulate()
{
    // For delay animation
    record += 1;
    if (record < 200) return;

    simulator->simulate();
    for (int i = 0; i < simulator->particles.size(); ++i)
    {
        glparticles->myupdate(i, simulator->particles[i].predPosition);
    }
    glparticles->hack();
    this->update();
}
