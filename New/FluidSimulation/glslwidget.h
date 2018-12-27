#ifndef GLSLWIDGET_H
#define GLSLWIDGET_H

#include "openglwidget.h"

#include <QWidget>

namespace Ui {
class GLSLWidget;
}

class GLSLWidget : public QWidget
{
    Q_OBJECT

public:
    explicit GLSLWidget(QWidget *parent = 0);
    ~GLSLWidget();

public:
    OpenGLWidget* widget;

private:
    Ui::GLSLWidget *ui;
};

#endif // GLSLWIDGET_H
