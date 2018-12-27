#include "glslwidget.h"
#include "ui_glslwidget.h"

#include <QGridLayout>

GLSLWidget::GLSLWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::GLSLWidget)
{
    ui->setupUi(this);
    this->setFixedSize(800, 800);

    QGridLayout* layout = new QGridLayout();
    widget = new OpenGLWidget();
    layout->addWidget(widget);
    this->setLayout(layout);
}

GLSLWidget::~GLSLWidget()
{
    delete ui;
}
