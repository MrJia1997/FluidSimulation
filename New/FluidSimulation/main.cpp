#include "glslwidget.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    GLSLWidget w;
    w.show();

    return a.exec();
}
