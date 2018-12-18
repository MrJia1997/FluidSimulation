#include "fluidsimulation.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    FluidSimulation w;
    w.show();
    return a.exec();
}
