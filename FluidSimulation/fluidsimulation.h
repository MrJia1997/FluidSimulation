#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_fluidsimulation.h"

class FluidSimulation : public QMainWindow
{
    Q_OBJECT

public:
    FluidSimulation(QWidget *parent = Q_NULLPTR);

private:
    Ui::FluidSimulationClass ui;
};
