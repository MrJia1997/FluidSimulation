#-------------------------------------------------
#
# Project created by QtCreator 2018-12-25T19:57:59
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = FluidSimulation
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.

DEFINES += QT_DEPRECATED_WARNINGS
CONFIG += c++11

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


SOURCES += \
        main.cpp \
        glslwidget.cpp \
    openglwidget.cpp \
    logic/geometry.cpp \
    logic/simulator.cpp \
    glslwidget.cpp \
    main.cpp \
    openglwidget.cpp \
    glwalls.cpp \
    glparticles.cpp

HEADERS += \
        glslwidget.h \
    openglwidget.h \
    logic/constants.h \
    logic/geometry.h \
    logic/kernel.h \
    logic/particle.h \
    logic/simulator.h \
    glslwidget.h \
    openglwidget.h \
    glwalls.h \
    glparticles.h

FORMS += \
        glslwidget.ui

LIBS += -lopengl32 -lGLU32
INCLUDEPATH += C:/dev/Eigen-3.3.5

DISTFILES += \
    scene/wall.frag \
    scene/wall.vert \
    scene/particle.frag \
    scene/particle.vert
