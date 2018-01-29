TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    CPoint.cpp \
    CVector.cpp \
    TP_OPENGL.cpp

HEADERS += \
    CPoint.h \
    CVector.h

LIBS += \
    -lglut \
    -lGL \
    -lGLU \
