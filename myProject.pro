TEMPLATE = app
TARGET = myProject
DEPENDPATH += .
INCLUDEPATH += .
INCLUDEPATH += ./src
QT += opengl xml
CONFIG += qt \
    release
MOC_DIR = ./tmp/moc
OBJECTS_DIR = ./tmp/obj

# Input
HEADERS += ./src/point3.h \
    ./src/Mesh.h \
    ./src/MyViewer.h \
    ./src/gl/GLUtilityMethods.h \
    ./src/gl/BasicColors.h
SOURCES += ./src/main.cpp \
    ./src/gl/GLUtilityMethods.cpp\
    ./src/gl/BasicColors.cpp



EXT_DIR = /home/gabriel/Documents/IGR/IGR205/basicQGLProject/basicQGLProject/extern



#QGLViewer
{
 INCLUDEPATH += $${EXT_DIR}/libQGLViewer-2.6.1
 LIBS +=    -L$${EXT_DIR}/libQGLViewer-2.6.1/QGLViewer/build-QGLViewerLibs \
            -lQGLViewer
}



LIBS += -lglut \
    -lGLU
LIBS += -lgsl \
    -lgomp
LIBS += -lgslcblas \
    -lgomp
#LIBS += -lblas \
#    -lgomp
release:QMAKE_CXXFLAGS_RELEASE += -O3 \
    -fopenmp
release:QMAKE_CFLAGS_RELEASE += -O3 \
    -fopenmp
