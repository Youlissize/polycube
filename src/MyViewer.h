#ifndef MYVIEWER_H
#define MYVIEWER_H

// Mesh stuff:
#include "Mesh.h"

// Parsing:
#include "BasicIO.h"
#include "parsing.h"

// opengl and basic gl utilities:
#define GL_GLEXT_PROTOTYPES
#include <gl/openglincludeQtComp.h>
#include <GL/glext.h>
#include <QOpenGLFunctions_3_0>
#include <QOpenGLFunctions>
#include <QGLViewer/qglviewer.h>

#include <gl/GLUtilityMethods.h>

// Qt stuff:
#include <QFormLayout>
#include <QToolBar>
#include <QColorDialog>
#include <QFileDialog>
#include <QKeyEvent>
#include <QInputDialog>
#include <QLineEdit>


#include "qt/QSmartAction.h"

// To get the right path
#include "fullpath.h"

// eigen
#include "extern/eigen3/Eigen/Sparse"

class MySparseMatrix {
    std::vector< std::map< unsigned int , double > > _ASparse;

    unsigned int _rows , _columns;

public:
    MySparseMatrix() {
        _rows = _columns = 0;
    }
    MySparseMatrix( int rows , int columns ) {
        setDimensions(rows , columns);
    }
    ~MySparseMatrix() {
    }

    void setDimensions( int rows , int columns ) {
        _rows = rows; _columns = columns;
        _ASparse.clear();
        _ASparse.resize(_rows);
    }

    double & operator() (unsigned int row , unsigned int column) {
        return _ASparse[row][column];
    }

    void convertToEigenFormat(Eigen::SparseMatrix<double> & _A) {
        // convert ad-hoc matrix to Eigen sparse format:
        {
            _A.resize(_rows , _columns);
            std::vector< Eigen::Triplet< double > > triplets;
            for( unsigned int r = 0 ; r < _rows ; ++r ) {
                for( std::map< unsigned int , double >::const_iterator it = _ASparse[r].begin() ; it != _ASparse[r].end() ; ++it ) {
                    unsigned int c = it->first;
                    double val = it->second;
                    triplets.push_back( Eigen::Triplet< double >(r,c,val) );
                }
            }
            _A.setFromTriplets( triplets.begin() , triplets.end() );
        }
    }
};

class MyViewer : public QGLViewer , public QOpenGLFunctions_3_0
{
    Q_OBJECT

    Mesh mesh;

    QWidget * controls;

public :

    MyViewer(QGLWidget * parent = NULL) : QGLViewer(parent) , QOpenGLFunctions_3_0() {
    }



    void add_actions_to_toolBar(QToolBar *toolBar)
    {
        // Specify the actions :
        DetailedAction * open_mesh = new DetailedAction( QIcon(fullpath+"/icons/open.png") , "Open Mesh" , "Open Mesh" , this , this , SLOT(open_mesh()) );
        DetailedAction * save_mesh = new DetailedAction( QIcon(fullpath+"/icons/save.png") , "Save model" , "Save model" , this , this , SLOT(save_mesh()) );
        DetailedAction * help = new DetailedAction( QIcon(fullpath+"/icons/help.png") , "HELP" , "HELP" , this , this , SLOT(help()) );
        DetailedAction * saveCamera = new DetailedAction( QIcon(fullpath+"/icons/camera.png") , "Save camera" , "Save camera" , this , this , SLOT(saveCamera()) );
        DetailedAction * openCamera = new DetailedAction( QIcon(fullpath+"/icons/open_camera.png") , "Open camera" , "Open camera" , this , this , SLOT(openCamera()) );
        DetailedAction * saveSnapShotPlusPlus = new DetailedAction( QIcon(fullpath+"/icons/save_snapshot.png") , "Save snapshot" , "Save snapshot" , this , this , SLOT(saveSnapShotPlusPlus()) );
        DetailedAction * work = new DetailedAction( QIcon(fullpath+"/icons/work.png") , "Work" , "Work" , this , this , SLOT(work()) );

        // Add them :
        toolBar->addAction( open_mesh );
        toolBar->addAction( save_mesh );
        toolBar->addAction( help );
        toolBar->addAction( saveCamera );
        toolBar->addAction( openCamera );
        toolBar->addAction( saveSnapShotPlusPlus );
        toolBar->addAction( work );
    }


    void draw() {
        glEnable(GL_DEPTH_TEST);
        glEnable( GL_LIGHTING );
        glColor3f(0.5,0.5,0.8);
        glBegin(GL_TRIANGLES);
        for( unsigned int t = 0 ; t < mesh.triangles.size() ; ++t ) {
            point3d const & p0 = mesh.vertices[ mesh.triangles[t][0] ].p;
            point3d const & p1 = mesh.vertices[ mesh.triangles[t][1] ].p;
            point3d const & p2 = mesh.vertices[ mesh.triangles[t][2] ].p;
            point3d const & n = point3d::cross( p1-p0 , p2-p0 ).direction();
            glNormal3f(n[0],n[1],n[2]);
            glVertex3f(p0[0],p0[1],p0[2]);
            glVertex3f(p1[0],p1[1],p1[2]);
            glVertex3f(p2[0],p2[1],p2[2]);
        }
        glEnd();
    }

    void pickBackgroundColor() {
        QColor _bc = QColorDialog::getColor( this->backgroundColor(), this);
        if( _bc.isValid() ) {
            this->setBackgroundColor( _bc );
            this->update();
        }
    }

    void adjustCamera( point3d const & bb , point3d const & BB ) {
        point3d const & center = ( bb + BB )/2.f;
        setSceneCenter( qglviewer::Vec( center[0] , center[1] , center[2] ) );
        setSceneRadius( 1.5f * ( BB - bb ).norm() );
        showEntireScene();
    }


    void init() {
        makeCurrent();
        initializeOpenGLFunctions();

        setMouseTracking(true);// Needed for MouseGrabber.

        setBackgroundColor(QColor(255,255,255));

        // Lights:
        GLTools::initLights();
        GLTools::setSunsetLight();
        GLTools::setDefaultMaterial();

        //
        glShadeModel(GL_SMOOTH);
        glFrontFace(GL_CCW); // CCW ou CW

        glEnable(GL_DEPTH);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);

        glEnable(GL_CLIP_PLANE0);

        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        glEnable(GL_COLOR_MATERIAL);

        //
        setSceneCenter( qglviewer::Vec( 0 , 0 , 0 ) );
        setSceneRadius( 10.f );
        showEntireScene();
    }

    QString helpString() const {
        QString text("<h2>Our cool project</h2>");
        text += "<p>";
        text += "This is a research application, it can explode.";
        text += "<h3>Participants</h3>";
        text += "<ul>";
        text += "<li>jmt</li>";
        text += "<li>...</li>";
        text += "</ul>";
        text += "<h3>Basics</h3>";
        text += "<p>";
        text += "<ul>";
        text += "<li>H   :   make this help appear</li>";
        text += "<li>Ctrl + mouse right button double click   :   choose background color</li>";
        text += "<li>Ctrl + T   :   change window title</li>";
        text += "</ul>";
        return text;
    }

    void updateTitle( QString text ) {
        this->setWindowTitle( text );
        emit windowTitleUpdated(text);
    }

    void keyPressEvent( QKeyEvent * event ) {
        if( event->key() == Qt::Key_H ) {
            help();
        }
        else if( event->key() == Qt::Key_T ) {
            if( event->modifiers() & Qt::CTRL )
            {
                bool ok;
                QString text = QInputDialog::getText(this, tr(""), tr("title:"), QLineEdit::Normal,this->windowTitle(), &ok);
                if (ok && !text.isEmpty())
                {
                    updateTitle(text);
                }
            }
        }
    }

    void mouseDoubleClickEvent( QMouseEvent * e )
    {
        if( (e->modifiers() & Qt::ControlModifier)  &&  (e->button() == Qt::RightButton) )
        {
            pickBackgroundColor();
            return;
        }

        if( (e->modifiers() & Qt::ControlModifier)  &&  (e->button() == Qt::LeftButton) )
        {
            showControls();
            return;
        }

        QGLViewer::mouseDoubleClickEvent( e );
    }

    void mousePressEvent(QMouseEvent* e ) {
        QGLViewer::mousePressEvent(e);
    }

    void mouseMoveEvent(QMouseEvent* e  ){
        QGLViewer::mouseMoveEvent(e);
    }

    void mouseReleaseEvent(QMouseEvent* e  ) {
        QGLViewer::mouseReleaseEvent(e);
    }

signals:
    void windowTitleUpdated( const QString & );

public slots:
    void open_mesh() {
        bool success = false;
        QString fileName = QFileDialog::getOpenFileName(NULL,"","");
        if ( !fileName.isNull() ) { // got a file name
            if(fileName.endsWith(QString(".off")))
                success = OFFIO::openTriMesh(fileName.toStdString() , mesh.vertices , mesh.triangles );
            else if(fileName.endsWith(QString(".obj")))
                success = OBJIO::openTriMesh(fileName.toStdString() , mesh.vertices , mesh.triangles );
            else if(fileName.endsWith(QString(".mesh")))
                success =  MeshIO::openTrisAndTets(fileName.toStdString() , mesh.vertices , mesh.triangles, mesh.tetras );

            if(success) {
                std::cout << fileName.toStdString() << " was opened successfully" << std::endl;
                point3d bb(FLT_MAX,FLT_MAX,FLT_MAX) , BB(-FLT_MAX,-FLT_MAX,-FLT_MAX);
                for( unsigned int v = 0 ; v < mesh.vertices.size() ; ++v ) {
                    bb = point3d::min(bb , mesh.vertices[v]);
                    BB = point3d::max(BB , mesh.vertices[v]);
                }
                adjustCamera(bb,BB);
                update();
            }
            else
                std::cout << fileName.toStdString() << " could not be opened" << std::endl;
        }
    }

    void save_mesh() {
        bool success = false;
        QString fileName = QFileDialog::getOpenFileName(NULL,"","");
        if ( !fileName.isNull() ) { // got a file name
            if(fileName.endsWith(QString(".off")))
                success = OFFIO::save(fileName.toStdString() , mesh.vertices , mesh.triangles );
            else if(fileName.endsWith(QString(".obj")))
                success = OBJIO::save(fileName.toStdString() , mesh.vertices , mesh.triangles );
            if(success)
                std::cout << fileName.toStdString() << " was saved" << std::endl;
            else
                std::cout << fileName.toStdString() << " could not be saved" << std::endl;
        }
    }

    void showControls()
    {
        // Show controls :
        controls->close();
        controls->show();
    }

    void saveCameraInFile(const QString &filename){
        std::ofstream out (filename.toUtf8());
        if (!out)
            exit (EXIT_FAILURE);
        // << operator for point3 causes linking problem on windows
        out << camera()->position()[0] << " \t" << camera()->position()[1] << " \t" << camera()->position()[2] << " \t" " " <<
                                          camera()->viewDirection()[0] << " \t" << camera()->viewDirection()[1] << " \t" << camera()->viewDirection()[2] << " \t" << " " <<
                                          camera()->upVector()[0] << " \t" << camera()->upVector()[1] << " \t" <<camera()->upVector()[2] << " \t" <<" " <<
                                          camera()->fieldOfView();
        out << std::endl;

        out.close ();
    }

    void openCameraFromFile(const QString &filename){

        std::ifstream file;
        file.open(filename.toStdString().c_str());

        qglviewer::Vec pos;
        qglviewer::Vec view;
        qglviewer::Vec up;
        float fov;

        file >> (pos[0]) >> (pos[1]) >> (pos[2]) >>
                                                    (view[0]) >> (view[1]) >> (view[2]) >>
                                                                                           (up[0]) >> (up[1]) >> (up[2]) >>
                                                                                                                            fov;

        camera()->setPosition(pos);
        camera()->setViewDirection(view);
        camera()->setUpVector(up);
        camera()->setFieldOfView(fov);

        camera()->computeModelViewMatrix();
        camera()->computeProjectionMatrix();

        update();
    }


    void openCamera(){
        QString fileName = QFileDialog::getOpenFileName(NULL,"","*.cam");
        if ( !fileName.isNull() ) {                 // got a file name
            openCameraFromFile(fileName);
        }
    }
    void saveCamera(){
        QString fileName = QFileDialog::getSaveFileName(NULL,"","*.cam");
        if ( !fileName.isNull() ) {                 // got a file name
            saveCameraInFile(fileName);
        }
    }

    void saveSnapShotPlusPlus(){
        QString fileName = QFileDialog::getSaveFileName(NULL,"*.png","");
        if ( !fileName.isNull() ) {                 // got a file name
            setSnapshotFormat("PNG");
            setSnapshotQuality(100);
            saveSnapshot( fileName );
            saveCameraInFile( fileName+QString(".cam") );
        }
    }

    void work(){
        float alpha = 0.0f;
        const float epsilon = 0.2f;

        Eigen::VectorXd gradient(3*mesh.vertices.size());

        Eigen::VectorXd pb(3*mesh.vertices.size());
        for( unsigned int t = 0 ; t < mesh.vertices.size() ; ++t ) {
            point3d p = mesh.vertices[ t ].p;
            pb(3*t) = p[0];
            pb(3*t+1) = p[1];
            pb(3*t+2) = p[2];
            gradient[3*t] = 0.0;
            gradient[3*t+1] = 0.0;
            gradient[3*t+2] = 0.0;
        }

        Eigen::SparseMatrix<double> A(3*mesh.vertices.size(), 3*mesh.vertices.size());
        MySparseMatrix A_mine( 3*mesh.vertices.size() , 3*mesh.vertices.size() );
        Eigen::VectorXd b(3*mesh.vertices.size());
        for(unsigned int t = 0; t < 3*mesh.vertices.size(); ++t) {
            b[t] = 0.0;
        }
        for( unsigned int t = 0 ; t < mesh.tetras.size() ; ++t ) {
            for(unsigned int it = 0; it < mesh.tetras[t].size(); ++it) {
                for(unsigned int jt = 0; jt < mesh.tetras[t].size(); ++jt) {
                    if(it != jt) {
                        int i = mesh.tetras[t][it];
                        int j = mesh.tetras[t][jt];

                        A_mine(3*i, 3*i) += 1;
                        A_mine(3*i+1, 3*i+1) += 1;
                        A_mine(3*i+2, 3*i+2) += 1;

                        A_mine(3*j, 3*j) += 1;
                        A_mine(3*j+1, 3*j+1) += 1;
                        A_mine(3*j+2, 3*j+2) += 1;

                        A_mine(3*i, 3*j) += -1;
                        A_mine(3*i+1, 3*j+1) += -1;
                        A_mine(3*i+2, 3*j+2) += -1;

                        A_mine(3*j, 3*i) += -1;
                        A_mine(3*j+1, 3*i+1) += -1;
                        A_mine(3*j+2, 3*i+2) += -1;

                        point3d pi = mesh.vertices[ i ].p;
                        point3d pj = mesh.vertices[ j ].p;
                        b(3*i)   += pi[0]-pj[0];
                        b(3*i+1) += pi[1]-pj[1];
                        b(3*i+2) += pi[2]-pj[2];
                        b(3*j)   += pj[0]-pi[0];
                        b(3*j+1) += pj[1]-pi[1];
                        b(3*j+2) += pj[2]-pi[2];
                    }
                }
            }
        }

        A_mine.convertToEigenFormat(A);
        for (int k = 0; k < 1; k++) {
            gradient = 2*A*pb - 2*b;

            for( unsigned int t = 0 ; t < mesh.triangles.size() ; ++t ) {
                int i0 = mesh.triangles[t][0];
                int i1 = mesh.triangles[t][1];
                int i2 = mesh.triangles[t][2];
                point3d p0 = point3d(pb(3*i0), pb(3*i0+1), pb(3*i0+2));
                point3d p1 = point3d(pb(3*i1), pb(3*i1+1), pb(3*i1+2));
                point3d p2 = point3d(pb(3*i2), pb(3*i2+1), pb(3*i2+2));
                point3d n = point3d::cross( p1-p0 , p2-p0 );
                double c0 = n[0];
                double c1 = n[1];
                double c2 = n[2];
                double c0b = std::sqrt(c0*c0 + epsilon);
                double c1b = std::sqrt(c1*c1 + epsilon);
                double c2b = std::sqrt(c2*c2 + epsilon);

                gradient(3*i0) += alpha * (c0/c0b * (p1[1] - p2[1]) + c2/c2b * (p2[2] - p1[2]));
                gradient(3*i1) += alpha * (c0/c0b * (p2[1] - p0[1]) + c2/c2b * (p0[2] - p2[2]));
                gradient(3*i2) += alpha * (c0/c0b * (p0[1] - p1[1]) + c2/c2b * (p1[2] - p0[2]));

                gradient(3*i0+1) += alpha * (c0/c0b * (p2[0] - p1[0]) + c1/c1b * (p1[2] - p2[2]));
                gradient(3*i1+1) += alpha * (c0/c0b * (p0[0] - p2[0]) + c1/c1b * (p2[2] - p0[2]));
                gradient(3*i2+1) += alpha * (c0/c0b * (p1[0] - p0[0]) + c1/c1b * (p0[2] - p1[2]));

                gradient(3*i0+2) += alpha * (c1/c1b * (p2[1] - p1[1]) + c2/c2b * (p1[0] - p2[0]));
                gradient(3*i1+2) += alpha * (c1/c1b * (p0[1] - p2[1]) + c2/c2b * (p2[0] - p0[0]));
                gradient(3*i2+2) += alpha * (c1/c1b * (p1[1] - p0[1]) + c2/c2b * (p0[0] - p1[0]));
            }
            pb = pb - 0.0001*gradient;
        }
        for( unsigned int t = 0 ; t < mesh.vertices.size() ; ++t ) {
            point3d & p = mesh.vertices[ t ].p;
            if(abs(p[0]-pb(3*t)) > 0.000000001)
                std::cout << "notoke" << std::endl;
            if(abs(p[1]-pb(3*t+1)) > 0.000000001)
                std::cout << "notoke" << std::endl;
            if(abs(p[2]-pb(3*t+2)) > 0.000000001)
                std::cout << "notoke" << std::endl;
            p[0] = pb(3*t);
            p[1] = pb(3*t+1);
            p[2] = pb(3*t+2);
        }
        draw();
    }
};

#endif // MYVIEWER_H
