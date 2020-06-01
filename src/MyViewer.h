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
#include "extern/eigen3/Eigen/Dense"
#include "extern/eigen3/Eigen/SVD"

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
                    mesh.vertices[v].pInit = mesh.vertices[v].p;
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

    float area(int triangle) {
        int i0 = mesh.triangles[triangle][0];
        int i1 = mesh.triangles[triangle][1];
        int i2 = mesh.triangles[triangle][2];
        point3d p0 = mesh.vertices[i0].pInit;
        point3d p1 = mesh.vertices[i1].pInit;
        point3d p2 = mesh.vertices[i2].pInit;
        point3d n = point3d::cross( p1-p0 , p2-p0 );
        return n.norm()*0.5;
    }

    void work(){
        float alpha = 0.5f;
        float beta = 0.5f; // for shape complexity
        float h = 0.001f; // only for gradient descent
        float epsilon = 0.6f;

        // Calcul once the mesh total surface
        float totArea = 0;
        for (unsigned int t=0; t<mesh.triangles.size(); t++){
            totArea += area(t);
        }

        // Initialise Edges
        // link edges with their 2 adjacent Triangles
        // Edges = map ( (point1,point2) , (triangle1,triangle2) )
        std::map<std::pair<int,int>,std::pair<int,int>> edges = std::map<std::pair<int,int>,std::pair<int,int>>();
        for (unsigned int t=0; t<mesh.triangles.size(); t++){
            for (int i=0; i<3; i++){
                int tmp1 = mesh.triangles[t][i];
                int tmp2 = mesh.triangles[t][((i+1)%3)];
                int p1 = std::min(tmp1,tmp2);
                int p2 = std::max(tmp1,tmp2);

                std::pair<int,int> e = std::pair<int,int>(p1,p2);
                auto it = edges.find(e);
                if( it == edges.end() ) { // edge (p1,p2) never seen
                    edges.insert(std::pair<std::pair<int,int>,std::pair<int,int>>(e, std::pair<int,int>(t,-1)));
                }
                else {
                    it->second.second = t;
                }
            }
        }



        for(int alphait = 0; alphait < 1; alphait++) {  // Main loop
        std::cout << alphait*100/20 << "%" << std::endl;
        /*if(alphait != 0 && alphait % 30 == 0) {
            alpha = alpha * 2;
            if (epsilon > 0.01f)
                epsilon = epsilon / 2;
            if (epsilon < 0.01f)
                epsilon = 0.01f;
        }*/

        double energy = 0;

        for( unsigned int t = 0 ; t < mesh.triangles.size() ; ++t ) {
            int i0 = mesh.triangles[t][0];
            int i1 = mesh.triangles[t][1];
            int i2 = mesh.triangles[t][2];
            point3d p0 = mesh.vertices[i0].p;
            point3d p1 = mesh.vertices[i1].p;
            point3d p2 = mesh.vertices[i2].p;
            point3d n = point3d::cross( p1-p0 , p2-p0 );
            energy += fabs(n[0]) + fabs(n[1]) + fabs(n[2]);
        }
        std::cout << "Energy : " << energy << std::endl;

        std::vector< mat33d > tetrahedron_rotation_matrix( mesh.tetras.size() );
        for( unsigned int t = 0 ; t < mesh.tetras.size() ; ++t ) {
            tetrahedron_rotation_matrix[ t ].setIdentity();
        }
        Eigen::VectorXd pb(3*mesh.vertices.size());


        for(unsigned int rotationIt = 0; rotationIt < 4; ++rotationIt) {
            Eigen::VectorXd gradient(3*mesh.vertices.size());

            // Initialize Values
            for( unsigned int t = 0 ; t < mesh.vertices.size() ; ++t ) {
                point3d p = mesh.vertices[ t ].p;
                pb(3*t) = p[0];
                pb(3*t+1) = p[1];
                pb(3*t+2) = p[2];
                gradient[3*t] = 0.0;
                gradient[3*t+1] = 0.0;
                gradient[3*t+2] = 0.0;
            }

            // Compute ARAP Gradient
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

                            point3d pi = tetrahedron_rotation_matrix[ t ] * mesh.vertices[ i ].pInit; //not current position, but initial position
                            point3d pj = tetrahedron_rotation_matrix[ t ] * mesh.vertices[ j ].pInit; //not current position, but initial position
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
            A_mine(0,0) += 1;
            A_mine(1,1) += 1;
            A_mine(2,2) += 1;

            A_mine.convertToEigenFormat(A);
            gradient = 2*A*pb - 2*b;

            Eigen::SparseMatrix<double> polycubeHessianSparce(3*mesh.vertices.size(), 3*mesh.vertices.size());
            MySparseMatrix polycubeHessian( 3*mesh.vertices.size() , 3*mesh.vertices.size() );
            for( unsigned int t = 0 ; t < mesh.triangles.size() ; ++t ) {
                int i0 = mesh.triangles[t][0];
                int i1 = mesh.triangles[t][1];
                int i2 = mesh.triangles[t][2];
                std::vector<int> indexes = {3*i0, 3*i0+1, 3*i0+2, 3*i1, 3*i1+1, 3*i1+2, 3*i2, 3*i2+1, 3*i2+1}; // x0, y0, z0 , x1, y1, z1, x2, y2, z2
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

                Eigen::VectorXd grad_c0(9);
                Eigen::VectorXd grad_c1(9);
                Eigen::VectorXd grad_c2(9);
                for( unsigned int tgrad = 0 ; tgrad < 9 ; ++tgrad ) {
                    grad_c0(tgrad) = 0;
                    grad_c1(tgrad) = 0;
                    grad_c2(tgrad) = 0;
                }
                grad_c0(0) = p1[1] - p2[1]; //y1-y2
                grad_c0(1) = p2[1] - p0[1]; //y2-y0
                grad_c0(2) = p0[1] - p1[1]; //y0-y2
                grad_c0(3) = p2[0] - p1[0]; //x2-x1
                grad_c0(4) = p0[0] - p2[0]; //x0-x2
                grad_c0(5) = p1[0] - p0[0]; //x1-x0

                grad_c1(3) = p1[2] - p2[2]; //z1-z2
                grad_c1(4) = p2[2] - p0[2]; //z2-z0
                grad_c1(5) = p0[2] - p1[2]; //z0-z1
                grad_c1(6) = p2[1] - p1[1]; //y2-y1
                grad_c1(7) = p0[1] - p2[1]; //y0-y2
                grad_c1(8) = p1[1] - p0[1]; //y1-y0

                grad_c2(0) = p2[2] - p1[2]; //z2-z1
                grad_c2(1) = p0[2] - p2[2]; //z0-z2
                grad_c2(2) = p1[2] - p0[2]; //z1-z0
                grad_c2(6) = p1[0] - p2[0]; //x1-x2
                grad_c2(7) = p2[0] - p0[0]; //x2-x0
                grad_c2(8) = p0[0] - p1[0]; //x0-x1

                for(int i = 0; i < 9; ++i) {
                    gradient(indexes[i]) += alpha * (c0/c0b * grad_c0(i) + c1/c1b * grad_c1(i) + c2/c2b * grad_c2(i));
                }

                Eigen::MatrixXd H_small = Eigen::MatrixXd(9,9);
                //Compute c0, c1 and c2 hessian's first part
                H_small += epsilon * (1/(c0b*c0b*c0b) * grad_c0 * grad_c0.transpose() + 1/(c1b*c1b*c1b) * grad_c1 * grad_c1.transpose() + 1/(c2b*c2b*c2b) * grad_c2 * grad_c2.transpose());
                //Compute c0 hessian's second part
                H_small(0,4) += c0/c0b;
                H_small(4,0) += c0/c0b;
                H_small(0,7) += -c0/c0b;
                H_small(7,0) += -c0/c0b;

                H_small(3,1) += -c0/c0b;
                H_small(1,3) += -c0/c0b;
                H_small(3,7) += c0/c0b;
                H_small(7,3) += c0/c0b;

                H_small(6,1) += c0/c0b;
                H_small(1,6) += c0/c0b;
                H_small(6,4) += -c0/c0b;
                H_small(4,6) += -c0/c0b;

                //Compute c1 hessian's second part
                H_small(1,5) += c1/c1b;
                H_small(5,1) += c1/c1b;
                H_small(1,8) += -c1/c1b;
                H_small(8,1) += -c1/c1b;

                H_small(4,2) += -c1/c1b;
                H_small(2,4) += -c1/c1b;
                H_small(4,8) += c1/c1b;
                H_small(8,4) += c1/c1b;

                H_small(7,3) += c1/c1b;
                H_small(3,7) += c1/c1b;
                H_small(7,5) += -c1/c1b;
                H_small(5,7) += -c1/c1b;

                //Compute c2 hessian's second part
                H_small(2,3) += c2/c2b;
                H_small(3,2) += c2/c2b;
                H_small(2,6) += -c2/c2b;
                H_small(6,2) += -c2/c2b;

                H_small(5,0) += -c2/c2b;
                H_small(0,5) += -c2/c2b;
                H_small(5,6) += c2/c2b;
                H_small(6,5) += c2/c2b;

                H_small(8,0) += c2/c2b;
                H_small(0,8) += c2/c2b;
                H_small(8,3) += -c2/c2b;
                H_small(3,8) += -c2/c2b;
                for(int i = 0; i < 9; ++i) {
                    for(int j = 0; j < 9; ++j) {
                        polycubeHessian(indexes[i],indexes[j]) += H_small(i,j);
                    }
                }
            }


            // Compute the Aeij (for the normal alignment gradient)
            Eigen::VectorXd edgeWeights(edges.size());  // A(eij,X)
            Eigen::VectorXd neighboorsArea(mesh.triangles.size());
            for (unsigned int i=0; i< mesh.triangles.size() ; i++ ){
                neighboorsArea[i]=0;
            }
            for (auto it=edges.begin(); it !=edges.end(); it++) {
                int i = it->second.first;
                int j = it->second.second;
                neighboorsArea[i] += area(j);
                neighboorsArea[j] += area(i);
            }
            int ind=0;
            for (auto it=edges.begin(); it !=edges.end(); it++) {
                int i = it->second.first;
                int j = it->second.second;

                float gammaij = area(j)/neighboorsArea(i);
                float gammaji = area(i)/neighboorsArea(j);
                edgeWeights[ind]= gammaij*area(i) + gammaji*area(j);
                ind ++;
            }
            // Compute normal alignment gradient
            //TODO


            if(false){  //Gradient Descent
                pb = pb - h*gradient;
            }
            else //or Newton Descent
            {
                polycubeHessian.convertToEigenFormat(polycubeHessianSparce);
                Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > _Hessian_LDLT;

                _Hessian_LDLT.analyzePattern( 2*A + alpha * polycubeHessianSparce);
                _Hessian_LDLT.compute( 2*A + alpha * polycubeHessianSparce);
                pb = pb - _Hessian_LDLT.solve( gradient );
            }

            //Optimize Rotation matrices
            for( unsigned int t = 0 ; t < mesh.tetras.size() ; ++t ) {
                Eigen::Matrix3d S = Eigen::Matrix3d(3,3);
                for(unsigned int it = 0; it < mesh.tetras[t].size(); ++it) {
                    for(unsigned int jt = it+1; jt < mesh.tetras[t].size(); ++jt) {
                        int i = mesh.tetras[t][it];
                        int j = mesh.tetras[t][jt];
                        Eigen::Vector3d eip (pb(3*i)-pb(3*j), pb(3*i+1)-pb(3*j+1), pb(3*i+2)-pb(3*j+2));
                        point3d pei = tetrahedron_rotation_matrix[ t ] * mesh.vertices[ i ].pInit - tetrahedron_rotation_matrix[ t ] * mesh.vertices[ j ].pInit;
                        Eigen::Vector3d ei (pei[0], pei[1], pei[2]);
                        S += eip * ei.transpose();
                    }
                }
                Eigen::JacobiSVD<Eigen::Matrix3d> svd(S,  Eigen::ComputeFullU | Eigen::ComputeFullV);
                Eigen::Matrix3d U = svd.matrixU();
                Eigen::Matrix3d V = svd.matrixV();
                Eigen::Matrix3d product = U*V.transpose();
                Eigen::Matrix3d R = U * Eigen::Vector3d(1, 1, product.determinant()).asDiagonal() * V.transpose();
                for(unsigned int i = 0; i < 3; i++) {
                    for(unsigned int j = 0; j < 3; j++) {
                        tetrahedron_rotation_matrix[ t ](i, j) = R(i, j);
                    }
                }
            }

        }

        // Update Positions
        for( unsigned int t = 0 ; t < mesh.vertices.size() ; ++t ) {
            point3d & p = mesh.vertices[ t ].p;
            p[0] = pb(3*t);
            p[1] = pb(3*t+1);
            p[2] = pb(3*t+2);
        }
        }
        update();
    }
};

#endif // MYVIEWER_H
