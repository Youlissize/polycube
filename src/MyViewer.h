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

#include <set>

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

    std::map< unsigned int , double > const & getRow( unsigned int r ) const { return _ASparse[r]; }

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

    float alphaInit = 0.05f;
    float alpha = alphaInit;
    float beta = 0.02f; // for shape complexity
    float h = 0.001f; // only for gradient descent
    float epsilonInit = 0.6f;
    float epsilon = epsilonInit;

    bool decrease_epsilon = true;

    bool increase_alpha = false;
    unsigned int step_by_alpha_increase = 50;

    unsigned int rotation_opti_by_step = 1;
    unsigned int total_steps = 4;

    bool useNormalAlignment = true;
    bool use_triangle_area_constraints = true;

    bool useStepLimiter = true;
    float maxStepDistance = 2.f;

    bool postProcessWork = true;

    std::vector< mat33d > tetrahedron_rotation_matrix;

    double meshLength;

    QWidget * controls;

    std::vector<int> orientation = std::vector<int>();
    bool displayOrientation = false; //Shows axis-aligned normals patches found
    bool isOrientationComputed = false;

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
        DetailedAction * showPatches = new DetailedAction( QIcon(fullpath+"/icons/work.png") , "showPatches" , "showPatches" , this , this , SLOT(showPatches()) );


        // Add them :
        toolBar->addAction( open_mesh );
        toolBar->addAction( save_mesh );
        toolBar->addAction( help );
        toolBar->addAction( saveCamera );
        toolBar->addAction( openCamera );
        toolBar->addAction( saveSnapShotPlusPlus );
        toolBar->addAction( work );
        toolBar->addAction( showPatches );

    }


    void draw() {
        glEnable(GL_DEPTH_TEST);
        //glEnable( GL_LIGHTING );
        //glColor3f(0.5,0.5,0.8);
        glDisable (GL_LIGHTING);

        glBegin(GL_TRIANGLES);
        for( unsigned int t = 0 ; t < mesh.triangles.size() ; ++t ) {
            point3d const & p0 = mesh.vertices[ mesh.triangles[t][0] ].p;
            point3d const & p1 = mesh.vertices[ mesh.triangles[t][1] ].p;
            point3d const & p2 = mesh.vertices[ mesh.triangles[t][2] ].p;
            point3d const & n = - point3d::cross( p1-p0 , p2-p0 ).direction();
            //glNormal3f(n[0],n[1],n[2]);

            if(displayOrientation){
                if(isOrientationComputed){
                    if(abs(orientation[t])==1){
                        glColor3f(1,0,0);
                    }
                    else if(abs(orientation[t])==2){
                        glColor3f(0,1,0);
                    }
                    else if(abs(orientation[t])==3){
                        glColor3f(0,0,1);
                    }
                }
                else{glColor3f(fabs(n[0]),fabs(n[1]),fabs(n[2]));}
            }
            else{glColor3f(fabs(n[0]),fabs(n[1]),fabs(n[2]));}
            //glColor3f(fabs(n[0]),fabs(n[1]),fabs(n[2]));
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
                point3d bb(FLT_MAX,FLT_MAX,FLT_MAX) , BB(-FLT_MAX,-FLT_MAX,-FLT_MAX);
                for( unsigned int v = 0 ; v < mesh.vertices.size() ; ++v ) {
                    mesh.vertices[v].pInit = mesh.vertices[v].p;
                    bb = point3d::min(bb , mesh.vertices[v]);
                    BB = point3d::max(BB , mesh.vertices[v]);
                }
                meshLength = (BB-bb).norm();
                // Adjust mesh size
                for( unsigned int v = 0 ; v < mesh.vertices.size() ; ++v ) {
                    mesh.vertices[v].pInit *= (500/meshLength);
                    mesh.vertices[v].p *= (500/meshLength);
                }
                bb *= (500/meshLength);
                BB *= (500/meshLength);
                meshLength = 500; //500 is the size of the hand model, which works

                adjustCamera(bb,BB);

                tetrahedron_rotation_matrix.resize( mesh.tetras.size() );
                for( unsigned int t = 0 ; t < mesh.tetras.size() ; ++t ) {
                    tetrahedron_rotation_matrix[ t ].setIdentity();
                }

                update();
                std::cout << "Opened " << fileName.toStdString() << " that contains " << mesh.tetras.size() << " tetrahedra, " << mesh.triangles.size()
                          << " boundary triangles, and " << mesh.vertices.size() << " vertices." << std::endl;
                std::cout << "Bounding box  :  " << bb << "    ----   "  << BB << std::endl << std::endl;
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

    bool isInTriangle(unsigned int i, Triangle T){
        for (int j=0; j<3; j++) {
            if (T.corners[j]==i) {
                return true;
            }
        }
        return false;
    }

    void work(){

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



        for(unsigned int alphait = 0; alphait < total_steps; alphait++) {  // Main loop
            std::cout << alphait*100/total_steps << "%" << std::endl;
            std::cout << "Alpha = "<<alpha<<" , Epsilon = "<<epsilon<<std::endl;
            if(increase_alpha && alphait != 0 && alphait % step_by_alpha_increase == 0) {
                alpha = alpha * 2;
                /*if (epsilon > 0.01f)
                    epsilon = epsilon / 2;
                if (epsilon < 0.01f)
                    epsilon = 0.01f;*/
            }

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

            Eigen::VectorXd pb(3*mesh.vertices.size());


            for(unsigned int rotationIt = 0; rotationIt < rotation_opti_by_step; ++rotationIt) {
                Eigen::VectorXd gradient(3*mesh.vertices.size());

                // Initialize Values
                for( unsigned int v = 0 ; v < mesh.vertices.size() ; ++v ) {
                    point3d p = mesh.vertices[ v ].p;
                    pb(3*v) = p[0];
                    pb(3*v+1) = p[1];
                    pb(3*v+2) = p[2];
                    gradient[3*v] = 0.0;
                    gradient[3*v+1] = 0.0;
                    gradient[3*v+2] = 0.0;
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
                            //                            for(unsigned int it = 0; it < 1; ++it) { // THIS GIVES WORSE RESULTS
                            //                                for(unsigned int jt = it+1; jt < mesh.tetras[t].size(); ++jt) { // THIS GIVES WORSE RESULTS
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

                // fix the first vertex:
                {
                    A_mine(0,0) += 1;
                    A_mine(1,1) += 1;
                    A_mine(2,2) += 1;
                    b[0] += mesh.vertices[0].pInit[0];
                    b[1] += mesh.vertices[0].pInit[1];
                    b[2] += mesh.vertices[0].pInit[2];
                }

                A_mine.convertToEigenFormat(A);
                gradient = 2*A*pb - 2*b;


                std::cout << "\t finished computing the ARAP Hessian and gradient" << std::endl;

                Eigen::SparseMatrix<double> polycubeHessianSparse(3*mesh.vertices.size(), 3*mesh.vertices.size());
                MySparseMatrix polycubeHessian( 3*mesh.vertices.size() , 3*mesh.vertices.size() );
                for( unsigned int t = 0 ; t < mesh.triangles.size() ; ++t ) {
                    int i0 = mesh.triangles[t][0];
                    int i1 = mesh.triangles[t][1];
                    int i2 = mesh.triangles[t][2];
                    std::vector<int> indexes = {3*i0, 3*i0+1, 3*i0+2, 3*i1, 3*i1+1, 3*i1+2, 3*i2, 3*i2+1, 3*i2+2}; // x0, y0, z0 , x1, y1, z1, x2, y2, z2
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
                    grad_c2(0) = p1[1] - p2[1]; //y1-y2
                    grad_c2(1) = p2[0] - p1[0]; //x2-x1
                    grad_c2(3) = p2[1] - p0[1]; //y2-y0
                    grad_c2(4) = p0[0] - p2[0]; //x0-x2
                    grad_c2(6) = p0[1] - p1[1]; //y0-y2
                    grad_c2(7) = p1[0] - p0[0]; //x1-x0

                    grad_c0(1) = p1[2] - p2[2]; //z1-z2
                    grad_c0(2) = p2[1] - p1[1]; //y2-y1
                    grad_c0(4) = p2[2] - p0[2]; //z2-z0
                    grad_c0(5) = p0[1] - p2[1]; //y0-y2
                    grad_c0(7) = p0[2] - p1[2]; //z0-z1
                    grad_c0(8) = p1[1] - p0[1]; //y1-y0

                    grad_c1(0) = p2[2] - p1[2]; //z2-z1
                    grad_c1(2) = p1[0] - p2[0]; //x1-x2
                    grad_c1(3) = p0[2] - p2[2]; //z0-z2
                    grad_c1(5) = p2[0] - p0[0]; //x2-x0
                    grad_c1(6) = p1[2] - p0[2]; //z1-z0
                    grad_c1(8) = p0[0] - p1[0]; //x0-x1

                    for(int i = 0; i < 9; ++i) {
                        gradient(indexes[i]) += alpha * meshLength * (c0/c0b * grad_c0(i) + c1/c1b * grad_c1(i) + c2/c2b * grad_c2(i));
                    }

                    Eigen::MatrixXd H_small = Eigen::MatrixXd(9,9);
                    for(int i = 0; i < 9; ++i) {
                        for(int j = 0; j < 9; ++j) {
                            H_small(i,j) = 0.0;
                        }
                    }
                    //Compute c0, c1 and c2 hessian's first part
                    H_small += epsilon * (1/(c0b*c0b*c0b) * grad_c0 * grad_c0.transpose() + 1/(c1b*c1b*c1b) * grad_c1 * grad_c1.transpose() + 1/(c2b*c2b*c2b) * grad_c2 * grad_c2.transpose());
                    //Compute c0 hessian's second part
                    H_small(1,5) += c0/c0b;
                    H_small(5,1) += c0/c0b;
                    H_small(1,8) += -c0/c0b;
                    H_small(8,1) += -c0/c0b;

                    H_small(4,2) += -c0/c0b;
                    H_small(2,4) += -c0/c0b;
                    H_small(4,8) += c0/c0b;
                    H_small(8,4) += c0/c0b;

                    H_small(7,2) += c0/c0b;
                    H_small(2,7) += c0/c0b;
                    H_small(7,5) += -c0/c0b;
                    H_small(5,7) += -c0/c0b;

                    //Compute c1 hessian's second part
                    H_small(2,3) += c1/c1b;
                    H_small(3,2) += c1/c1b;
                    H_small(2,6) += -c1/c1b;
                    H_small(6,2) += -c1/c1b;

                    H_small(5,0) += -c1/c1b;
                    H_small(0,5) += -c1/c1b;
                    H_small(5,6) += c1/c1b;
                    H_small(6,5) += c1/c1b;

                    H_small(8,0) += c1/c1b;
                    H_small(0,8) += c1/c1b;
                    H_small(8,3) += -c1/c1b;
                    H_small(3,8) += -c1/c1b;

                    //Compute c2 hessian's second part
                    H_small(0,4) += c2/c2b;
                    H_small(4,0) += c2/c2b;
                    H_small(0,7) += -c2/c2b;
                    H_small(7,0) += -c2/c2b;

                    H_small(3,1) += -c2/c2b;
                    H_small(1,3) += -c2/c2b;
                    H_small(3,7) += c2/c2b;
                    H_small(7,3) += c2/c2b;

                    H_small(6,1) += c2/c2b;
                    H_small(1,6) += c2/c2b;
                    H_small(6,4) += -c2/c2b;
                    H_small(4,6) += -c2/c2b;

                    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H_small);
                    Eigen::VectorXd D = es.eigenvalues();
                    Eigen::MatrixXd Q = es.eigenvectors();
                    for(int i = 0; i < 9; i++) {
                        if(D(i) < 0)
                            D(i) = 0.0001;
                    }
                    H_small = Q * D.asDiagonal() * Q.transpose();

                   /* Eigen::JacobiSVD<Eigen::MatrixXd> svd(H_small,  Eigen::ComputeFullU | Eigen::ComputeFullV);
                    Eigen::MatrixXd U = svd.matrixU();
                    Eigen::MatrixXd V = svd.matrixV();
                    Eigen::VectorXd S = svd.singularValues();
                    for(int i = 0; i < 9; i++) {
                        if(S(i) < 0)
                            S(i) = 0.000001;
                    }
                    H_small = U * S.asDiagonal() * U.transpose();*/

                    for(int i = 0; i < 9; ++i) {
                        for(int j = 0; j < 9; ++j) {
                            polycubeHessian(indexes[i],indexes[j]) += meshLength * H_small(i,j);
                        }
                    }
                }

                std::cout << "\t finished computing the normal L1 norm Hessian and gradient" << std::endl;

                //Eigen::SparseMatrix<double> alignementHessianSparse(3*mesh.vertices.size(), 3*mesh.vertices.size());
                MySparseMatrix alignementHessian( 3*mesh.vertices.size() , 3*mesh.vertices.size() );
                if (useNormalAlignment){
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
                    ind=0;
                    for (auto it=edges.begin(); it !=edges.end(); it++) {
                        int i = it->second.first;
                        int j = it->second.second;
                        float airei = area(i);
                        float airej=area(j);
                        float Aeij = edgeWeights[ind];

                        unsigned int indA,indB,indC,indD;
                        Triangle Ti = mesh.triangles[i];
                        Triangle Tj = mesh.triangles[j];
                        int temp;
                        for(int k=0; k<3; k++) {
                            if(not(isInTriangle(Ti.corners[k], Tj))){
                                indC=Ti.corners[k];
                                temp=k;
                            }
                        }
                        indA = Ti.corners[((temp+1)%3)];
                        indB = Ti.corners[((temp+2)%3)];
                        for(int k=0; k<3; k++) {
                            if(not(isInTriangle(Tj.corners[k], Ti))){
                                indD=Ti.corners[k];
                            }
                        }
                        point3d pA = mesh.vertices[indA];
                        point3d pB = mesh.vertices[indB];
                        point3d pC = mesh.vertices[indC];
                        point3d pD = mesh.vertices[indD];


                        point3d ni = point3d::cross(pB-pA , pC-pA);
                        point3d nj = point3d::cross(pD-pA , pB-pA);

                        gradient[indA*3+0] += beta*2*(Aeij/totArea)* (
                                      (ni.y()-nj.y())*( (pC.z()-pB.z())/airei +  (pD.z()-pB.z())/airej)
                                    + (ni.z()-nj.z())*( (-pC.y()+pB.y())/airei +  (-pD.y()+pB.y())/airej)   );
                        gradient[indA*3+1] += beta*2*(Aeij/totArea)* (
                                      (ni.x()-nj.x())*( (-pC.z()+pB.z())/airei +  (-pD.z()+pB.z())/airej)
                                    + (ni.z()-nj.z())*( (pC.x()-pB.x())/airei +  (pD.x()-pB.x())/airej)   );
                        gradient[indA*3+2] += beta*2*(Aeij/totArea)* (
                                      (ni.x()-nj.x())*( (pC.y()-pB.y())/airei +  (pD.y()-pB.y())/airej)
                                    + (ni.y()-nj.y())*( (-pC.x()+pB.x())/airei +  (-pD.x()+pB.x())/airej)   );

                        gradient[indB*3+0] += -beta*2*(Aeij/totArea)* (
                                      (ni.y()-nj.y())*( (pC.z()-pA.z())/airei +  (pD.z()-pA.z())/airej)
                                    + (ni.z()-nj.z())*( (-pC.y()+pA.y())/airei +  (-pD.y()+pA.y())/airej)   );
                        gradient[indB*3+1] += -beta*2*(Aeij/totArea)* (
                                      (ni.x()-nj.x())*( (-pC.z()+pA.z())/airei +  (-pD.z()+pA.z())/airej)
                                    + (ni.z()-nj.z())*( (pC.x()-pA.x())/airei +  (pD.x()-pA.x())/airej)   );
                        gradient[indB*3+2] += -beta*2*(Aeij/totArea)* (
                                      (ni.x()-nj.x())*( (pC.y()-pA.y())/airei +  (pD.y()-pA.y())/airej)
                                    + (ni.y()-nj.y())*( (-pC.x()+pA.x())/airei +  (-pD.x()+pA.x())/airej)   );

                        gradient[indC*3+0] += beta*2*(Aeij/totArea)* (
                                      (ni.y()-nj.y())*( (pB.z()-pA.z())/airei )
                                    + (ni.z()-nj.z())*( (pA.y()-pB.y())/airei )   );
                        gradient[indC*3+1] += beta*2*(Aeij/totArea)* (
                                      (ni.x()-nj.x())*( (pA.z()-pB.z())/airei )
                                    + (ni.z()-nj.z())*( (pB.x()-pA.x())/airei )   );
                        gradient[indC*3+2] += beta*2*(Aeij/totArea)* (
                                      (ni.x()-nj.x())*( (pB.y()-pA.y())/airei )
                                    + (ni.y()-nj.y())*( (pA.x()-pB.x())/airei )   );

                        gradient[indD*3+0] += beta*2*(Aeij/totArea)* (
                                      (ni.y()-nj.y())*( (pB.z()-pA.z())/airej )
                                    + (ni.z()-nj.z())*( (pA.y()-pB.y())/airej )   );
                        gradient[indD*3+1] += beta*2*(Aeij/totArea)* (
                                      (ni.x()-nj.x())*( (pA.z()-pB.z())/airej )
                                    + (ni.z()-nj.z())*( (pB.x()-pA.x())/airej )   );
                        gradient[indD*3+2] += beta*2*(Aeij/totArea)* (
                                      (ni.x()-nj.x())*( (pB.y()-pA.y())/airej )
                                    + (ni.y()-nj.y())*( (pA.x()-pB.x())/airej )   );



                        ind ++;
                    }


                    for (unsigned int i=0; i<3*mesh.vertices.size(); i++) {
                        alignementHessian(i,i) += 1;
                    }

                    std::cout << "\t finished computing the normal alignement gradient" << std::endl;
                }

                if(false){  //Gradient Descent
                    pb = pb - h*gradient;
                }
                else //or Newton Descent
                {
                    if( use_triangle_area_constraints ) {
                        MySparseMatrix LHSsystem_mine( 3*mesh.vertices.size() + mesh.triangles.size() , 3*mesh.vertices.size() + mesh.triangles.size() );
                        Eigen::VectorXd RHSsystem( 3*mesh.vertices.size() + mesh.triangles.size() );
                        for( unsigned int coord = 0 ; coord < 3*mesh.vertices.size() ; ++coord ) RHSsystem[coord] = - gradient[coord];

                        // set energy Hessian to 2*A + alpha*polycubeHessian + beta*alignementHessian:
                        for( unsigned int r = 0 ; r < 3*mesh.vertices.size() ; ++r ) {
                            std::map< unsigned int , double > const & row_r = A_mine.getRow(r);
                            for( std::map< unsigned int , double >::const_iterator it = row_r.begin() ; it != row_r.end() ; ++it ) {
                                unsigned int c = it->first;
                                double val_rc = it->second;
                                LHSsystem_mine( r , c ) += 2 * val_rc;
                            }
                        }
                        for( unsigned int r = 0 ; r < 3*mesh.vertices.size() ; ++r ) {
                            std::map< unsigned int , double > const & row_r = polycubeHessian.getRow(r);
                            for( std::map< unsigned int , double >::const_iterator it = row_r.begin() ; it != row_r.end() ; ++it ) {
                                unsigned int c = it->first;
                                double val_rc = it->second;
                                LHSsystem_mine( r , c ) += alpha * val_rc;
                            }
                        }
                        for( unsigned int r = 0 ; r < 3*mesh.vertices.size() ; ++r ) {
                            std::map< unsigned int , double > const & row_r = alignementHessian.getRow(r);
                            for( std::map< unsigned int , double >::const_iterator it = row_r.begin() ; it != row_r.end() ; ++it ) {
                                unsigned int c = it->first;
                                double val_rc = it->second;
                                LHSsystem_mine( r , c ) += beta * val_rc;
                            }
                        }


                        // compute Jacobian of the constraints, and set appropriate values:
                        for( unsigned int t = 0 ; t < mesh.triangles.size() ; ++t ) {
                            int i0 = mesh.triangles[t][0];
                            int i1 = mesh.triangles[t][1];
                            int i2 = mesh.triangles[t][2];
                            std::vector<int> indexes = {3*i0, 3*i0+1, 3*i0+2, 3*i1, 3*i1+1, 3*i1+2, 3*i2, 3*i2+1, 3*i2+2}; // x0, y0, z0 , x1, y1, z1, x2, y2, z2
                            point3d p0 = point3d(pb(3*i0), pb(3*i0+1), pb(3*i0+2));
                            point3d p1 = point3d(pb(3*i1), pb(3*i1+1), pb(3*i1+2));
                            point3d p2 = point3d(pb(3*i2), pb(3*i2+1), pb(3*i2+2));
                            point3d n = point3d::cross( p1-p0 , p2-p0 );
                            point3d nOriginal = point3d::cross( mesh.vertices[i1].pInit-mesh.vertices[i0].pInit , mesh.vertices[i2].pInit-mesh.vertices[i0].pInit );
                            double c0 = n[0];
                            double c1 = n[1];
                            double c2 = n[2];

                            Eigen::VectorXd grad_c0(9);
                            Eigen::VectorXd grad_c1(9);
                            Eigen::VectorXd grad_c2(9);
                            for( unsigned int tgrad = 0 ; tgrad < 9 ; ++tgrad ) {
                                grad_c0(tgrad) = 0;
                                grad_c1(tgrad) = 0;
                                grad_c2(tgrad) = 0;
                            }
                            grad_c2(0) = p1[1] - p2[1]; //y1-y2
                            grad_c2(1) = p2[0] - p1[0]; //x2-x1
                            grad_c2(3) = p2[1] - p0[1]; //y2-y0
                            grad_c2(4) = p0[0] - p2[0]; //x0-x2
                            grad_c2(6) = p0[1] - p1[1]; //y0-y2
                            grad_c2(7) = p1[0] - p0[0]; //x1-x0

                            grad_c0(1) = p1[2] - p2[2]; //z1-z2
                            grad_c0(2) = p2[1] - p1[1]; //y2-y1
                            grad_c0(4) = p2[2] - p0[2]; //z2-z0
                            grad_c0(5) = p0[1] - p2[1]; //y0-y2
                            grad_c0(7) = p0[2] - p1[2]; //z0-z1
                            grad_c0(8) = p1[1] - p0[1]; //y1-y0

                            grad_c1(0) = p2[2] - p1[2]; //z2-z1
                            grad_c1(2) = p1[0] - p2[0]; //x1-x2
                            grad_c1(3) = p0[2] - p2[2]; //z0-z2
                            grad_c1(5) = p2[0] - p0[0]; //x2-x0
                            grad_c1(6) = p1[2] - p0[2]; //z1-z0
                            grad_c1(8) = p0[0] - p1[0]; //x0-x1

                            for(int i = 0; i < 9; ++i) {
                                LHSsystem_mine( 3*mesh.vertices.size() + t , indexes[i] ) += 2*(c0 * grad_c0(i) + c1 * grad_c1(i) + c2 * grad_c2(i));
                                LHSsystem_mine( indexes[i] , 3*mesh.vertices.size() + t ) += 2*(c0 * grad_c0(i) + c1 * grad_c1(i) + c2 * grad_c2(i));
                            }

                            RHSsystem[ 3*mesh.vertices.size() + t ] = - ( n.sqrnorm() - nOriginal.sqrnorm() );
                        }

                        // add an epsilon regularization on the diagonal:
                        {
                            double epsilonReg = 0.000000001;
                            for( unsigned int r = 0 ; r < 3*mesh.vertices.size() + mesh.triangles.size() ; ++r )
                                LHSsystem_mine( r , r ) += epsilonReg;
                        }

                        std::cout << "\t finished computing the LHS and RHS" << std::endl;

                        Eigen::SparseMatrix< double > LHSsystem;
                        LHSsystem_mine.convertToEigenFormat( LHSsystem );

                        std::cout << "\t finished converting the LHS to Eigen's format" << std::endl;

                        //                        Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > cg;
                        //                        std::cout << "\t\t cg.compute(LHSsystem)" << std::endl;
                        //                        cg.compute(LHSsystem);
                        //                        std::cout << "\t\t cg.compute(LHSsystem) OK" << std::endl;
                        //                        std::cout << "\t\t deltaX = cg.solve(RHSsystem)" << std::endl;
                        //                        Eigen::VectorXd deltaX = cg.solve(RHSsystem);
                        //                        std::cout << "\t\t deltaX = cg.solve(RHSsystem) OK" << std::endl;


                        Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > _LHS_LDLT;
                        _LHS_LDLT.analyzePattern( LHSsystem );
                        _LHS_LDLT.compute( LHSsystem );
                        Eigen::VectorXd deltaX = _LHS_LDLT.solve( RHSsystem );

                        std::cout << "\t finished solving the system" << std::endl;

                        for( unsigned int coord = 0 ; coord < 3 * mesh.vertices.size() ; ++coord ) {
                            pb[coord] += deltaX[coord];
                        }
                        std::cout << "\t finished update" << std::endl << std::endl;
                    }
                    else {
                        polycubeHessian.convertToEigenFormat(polycubeHessianSparse);
                        Eigen::SparseMatrix<double> Hessian = 2*A + alpha * polycubeHessianSparse;
                        Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > _Hessian_LDLT;

                        _Hessian_LDLT.analyzePattern( Hessian );
                        _Hessian_LDLT.compute( Hessian );
                        pb = pb - _Hessian_LDLT.solve( gradient );
                    }
                }

                //Optimize Rotation matrices
                for( unsigned int t = 0 ; t < mesh.tetras.size() ; ++t ) {
                    Eigen::Matrix3d S = Eigen::Matrix3d(3,3);
                    for( int i = 0 ; i < 3 ; ++i )
                        for( int j = 0 ; j < 3 ; ++j )
                            S(i,j) = 0.0;

                    //                    for(unsigned int it = 0; it < 1; ++it) {
                    //                        for(unsigned int jt = it+1; jt < mesh.tetras[t].size(); ++jt) {
                    for(unsigned int it = 0; it < mesh.tetras[t].size(); ++it) {
                        for(unsigned int jt = 0; jt < mesh.tetras[t].size(); ++jt) {
                            if( it != jt ) {
                                int i = mesh.tetras[t][it];
                                int j = mesh.tetras[t][jt];
                                Eigen::Vector3d ei_current (pb(3*i)-pb(3*j), pb(3*i+1)-pb(3*j+1), pb(3*i+2)-pb(3*j+2));
                                //   point3d pei = tetrahedron_rotation_matrix[ t ] * mesh.vertices[ i ].pInit - tetrahedron_rotation_matrix[ t ] * mesh.vertices[ j ].pInit;
                                point3d pei = mesh.vertices[ i ].pInit - mesh.vertices[ j ].pInit;
                                Eigen::Vector3d ei_init (pei[0], pei[1], pei[2]);
                                S += ei_current * ei_init.transpose();
                            }
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
            for( unsigned int v = 0 ; v < mesh.vertices.size() ; ++v ) {
                point3d & p = mesh.vertices[ v ].p;
                point3d target = point3d(pb(3*v), pb(3*v+1), pb(3*v+2));
                double distance =(target-p).norm();
                if(useStepLimiter && distance > maxStepDistance) {
                    target = maxStepDistance * (target-p) / distance + p;
                }
                p[0] = target[0];
                p[1] = target[1];
                p[2] = target[2];
            }
        }

        if( increase_alpha )
            alpha *= 1.1;
        if (decrease_epsilon){
            epsilon *= 0.95;
        }

        // Clean the result and extract an hexaedrale mesh
        if(postProcessWork) {

            // Compute orientation of each triangle
            orientation = std::vector<int>(mesh.triangles.size(),0);
            for (unsigned int t=0; t<mesh.triangles.size(); t++) {
                int i0 = mesh.triangles[t][0];
                int i1 = mesh.triangles[t][1];
                int i2 = mesh.triangles[t][2];
                std::vector<int> indexes = {3*i0, 3*i0+1, 3*i0+2, 3*i1, 3*i1+1, 3*i1+2, 3*i2, 3*i2+1, 3*i2+2}; // x0, y0, z0 , x1, y1, z1, x2, y2, z2
                point3d p0 = mesh.vertices[i0].p;
                point3d p1 = mesh.vertices[i1].p;
                point3d p2 = mesh.vertices[i2].p;
                point3d n = point3d::cross( p1-p0 , p2-p0 );
                n.normalize();
                orientation[t] = closestOrientation(n);
                // 1,2,3 = X,Y,Z , et + ou - selon la direction
            }

            // Check orientation validity and fill triangleNeighboors
            int nb_oppositeTriangles =0;
            std::vector<std::vector<unsigned int>> triangleNeighboors = std::vector<std::vector<unsigned int>>(mesh.triangles.size(),std::vector<unsigned int>());
            for (auto it=edges.begin(); it !=edges.end(); it++) {

                int ti = it->second.first;
                int tj = it->second.second;
                if(orientation[ti]+orientation[tj]==0){
                    nb_oppositeTriangles ++;
                }
                triangleNeighboors[ti].push_back(tj);
                triangleNeighboors[tj].push_back(ti);
            }
            std::cout<<"\tNumbers of triangles with opposite direction found : "<<nb_oppositeTriangles<<std::endl;

            // Relabelling some triangles
            for (unsigned int t=0; t<mesh.triangles.size(); t++) {
                if(triangleNeighboors[t].size()==3){
                    int o1 = orientation[triangleNeighboors[t][0]];
                    int o2 = orientation[triangleNeighboors[t][1]];
                    int o3 = orientation[triangleNeighboors[t][2]];
                    if(o1==o2 && orientation[t]==o3){
                        orientation[t] = o1;
                    }
                    if(o1==o3 && orientation[t]==o2){
                        orientation[t] = o1;
                    }
                    if(o3==o2 && orientation[t]==o1){
                        orientation[t] = o3;
                    }
                }
            }
            /*
            //Find Patches, and their neighboors
            std::vector<std::vector<unsigned int>> patches = std::vector<std::vector<unsigned int>>();
            std::vector<std::set<int>> patchesNieghboors = std::vector<std::set<int>>(); //only store the neighboor's orientation
            std::vector<bool> visitedTriangles = std::vector<bool>(mesh.triangles.size(),false);
            unsigned int ind = 0;
            int nb_patches = 0;
            while(ind<mesh.triangles.size()) {
                if (visitedTriangles[ind]){
                    ind++;
                }
                else{ // we have found a triangle from another patch
                    nb_patches ++;
                    int current_orientation = orientation[ind];
                    patches.push_back(std::vector<unsigned int>());
                    patches[nb_patches].push_back(ind);

                    //BFS to find all patch triangles
                    std::vector<int> stack = std::vector<int>();
                    stack.push_back(ind);
                    while(stack.size()>0) {
                        int t = stack[stack.size()-1];
                        stack.pop_back();

                        if(not(visitedTriangles[t])){
                            visitedTriangles[t]=true;
                            for(int k=0; k<triangleNeighboors[t].size();k ++){
                                int neigh = triangleNeighboors[t][k];
                                if(orientation[neigh]==current_orientation){
                                    stack.push_back(neigh);
                                    patches[nb_patches].push_back(neigh);
                                }
                                else{ //we found an edge
                                    patchesNieghboors[nb_patches].insert(orientation[neigh]);
                                }
                            }
                        }
                    } //BFS end
                    //nb_patches++;

                }

            }

            //Clean Patches issues
            for(unsigned int p=0; p<patches.size(); p++) {
                if(patchesNieghboors[p].size()==1){
                    int new_orientation = *patchesNieghboors[p].begin();
                    for(unsigned int i=0; i<patches[p].size(); i++){
                        orientation[patches[p][i]]= new_orientation;
                    }
                }
            }
            for(unsigned int p=0; p<patches.size(); p++) {
                if(patchesNieghboors[p].size()==2){
                    // should take the longuest border to choose the new orientation
                    int new_orientation = *patchesNieghboors[p].begin();
                    for(unsigned int i=0; i<patches[p].size(); i++){
                        orientation[patches[p][i]]= new_orientation;
                    }
                }
            }
            std::cout<<"\tNumber of patches found : "<<nb_patches<<std::endl;
            */
            isOrientationComputed = true;


        }

        update();
    }

    void showPatches(){
        displayOrientation = not(displayOrientation);
        update();
    }

    int closestOrientation(point3d n){
        float nx = fabs(n.x());
        float ny = fabs(n.y());
        float nz = fabs(n.z());

        if(nx>ny){
            if(nx>nz){
                if(n.x()>0){
                    return 1;
                }
                else{
                    return -1;
                }
            }
            else{
                if(n.z()>0){
                    return 3;
                }
                else{
                    return -3;
                }
            }
        }
        else if (ny>nz){
            if(n.y()>0){
                return 2;
            }
            else{
                return -2;
            }
        }
        else {
            if(n.z()>0){
                return 3;
            }
            else{
                return -3;
            }
        }
    }

};

#endif // MYVIEWER_H
