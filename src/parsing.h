#ifndef PARSING_H
#define PARSING_H

#endif // PARSING_H


#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include"Mesh.h"




namespace MeshIO{
template< class point_t > bool open( const std::string & filename , std::vector< point_t > & vertices , std::vector< int > & tetras )
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    std::string dummy_s;

    unsigned int meshVersion , dataDimension , nVertices , nTriangles , nTetras;

    myfile >> dummy_s >> meshVersion;
    myfile >> dummy_s >> dataDimension;

    if( dataDimension != 3 )
    {
        std::cout << "dataDimension != 3 , " << filename << " is not a tetrahedral mesh file !" << std::endl;
        myfile.close();
        return false;
    }

    myfile >> dummy_s >> nVertices;
    vertices.resize( nVertices );
    for( unsigned int v = 0 ; v < nVertices ; ++v )
    {
        double x , y , z , dummy_i;
        myfile >> x >> y >> z >> dummy_i;
        vertices[v] = point_t(x,y,z);
    }

    myfile >> dummy_s >> nTriangles;
    for( unsigned int v = 0 ; v < nTriangles ; ++v )
    {
        int x , y , z , dummy_i;
        myfile >> x >> y >> z >> dummy_i;
    }

    myfile >> dummy_s >> nTetras;
    tetras.resize(4 * nTetras );
    for( unsigned int t = 0 ; t < nTetras ; ++t )
    {
        int w , x , y , z , dummy_i;
        myfile >> w >> x >> y >> z >> dummy_i;
        tetras[ 4*t ] = w-1;
        tetras[ 4*t + 1 ] = x-1;
        tetras[ 4*t + 2 ] = y-1;
        tetras[ 4*t + 3 ] = z-1;
    }
    myfile.close();

    return true;
}


template< class point_t , class int_t > bool openTrisAndTets( const std::string & filename , std::vector< point_t > & vertices , std::vector< Triangle > & tris , std::vector< std::vector< int_t > > & tetras ,
                                                              bool parseBoundaryTrisOnly = true , bool reverseTriangleOrientation = false )
{

    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    std::string dummy_s;

    unsigned int meshVersion , dataDimension , nVertices , nTriangles , nTetras;

    myfile >> dummy_s >> meshVersion;
    myfile >> dummy_s >> dataDimension;

    if( dataDimension != 3 )
    {
        std::cout << "dataDimension != 3 , " << filename << " is not a tetrahedral mesh file !" << std::endl;
        myfile.close();
        return false;
    }

    myfile >> dummy_s >> nVertices;
    vertices.resize( nVertices );
    for( unsigned int v = 0 ; v < nVertices ; ++v )
    {
        double x , y , z , dummy_i;
        myfile >> x >> y >> z >> dummy_i;
        vertices[v] = point_t(x,y,z);
    }
    std::cout << "\t parsed vertices" << std::endl;

    myfile >> dummy_s >> nTriangles;
    tris.clear();
    if(!parseBoundaryTrisOnly) tris.resize(nTriangles );
    for( unsigned int t = 0 ; t < nTriangles ; ++t )
    {
        int x , y , z , dummy_i;
        if(reverseTriangleOrientation)
            myfile >> x >> z >> y >> dummy_i;
            else
        myfile >> x >> y >> z >> dummy_i;
        if(!parseBoundaryTrisOnly) {
            //tris[t].resize(3);
            tris[t].corners[ 0 ] = x-1;
            tris[t].corners[ 1 ] = y-1;
            tris[t].corners[ 2 ] = z-1;
        }
        else {
            if( dummy_i != 0 ) {
                Triangle newTri;
                newTri.corners[ 0 ] = x-1;
                newTri.corners[ 1 ] = y-1;
                newTri.corners[ 2 ] = z-1;
                tris.push_back(newTri);
            }
        }
    }
    std::cout << "\t parsed tris : " << tris.size() << " tris out of " << nTriangles << " in the file" << std::endl;

    myfile >> dummy_s >> nTetras;
    tetras.resize(nTetras );
    for( unsigned int t = 0 ; t < nTetras ; ++t )
    {
        int w , x , y , z , dummy_i;
        myfile >> w >> x >> y >> z >> dummy_i;
        tetras[ t ].resize(4);
        tetras[ t ][0] = w-1;
        tetras[ t ][1] = x-1;
        tetras[ t ][2] = y-1;
        tetras[ t ][3] = z-1;
    }
    std::cout << "\t parsed tets" << std::endl;
    myfile.close();

    return true;
}


template< class point_t , class int_t > bool openTris( const std::string & filename , std::vector< point_t > & vertices , std::vector< std::vector< int_t > > & tris )
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    std::string dummy_s;

    unsigned int meshVersion , dataDimension , nVertices , nTriangles;

    myfile >> dummy_s >> meshVersion;
    myfile >> dummy_s >> dataDimension;

    if( dataDimension != 3 )
    {
        std::cout << "dataDimension != 3 , " << filename << " is not a tetrahedral mesh file !" << std::endl;
        myfile.close();
        return false;
    }

    myfile >> dummy_s >> nVertices;
    std::cout << "allocate " << nVertices << " vertices" << std::endl;
    vertices.resize( nVertices );
    for( unsigned int v = 0 ; v < nVertices ; ++v )
    {
        double x , y , z , dummy_i;
        myfile >> x >> y >> z >> dummy_i;
        vertices[v] = point_t(x,y,z);
    }

    myfile >> dummy_s >> nTriangles;
    std::cout << "allocate " << nTriangles << " tris" << std::endl;
    tris.resize( nTriangles );
    for( unsigned int t = 0 ; t < nTriangles ; ++t )
    {
        tris[t].resize(3);
        int_t x , y , z , dummy_i;
        myfile >> x >> y >> z >> dummy_i;
        tris[ t ][ 0 ] = x-1;
        tris[ t ][ 1 ] = y-1;
        tris[ t ][ 2 ] = z-1;
    }

    myfile.close();

    return true;
}

}
