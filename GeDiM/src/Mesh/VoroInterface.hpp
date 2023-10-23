#ifndef __VoroInterface_H
#define __VoroInterface_H

#include "GeometryUtilities.hpp"

#include "IMeshDAO.hpp"
using namespace std;

#if ENABLE_VORO == 1
#include "voro++.hh"
#endif

namespace Gedim
{
class VoroInterface final
{
public:
    struct Cell0D
    {
        const double x;
        const double y;
        const double z;

        unsigned int id;
        unsigned int marker;

        Cell0D(const double x,
               const double y,
               const double z): x(x), y(y), z(z) {};
    };

    struct Cell2D
    {
        unsigned int id;
        int marker;
        vector<unsigned int> vertices;
        vector<unsigned int> edges;
    };

private:

    const Gedim::GeometryUtilities& geometryUtilities;

public:
    VoroInterface(const Gedim::GeometryUtilities& geometryUtilities);

    void GenerateVoronoiTassellations3D(const Eigen::MatrixXd& polyhedronVertices,
                                        const Eigen::MatrixXi& polyhedronEdges,
                                        const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                        const unsigned int &numPoints,
                                        const unsigned int& numIterations,
                                        Gedim::IMeshDAO& mesh);

    void GenerateVoronoiTassellations2D(const Eigen::MatrixXd& polygonVertices,
                                        const unsigned int& numPoints,
                                        const unsigned int& numIterations,
                                        Gedim::IMeshDAO& mesh);

private:

#if ENABLE_VORO == 1
    bool InsertNewPoints(Cell0D& cell0D, list<Cell0D>& cell0Ds);
    inline double rnd() {return double(rand())/RAND_MAX;}

    void GenerateRandomPoints(const Eigen::MatrixXd& domainVertices,
                              const unsigned int& numPoints,
                              voro::container& con);

    void GenerateCartesianPoints3D(const Eigen::MatrixXd &polyhedronVertices,
                                   const unsigned int &numPoints,
                                   voro::container &con);

    void GenerateRandomPoints(const Eigen::MatrixXd& domainVertices,
                              const unsigned int& numPoints,
                              Eigen::MatrixXd& VoronoiPoints);
#endif
};

}

#endif

