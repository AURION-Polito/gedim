// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#ifndef __VoroInterface_H
#define __VoroInterface_H

#include "Gedim_Macro.hpp"
#include "GeometryUtilities.hpp"

#include "IMeshDAO.hpp"

#if ENABLE_VORO == 1
#include "voro++.hh"
#endif

#define TEST_VORO

namespace Gedim
{
class VoroInterface final
{
  public:
    struct Cell0D
    {
        const Eigen::VectorXd coordinates;
        unsigned int marker;
        bool active = true;
        std::vector<unsigned int> neighbors_1D;
        std::vector<unsigned int> neighbors_2D;

        Cell0D(const Eigen::VectorXd &coordinates) : coordinates(coordinates){};
    };

    struct Cell1D
    {
        unsigned int marker;
        unsigned int origin;
        unsigned int end;
        bool active = true;
        std::vector<unsigned int> neighbors_2D;
    };

    struct Cell2D
    {
        unsigned int marker;
        std::vector<unsigned int> vertices;
        std::vector<unsigned int> edges;
        bool active = true;
        std::vector<int> neighbors_of_related_3D_cells; // only positive ones
    };

    struct Cell3D
    {
        unsigned int marker;
        std::vector<unsigned int> vertices;
        std::set<unsigned int> edges;
        std::vector<unsigned int> faces;
        std::vector<int> neighbors;
        bool active = true;
    };

  private:
    const Gedim::GeometryUtilities &geometryUtilities;

#if ENABLE_VORO == 1
    unsigned int InsertNewPoint(const Cell0D &cell0D, std::map<unsigned int, Cell0D> &cell0Ds);

    inline double rnd()
    {
        return double(rand()) / RAND_MAX;
    }

    void GenerateCartesianPoints3D(const Eigen::MatrixXd &polyhedronVertices, const unsigned int &numPoints, voro::container &con);
#endif
  public:
    VoroInterface(const Gedim::GeometryUtilities &geometryUtilities);

    Eigen::MatrixXd GenerateRandomPoints(const Eigen::MatrixXd &domainVertices,
                                         const unsigned int &numPoints,
                                         const unsigned int random_seed = static_cast<unsigned int>(time(nullptr)));

    void GenerateVoronoiTassellations3D(const Eigen::MatrixXd &polyhedronVertices,
                                        const Eigen::MatrixXi &polyhedronEdges,
                                        const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                                        const unsigned int &numPoints,
                                        const unsigned int &numIterations,
                                        Gedim::IMeshDAO &mesh,
                                        const unsigned int random_seed = static_cast<unsigned int>(time(nullptr)));

    void GenerateVoronoiTassellations2D(const Eigen::MatrixXd &polygonVertices,
                                        const unsigned int &numPoints,
                                        const unsigned int &numIterations,
                                        Gedim::IMeshDAO &mesh,
                                        const unsigned int random_seed = static_cast<unsigned int>(time(nullptr)));

    void GenerateVoronoiTassellations2D(const Eigen::MatrixXd &polygonVertices,
                                        const unsigned int &numIterations,
                                        Eigen::MatrixXd &VoronoiPoints,
                                        Gedim::IMeshDAO &mesh);

    void GenerateVoronoiTassellations3D(const Eigen::MatrixXd &domain_vertices,
                                        const Eigen::MatrixXi &domain_edges,
                                        const std::vector<Eigen::MatrixXi> &domain_faces,
                                        const unsigned int &num_iterations,
                                        Eigen::MatrixXd &VoronoiPoints,
                                        Gedim::IMeshDAO &mesh);
};

} // namespace Gedim

#endif
