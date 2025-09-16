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

#ifndef __MapHexahedron_H
#define __MapHexahedron_H

#include "Eigen/Eigen"

namespace Gedim
{
class MapHexahedron
{
  public:
    struct MapHexahedronData final
    {
        Eigen::MatrixXd ReferenceVertices;
        Eigen::MatrixXd A;
        Eigen::MatrixXd Coefficients;
        Eigen::MatrixXd Vertices;
        std::map<unsigned int, unsigned int> VertexOrder;
    };

  private:
    static inline Eigen::MatrixXd Psi(const Eigen::MatrixXd &points)
    {
        const Eigen::ArrayXd x = points.row(0);
        const Eigen::ArrayXd y = points.row(1);
        const Eigen::ArrayXd z = points.row(2);

        Eigen::MatrixXd evalPsi = Eigen::MatrixXd::Zero(8, points.cols());
        evalPsi.row(0) = (1.0 - x) * (1.0 - y) * (1.0 - z);
        evalPsi.row(1) = x * (1.0 - y) * (1.0 - z);
        evalPsi.row(2) = x * y * (1.0 - z);
        evalPsi.row(3) = y * (1.0 - x) * (1.0 - z);
        evalPsi.row(4) = (1.0 - x) * (1.0 - y) * z;
        evalPsi.row(5) = x * (1.0 - y) * z;
        evalPsi.row(6) = x * y * z;
        evalPsi.row(7) = y * (1.0 - x) * z;

        return evalPsi;
    }

    static inline std::vector<Eigen::MatrixXd> dPsi(const Eigen::MatrixXd &points)
    {
        const Eigen::ArrayXd x = points.row(0);
        const Eigen::ArrayXd y = points.row(1);
        const Eigen::ArrayXd z = points.row(2);

        std::vector<Eigen::MatrixXd> evaldPsi(3, Eigen::MatrixXd(8, points.cols()));

        evaldPsi[0].row(0) = -(1.0 - y) * (1.0 - z);
        evaldPsi[0].row(1) = (1.0 - y) * (1.0 - z);
        evaldPsi[0].row(2) = y * (1.0 - z);
        evaldPsi[0].row(3) = -y * (1.0 - z);
        evaldPsi[0].row(4) = -(1.0 - y) * z;
        evaldPsi[0].row(5) = (1.0 - y) * z;
        evaldPsi[0].row(6) = y * z;
        evaldPsi[0].row(7) = -y * z;

        evaldPsi[1].row(0) = -(1.0 - x) * (1.0 - z);
        evaldPsi[1].row(1) = -x * (1.0 - z);
        evaldPsi[1].row(2) = x * (1.0 - z);
        evaldPsi[1].row(3) = (1.0 - x) * (1.0 - z);
        evaldPsi[1].row(4) = -(1.0 - x) * z;
        evaldPsi[1].row(5) = -x * z;
        evaldPsi[1].row(6) = x * z;
        evaldPsi[1].row(7) = (1.0 - x) * z;

        evaldPsi[2].row(0) = -(1.0 - x) * (1.0 - y);
        evaldPsi[2].row(1) = -x * (1.0 - y);
        evaldPsi[2].row(2) = -x * y;
        evaldPsi[2].row(3) = -y * (1.0 - x);
        evaldPsi[2].row(4) = (1.0 - x) * (1.0 - y);
        evaldPsi[2].row(5) = x * (1.0 - y);
        evaldPsi[2].row(6) = x * y;
        evaldPsi[2].row(7) = y * (1.0 - x);

        return evaldPsi;
    }

  public:
    MapHexahedron()
    {
    }
    ~MapHexahedron()
    {
    }

    static MapHexahedronData Compute(const Eigen::MatrixXd &vertices, const std::vector<Eigen::MatrixXi> &faces);

    static Eigen::MatrixXd F(const MapHexahedron::MapHexahedronData &mapData, const Eigen::MatrixXd &referencePoints);

    static Eigen::MatrixXd FInv(const MapHexahedron::MapHexahedronData &mapData, const Eigen::MatrixXd &points);

    static Eigen::MatrixXd J(const MapHexahedron::MapHexahedronData &mapData, const Eigen::MatrixXd &referencePoints);

    static Eigen::MatrixXd JInv(const MapHexahedron::MapHexahedronData &mapData, const Eigen::MatrixXd &referencePoints);

    static Eigen::VectorXd DetJ(const MapHexahedron::MapHexahedronData &mapData, const Eigen::MatrixXd &referencePoints);
};
} // namespace Gedim

#endif // __MapHexahedron_H
