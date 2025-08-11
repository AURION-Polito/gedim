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

#ifndef __MapQuadrilateral_H
#define __MapQuadrilateral_H

#include "Eigen/Eigen"

namespace Gedim
{
class MapQuadrilateral
{
  public:
    Eigen::MatrixXd ReferencePoints;

  private:
    inline Eigen::MatrixXd Psi(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd x = points.row(0);
        const Eigen::ArrayXd y = points.row(1);

        Eigen::MatrixXd evalPsi = Eigen::MatrixXd::Zero(4, points.cols());
        evalPsi.row(0) = (1.0 - x) * (1.0 - y);
        evalPsi.row(1) = x * (1.0 - y);
        evalPsi.row(2) = x * y;
        evalPsi.row(3) = y * (1.0 - x);

        return evalPsi;
    }

    inline std::vector<Eigen::MatrixXd> dPsi(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd x = points.row(0);
        const Eigen::ArrayXd y = points.row(1);

        std::vector<Eigen::MatrixXd> evaldPsi(2, Eigen::MatrixXd(4, points.cols()));

        evaldPsi[0].row(0) = -(1.0 - y);
        evaldPsi[0].row(1) = (1.0 - y);
        evaldPsi[0].row(2) = y;
        evaldPsi[0].row(3) = -y;
        evaldPsi[1].row(0) = -(1.0 - x);
        evaldPsi[1].row(1) = -x;
        evaldPsi[1].row(2) = x;
        evaldPsi[1].row(3) = (1.0 - x);

        return evaldPsi;
    }

  public:
    MapQuadrilateral()
    {
        ReferencePoints.resize(3, 4);
        ReferencePoints.col(0) << 0.0, 0.0, 0.0;
        ReferencePoints.col(1) << 1.0, 0.0, 0.0;
        ReferencePoints.col(2) << 1.0, 1.0, 0.0;
        ReferencePoints.col(3) << 0.0, 1.0, 0.0;
    }
    ~MapQuadrilateral()
    {
    }

    Eigen::MatrixXd F(const Eigen::MatrixXd &vertices, const Eigen::MatrixXd &referencePoints) const;

    Eigen::MatrixXd FInv(const Eigen::MatrixXd &vertices, const Eigen::MatrixXd &points) const;

    Eigen::MatrixXd J(const Eigen::MatrixXd &vertices, const Eigen::MatrixXd &referencePoints) const;

    Eigen::MatrixXd JInv(const Eigen::MatrixXd &vertices, const Eigen::MatrixXd &referencePoints) const;

    Eigen::VectorXd DetJ(const Eigen::MatrixXd &vertices, const Eigen::MatrixXd &referencePoints) const;
};
} // namespace Gedim

#endif // __MappingQuadrilateral_H
