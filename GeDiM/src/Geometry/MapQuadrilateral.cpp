#include "MapQuadrilateral.hpp"

using namespace Eigen;
using namespace std;

namespace Gedim
{
// ***************************************************************************
MatrixXd MapQuadrilateral::F(const MatrixXd &vertices, const MatrixXd &x) const
{
    MatrixXd shiftedPoints;
    ShiftSquareQuadraturePoints(x, shiftedPoints);

    MatrixXd evalPsi(4, x.cols());
    evalPsi.row(0) << Psi1(shiftedPoints).transpose();
    evalPsi.row(1) << Psi2(shiftedPoints).transpose();
    evalPsi.row(2) << Psi3(shiftedPoints).transpose();
    evalPsi.row(3) << Psi4(shiftedPoints).transpose();

    return vertices * evalPsi;
}
// ***************************************************************************
MatrixXd MapQuadrilateral::J(const MatrixXd &vertices, const MatrixXd &x) const
{
    MatrixXd shiftedPoints;
    ShiftSquareQuadraturePoints(x, shiftedPoints);

    const unsigned int numPoints = x.cols();
    MatrixXd jacb;
    jacb.setZero(3, 3 * numPoints);
    MatrixXd evaldPsi(8, numPoints);

    evaldPsi.row(0) << dPsi11(shiftedPoints).transpose();
    evaldPsi.row(1) << dPsi21(shiftedPoints).transpose();
    evaldPsi.row(2) << dPsi31(shiftedPoints).transpose();
    evaldPsi.row(3) << dPsi41(shiftedPoints).transpose();
    evaldPsi.row(4) << dPsi12(shiftedPoints).transpose();
    evaldPsi.row(5) << dPsi22(shiftedPoints).transpose();
    evaldPsi.row(6) << dPsi32(shiftedPoints).transpose();
    evaldPsi.row(7) << dPsi42(shiftedPoints).transpose();

    MatrixXd pointJacb;
    for (unsigned int p = 0; p < numPoints; p++)
    {
        pointJacb.setZero(2, 4);
        pointJacb.row(0) << evaldPsi(0, p), evaldPsi(1, p), evaldPsi(2, p), evaldPsi(3, p);
        pointJacb.row(1) << evaldPsi(4, p), evaldPsi(5, p), evaldPsi(6, p), evaldPsi(7, p);

        jacb.block(0, 3 * p, 2, 2) = vertices.block(0, 0, 2, 4) * pointJacb.transpose();
        jacb.block(2, 3 * p, 1, 3) << 0.0, 0.0, 0.5;
    }

    return 2.0 * jacb;
}
// ***************************************************************************
VectorXd MapQuadrilateral::DetJ(const MatrixXd &vertices, const MatrixXd &x) const
{
    MatrixXd jacb = J(vertices, x);

    const unsigned int numPoints = x.cols();
    VectorXd detJacb(numPoints);
    for (unsigned int p = 0; p < numPoints; p++)
        detJacb[p] = jacb.block(0, 3 * p, 3, 3).determinant();

    return detJacb;
}
// ***************************************************************************
} // namespace Gedim
