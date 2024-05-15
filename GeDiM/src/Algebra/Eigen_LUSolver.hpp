#ifndef __EIGENLUSOLVER_H
#define __EIGENLUSOLVER_H

#include "ILinearSolver.hpp"
#include "Eigen/Eigen"

namespace Gedim
{
/// \brief Eigen LU Linear solver
template<typename Eigen_ArrayType = Eigen::VectorXd,
         typename Eigen_SparseArrayType = Eigen::SparseMatrix<double>>
class Eigen_LUSolver final : public ILinearSolver
{
private:
    Eigen::SparseLU<Eigen_SparseArrayType> linearSolver; ///< The solver
    const IArray* _rightHandSide; ///< The rightHandSide of the linear system
    IArray* _solution; ///< The solution of the linear syste

public:
    Eigen_LUSolver();
    ~Eigen_LUSolver();

    void Initialize(const ISparseArray& matrix,
                    const IArray& rightHandSide,
                    IArray& solution);

    void Initialize(const IArray& rightHandSide,
                    IArray& solution);

    void Solve() const;

    void Initialize(const ISparseArray& matrix);

    void Solve(const IArray& rightHandSide,
               IArray& solution) const;
};
}

#endif // __EIGENLUSOLVER_H
