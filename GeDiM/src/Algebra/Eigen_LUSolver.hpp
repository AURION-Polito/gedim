#ifndef __EIGENLUSOLVER_H
#define __EIGENLUSOLVER_H

#include "Eigen/Eigen"
#include "ILinearSolver.hpp"

namespace Gedim
{
/// \brief Eigen LU Linear solver
template <typename Eigen_ArrayType = Eigen::VectorXd, typename Eigen_SparseArrayType = Eigen::SparseMatrix<double>>
class Eigen_LUSolver final : public ILinearSolver
{
  private:
    Eigen::SparseLU<Eigen_SparseArrayType> linearSolver; ///< The solver
    const IArray *_rightHandSide;                        ///< The rightHandSide of the linear system
    IArray *_solution;                                   ///< The solution of the linear syste

  public:
    Eigen_LUSolver()
    {
        _rightHandSide = nullptr;
        _solution = nullptr;
    }
    ~Eigen_LUSolver()
    {
        _rightHandSide = nullptr;
        _solution = nullptr;
    }

    void Initialize(const ISparseArray &matrix, const IArray &rightHandSide, IArray &solution, const Configuration &config = Configuration());

    void Initialize(const IArray &rightHandSide, IArray &solution, const Configuration &config = Configuration());

    SolutionInfo Solve() const;

    void Initialize(const ISparseArray &matrix, const Configuration &config = Configuration());

    SolutionInfo Solve(const IArray &rightHandSide, IArray &solution) const;
};
} // namespace Gedim

#endif // __EIGENLUSOLVER_H
