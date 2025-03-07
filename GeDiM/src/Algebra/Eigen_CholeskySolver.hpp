#ifndef __EIGENCHOLESKYSOLVER_H
#define __EIGENCHOLESKYSOLVER_H

#include "Eigen/Eigen"
#include "ILinearSolver.hpp"

namespace Gedim
{
/// \brief Eigen Cholesky Linear solver
template <typename Eigen_ArrayType = Eigen::VectorXd,
          typename Eigen_SparseArrayType = Eigen::SparseMatrix<double>,
          typename Eigen_SolverType = Eigen::SimplicialLDLT<Eigen_SparseArrayType, Eigen::Lower, Eigen::AMDOrdering<int>>>
class Eigen_CholeskySolver final : public ILinearSolver
{
  private:
    Eigen_SolverType linearSolver; ///< The solver
    const IArray *_rightHandSide;  ///< The rightHandSide of the linear syste
    IArray *_solution;             ///< The solution of the linear syste

  public:
    Eigen_CholeskySolver()
    {
        _rightHandSide = nullptr;
        _solution = nullptr;
    }
    ~Eigen_CholeskySolver()
    {
        _rightHandSide = nullptr;
        _solution = nullptr;
    }

    void Initialize(const ISparseArray &matrix, const IArray &rightHandSide, IArray &solution, const Configuration &config = Configuration());

    void Initialize(const IArray &rightHandSide, IArray &solution, const Configuration &config = Configuration());

    ILinearSolver::SolutionInfo Solve() const;

    void Initialize(const ISparseArray &matrix, const Configuration &config = Configuration());
    ILinearSolver::SolutionInfo Solve(const IArray &rightHandSide, IArray &solution) const;
};
} // namespace Gedim

#endif // __EIGENCHOLESKYSOLVER_H
