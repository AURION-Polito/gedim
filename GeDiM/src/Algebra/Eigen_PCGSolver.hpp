#ifndef __EIGENPCGSOLVER_H
#define __EIGENPCGSOLVER_H

#include "ILinearSolver.hpp"
#include "Eigen/Eigen"

namespace Gedim
{
  /// \brief Eigen PCG Linear solver
  template<typename Eigen_ArrayType = Eigen::VectorXd,
           typename Eigen_SparseArrayType = Eigen::SparseMatrix<double>,
           typename Eigen_SolverType = Eigen::ConjugateGradient<Eigen_SparseArrayType, Eigen::Lower, Eigen::DiagonalPreconditioner<double>>>
  class Eigen_PCGSolver final : public ILinearSolver
  {
    private:
      Eigen_SolverType linearSolver; ///< The solver
      const IArray* _rightHandSide; ///< The rightHandSide of the linear syste
      IArray* _solution; ///< The solution of the linear syste
      Configuration _config;

    public:
      Eigen_PCGSolver();
      ~Eigen_PCGSolver();

      void Initialize(const ISparseArray& matrix,
                      const IArray& rightHandSide,
                      IArray& solution,
                      const Configuration& config = { 100, 1e-6 });

      void Initialize(const IArray& rightHandSide,
                      IArray& solution,
                      const Configuration& config = { 100, 1e-6 });

      ILinearSolver::SolutionInfo Solve() const;

      void Initialize(const ISparseArray& matrix,
                      const Configuration& config = { 100, 1e-6 });
      ILinearSolver::SolutionInfo Solve(const IArray& rightHandSide,
                                        IArray& solution) const;
  };
}

#endif // __EIGENPCGSOLVER_H
