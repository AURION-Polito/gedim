#ifndef __Pardiso_LUSolver_H
#define __Pardiso_LUSolver_H

#include "ILinearSolver.hpp"
#include "Eigen/Eigen"
#include "Macro.hpp"

#if ENABLE_MKL
#include <Eigen/PardisoSupport>
#endif

namespace Gedim
{
  /// \brief Pardiso LU Linear solver
  template<typename Eigen_ArrayType = Eigen::VectorXd,
           typename Eigen_SparseArrayType = Eigen::SparseMatrix<double>>
  class Pardiso_LUSolver final : public ILinearSolver
  {
    private:
#if ENABLE_MKL
      Eigen::PardisoLU<Eigen_SparseArrayType> linearSolver; ///< The solver
#endif
      const IArray* _rightHandSide; ///< The rightHandSide of the linear system
      IArray* _solution; ///< The solution of the linear syste

    public:
      Pardiso_LUSolver();
      ~Pardiso_LUSolver();

      void Initialize(const ISparseArray& matrix,
                      const IArray& rightHandSide,
                      IArray& solution,
                      const ILinearSolver::Configuration& config = Configuration());

      ILinearSolver::SolutionInfo Solve() const;

      void Initialize(const ISparseArray& matrix,
                      const ILinearSolver::Configuration& config = Configuration());

      ILinearSolver::SolutionInfo Solve(const IArray& rightHandSide,
                                        IArray& solution) const;
  };
}

#endif // __Pardiso_LUSolver_H
