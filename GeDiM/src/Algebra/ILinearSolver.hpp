#ifndef __ILINEARSOLVER_H
#define __ILINEARSOLVER_H

#include "ISparseArray.hpp"
#include "IArray.hpp"

namespace Gedim
{
  /// \brief Interface used for algebra linear solvers
  class ILinearSolver
  {
    public:
      struct Configuration final
      {
          unsigned int MaxIterations;
          double Tolerance;
      };

      struct SolutionInfo final
      {
          unsigned int Iterations;
          double Residual;
      };

    public:
      virtual ~ILinearSolver() { }

      /// \brief Initialize the linear solver Ax = b
      /// \param matrix The matrix A
      /// \param rightHandSide The right-hand side b
      /// \param solution The solution x
      virtual void Initialize(const ISparseArray& matrix,
                              const IArray& rightHandSide,
                              IArray& solution,
                              const Configuration& config = { 100, 1e-6 }) = 0;

      /// \brief Compute the solution
      virtual SolutionInfo Solve() const = 0;

      /// \brief Initialize the linear solver for system Ax = b
      /// \param matrix The matrix A
      virtual void Initialize(const ISparseArray& matrix,
                              const Configuration& config = { 100, 1e-6 }) = 0;

      /// \brief Compute the solution for system Ax = b
      /// \param rightHandSide The right-hand side b
      /// \param solution The solution x
      virtual SolutionInfo Solve(const IArray& rightHandSide,
                                 IArray& solution) const = 0;
  };
}

#endif // __ILINEARSOLVER_H
