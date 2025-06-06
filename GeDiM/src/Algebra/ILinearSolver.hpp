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

#ifndef __ILINEARSOLVER_H
#define __ILINEARSOLVER_H

#include "IArray.hpp"
#include "ISparseArray.hpp"

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
    virtual ~ILinearSolver()
    {
    }

    /// \brief Initialize the linear solver Ax = b
    /// \param matrix The matrix A
    /// \param rightHandSide The right-hand side b
    /// \param solution The solution x
    virtual void Initialize(const ISparseArray &matrix,
                            const IArray &rightHandSide,
                            IArray &solution,
                            const Configuration &config = {100, 1e-6}) = 0;

    /// \brief Compute the solution
    virtual SolutionInfo Solve() const = 0;

    /// \brief Initialize the linear solver for system Ax = b
    /// \param matrix The matrix A
    virtual void Initialize(const ISparseArray &matrix, const Configuration &config = {100, 1e-6}) = 0;

    /// \brief Compute the solution for system Ax = b
    /// \param rightHandSide The right-hand side b
    /// \param solution The solution x
    virtual SolutionInfo Solve(const IArray &rightHandSide, IArray &solution) const = 0;
};
} // namespace Gedim

#endif // __ILINEARSOLVER_H
