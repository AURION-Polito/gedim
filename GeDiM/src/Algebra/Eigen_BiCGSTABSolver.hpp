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

#ifndef __EIGENBiCGSTABSOLVER_H
#define __EIGENBiCGSTABSOLVER_H

#include "Eigen/Eigen"
#include "ILinearSolver.hpp"

namespace Gedim
{
/// \brief Eigen BiCGSTAB Linear solver
template <typename Eigen_ArrayType = Eigen::VectorXd,
          typename Eigen_SparseArrayType = Eigen::SparseMatrix<double>,
          typename Eigen_SolverType = Eigen::BiCGSTAB<Eigen_SparseArrayType, Eigen::DiagonalPreconditioner<double>>>
class Eigen_BiCGSTABSolver final : public ILinearSolver
{
  public:
    Eigen_SolverType linearSolver; ///< The solver
    const IArray *_rightHandSide;  ///< The rightHandSide of the linear syste
    IArray *_solution;             ///< The solution of the linear syste
    Configuration _config;

  public:
    Eigen_BiCGSTABSolver()
    {
        _rightHandSide = nullptr;
        _solution = nullptr;
    }
    ~Eigen_BiCGSTABSolver()
    {
        _rightHandSide = nullptr;
        _solution = nullptr;
    }

    void Initialize(const ISparseArray &matrix, const IArray &rightHandSide, IArray &solution, const Configuration &config = {100, 1e-6});

    void Initialize(const IArray &rightHandSide, IArray &solution, const Configuration &config = {100, 1e-6});

    ILinearSolver::SolutionInfo Solve() const;

    void Initialize(const ISparseArray &matrix, const Configuration &config = {100, 1e-6});
    ILinearSolver::SolutionInfo Solve(const IArray &rightHandSide, IArray &solution) const;
};
} // namespace Gedim

#endif // __EIGENBiCGSTABSOLVER_H
