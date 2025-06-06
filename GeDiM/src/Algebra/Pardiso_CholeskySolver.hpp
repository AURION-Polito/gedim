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

#ifndef __Pardiso_CholeskySolver_H
#define __Pardiso_CholeskySolver_H

#include "Eigen/Eigen"
#include "ILinearSolver.hpp"

#include "Gedim_Macro.hpp"

#if ENABLE_MKL
#include <Eigen/PardisoSupport>
#endif

namespace Gedim
{
/// \brief Pardiso Cholesky Linear solver
template <typename Eigen_ArrayType = Eigen::VectorXd, typename Eigen_SparseArrayType = Eigen::SparseMatrix<double>>
class Pardiso_CholeskySolver final : public ILinearSolver
{
  private:
#if ENABLE_MKL
    Eigen::PardisoLDLT<Eigen_SparseArrayType, Eigen::Lower> linearSolver; ///< The solver
#endif
    const IArray *_rightHandSide; ///< The rightHandSide of the linear syste
    IArray *_solution;            ///< The solution of the linear syste

  public:
    Pardiso_CholeskySolver();
    ~Pardiso_CholeskySolver();

    void Initialize(const ISparseArray &matrix,
                    const IArray &rightHandSide,
                    IArray &solution,
                    const ILinearSolver::Configuration &config = ILinearSolver::Configuration());

    ILinearSolver::SolutionInfo Solve() const;

    void Initialize(const ISparseArray &matrix, const ILinearSolver::Configuration &config = Configuration());

    ILinearSolver::SolutionInfo Solve(const IArray &rightHandSide, IArray &solution) const;
};
} // namespace Gedim

#endif // __Pardiso_CholeskySolver_H
