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

#ifndef __PETSCKSPSOLVER_H
#define __PETSCKSPSOLVER_H

#include "Gedim_Macro.hpp"

#if ENABLE_PETSC == 1

#include "petscksp.h"
#include "petscmat.h"
#include "petscvec.h"

#include "ILinearSolver.hpp"

namespace Gedim
{
enum struct PETSc_SolverTypes
{
    PETSc_KSPCG = 0,
    PETSc_KSPBICG = 1,
    PETSc_KSPGMRES = 2,
    PETSc_KSPBICGS = 3
};

enum struct PETSc_Preconditioners
{
    PETSc_DEFAULT = -1,
    PETSc_PCNONE = 0,
    PETSc_PCILU = 1,
    PETSc_PCJACOBI = 2,
    PETSc_PCFIELDSPLIT = 3
};

/// \brief PETSc PCG Linear solver
template <typename PETSc_ArrayType = Vec, typename PETSc_SparseArrayType = Mat, PETSc_SolverTypes PETSc_SolverType = PETSc_SolverTypes::PETSc_KSPGMRES, PETSc_Preconditioners PETSc_Preconditioner = PETSc_Preconditioners::PETSc_DEFAULT>
class PETSc_KSPSolver final : public ILinearSolver
{
  private:
    PC preconditioner;
    KSP linearSolver;             ///< The solver
    const IArray *_rightHandSide; ///< The rightHandSide of the linear syste
    IArray *_solution;            ///< The solution of the linear syste
    Configuration _config;

  public:
    PETSc_KSPSolver();
    ~PETSc_KSPSolver();

    void Initialize(const ISparseArray &matrix, const IArray &rightHandSide, IArray &solution, const Configuration &config = {100, 1e-6});

    ILinearSolver::SolutionInfo Solve() const;

    void Initialize(const ISparseArray &matrix, const Configuration &config = {100, 1e-6});
    ILinearSolver::SolutionInfo Solve(const IArray &rightHandSide, IArray &solution) const;
};
} // namespace Gedim

#endif // ENABLE_PETSC

#endif // __PETScPCGSOLVER_H
