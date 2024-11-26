#ifndef __PETSCKSPSOLVER_H
#define __PETSCKSPSOLVER_H

#include "Macro.hpp"

#if ENABLE_PETSC == 1

#include "petscmat.h"
#include "petscvec.h"
#include "petscksp.h"

#include "ILinearSolver.hpp"

namespace Gedim
{
  enum struct PETSc_SolverTypes
  {
    PETSc_KSPCG = 0,
    PETSc_KSPBICG = 1,
    PETSc_KSPGMRES = 2
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
  template<typename PETSc_ArrayType = Vec,
           typename PETSc_SparseArrayType = Mat,
           PETSc_SolverTypes PETSc_SolverType = PETSc_SolverTypes::PETSc_KSPGMRES,
           PETSc_Preconditioners PETSc_Preconditioner = PETSc_Preconditioners::PETSc_DEFAULT>
  class PETSc_KSPSolver final : public ILinearSolver
  {
    private:
      PC preconditioner;
      KSP linearSolver; ///< The solver
      const IArray* _rightHandSide; ///< The rightHandSide of the linear syste
      IArray* _solution; ///< The solution of the linear syste
      Configuration _config;

    public:
      PETSc_KSPSolver();
      ~PETSc_KSPSolver();

      void Initialize(const ISparseArray& matrix,
                      const IArray& rightHandSide,
                      IArray& solution,
                      const Configuration& config = { 100, 1e-6 });

      ILinearSolver::SolutionInfo Solve() const;

      void Initialize(const ISparseArray& matrix,
                      const Configuration& config = { 100, 1e-6 });
      ILinearSolver::SolutionInfo Solve(const IArray& rightHandSide,
                                        IArray& solution) const;
  };
}

#endif // ENABLE_PETSC

#endif // __PETScPCGSOLVER_H
