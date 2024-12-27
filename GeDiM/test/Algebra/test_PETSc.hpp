#ifndef __TEST_PETSc_H
#define __TEST_PETSc_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "Gedim_Macro.hpp"

#if ENABLE_PETSC == 1

#include "PETSc_Array.hpp"
#include "PETSc_KSPSolver.hpp"
#include "PETSc_SparseArray.hpp"

namespace UnitTesting
{
TEST(TestPETSc, TestPETSc_Solvers)
{
    PetscInitialize(nullptr, nullptr, nullptr, nullptr);

    // generic matrix
    Gedim::PETSc_SparseArray<> A;
    A.SetSize(2, 2);
    A.Triplet(0, 0, 17);
    A.Triplet(1, 1, 85);
    A.Triplets({0, 1}, {1, 0}, {38.0, 38.0});
    A.Create();

    Gedim::PETSc_Array<> b;
    b.SetSize(2);
    b.SetValues({0}, {50.0});
    b.AddValues({0}, {5.0});
    b.SetValue(1, 120.0);
    b.AddValue(1, 3.0);
    b.Create();

    Gedim::PETSc_Array<> x;
    x.SetSize(2);
    x.Create();

    // PCG
    {
        Gedim::PETSc_KSPSolver<Vec, Mat, Gedim::PETSc_SolverTypes::PETSc_KSPCG> solver;
        solver.Initialize(A, b, x, {1000, 1.0e-12});
        const auto solver_result = solver.Solve();

        const auto solution = x.GetValues();

        ASSERT_TRUE(solver_result.Iterations < 1000);
        ASSERT_TRUE(solver_result.Residual < 1.0e-12);

        Gedim::PETSc_Array<> exact_x;
        exact_x.SetSize(2);
        exact_x.Create();
        exact_x.Ones();

        exact_x -= x;
        ASSERT_TRUE(exact_x.Norm() <= 1.0e-11 * x.Norm());
    }

    // BICG
    {
        Gedim::PETSc_KSPSolver<Vec, Mat, Gedim::PETSc_SolverTypes::PETSc_KSPBICG> solver;
        solver.Initialize(A, b, x, {1000, 1.0e-12});
        const auto solver_result = solver.Solve();

        const auto solution = x.GetValues();

        ASSERT_TRUE(solver_result.Iterations < 1000);
        ASSERT_TRUE(solver_result.Residual < 1.0e-12);

        Gedim::PETSc_Array<> exact_x;
        exact_x.SetSize(2);
        exact_x.Create();
        exact_x.Ones();

        exact_x -= x;
        ASSERT_TRUE(exact_x.Norm() <= 1.0e-11 * x.Norm());
    }

    // GMRES
    {
        Gedim::PETSc_KSPSolver<Vec, Mat, Gedim::PETSc_SolverTypes::PETSc_KSPGMRES> solver;
        solver.Initialize(A, b, x, {1000, 1.0e-12});
        const auto solver_result = solver.Solve();

        const auto solution = x.GetValues();

        ASSERT_TRUE(solver_result.Iterations < 1000);
        ASSERT_TRUE(solver_result.Residual < 1.0e-12);

        Gedim::PETSc_Array<> exact_x;
        exact_x.SetSize(2);
        exact_x.Create();
        exact_x.Ones();

        exact_x -= x;
        ASSERT_TRUE(exact_x.Norm() <= 1.0e-11 * x.Norm());
    }

    PetscFinalize();
}
} // namespace UnitTesting

#endif // ENABLE_PETSC

#endif // __TEST_Eigen_H
