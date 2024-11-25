#ifndef __TEST_PETSc_H
#define __TEST_PETSc_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "PETSc_Array.hpp"

#if ENABLE_PETSC == 1

namespace UnitTesting
{
  TEST(TestPETSc, TestLU_GenericMatrix)
  {
    PetscInitialize(nullptr, nullptr, nullptr, nullptr);

    // generic matrix
    // Gedim::Eigen_SparseArray<Eigen::VectorXd, Eigen::SparseMatrix<double>> A;
    // A.SetSize(2, 2);
    // A.Triplet(0, 0, 17);
    // A.Triplet(0, 1, 38);
    // A.Triplet(1, 0, 38);
    // A.Triplet(1, 1, 85);
    // A.Create();

    Gedim::PETSc_Array b;
    b.SetSize(2);
    b.SetValues({0, 1}, { 55, 123 });

    Gedim::PETSc_Array x;

    //Gedim::Eigen_LUSolver<Eigen::VectorXd, Eigen::SparseMatrix<double>> solver;
    //solver.Initialize(A, b, x);
    //solver.Solve();

    //ASSERT_DOUBLE_EQ(1.0, x[0]);
    //ASSERT_DOUBLE_EQ(1.0, x[1]);

    b.Destroy();

    PetscFinalize();
  }
}

#endif // ENABLE_PETSC

#endif // __TEST_Eigen_H
