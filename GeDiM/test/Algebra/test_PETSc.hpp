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
    try
    {
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
      b[0] = 55;
      b[1] = 123;

      Gedim::PETSc_Array x;

      std::cout<< b<< std::endl;

      //Gedim::Eigen_LUSolver<Eigen::VectorXd, Eigen::SparseMatrix<double>> solver;
      //solver.Initialize(A, b, x);
      //solver.Solve();

      //ASSERT_DOUBLE_EQ(1.0, x[0]);
      //ASSERT_DOUBLE_EQ(1.0, x[1]);
    }
    catch (const std::exception& exception)
    {
      std::cerr<< exception.what()<< std::endl;
      FAIL();
    }
  }


}

#endif // ENABLE_PETSC

#endif // __TEST_Eigen_H
