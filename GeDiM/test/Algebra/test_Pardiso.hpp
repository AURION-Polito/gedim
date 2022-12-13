#ifndef __TEST_Pardiso_H
#define __TEST_Pardiso_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "Macro.hpp"

#if ENABLE_MKL

#include "Pardiso_LUSolver.hpp"
#include "Pardiso_CholeskySolver.hpp"
#include "Eigen_SparseArray.hpp"
#include "Eigen_Array.hpp"

namespace UnitTesting
{
  TEST(TestPardiso, TestLU_GenericMatrix)
  {
    try
    {
      // generic matrix
      Gedim::Eigen_SparseArray<Eigen::VectorXd, Eigen::SparseMatrix<double>> A;
      A.SetSize(2, 2);
      A.Triplet(0, 0, 17);
      A.Triplet(0, 1, 38);
      A.Triplet(1, 0, 38);
      A.Triplet(1, 1, 85);
      A.Create();

      Gedim::Eigen_Array<Eigen::VectorXd, Eigen::SparseMatrix<double>> b;
      b.SetSize(2);
      b[0] = 55;
      b[1] = 123;

      Gedim::Eigen_Array<Eigen::VectorXd, Eigen::SparseMatrix<double>> x;

      Gedim::Pardiso_LUSolver<Eigen::VectorXd, Eigen::SparseMatrix<double>> solver;
      solver.Initialize(A, b, x);
      solver.Solve();

      ASSERT_TRUE(abs(x[0] - 1.0) < 1e-12);
      ASSERT_TRUE(abs(x[1] - 1.0) < 1e-12);
    }
    catch (const std::exception& exception)
    {
      std::cerr<< exception.what()<< std::endl;
      FAIL();
    }
  }

  TEST(TestPardiso, TestLU_SymmetricMatrix)
  {
    try
    {
      // symmetric matrix
      Gedim::Eigen_SparseArray<Eigen::VectorXd, Eigen::SparseMatrix<double>> A;
      A.SetSize(2, 2, Gedim::ISparseArray::SparseArrayTypes::Symmetric);
      A.Triplet(0, 0, 17);
      A.Triplet(1, 0, 38);
      A.Triplet(1, 1, 85);
      A.Create();

      Gedim::Eigen_Array<Eigen::VectorXd, Eigen::SparseMatrix<double>> b;
      b.SetSize(2);
      b[0] = 55;
      b[1] = 123;

      Gedim::Eigen_Array<Eigen::VectorXd, Eigen::SparseMatrix<double>> x;

      Gedim::Pardiso_LUSolver<Eigen::VectorXd, Eigen::SparseMatrix<double>> solver;
      solver.Initialize(A, b, x);
      solver.Solve();

      ASSERT_TRUE(abs(x[0] - 1.0) < 1e-12);
      ASSERT_TRUE(abs(x[1] - 1.0) < 1e-12);
    }
    catch (const std::exception& exception)
    {
      std::cerr<< exception.what()<< std::endl;
      FAIL();
    }
  }

  TEST(TestPardiso, TestCholesky_SymmetricMatrix)
  {
    try
    {
      // symmetric matrix
      Gedim::Eigen_SparseArray<Eigen::VectorXd, Eigen::SparseMatrix<double>> A;
      A.SetSize(2, 2, Gedim::ISparseArray::SparseArrayTypes::Symmetric);
      A.Triplet(0, 0, 17);
      A.Triplet(1, 0, 38);
      A.Triplet(1, 1, 85);
      A.Create();

      Gedim::Eigen_Array<Eigen::VectorXd, Eigen::SparseMatrix<double>> b;
      b.SetSize(2);
      b[0] = 55;
      b[1] = 123;

      Gedim::Eigen_Array<Eigen::VectorXd, Eigen::SparseMatrix<double>> x;

      Gedim::Pardiso_CholeskySolver<Eigen::VectorXd, Eigen::SparseMatrix<double>> solver;
      solver.Initialize(A, b, x);
      solver.Solve();

      ASSERT_TRUE(abs(x[0] - 1.0) < 1e-12);
      ASSERT_TRUE(abs(x[1] - 1.0) < 1e-12);
    }
    catch (const std::exception& exception)
    {
      std::cerr<< exception.what()<< std::endl;
      FAIL();
    }
  }
}

#endif // ENABLE_MKL

#endif // __TEST_Pardiso_H
