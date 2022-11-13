#ifndef __TEST_Eigen_H
#define __TEST_Eigen_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "Eigen_LUSolver.hpp"
#include "Eigen_CholeskySolver.hpp"
#include "Eigen_SparseArray.hpp"
#include "Eigen_Array.hpp"

namespace UnitTesting
{
  TEST(TestEigen, TestLU_GenericMatrix)
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

      Gedim::Eigen_LUSolver<Eigen::VectorXd, Eigen::SparseMatrix<double>> solver;
      solver.Initialize(A, b, x);
      solver.Solve();

      ASSERT_DOUBLE_EQ(1.0, x[0]);
      ASSERT_DOUBLE_EQ(1.0, x[1]);
    }
    catch (const std::exception& exception)
    {
      std::cerr<< exception.what()<< std::endl;
      FAIL();
    }
  }

  TEST(TestEigen, TestLU_SymmetricMatrix)
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

      Gedim::Eigen_LUSolver<Eigen::VectorXd, Eigen::SparseMatrix<double>> solver;
      solver.Initialize(A, b, x);
      solver.Solve();

      ASSERT_DOUBLE_EQ(1.0, x[0]);
      ASSERT_DOUBLE_EQ(1.0, x[1]);
    }
    catch (const std::exception& exception)
    {
      std::cerr<< exception.what()<< std::endl;
      FAIL();
    }
  }

  TEST(TestEigen, TestCholesky_SymmetricMatrix)
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

      Gedim::Eigen_CholeskySolver<Eigen::VectorXd, Eigen::SparseMatrix<double>> solver;
      solver.Initialize(A, b, x);
      solver.Solve();

      ASSERT_DOUBLE_EQ(1.0, x[0]);
      ASSERT_DOUBLE_EQ(1.0, x[1]);
    }
    catch (const std::exception& exception)
    {
      std::cerr<< exception.what()<< std::endl;
      FAIL();
    }
  }
}

#endif // __TEST_Eigen_H
