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

#ifndef __TEST_LP_UTILITIES_H
#define __TEST_LP_UTILITIES_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "LPUtilities.hpp"

namespace UnitTesting
{
TEST(TestLPUtilities, TestLPUtilities_StandardForm_1)
{
    const unsigned int n = 3;
    const unsigned int m = 3;
    Eigen::VectorXd cost_vector_data = Eigen::VectorXd::Zero(n);
    cost_vector_data << -4.0, -5.0, 1.0;
    Eigen::MatrixXd A_data = Eigen::MatrixXd::Zero(m, n);
    A_data << 2.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 2.0, 0.0;
    Eigen::VectorXd b_data = Eigen::VectorXd::Zero(m);
    b_data << 7.0, 16.0, 8.0;
    std::vector<Gedim::LPUtilties::Simplex::ConstraintType> constraint_types_data = {
        Gedim::LPUtilties::Simplex::ConstraintType::GE,
        Gedim::LPUtilties::Simplex::ConstraintType::LE,
        Gedim::LPUtilties::Simplex::ConstraintType::EQ};

    Eigen::VectorXd LB_data = Eigen::VectorXd::Zero(n);
    LB_data << 0.0, std::nan(""), std::nan("");
    Eigen::VectorXd UB_data = Eigen::VectorXd::Zero(n);
    UB_data << std::nan(""), 0.0, std::nan("");
    Gedim::LPUtilties::Simplex simplex(n, cost_vector_data, A_data, b_data, constraint_types_data, LB_data, UB_data);

    Eigen::VectorXd cost_vector_final = Eigen::VectorXd::Zero(6);
    cost_vector_final << -4.0, 5.0, 1.0, -1.0, 0.0, 0.0;

    Eigen::VectorXd b_final = Eigen::VectorXd::Zero(3);
    b_final << -7.0, 16.0, 8.0;

    Eigen::MatrixXd A_final = Eigen::MatrixXd::Zero(3, 6);
    A_final << -2.0, 0.0, -1.0, 1.0, 1.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 1.0, 1.0, -2.0, 0.0, 0.0, 0.0, 0.0;

    Eigen::VectorXd trasl_sol = Eigen::VectorXd::Zero(n);
    trasl_sol << 0.0, 0.0, 0.0;

    Eigen::MatrixXd map_sol = Eigen::MatrixXd::Zero(3, 6);
    map_sol << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0;

    ASSERT_TRUE(simplex.n_final == 6);
    ASSERT_TRUE(simplex.m_final == 3);

    ASSERT_TRUE((simplex.cost_vector_final - cost_vector_final).norm() <= 1.e-12);
    ASSERT_TRUE((simplex.b_final - b_final).norm() <= 1.e-12);
    ASSERT_TRUE((simplex.A_final - A_final).norm() <= 1.e-12);
    ASSERT_TRUE(abs(simplex.trasl_cost) <= 1.e-12);
    ASSERT_TRUE((simplex.trasl_sol - trasl_sol).norm() <= 1.e-12);
    ASSERT_TRUE((simplex.map_sol - map_sol).norm() <= 1.e-12);

    std::vector<unsigned int> slack_variable = {4, 5};
    std::vector<unsigned int> no_slack_variable = {0, 1, 2, 3};

    for (unsigned int i = 0; i < no_slack_variable.size(); i++)
        ASSERT_TRUE(simplex.no_slack_indices[i] == no_slack_variable[i]);

    for (unsigned int i = 0; i < slack_variable.size(); i++)
        ASSERT_TRUE(simplex.slack_indices[i] == slack_variable[i]);
}

TEST(TestLPUtilities, TestLPUtilities_StandardForm_2)
{
    const unsigned int n = 2;
    const unsigned int m = 2;
    Eigen::VectorXd cost_vector_data = Eigen::VectorXd::Zero(n);
    cost_vector_data << -8.0, -3.0;
    Eigen::MatrixXd A_data = Eigen::MatrixXd::Zero(m, n);
    A_data << 4.0, 5.0, 4.0, 10.0;
    Eigen::VectorXd b_data = Eigen::VectorXd::Zero(m);
    b_data << 10.0, 15.0;
    std::vector<Gedim::LPUtilties::Simplex::ConstraintType> constraint_types_data = {
        Gedim::LPUtilties::Simplex::ConstraintType::LE,
        Gedim::LPUtilties::Simplex::ConstraintType::LE};

    Eigen::VectorXd LB_data = Eigen::VectorXd::Zero(n);
    LB_data << 0.0, 0.0;
    Eigen::VectorXd UB_data = Eigen::VectorXd::Zero(n);
    UB_data << std::nan(""), 1.0;
    Gedim::LPUtilties::Simplex simplex(n, cost_vector_data, A_data, b_data, constraint_types_data, LB_data, UB_data);

    Eigen::VectorXd cost_vector_final = Eigen::VectorXd::Zero(5);
    cost_vector_final << -8.0, -3.0, 0.0, 0.0, 0.0;

    Eigen::VectorXd b_final = Eigen::VectorXd::Zero(3);
    b_final << 10.0, 15.0, 1.0;

    Eigen::MatrixXd A_final = Eigen::MatrixXd::Zero(3, 5);
    A_final << 4.0, 5.0, 0.0, 1.0, 0.0, 4.0, 10.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0;

    Eigen::VectorXd trasl_sol = Eigen::VectorXd::Zero(n);
    trasl_sol << 0.0, 0.0;

    Eigen::MatrixXd map_sol = Eigen::MatrixXd::Zero(n, 5);
    map_sol << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0;

    ASSERT_TRUE(simplex.n_final == 5);
    ASSERT_TRUE(simplex.m_final == 3);

    ASSERT_TRUE((simplex.cost_vector_final - cost_vector_final).norm() <= 1.e-12);
    ASSERT_TRUE((simplex.b_final - b_final).norm() <= 1.e-12);
    ASSERT_TRUE((simplex.A_final - A_final).norm() <= 1.e-12);
    ASSERT_TRUE(abs(simplex.trasl_cost) <= 1.e-12);
    ASSERT_TRUE((simplex.trasl_sol - trasl_sol).norm() <= 1.e-12);
    ASSERT_TRUE((simplex.map_sol - map_sol).norm() <= 1.e-12);

    std::vector<unsigned int> slack_variable = {2, 3, 4};
    std::vector<unsigned int> no_slack_variable = {0, 1};

    for (unsigned int i = 0; i < no_slack_variable.size(); i++)
        ASSERT_TRUE(simplex.no_slack_indices[i] == no_slack_variable[i]);

    for (unsigned int i = 0; i < slack_variable.size(); i++)
        ASSERT_TRUE(simplex.slack_indices[i] == slack_variable[i]);
}

TEST(TestLPUtilities, TestLPUtilities_SolvePrimal_1)
{
    const unsigned int n = 2;
    const unsigned int m = 3;
    Eigen::VectorXd cost_vector_data = Eigen::VectorXd::Zero(n);
    cost_vector_data << 2.0, 1.0;
    Eigen::MatrixXd A_data = Eigen::MatrixXd::Zero(m, n);
    A_data << 2.0, 1.0, 1.0, -1.0, 0.0, 1.0;
    Eigen::VectorXd b_data = Eigen::VectorXd::Zero(m);
    b_data << 5.0, 1.0, 4.0;
    std::vector<Gedim::LPUtilties::Simplex::ConstraintType> constraint_types_data = {
        Gedim::LPUtilties::Simplex::ConstraintType::LE,
        Gedim::LPUtilties::Simplex::ConstraintType::LE,
        Gedim::LPUtilties::Simplex::ConstraintType::LE};

    Eigen::VectorXd LB_data = Eigen::VectorXd::Zero(n);
    LB_data << 0.0, 0.0;
    Eigen::VectorXd UB_data = Eigen::VectorXd::Zero(n);
    UB_data << std::nan(""), std::nan("");
    Gedim::LPUtilties::Simplex simplex(n, cost_vector_data, A_data, b_data, constraint_types_data, LB_data, UB_data);

    Eigen::VectorXd cost_vector_final = Eigen::VectorXd::Zero(5);
    cost_vector_final << 2.0, 1.0, 0.0, 0.0, 0.0;

    Eigen::VectorXd b_final = Eigen::VectorXd::Zero(3);
    b_final << 5.0, 1.0, 4.0;

    Eigen::MatrixXd A_final = Eigen::MatrixXd::Zero(3, 5);
    A_final << 2.0, 1.0, 1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0;

    Eigen::VectorXd trasl_sol = Eigen::VectorXd::Zero(n);
    trasl_sol << 0.0, 0.0;

    Eigen::MatrixXd map_sol = Eigen::MatrixXd::Zero(n, 5);
    map_sol << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0;

    ASSERT_TRUE(simplex.n_final == 5);
    ASSERT_TRUE(simplex.m_final == 3);

    ASSERT_TRUE((simplex.cost_vector_final - cost_vector_final).norm() <= 1.e-12);
    ASSERT_TRUE((simplex.b_final - b_final).norm() <= 1.e-12);
    ASSERT_TRUE((simplex.A_final - A_final).norm() <= 1.e-12);
    ASSERT_TRUE(abs(simplex.trasl_cost) <= 1.e-12);
    ASSERT_TRUE((simplex.trasl_sol - trasl_sol).norm() <= 1.e-12);
    ASSERT_TRUE((simplex.map_sol - map_sol).norm() <= 1.e-12);

    std::vector<unsigned int> slack_variable = {2, 3, 4};
    std::vector<unsigned int> no_slack_variable = {0, 1};

    for (unsigned int i = 0; i < no_slack_variable.size(); i++)
        ASSERT_TRUE(simplex.no_slack_indices[i] == no_slack_variable[i]);

    for (unsigned int i = 0; i < slack_variable.size(); i++)
        ASSERT_TRUE(simplex.slack_indices[i] == slack_variable[i]);

    simplex.solve_primal_simplex(slack_variable, no_slack_variable);

    std::vector<unsigned int> basis = {1, 0, 4};
    std::vector<unsigned int> no_basis = {3, 2};

    for (unsigned int i = 0; i < no_slack_variable.size(); i++)
        ASSERT_TRUE(no_basis[i] == no_slack_variable[i]);

    for (unsigned int i = 0; i < slack_variable.size(); i++)
        ASSERT_TRUE(basis[i] == slack_variable[i]);

    ASSERT_TRUE(abs(simplex.result.solution(0) - 2.0) < 1.0e-12);

    ASSERT_TRUE(abs(simplex.result.solution(1) - 1.0) < 1.0e-12);
}

} // namespace UnitTesting

#endif
