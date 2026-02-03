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

#ifndef __LPUtilties_HPP
#define __LPUtilties_HPP

#include "Eigen/Eigen"
#include <stdlib.h>

namespace Gedim
{
namespace LPUtilties
{

/// primal simplex algorithm to solve
/// max  c' x
/// s.t. Ax = b (slack variables are added in case <= or >=)
///      ub >= x >= lb
struct Simplex
{
    struct PrimalResult
    {
        enum struct ExitValue
        {
            SUCCESS = 1,   // success
            UNBOUNDED = 2, // primal is unbounded
            NOSLACK = 3,   // no starting basis
            MAXIT = 4
        };

        Eigen::VectorXd solution;
        double cost_value;

        Eigen::VectorXd solution_final;

        ExitValue exit_value = ExitValue::UNBOUNDED;
    };

    enum struct ConstraintType
    {
        GE = 1, // greater or equal
        EQ = 2, // equality
        LE = 3, // lower or equal
    };
    unsigned int max_iteration = 1000;

    unsigned int m; // num constraints
    unsigned int n; // num variables

    Eigen::VectorXd cost_vector; // cost functional
    Eigen::MatrixXd A;           // matrix constraints
    Eigen::VectorXd b;           // rhs constraints
    std::vector<ConstraintType> constraint_types;
    Eigen::VectorXd LB; // lower bound, std::nan("") means free
    Eigen::VectorXd UB; // upper bound, std::nan("") means free

    double trasl_cost;
    Eigen::VectorXd trasl_sol;
    Eigen::VectorXd sign_sol;
    Eigen::MatrixXd map_sol;

    unsigned int m_final;              // num final equality constraint after preprocessing
    unsigned int n_final;              // num final variables after preprocessing
    Eigen::MatrixXd A_final;           // matrix constraints after preprocessing
    Eigen::VectorXd b_final;           // rhs constraints after preprocessing
    Eigen::VectorXd cost_vector_final; // cost functional after preprocessing

    std::vector<unsigned int> slack_indices;
    std::vector<unsigned int> no_slack_indices;

    PrimalResult result;

    Simplex(const unsigned int n_data,
            const Eigen::VectorXd &cost_vector_data,
            const Eigen::MatrixXd &A_data,
            const Eigen::VectorXd &b_data,
            const std::vector<ConstraintType> &constraint_types_data,
            const Eigen::VectorXd &LB_data,
            const Eigen::VectorXd &UB_data);

    void solve_primal_simplex(std::vector<unsigned int> &basis, std::vector<unsigned int> &no_basis);
};

} // namespace LPUtilties
} // namespace Gedim

#endif
