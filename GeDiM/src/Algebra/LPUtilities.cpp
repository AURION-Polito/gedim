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

#include "LPUtilities.hpp"
#include <numeric>

using namespace std;
using namespace Eigen;

namespace Gedim
{
namespace LPUtilties
{
//**********************************************************
Simplex::Simplex(const unsigned int n_data,
                 const Eigen::VectorXd &cost_vector_data,
                 const Eigen::MatrixXd &A_data,
                 const Eigen::VectorXd &b_data,
                 const std::vector<ConstraintType> &constraint_types_data,
                 const Eigen::VectorXd &LB_data,
                 const Eigen::VectorXd &UB_data)
    : n(n_data), cost_vector(cost_vector_data), A(A_data), b(b_data), constraint_types(constraint_types_data),
      LB(LB_data), UB(UB_data)
{
    m = A.rows();
    trasl_cost = 0.0;
    trasl_sol = Eigen::VectorXd::Zero(n);
    sign_sol = Eigen::VectorXd::Ones(n);
    n_final = n;
    m_final = m;

    no_slack_indices.resize(n);
    std::iota(no_slack_indices.begin(), no_slack_indices.end(), 0);

    for (unsigned int i = 0; i < n; i++)
    {
        if (UB[i] <= 0.0 && std::isnan(LB[i]))
        {
            sign_sol(i) = -1.0;
            A.col(i) = -A.col(i);
            cost_vector(i) = -cost_vector(i);
            LB[i] = -UB[i];
            UB[i] = std::nan("");
        }

        if (!std::isnan(LB[i]))
        {
            trasl_sol(i) = this->LB(i);
            trasl_cost += this->cost_vector(i) * this->LB(i);
            this->b = this->b - this->A.col(i) * this->LB(i);
            this->UB(i) = this->UB(i) - this->LB(i);
        }
        else
            n_final++;

        if (!std::isnan(UB[i]))
        {
            m_final++;
            n_final++;
        }
    }

    for (unsigned int i = 0; i < m; i++)
    {
        switch (constraint_types[i])
        {
        case ConstraintType::GE: {
            A.row(i) = -A.row(i);
            b(i) = -b(i);
            n_final++;
        }
        break;
        case ConstraintType::EQ:
            break;
        case ConstraintType::LE:
            n_final++;
            break;
        }
    }

    b_final = Eigen::VectorXd::Zero(m_final);
    b_final.segment(0, m) = b;
    A_final = Eigen::MatrixXd::Zero(m_final, n_final);
    A_final.topLeftCorner(m, n) = A;
    cost_vector_final = Eigen::VectorXd::Zero(n_final);
    cost_vector_final.segment(0, n) = cost_vector;
    map_sol = Eigen::MatrixXd::Zero(n, n_final);

    unsigned int count_add_var = n;
    unsigned int count_add_eq = m;
    for (unsigned int i = 0; i < n; i++)
    {
        if (!std::isnan(UB[i]))
        {
            A_final(count_add_eq, i) = 1.0;
            A_final(count_add_eq, count_add_var) = 1.0;
            b_final(count_add_eq) = UB[i];

            slack_indices.push_back(count_add_var);

            count_add_var++;
            count_add_eq++;
        }

        if (std::isnan(LB[i]))
        {
            map_sol(i, count_add_var) = -sign_sol(i);
            A_final.col(count_add_var) = -A_final.col(i);
            cost_vector_final(count_add_var) = -cost_vector(i);
            no_slack_indices.push_back(count_add_var);
            count_add_var++;
        }

        map_sol(i, i) = sign_sol(i);
    }

    for (unsigned int i = 0; i < m; i++)
    {
        switch (constraint_types[i])
        {
        case ConstraintType::LE:
        case ConstraintType::GE: {
            A_final(i, count_add_var) = 1.0;
            slack_indices.push_back(count_add_var);
            count_add_var++;
        }
        break;
        case ConstraintType::EQ:
            break;
        }
    }

    if (slack_indices.size() != m_final)
        result.exit_value = PrimalResult::ExitValue::NOSLACK;
}
//**********************************************************
void Simplex::solve_primal_simplex(std::vector<unsigned int> &basis,
                                   std::vector<unsigned int> &no_basis,
                                   const unsigned int max_iteration,
                                   const double &tolerance)
{
    if (basis.size() != m_final && basis.size() + no_basis.size() != n_final)
        throw std::runtime_error("not valid basis");

    Eigen::MatrixXd AB = A_final(Eigen::all, basis);
    Eigen::MatrixXd ANB = A_final(Eigen::all, no_basis);

    auto AB_lu = AB.lu();
    Eigen::MatrixXd alpha = -AB_lu.solve(ANB);
    Eigen::VectorXd beta = AB_lu.solve(b_final);
    Eigen::VectorXd r = cost_vector_final(no_basis) + alpha.transpose() * cost_vector_final(basis);

    unsigned int q = std::numeric_limits<unsigned int>::max();
    bool new_iteration = false;
    for (unsigned int i = 0; i < no_basis.size(); i++)
    {
        if (r(i) > tolerance)
        {
            q = i;
            new_iteration = true;
            break;
        }
    }

    unsigned int it = 0;
    while (new_iteration)
    {
        if (it >= max_iteration)
        {
            result.exit_value = PrimalResult::ExitValue::MAXIT;
            return;
        }

        unsigned int p = std::numeric_limits<unsigned int>::max();
        double min_beta_alpha = std::numeric_limits<double>::max();

        for (unsigned int i = 0; i < basis.size(); i++)
        {
            if (alpha(i, q) >= -tolerance)
                continue;

            const double ratio_alpha_beta = -beta(i) / alpha(i, q);
            if (min_beta_alpha > ratio_alpha_beta)
            {
                min_beta_alpha = ratio_alpha_beta;
                p = i;
            }
        }

        if (p == std::numeric_limits<unsigned int>::max())
        {
            result.exit_value = PrimalResult::ExitValue::UNBOUNDED;
            return;
        }

        const unsigned int tmp = basis[p];
        basis[p] = no_basis[q];
        no_basis[q] = tmp;

        AB = A_final(Eigen::all, basis);
        ANB = A_final(Eigen::all, no_basis);

        AB_lu = AB.lu();
        alpha = -AB_lu.solve(ANB);
        beta = AB_lu.solve(b_final);
        r = cost_vector_final(no_basis) + alpha.transpose() * cost_vector_final(basis);

        q = std::numeric_limits<unsigned int>::max();
        new_iteration = false;
        for (unsigned int i = 0; i < no_basis.size(); i++)
        {
            if (r(i) > tolerance)
            {
                q = i;
                new_iteration = true;
                break;
            }
        }

        it++;
    }

    result.solution_final = Eigen::VectorXd::Zero(n_final);
    result.solution_final(basis) = beta;

    result.solution = map_sol * result.solution_final + trasl_sol;
    result.cost_value = cost_vector_final.transpose() * result.solution_final;

    result.exit_value = PrimalResult::ExitValue::SUCCESS;
}
//**********************************************************
} // namespace LPUtilties
} // namespace Gedim
