#ifndef __SUITESPARSE_UTILITIES_H
#define __SUITESPARSE_UTILITIES_H

#include <iostream>
#include <string>

#include "Eigen/Eigen"
#include "Eigen_Utilities.hpp"
#include "Gedim_Macro.hpp"

#if ENABLE_SUITESPARSE == 1
#include <cholmod.h>
#include <klu.h>
#endif

namespace Gedim
{
#if ENABLE_SUITESPARSE == 1
struct SuiteSparse_Utilities
{
    template <typename Derived> static double condest(const Eigen::SparseMatrix<Derived> &sparseA)
    {

        const std::vector<Eigen::Triplet<Derived>> Aentries = Gedim::Eigen_Utilities::SparseMatrixToTriplets<Derived>(sparseA);
        const unsigned int N = sparseA.rows();

        cholmod_common CholCommon, *cc;
        cc = &CholCommon;
        cholmod_l_start(cc);

        cholmod_triplet *Gct = cholmod_l_allocate_triplet(N,
                                                          N,
                                                          Aentries.size(),
                                                          0, // stype: both upper and lower are stored
                                                          CHOLMOD_REAL,
                                                          cc);
        for (size_t i = 0; i < Aentries.size(); ++i)
        {
            reinterpret_cast<long *>(Gct->i)[i] = Aentries[i].row();
            reinterpret_cast<long *>(Gct->j)[i] = Aentries[i].col();
            reinterpret_cast<double *>(Gct->x)[i] = Aentries[i].value();
        }
        Gct->nnz = Aentries.size();

        // convert triplet matrix to sparse
        cholmod_sparse *G = cholmod_l_triplet_to_sparse(Gct, Aentries.size(), cc);

        klu_l_common Common;
        klu_l_defaults(&Common);
        long *Ap = reinterpret_cast<long *>(G->p);
        long *Ai = reinterpret_cast<long *>(G->i);
        double *Ax = reinterpret_cast<double *>(G->x);
        klu_l_symbolic *Symbolic = klu_l_analyze(N, Ap, Ai, &Common);
        klu_l_numeric *Numeric = klu_l_factor(Ap, Ai, Ax, Symbolic, &Common);

        /* Computes a reasonably accurate estimate of the 1-norm condition number, using
         * Hager’s method, as modified by Higham and Tisseur (same method as used in
         * MATLAB’s condest */
        int info = klu_l_condest(Ap, Ax, Symbolic, Numeric, &Common);

        if (!info)
        {
            cholmod_l_free_sparse(&G, &CholCommon);
            cholmod_l_free_triplet(&Gct, &CholCommon);
            klu_l_free_numeric(&Numeric, &Common);
            klu_l_free_symbolic(&Symbolic, &Common);
            throw std::runtime_error("Suitesparse fails.");
        }
        else
        {
            double condest = Common.condest;

            cholmod_free_sparse(&G, &CholCommon);
            cholmod_l_free_triplet(&Gct, &CholCommon);
            klu_l_free_numeric(&Numeric, &Common);
            klu_l_free_symbolic(&Symbolic, &Common);
            return condest;
        }
    }
};
#endif
} // namespace Gedim
#endif // __SUITESPARSE_UTILITIES_H
