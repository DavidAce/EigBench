//
// Created by david on 2018-10-30.
//

#include <algorithm>
#include "arpack_base.h"
#include "matrix_product_dense.h"
#include "matrix_product_sparse.h"




//template<typename Scalar>
//void arpack_base<Scalar>::shift_invert_eigvals(Scalar sigma) {
//    if (solution.meta.eigvals_found){
//        std::transform(solution.eigvals.begin(), solution.eigvals.end(), solution.eigvals.begin(),
//                       [sigma](Scalar num) -> Scalar
//                       { return 1.0/(num - sigma); });
//    }else{
//        std::cerr << "Eigenvalues haven't been computed yet. Can't invert. Exiting " << std::endl;
//    }
//
//}



template<typename MatrixType>
arpack_base<MatrixType>::arpack_base(
        MatrixType               &matrix_,
        arpack_extra::SolverConf &solverConf_,
        arpack_extra::Solution   &solution_)
        :
        matrix(matrix_),
        solverConf(solverConf_),
        solution  (solution_)
{

    t_sol.set_properties(profile_arpack, 10,"Time iterating  ");
    t_get.set_properties(profile_arpack, 10,"Time getting sol");
    t_sub.set_properties(profile_arpack, 10,"Time subtracting");
    t_all.set_properties(profile_arpack, 10,"Time doing all  ");
}




template<typename MatrixType>
void arpack_base<MatrixType>::subtract_phase() {

    if constexpr (std::is_same<Scalar, std::complex<double>>::value) {
        if (solution.meta.eigvecs_found){
            using namespace std::complex_literals;
            for (int i = 0; i < solution.meta.nev_found; i++) {
                auto begin = solution.eigvecs.begin() + i * solution.meta.rows;
                auto end = begin + solution.meta.rows;
                Scalar inv_phase = -1.0i * std::arg(solution.eigvecs[i * solution.meta.rows]);
                Scalar exp_inv_phase = std::exp(inv_phase);
                std::transform(begin, end, begin,
                               [exp_inv_phase](std::complex<double> num) -> std::complex<double>
                               { return (num * exp_inv_phase); });
            }
        }else{
            std::cerr << "Eigenvalues haven't been computed yet. Can't subtract phase. Exiting " << std::endl;
        }

    }
}




template class arpack_base <arpack_extra::DenseType<double>>;
template class arpack_base <arpack_extra::SparseType<double>>;
template class arpack_base <arpack_extra::DenseType<std::complex<double>>>;
template class arpack_base <arpack_extra::SparseType<std::complex<double>>>;
template class arpack_base <DenseMatrixProduct<double>>;
template class arpack_base <DenseMatrixProduct<std::complex<double>>>;
template class arpack_base <SparseMatrixProduct<double>>;
template class arpack_base <SparseMatrixProduct<std::complex<double>>>;
