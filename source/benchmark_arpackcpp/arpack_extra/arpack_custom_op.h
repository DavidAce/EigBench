//
// Created by david on 2018-11-02.
//

#ifndef ARPACK_CUSTOM_OP_H
#define ARPACK_CUSTOM_OP_H

#include "arpack_base.h"
#include "matrix_product_dense.h"
#include "matrix_product_sparse.h"
#include "nmspc_arpack_extra.h"

template<typename MatrixType>
class arpack_custom_op : public arpack_base<MatrixType> {
private:
    int nev_internal;
    int ncv_internal;


    void eigs_sym();
    void eigs_nsym();
    void eigs_comp();


public:
//    using MatrixType = DenseMatrixProduct<Scalar>;
    using Scalar = typename arpack_base<MatrixType>::Scalar;
    using arpack_base<MatrixType>::arpack_base;
    using arpack_base<MatrixType>::find_solution;
    using arpack_base<MatrixType>::copy_solution_symm;
    using arpack_base<MatrixType>::copy_solution_nsym;
    using arpack_base<MatrixType>::subtract_phase;
    using arpack_base<MatrixType>::solverConf;
    using arpack_base<MatrixType>::solution;
    using arpack_base<MatrixType>::matrix;
    using arpack_base<MatrixType>::t_all      ;
    using arpack_base<MatrixType>::t_get      ;
    using arpack_base<MatrixType>::t_sol      ;
    using arpack_base<MatrixType>::t_sub      ;
    void eigs();
};



#endif //EIGBENCH_ARPACK_DENSE_ALT_H
