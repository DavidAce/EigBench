//
// Created by david on 2018-11-02.
//

#ifndef ARPACK_SPARSE_MATPROD_H
#define ARPACK_SPARSE_MATPROD_H

#include "arpack_base.h"
#include "matrix_product_sparse.h"
#include "nmspc_arpack_extra.h"

template<typename Scalar>
class arpack_sparse_matProd : public arpack_base<SparseMatrixProduct<Scalar>> {
private:
    int nev_internal;
    int ncv_internal;


    void eigs_sym();
    void eigs_nsym();
    void eigs_comp();


public:
    using MatrixType = SparseMatrixProduct<Scalar>;
    using arpack_base<MatrixType>::arpack_base;
    using arpack_base<MatrixType>::find_solution;
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



#endif //EIGBENCH_ARPACK_SPARSE_ALT_H
