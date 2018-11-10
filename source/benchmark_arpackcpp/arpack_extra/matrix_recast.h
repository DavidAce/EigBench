//
// Created by david on 2018-10-29.
//

#ifndef MATRIX_RECAST_H
#define MATRIX_RECAST_H

#include <complex>
#include <vector>
#include "nmspc_arpack_extra.h"
#include "matrix_product_dense.h"
#include "matrix_product_sparse.h"


template<typename Scalar>
class matrix_recast {
private:


    bool isSparse;
    bool isReal;
    bool isHermitian;
    double sparcity;

    const Scalar *matrix_ptr;
    int     L;

    void check_if_real();
    void check_if_sparse();
    void check_if_hermitian();

//    arpack_extra::DenseType<double>                matrix_real_dense;
//    arpack_extra::DenseType<std::complex<double>>  matrix_cplx_dense;
//    arpack_extra::SparseType<double>               matrix_real_sparse;
//    arpack_extra::SparseType<std::complex<double>> matrix_cplx_sparse;
//    DenseMatrixProduct<double>                matrix_real_dense;
//    DenseMatrixProduct<std::complex<double>>  matrix_cplx_dense;
//    SparseMatrixProduct<double>               matrix_real_sparse;
//    SparseMatrixProduct<std::complex<double>> matrix_cplx_sparse;


public:
    matrix_recast(const Scalar *matrix_ptr_, int L_);

    DenseMatrixProduct<double>               get_as_real_dense();
    DenseMatrixProduct<std::complex<double>> get_as_cplx_dense();
    SparseMatrixProduct<double>               get_as_real_sparse();
    SparseMatrixProduct<std::complex<double>> get_as_cplx_sparse();

//    void convert_to_real_dense();
//    void convert_to_cplx_dense();
//    void convert_to_real_sparse();
//    void convert_to_cplx_sparse();

    bool is_sparse      ()   {return isSparse;}
    bool is_real        ()   {return isReal;}
    bool is_symmetric()   {return isHermitian;}
};


#endif //MATRIX_RECAST_H
