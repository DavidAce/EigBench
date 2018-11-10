//
// Created by david on 2018-10-29.
//

#ifndef EIGBENCH_CLASS_EIGSOLVER_ARPACK_2_H
#define EIGBENCH_CLASS_EIGSOLVER_ARPACK_2_H

#include "arpack_extra/matrix_product_dense.h"
#include "arpack_extra/matrix_product_sparse.h"
#include "arpack_extra/arpack_dense_matProd.h"
#include "arpack_extra/arpack_sparse_matProd.h"
#include "arpack_extra/arpack_custom_op.h"
#include "arpack_extra/nmspc_arpack_extra.h"
#include "arpack_extra/matrix_recast.h"

class class_eigsolver {

//private:
//    template<typename Scalar>
//    void eigs

private:


public:

    arpack_extra::SolverConf solverConf;
    arpack_extra::Solution   solution;

    void conf_init(const int L,
                   const int nev,
                   const int ncv,
                   const std::complex<double> sigma            = std::numeric_limits<std::complex<double>>::quiet_NaN(),
                   const arpack_extra::modes::Type type        = arpack_extra::modes::Type::REAL,
                   const arpack_extra::modes::Form form        = arpack_extra::modes::Form::NONSYMMETRIC,
                   const arpack_extra::modes::Ritz ritz        = arpack_extra::modes::Ritz::SR,
                   const arpack_extra::modes::Side side        = arpack_extra::modes::Side::R,
                   const arpack_extra::modes::Storage storage  = arpack_extra::modes::Storage::DENSE,
                   const bool compute_eigvecs_           = false,
                   const bool remove_phase_              = false
                   );


    template<typename Scalar>
    void eigs_auto(const Scalar *matrix_data,
                   const int L,
                   const int nev,
                   const bool compute_eigvecs_          = false,
                   const arpack_extra::modes::Ritz ritz = arpack_extra::modes::Ritz::SR,
                   const std::complex<double> sigma     = std::numeric_limits<std::complex<double>>::quiet_NaN(),
                   const arpack_extra::modes::Side side = arpack_extra::modes::Side::R,
                   const bool remove_phase_             = false,
                   Scalar *residual_                    = nullptr);

    template<typename Scalar>
    void eigs_dense(const  Scalar *matrix,
                   const int L,
                   const int nev,
                   const int ncv,
                   const std::complex<double> sigma     = std::numeric_limits<std::complex<double>>::quiet_NaN(),
                   const arpack_extra::modes::Form form = arpack_extra::modes::Form::NONSYMMETRIC,
                   const arpack_extra::modes::Ritz ritz = arpack_extra::modes::Ritz::SR,
                   const arpack_extra::modes::Side side = arpack_extra::modes::Side::R,
                   const bool compute_eigvecs_          = false,
                   const bool remove_phase_             = false,
                   Scalar *residual_                    = nullptr);

    template<typename Scalar>
    void eigs_sparse(const Scalar *matrix,
                     const int L,
                     const int nev,
                     const int ncv,
                     const std::complex<double> sigma     = std::numeric_limits<std::complex<double>>::quiet_NaN(),
                     const arpack_extra::modes::Form form = arpack_extra::modes::Form::NONSYMMETRIC,
                     const arpack_extra::modes::Ritz ritz = arpack_extra::modes::Ritz::SR,
                     const arpack_extra::modes::Side side = arpack_extra::modes::Side::R,
                     const bool compute_eigvecs_          = false,
                     const bool remove_phase_             = false,
                     Scalar *residual_                    = nullptr);

};


template<typename Scalar>
void class_eigsolver::eigs_dense   (const Scalar *matrix,
                                   const int L,
                                   const int nev,
                                   const int ncv,
                                   const std::complex<double> sigma,
                                   const arpack_extra::modes::Form form,
                                   const arpack_extra::modes::Ritz ritz,
                                   const arpack_extra::modes::Side side,
                                   const bool compute_eigvecs_,
                                   const bool remove_phase_,
                                   Scalar *residual_)
{
    using namespace arpack_extra::modes;
    bool is_cplx = std::is_same<std::complex<double>,Scalar>::value;
    Type type = is_cplx ? Type::CPLX : Type::REAL;
    Storage storage = Storage::DENSE;
    conf_init(L, nev, ncv, sigma,type,form,ritz,side,storage,compute_eigvecs_,remove_phase_);
    auto matrix_dense = DenseMatrixProduct<Scalar> (matrix,L);
    arpack_dense_matProd<Scalar> arpacksolver(matrix_dense, solverConf, solution);
    arpacksolver.eigs();
}


template<typename Scalar>
void class_eigsolver::eigs_sparse   (const Scalar *matrix,
                                    const int L,
                                    const int nev,
                                    const int ncv,
                                    const std::complex<double> sigma,
                                    const arpack_extra::modes::Form form,
                                    const arpack_extra::modes::Ritz ritz,
                                    const arpack_extra::modes::Side side,
                                    const bool compute_eigvecs_,
                                    const bool remove_phase_,
                                    Scalar *residual_)
{
    using namespace arpack_extra::modes;
    bool is_cplx = std::is_same<std::complex<double>,Scalar>::value;
    Type type = is_cplx ? Type::CPLX : Type::REAL;
    Storage storage = Storage::SPARSE;
    conf_init(L, nev, ncv, sigma,type,form,ritz,side,storage,compute_eigvecs_,remove_phase_);
    auto matrix_sparse = SparseMatrixProduct<Scalar> (matrix,L);
    arpack_sparse_matProd<Scalar> arpacksolver(matrix_sparse, solverConf, solution);
    arpacksolver.eigs();
}




template<typename Scalar>
void class_eigsolver::eigs_auto   (const Scalar *matrix_data,
                                   const int L,
                                   const int nev,
                                   const bool compute_eigvecs_          ,
                                   const arpack_extra::modes::Ritz ritz ,
                                   const std::complex<double> sigma     ,
                                   const arpack_extra::modes::Side side ,
                                   const bool remove_phase_             ,
                                   Scalar *residual_
                                   )
{
    using namespace arpack_extra::modes;
    matrix_recast<Scalar> matRecast(matrix_data,L);
    bool is_sparse    = matRecast.is_sparse();
    bool is_real      = matRecast.is_real();
    bool is_symmetric = matRecast.is_symmetric();

    Form form        = is_symmetric ? Form::SYMMETRIC : Form::NONSYMMETRIC;
    Type type        = is_real      ? Type::REAL      : Type ::CPLX;
    Storage storage  = is_sparse    ? Storage::SPARSE : Storage::DENSE;

    conf_init(L, nev, -1, sigma,type, form,ritz,side,storage,compute_eigvecs_,remove_phase_);

    using namespace arpack_extra;

    if(is_real) {
        if(is_sparse) {
            auto matrix = matRecast.get_as_real_sparse();
            arpack_custom_op<SparseMatrixProduct<double>> arpacksolver(matrix, solverConf, solution);
            arpacksolver.eigs();
        }else {
            auto matrix = matRecast.get_as_real_dense();
            arpack_custom_op<DenseMatrixProduct<double>> arpacksolver(matrix, solverConf, solution);
            arpacksolver.eigs();
        }
    }else {
        if(is_sparse) {
            auto matrix = matRecast.get_as_cplx_sparse();
            arpack_custom_op<SparseMatrixProduct<std::complex<double>>> arpacksolver(matrix, solverConf, solution);
            arpacksolver.eigs();
        }else {
            auto matrix = matRecast.get_as_cplx_dense();
            arpack_custom_op<DenseMatrixProduct<std::complex<double>>> arpacksolver(matrix, solverConf, solution);
            arpacksolver.eigs();
        }
    }

}

















#endif //EIGBENCH_CLASS_EIGSOLVER_ARPACK_2_H
