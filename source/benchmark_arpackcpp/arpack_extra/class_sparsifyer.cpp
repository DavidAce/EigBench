//
// Created by david on 2018-10-25.
//

#include "class_sparsifyer.h"
#include <Eigen/Core>
#include <Eigen/Sparse>


template<typename Scalar>
class_sparsifyer<Scalar>::class_sparsifyer(const Scalar *data, int L, double sparcity_threshold, bool prune) {
    auto matrix_dense  = Eigen::Map<Eigen::Matrix<const Scalar,Eigen::Dynamic, Eigen::Dynamic>> (data, L,L);
    Eigen::SparseMatrix<Scalar> matrix_sparse = matrix_dense.sparseView().pruned(sparcity_threshold);
    matrix_sparse.makeCompressed();
    sparcity = (double) (matrix_dense.array().cwiseAbs() > sparcity_threshold )
            .select(Eigen::MatrixXd::Ones(L,L),0).sum() / matrix_dense.size();

    nx       = L;
    n        = L*L;
    nnz      = (int) matrix_sparse.nonZeros();
    irow     = std::vector<int>(matrix_sparse.innerNonZeroPtr() , matrix_sparse.innerNonZeroPtr() + matrix_sparse.innerSize());
    pcol     = std::vector<int>(matrix_sparse.outerIndexPtr()   , matrix_sparse.outerIndexPtr()   + matrix_sparse.outerSize());
    valA     = std::vector<Scalar>(matrix_sparse.valuePtr(), matrix_sparse.valuePtr() + matrix_sparse.nonZeros());

}

//template<>
//class_sparsifyer<double>::class_sparsifyer(double *data, int L, double sparcity_threshold, bool prune){
//    auto matrix_dense  = Eigen::Map<Eigen::MatrixXd> (data, L,L);
//    Eigen::SparseMatrix<double> matrix_sparse = matrix_dense.sparseView().pruned(sparcity_threshold);
//    matrix_sparse.makeCompressed();
//    sparcity = (double) (matrix_dense.array().cwiseAbs() > sparcity_threshold )
//            .select(Eigen::MatrixXd::Ones(L,L),0).sum() / matrix_dense.size();
//
//    L       = L;
//    L        = L*L;
//    nnz      = (int) matrix_sparse.nonZeros();
//    irow     = std::vector<int>(matrix_sparse.innerNonZeroPtr() , matrix_sparse.innerNonZeroPtr() + matrix_sparse.innerSize());
//    pcol     = std::vector<int>(matrix_sparse.outerIndexPtr()   , matrix_sparse.outerIndexPtr()   + matrix_sparse.outerSize());
//    valA     = std::vector<double>(matrix_sparse.valuePtr(), matrix_sparse.valuePtr() + matrix_sparse.nonZeros());
//}
//
//
//
//
//template<>
//class_sparsifyer<std::complex<double>>::class_sparsifyer(std::complex<double> *data, int L, double sparcity_threshold, bool prune){
//    auto matrix_dense  = Eigen::Map<Eigen::MatrixXcd> (data, L,L);
//    Eigen::SparseMatrix<std::complex<double>> matrix_sparse = matrix_dense.sparseView().pruned(sparcity_threshold);
//    matrix_sparse.makeCompressed();
//    sparcity = (double) (matrix_dense.array().cwiseAbs() > sparcity_threshold )
//            .select(Eigen::MatrixXd::Ones(L,L),0).sum() / matrix_dense.size();
//
//    L       = L;
//    L        = L*L;
//    nnz      = (int) matrix_sparse.nonZeros();
//    irow     = std::vector<int>(matrix_sparse.innerNonZeroPtr() , matrix_sparse.innerNonZeroPtr() + matrix_sparse.innerSize());
//    pcol     = std::vector<int>(matrix_sparse.outerIndexPtr()   , matrix_sparse.outerIndexPtr()   + matrix_sparse.outerSize());
//    valA     = std::vector<std::complex<double>>(matrix_sparse.valuePtr(), matrix_sparse.valuePtr() + matrix_sparse.nonZeros());
//}
//


