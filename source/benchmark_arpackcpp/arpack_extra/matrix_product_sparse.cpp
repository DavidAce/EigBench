//
// Created by david on 2018-05-08.
//

#include "matrix_product_sparse.h"
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using namespace arpack_extra::modes;

//Eigen::Map< SparseMatrixType >::Map	(
//       	Index 	rows,
//          Index 	cols,
//          Index 	nnz,
//          StorageIndex * 	outerIndexPtr,
//          StorageIndex * 	innerIndexPtr,
//          Scalar       * 	valuePtr,
//          StorageIndex * 	innerNonZerosPtr = 0
//)
//
//template<typename Scalar>
//SparseMatrixProduct<Scalar>::SparseMatrixProduct(const Scalar *const A_, const int L_,
//                                                 const arpack_extra::modes::Form form_,
//                                                 const arpack_extra::modes::Side side_)
//                                                 : L(L_), form(form_), side(side_)
//{
//
//}
//
//template<typename Scalar>
//SparseMatrixProduct<Scalar>::SparseMatrixProduct(const Eigen::EigenBase<Derived> &matrix_sparse,
//                                                 const arpack_extra::modes::Form form_,
//                                                 const arpack_extra::modes::Side side_)
//        : A_matrix(matrix_sparse), L(A_matrix.rows()), form(form_), side(side_)
//{
//    A_matrix.makeCompressed();
//}

//
//template<typename Scalar>
//void SparseMatrixProduct<Scalar>::clear() {
//}
//
//template<typename Scalar>
//void SparseMatrixProduct<Scalar>::print() const {
//    std::cout << "A_matrix: \n" << A_matrix << std::endl;
//}
//
//
//template<typename Scalar>
//void SparseMatrixProduct<Scalar>::FactorOP()
//
///*  Full pivot LU decomposition
// *  Factors (A-shift*L) = P^-1 LU Q^-1
// */
//
//{
//    assert(readyShift and "Shift value sigma has not been set.");
//
//
//
//
//    Scalar sigma;
//    if constexpr(std::is_same<Scalar,double>::value)
//    {sigma = sigmaR;}
//    else
//    {sigma = std::complex<double>(sigmaR,sigmaI);}
//    Eigen::SparseMatrix<Scalar> Id(L,L);
//    Id.setIdentity();
//    Id.makeCompressed();
//    Eigen::SparseMatrix<Scalar> A_shift = (A_matrix - sigma * Id);
//
//    lu.analyzePattern(A_shift);
//    lu.factorize(A_shift);
//    readyFactorOp = true;
//}
//
//template<typename Scalar>
//void SparseMatrixProduct<Scalar>::MultOPv(Scalar* x_in_ptr, Scalar* x_out_ptr) {
//    assert(readyFactorOp and "FactorOp() has not been run yet.");
//    using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
//
//    Eigen::Map<VectorType>       x_in    (x_in_ptr,L);
//    Eigen::Map<VectorType>       x_out   (x_out_ptr,L);
//
//    switch (side){
//        case Side::R: {
//            x_out = lu.solve(x_in);
//            break;
//        }
//        case Side::L: {
//            std::cerr << "Left sided shift invert hasn't been implemented yet.." << std::endl;
//            exit(1);
//            break;
//        }
//    }
//    counter++;
//
//
//}
//
//template<typename Scalar>
//void SparseMatrixProduct<Scalar>::MultAx(Scalar* x_in, Scalar* x_out) {
//    switch (form){
//        case Form::NONSYMMETRIC:
//            switch (side) {
//                case Side::R: {
//                    Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic,1>> x_vec_in (x_in,  L);
//                    Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic,1>> x_vec_out(x_out, L);
//                    x_vec_out.noalias() = A_matrix * x_vec_in ;
//
//                    break;
//                }
//                case Side::L: {
//                    Eigen::Map<Eigen::Matrix<Scalar, 1, Eigen::Dynamic>> x_vec_in(x_in, L);
//                    Eigen::Map<Eigen::Matrix<Scalar, 1, Eigen::Dynamic>> x_vec_out(x_out, L);
//                    x_vec_out.noalias() = x_vec_in * A_matrix;
//                    break;
//                }
//            }
//            break;
//        case Form::SYMMETRIC: {
//            Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> x_vec_in(x_in, L);
//            Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> x_vec_out(x_out, L);
//            x_vec_out.noalias() = A_matrix.template selfadjointView<Eigen::Upper>() * x_vec_in;
//            break;
//        }
//    }
//    counter++;
//}

//
//template class SparseMatrixProduct<std::complex<double>>;
//template class SparseMatrixProduct<double>;
//
//


