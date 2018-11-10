//
// Created by david on 2018-05-08.
//

#include "matrix_product_dense.h"
#define profile_matrix_product_dense 0

using namespace arpack_extra::modes;
//
//template<typename Scalar>
//void DenseMatrixProduct<Scalar>::print() const {
//    std::cout << "A_matrix: \n" << A_matrix << std::endl;
//}
//
//
//template<typename Scalar>
//void DenseMatrixProduct<Scalar>::FactorOP()
//
///*  Full pivot LU decomposition
// *  Factors (A-shift*L) = P^-1 LU Q^-1
// */
//{
//    assert(readyShift and "Shift value sigma has not been set.");
//    Scalar sigma;
//    if constexpr(std::is_same<Scalar,double>::value)
//    {sigma = sigmaR;}
//    else
//    {sigma = std::complex<double>(sigmaR,sigmaI);}
//    lu.compute(A_matrix - sigma * MatrixType::Identity(L,L));
//    readyFactorOp = true;
//}
//
//
//
//
//template<typename Scalar>
//void DenseMatrixProduct<Scalar>::MultOPv(Scalar* x_in_ptr, Scalar* x_out_ptr) {
//    assert(readyFactorOp and "FactorOp() has not been run yet.");
//    using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
//    Eigen::Map<VectorType>       x_in    (x_in_ptr,L);
//    Eigen::Map<VectorType>       x_out   (x_out_ptr,L);
//
//    switch (side){
//        case Side::R: {
//            x_out = lu.solve(x_in);
//            break;
//        }
//        case Side::L: {
//            x_out = x_in * lu.inverse();
//        }
//    }
//    counter++;
//}
//
//
//
//
//template<typename Scalar>
//void DenseMatrixProduct<Scalar>::MultAx(Scalar* x_in, Scalar* x_out) {
//    switch (form){
//        case Form::NONSYMMETRIC:
//          switch (side) {
//              case Side::R: {
//                  Eigen::Map<VectorType> x_vec_in (x_in,  L);
//                  Eigen::Map<VectorType> x_vec_out(x_out, L);
//                  x_vec_out.noalias() = A_matrix * x_vec_in ;
//                  break;
//              }
//              case Side::L: {
//                  Eigen::Map<VectorTypeT> x_vec_in(x_in, L);
//                  Eigen::Map<VectorTypeT> x_vec_out(x_out, L);
//                  x_vec_out.noalias() = x_vec_in * A_matrix;
//                  break;
//              }
//          }
//          break;
//        case Form::SYMMETRIC: {
//            Eigen::Map<VectorType> x_vec_in(x_in, L);
//            Eigen::Map<VectorType> x_vec_out(x_out, L);
//            x_vec_out.noalias() = A_matrix.template selfadjointView<Eigen::Upper>() * x_vec_in;
//            break;
//        }
//    }
//    counter++;
//}
//
//template class DenseMatrixProduct<std::complex<double>>;
//template class DenseMatrixProduct<double>;
//
//
//
//

