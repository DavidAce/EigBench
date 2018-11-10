//
// Created by david on 2018-05-08.
//

#ifndef MATRIX_PRODUCT_DENSE_H
#define MATRIX_PRODUCT_DENSE_H

#ifdef EIGEN_USE_BLAS
#define EIGEN_USE_BLAS_SUSPEND
#undef EIGEN_USE_BLAS
#endif


#include "nmspc_arpack_extra.h"
#include <iostream>
#include <Eigen/Core>
#include <Eigen/LU>


template <typename Scalar_>
class DenseMatrixProduct {
public:
    using Scalar      = Scalar_;
    using MatrixType  = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic, Eigen::ColMajor>;
    using VectorType  = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
    using VectorTypeT = Eigen::Matrix<Scalar,1,Eigen::Dynamic>;
private:

    const MatrixType A_matrix;           // The actual matrix. Given matrices will be copied into this one.
    const int L;                        // The linear matrix dimension
    arpack_extra::modes::Form form;     // Chooses SYMMETRIC / NONSYMMETRIC mode
    arpack_extra::modes::Side side;     // Chooses whether to find (R)ight or (L)eft eigenvectors

    // Shift-invert mode stuff
    Eigen::PartialPivLU<MatrixType> lu;                         // Object for dense LU decomposition used in shift-invert mode
    double sigmaR = std::numeric_limits<double>::quiet_NaN();   // The real part of the shift
    double sigmaI = std::numeric_limits<double>::quiet_NaN();   // The imag part of the shift
    bool readyFactorOp = false;                                 // Flag to make sure LU factorization has occurred
    bool readyShift = false;                                    // Flag to make sure

public:
    // Pointer to data constructor, copies the matrix into an internal Eigen matrix.
    DenseMatrixProduct(
            const Scalar * const A_,
            const int L_,
            const arpack_extra::modes::Form form_ = arpack_extra::modes::Form::NONSYMMETRIC,
            const arpack_extra::modes::Side side_ = arpack_extra::modes::Side::R

    ): A_matrix(Eigen::Map<const MatrixType>(A_,L_,L_)),
       L(L_), form(form_), side(side_) {}

    // Eigen type constructor. Pass any copy-assignable eigen type into an internal Eigen matrix.
    template<typename Derived>
    explicit DenseMatrixProduct(
            const Eigen::EigenBase<Derived> &matrix_,
            const arpack_extra::modes::Form form_ = arpack_extra::modes::Form::NONSYMMETRIC,
            const arpack_extra::modes::Side side_ = arpack_extra::modes::Side::R)
            : A_matrix(matrix_), L(A_matrix.rows()), form(form_), side(side_)
    {}

    // Functions used in in Arpack++ solver
    int rows() const {return L;};
    int cols() const {return L;};
    void FactorOP();                                      //  Factors (A-sigma*I) into PLU
    void MultOPv(Scalar* x_in_ptr, Scalar* x_out_ptr);    //   Computes the matrix-vector product x_out <- inv(A-sigma*I)*x_in.
    void MultAx (Scalar* x_in_ptr, Scalar* x_out_ptr);    //   Computes the matrix-vector multiplication x_out <- A*x_in.

    // Various utility functions
    int counter = 0;
    void print()const;
    void set_shift(std::complex<double> sigma_)   {sigmaR=std::real(sigma_);sigmaI=std::imag(sigma_) ;readyShift = true;}
    void set_shift(double               sigma_)   {sigmaR=sigma_, sigmaI = 0.0;readyShift = true;}
    void set_shift(double sigmaR_, double sigmaI_){sigmaR=sigmaR_;sigmaI = sigmaI_ ;readyShift = true;}
    void set_mode(const arpack_extra::modes::Form form_){form = form_;}
    void set_side(const arpack_extra::modes::Side side_){side = side_;}
    const MatrixType & get_matrix()const{return A_matrix;}
    const arpack_extra::modes::Form &get_form()const{return form;}
    const arpack_extra::modes::Side &get_side()const{return side;}
};




// Function definitions



template<typename Scalar>
void DenseMatrixProduct<Scalar>::print() const {
    std::cout << "A_matrix: \n" << A_matrix << std::endl;
}


template<typename Scalar>
void DenseMatrixProduct<Scalar>::FactorOP()

/*  Partial pivot LU decomposition
 *  Factors P(A-sigma*I) = LU
 */
{
    assert(readyShift and "Shift value sigma has not been set.");
    Scalar sigma;
    if constexpr(std::is_same<Scalar,double>::value)
    {sigma = sigmaR;}
    else
    {sigma = std::complex<double>(sigmaR,sigmaI);}
    lu.compute(A_matrix - sigma * MatrixType::Identity(L,L));
    readyFactorOp = true;
}




template<typename Scalar>
void DenseMatrixProduct<Scalar>::MultOPv(Scalar* x_in_ptr, Scalar* x_out_ptr) {
    using namespace arpack_extra::modes;
    assert(readyFactorOp and "FactorOp() has not been run yet.");
    switch (side){
        case Side::R: {
            Eigen::Map<VectorType>       x_in    (x_in_ptr,L);
            Eigen::Map<VectorType>       x_out   (x_out_ptr,L);
            x_out.noalias() = lu.solve(x_in);
            break;
        }
        case Side::L: {
            Eigen::Map<VectorTypeT>       x_in    (x_in_ptr,L);
            Eigen::Map<VectorTypeT>       x_out   (x_out_ptr,L);
            x_out.noalias() = x_in * lu.inverse();
            break;
        }
    }
    counter++;
}




template<typename Scalar>
void DenseMatrixProduct<Scalar>::MultAx(Scalar* x_in, Scalar* x_out) {
    using namespace arpack_extra::modes;
    switch (form){
        case Form::NONSYMMETRIC:
            switch (side) {
                case Side::R: {
                    Eigen::Map<VectorType> x_vec_in (x_in,  L);
                    Eigen::Map<VectorType> x_vec_out(x_out, L);
                    x_vec_out.noalias() = A_matrix * x_vec_in ;
                    break;
                }
                case Side::L: {
                    Eigen::Map<VectorTypeT> x_vec_in(x_in, L);
                    Eigen::Map<VectorTypeT> x_vec_out(x_out, L);
                    x_vec_out.noalias() = x_vec_in * A_matrix;
                    break;
                }
            }
            break;
        case Form::SYMMETRIC: {
            Eigen::Map<VectorType> x_vec_in(x_in, L);
            Eigen::Map<VectorType> x_vec_out(x_out, L);
            x_vec_out.noalias() = A_matrix.template selfadjointView<Eigen::Upper>() * x_vec_in;
            break;
        }
    }
    counter++;
}


#ifdef EIGEN_USE_BLAS_SUSPEND
#define EIGEN_USE_BLAS
#undef EIGEN_USE_BLAS_SUSPEND
#endif



#endif