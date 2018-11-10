//
// Created by david on 2018-05-08.
//

#ifndef MATRIX_PRODUCT_SPARSE_H
#define MATRIX_PRODUCT_SPARSE_H

#ifdef EIGEN_USE_BLAS
#define EIGEN_USE_BLAS_SUSPEND
#undef EIGEN_USE_BLAS
#endif


#include "nmspc_arpack_extra.h"
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#define profile_matrix_product_sparse 0


template <typename Scalar_>
class SparseMatrixProduct {
public:
    using Scalar      = Scalar_;
    using MatrixType  = Eigen::SparseMatrix<Scalar>;
    using DenseMatrixType = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
    using VectorType  = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
    using VectorTypeT = Eigen::Matrix<Scalar,1,Eigen::Dynamic>;
private:
    MatrixType A_matrix;          // The actual matrix. Given matrices will be copied into this one.
    const int L;                        // The linear matrix dimension
    arpack_extra::modes::Form form;     // Chooses SYMMETRIC / NONSYMMETRIC mode
    arpack_extra::modes::Side side;     // Chooses whether to find (R)ight or (L)eft eigenvectors

    // Shift-invert mode stuff
    Eigen::SparseLU<MatrixType> lu;                             // Object for sparse LU decomposition used in shift-invert mode
    double sigmaR = std::numeric_limits<double>::quiet_NaN();   // The real part of the shift
    double sigmaI = std::numeric_limits<double>::quiet_NaN();   // The imag part of the shift
    bool readyFactorOp = false;                                 // Flag to make sure LU factorization has occurred
    bool readyShift = false;

public:
    // Pointer to data constructor, copies the matrix into an internal Eigen matrix.
    SparseMatrixProduct(
            const Scalar * A_,
            const int L_,
            const arpack_extra::modes::Form form_ = arpack_extra::modes::Form::NONSYMMETRIC,
            const arpack_extra::modes::Side side_ = arpack_extra::modes::Side::R)
            : A_matrix(Eigen::Map<const DenseMatrixType>(A_,L_,L_).sparseView()), L(L_), form(form_), side(side_)
    {A_matrix.makeCompressed();}

    // Eigen type constructor. Pass any copy-assignable eigen type into an internal Eigen matrix.
    template<typename Derived>
    explicit SparseMatrixProduct(
            const Eigen::EigenBase<Derived> &matrix_,
            const arpack_extra::modes::Form form_ = arpack_extra::modes::Form::NONSYMMETRIC,
            const arpack_extra::modes::Side side_ = arpack_extra::modes::Side::R)
            : A_matrix(matrix_.derived().sparseView()), L(A_matrix.rows()), form(form_), side(side_)
    {A_matrix.makeCompressed();}

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
    const MatrixType & get_matrix()const {return A_matrix;}
    const arpack_extra::modes::Form &get_form()const{return form;}
    const arpack_extra::modes::Side &get_side()const{return side;}
};




// Function definitions



template<typename Scalar>
void SparseMatrixProduct<Scalar>::print() const {
    std::cout << "A_matrix: \n" << A_matrix << std::endl;
}


template<typename Scalar>
void SparseMatrixProduct<Scalar>::FactorOP()

/*  Sparse decomposition
 *  Factors P(A-sigma*I) = LU
 */

{
    assert(readyShift and "Shift value sigma has not been set.");
    Scalar sigma;
    if constexpr(std::is_same<Scalar,double>::value)
    {sigma = sigmaR;}
    else
    {sigma = std::complex<double>(sigmaR,sigmaI);}
    Eigen::SparseMatrix<Scalar> Id(L,L);
    Id.setIdentity();
    Id.makeCompressed();
    Eigen::SparseMatrix<Scalar> A_shift = (A_matrix - sigma * Id);
    lu.analyzePattern(A_shift);
    lu.factorize(A_shift);
    readyFactorOp = true;
}

template<typename Scalar>
void SparseMatrixProduct<Scalar>::MultOPv(Scalar* x_in_ptr, Scalar* x_out_ptr) {
    using namespace arpack_extra::modes;
    assert(readyFactorOp and "FactorOp() has not been run yet.");
    using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
    Eigen::Map<VectorType>       x_in    (x_in_ptr,L);
    Eigen::Map<VectorType>       x_out   (x_out_ptr,L);
    switch (side){
        case Side::R: {
            x_out.noalias() = lu.solve(x_in);
            break;
        }
        case Side::L: {
            std::cerr << "Left sided sparse shift invert hasn't been implemented yet..." << std::endl;
            exit(1);
            break;
        }
    }
    counter++;
}

template<typename Scalar>
void SparseMatrixProduct<Scalar>::MultAx(Scalar* x_in, Scalar* x_out) {
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
            Eigen::Map<VectorType> x_vec_in (x_in,  L);
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




#endif //
