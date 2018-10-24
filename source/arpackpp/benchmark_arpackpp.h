//
// Created by david on 2018-10-24.
//

#ifndef TRAINING_BENCHMARK_ARPACKPP_H
#define TRAINING_BENCHMARK_ARPACKPP_H
#include <Eigen/Core>
#include <complex>
#include <general/class_tic_toc.h>
#include <arpackpp/class_eigsolver_arpack.h>

class benchmark_arpackpp {
public:

    benchmark_arpackpp() = default;


    template<typename Scalar>
    Eigen::ArrayXXcd eigs_Arpackpp_dense(Eigen::MatrixXcd matrix, int nev, class_tic_toc &timer){
        class_eigsolver_arpack<Scalar,Form::SYMMETRIC> solver;
        int ncv = std::min((int)matrix.rows()/2, nev*8);
        Eigen::Matrix<Scalar,Eigen::Dynamic, Eigen::Dynamic> matrix_scalar;
        if constexpr (std::is_same<Scalar,std::complex<double>>::value){
            matrix_scalar = matrix;
        }
        else
        {
            matrix_scalar = matrix.real();
        }


        int L   = matrix_scalar.rows();
        timer.tic();
        solver.eig(matrix_scalar.data(), L , nev,ncv, Ritz::LM, Side::R, true, false);
        timer.toc();

        return Eigen::Map<const Eigen::Array<Scalar,Eigen::Dynamic,Eigen::Dynamic>>(solver.ref_eigvecs().data(),solver.Rows(), solver.Cols());
    }

    template<typename Scalar>
    Eigen::ArrayXXcd eigs_Arpackpp_sparse(Eigen::MatrixXcd matrix, int nev, class_tic_toc &timer){
        class_eigsolver_arpack<Scalar,Form::SYMMETRIC> solver;
        int ncv = std::min((int)matrix.rows()/2, nev*8);
        Eigen::Matrix<Scalar,Eigen::Dynamic, Eigen::Dynamic> matrix_scalar;
        if constexpr (std::is_same<Scalar,std::complex<double>>::value){
            matrix_scalar = matrix;
        }
        else
        {
            matrix_scalar = matrix.real();
        }


        int L   = matrix_scalar.rows();
        timer.tic();
        solver.eig(matrix_scalar.data(), L , nev,ncv, Ritz::LM, Side::R, true, false);
        timer.toc();

        return Eigen::Map<const Eigen::Array<Scalar,Eigen::Dynamic,Eigen::Dynamic>>(solver.ref_eigvecs().data(),solver.Rows(), solver.Cols());
    }


};


#endif //TRAINING_BENCHMARK_ARPACKPP_H
