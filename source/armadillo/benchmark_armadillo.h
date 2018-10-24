//
// Created by david on 2018-10-24.
//

#ifndef TRAINING_BENCHMARK_ARMADILLO_H
#define TRAINING_BENCHMARK_ARMADILLO_H
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <general/class_tic_toc.h>
#include <armadillo>


class benchmark_armadillo {
public:
    benchmark_armadillo() = default;

    template <typename Derived>
    Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>  eig_Armadillo_dense(Eigen::MatrixBase<Derived> matrix, class_tic_toc &timer){
        using namespace arma;
        int L = matrix.rows();
        Mat<typename Derived::Scalar> D(matrix.derived().data(), L,L,false, false);
        vec eigval;
        cx_mat eigvec;
        timer.tic();
        eig_sym(eigval, eigvec, D);
        timer.toc();
        return Eigen::Map<Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>> (eigvec.begin(), eigvec.n_rows, eigvec.n_cols);
    }


    Eigen::ArrayXXcd  eigs_Armadillo_sparse(Eigen::MatrixXcd matrix, int nev, class_tic_toc &timer){
        using namespace arma;
        long L = matrix.rows();
        Eigen::MatrixXd matrixreal = matrix.real();
        mat D(matrixreal.data(), L,L,false, false);
        sp_mat S(D);
        vec eigval;
        mat eigvec;
        timer.tic();
        eigs_sym(eigval,eigvec,S,nev);
        timer.toc();
        return Eigen::Map<Eigen::ArrayXXd> (eigvec.begin(), eigvec.n_rows, eigvec.n_cols);
    }

};


#endif //TRAINING_BENCHMARK_ARMADILLO_H
