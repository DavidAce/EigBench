//
// Created by david on 2018-10-24.
//

#ifndef TRAINING_BENCHMARK_LAPACKE_H
#define TRAINING_BENCHMARK_LAPACKE_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <general/class_tic_toc.h>
#include <lapacke.h>
class benchmark_lapacke {
public:

    benchmark_lapacke() = default;

    auto eig_sym_Lapacke_zheev(Eigen::MatrixXcd matrix, class_tic_toc &timer){
        //Complex full eigenvalue solver
        int L = matrix.rows();
        Eigen::ArrayXd  wl(L);
        timer.tic();
        LAPACKE_zheev(LAPACK_COL_MAJOR,'V','U', L, reinterpret_cast< __complex__ double*>(matrix.data()), L, wl.data());
        timer.toc();
        return
        std::make_pair(
                wl,
                Eigen::Map<Eigen::ArrayXXcd> (matrix.data(),L, L)
                );
//        return eigvals;
    }


    auto eig_Lapacke_dgeev(Eigen::MatrixXd matrix, class_tic_toc &timer){
        //Real full eigenvalue solver

        int L = matrix.rows();
        Eigen::ArrayXd  wr(L);
        Eigen::ArrayXd  wi(L);
        Eigen::ArrayXXd vl(L,L);
        Eigen::ArrayXXd vr(L,L);
        timer.tic();
        LAPACKE_dgeev(LAPACK_COL_MAJOR,'V','V', L ,matrix.data(), L, wr.data(), wi.data(), vl.data(),L, vr.data(),L);
        timer.toc();
        // Write back into a complex container
        Eigen::ArrayXcd eigvals(L);
        Eigen::ArrayXXd eigvecs(L,L);
        eigvals.real() = Eigen::Map<Eigen::ArrayXd>(wr.data(), L);
        eigvals.imag() = Eigen::Map<Eigen::ArrayXd>(wi.data(), L);
        eigvecs.real() = Eigen::Map<Eigen::ArrayXd>(vr.data(), L);
//        eigvecs.imag() = Eigen::Map<Eigen::ArrayXd>(vi.data(), L);
        return  std::make_pair(eigvals,eigvecs );
    }

    auto eig_Lapacke_zheevd(Eigen::MatrixXcd matrix, class_tic_toc &timer){
        int L = matrix.rows();
        Eigen::ArrayXd  wl(L);

        timer.tic();
        LAPACKE_zheevd(LAPACK_COL_MAJOR,'V','U', L,reinterpret_cast< __complex__ double*>(matrix.data()), L, wl.data());
        timer.toc();

//    timer.tic();
//    LAPACKE_zheevd(LAPACK_COL_MAJOR,'N','U', L,reinterpret_cast< __complex__ double*>(matrix), L, wl);
//    timer.toc();
//        Eigen::ArrayXd eigvals = Eigen::Map<Eigen::ArrayXd> (wl, L);
        return std::make_pair(wl,Eigen::Map<Eigen::ArrayXXcd> (matrix.data(), L,L));
    }

};


#endif //TRAINING_BENCHMARK_LAPACKE_H
