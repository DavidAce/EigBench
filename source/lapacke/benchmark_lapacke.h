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

    Eigen::ArrayXd eig_Lapacke_zheev(Eigen::MatrixXd &matrix, int L, class_tic_toc &timer){
        auto *wl=new double[L];//eigenvalues
        timer.tic();
        LAPACKE_zheev(LAPACK_COL_MAJOR,'V','U', L, reinterpret_cast< __complex__ double*>(matrix.data()), L, wl);
        timer.toc();
        Eigen::ArrayXd eigvals = Eigen::Map<Eigen::ArrayXd> (wl, L);
        delete[] wl;
        return eigvals;
    }


    Eigen::ArrayXcd eig_Lapacke_dgeev(Eigen::MatrixXd matrix, class_tic_toc &timer){
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
        eigvals.real() = Eigen::Map<Eigen::ArrayXd>(wr.data(), L);
        eigvals.imag() = Eigen::Map<Eigen::ArrayXd>(wi.data(), L);
        return eigvals;
    }

    Eigen::ArrayXd eig_Lapacke_zheevd(Eigen::MatrixXd matrix, class_tic_toc &timer){
        int L = matrix.rows();
        auto* wl=new double[L];//eigenvalues

        timer.tic();
        LAPACKE_zheevd(LAPACK_COL_MAJOR,'V','U', L,reinterpret_cast< __complex__ double*>(matrix.data()), L, wl);
        timer.toc();

//    timer.tic();
//    LAPACKE_zheevd(LAPACK_COL_MAJOR,'N','U', L,reinterpret_cast< __complex__ double*>(matrix), L, wl);
//    timer.toc();
        Eigen::ArrayXd eigvals = Eigen::Map<Eigen::ArrayXd> (wl, L);
        delete[] wl;
        return eigvals;
    }

};


#endif //TRAINING_BENCHMARK_LAPACKE_H
