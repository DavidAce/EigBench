//
// Created by david on 2018-10-24.
//

#ifndef TRAINING_BENCHMARK_EIGEN3_H
#define TRAINING_BENCHMARK_EIGEN3_H
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <general/class_tic_toc.h>

class benchmark_eigen3 {
public:

    benchmark_eigen3() = default;

    template <typename Derived>
    Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> eig_sym_Eigen(const Eigen::EigenBase<Derived> &matrix, class_tic_toc &timer){
        timer.tic();
        Eigen::SelfAdjointEigenSolver<Derived> es(matrix);
        timer.toc();
        return es.eigenvectors();
    }

    template <typename Derived>
    Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> eig_gen_Eigen(const Eigen::EigenBase<Derived> &matrix, class_tic_toc &timer){
        timer.tic();
        Eigen::EigenSolver<Derived> es(matrix);
        timer.toc();
        return es.eigenvectors();
    }

};


#endif //TRAINING_BENCHMARK_EIGEN3_H
