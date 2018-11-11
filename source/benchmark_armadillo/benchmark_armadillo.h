//
// Created by david on 2018-10-24.
//

#ifndef TRAINING_BENCHMARK_ARMADILLO_H
#define TRAINING_BENCHMARK_ARMADILLO_H
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <general/class_tic_toc.h>
#include <general/matrix_generator.h>
#include <armadillo>


class benchmark_armadillo {
public:
    benchmark_armadillo() = default;

    auto run_eigs_real_symm_sparse(int reps = 100, std::vector<int> nev_list = {1,8}, std::vector<int> L_list = {20,40}, std::vector<double> sparcity_list = {0.1,0.5}, int seed = 0 ){
        using namespace arma;
        rn::seed(seed);
        class_tic_toc timer(true,5,"");
        matrix_generator matgen;
        std::vector<std::tuple<int,double,int,double,double,double,double,double,double>>    results;
        for (auto L : L_list){
            for (auto s : sparcity_list) {
                for(auto nev : nev_list) {
                    Eigen::ArrayXd timings(reps);
//                    Eigen::ArrayXd iterations(R);
                    for (int r = 0; r < reps; r++) {
                        auto matrix = matgen.make_symmetric_dense<double>(L, s,false);
//                        int ncv = std::min((int)matrix.rows()/2, nev*8);
                        mat D(matrix.data(), L,L,false, false);
                        sp_mat S(D);
                        vec eigval;
                        mat eigvec;
                        timer.tic();
                        eigs_sym(eigval,eigvec,S,nev);
                        timer.toc();
                        timings(r) = timer.get_last_time_interval();
//                        iterations(r) =
                    }
                    double tavg = timings.mean();
                    double tstd = std::sqrt((timings - tavg).square().sum()/(timings.size()-1));
                    double nan = std::numeric_limits<double>::quiet_NaN();
//                    double iavg = iterations.mean();
//                    double istd = std::sqrt((iterations - iavg).square().sum()/(iterations.size()-1));
                    results.emplace_back(L,s,nev,tavg, tstd, nan, nan,nan,nan);
                }
            }
        }
        return results;
    }



    template <typename Derived>
    std::pair<Eigen::Array<double, Eigen::Dynamic, 1>,
              Eigen::Array<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>>
              eig_sym_Armadillo_dense(Eigen::MatrixBase <Derived> matrix, class_tic_toc &timer)
    {
        using namespace arma;
        long L = matrix.rows();
        Mat<typename Derived::Scalar> D(matrix.derived().data(), L,L,false, false);
        vec eigval;
        Mat<typename Derived::Scalar> eigvec;
        timer.tic();
        eig_sym(eigval, eigvec, D);
        timer.toc();

        return
        std::make_pair(
                Eigen::Map<Eigen::Array<double, Eigen::Dynamic, 1                           >> (eigval.begin(), eigval.n_rows),
                Eigen::Map<Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>> (eigvec.begin(), eigvec.n_rows, eigvec.n_cols)
        );
    }

    template <typename Derived>
    std::pair<Eigen::Array<std::complex<double>, Eigen::Dynamic, 1>,
              Eigen::Array<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>>
    eig_Armadillo_nsym_dense(Eigen::MatrixBase<Derived> matrix, class_tic_toc &timer)
    {
        using namespace arma;
        long L = matrix.rows();
        Mat<typename Derived::Scalar> D(matrix.derived().data(), L,L,false, false);
        cx_vec eigval;
        cx_mat eigvec;
        timer.tic();
        eig_gen(eigval,eigvec,  D);
        timer.toc();
        std::cout << eigval << std::endl;
        return
                std::make_pair(
                        Eigen::Map<Eigen::Array<std::complex<double>, Eigen::Dynamic, 1             >> (eigval.begin(), eigval.n_rows),
                        Eigen::Map<Eigen::Array<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> (eigvec.begin(), eigvec.n_rows, eigvec.n_cols)
                );
    }


    std::pair<Eigen::Array<double, Eigen::Dynamic, 1>,
              Eigen::Array<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>>
    eigs_Armadillo_real_symm_sparse(Eigen::MatrixXd matrix, int nev, class_tic_toc &timer)
    {
        using namespace arma;
        long L = matrix.rows();
        mat D(matrix.data(), L,L,false, false);
        sp_mat S(D);
        vec eigval;
        mat eigvec;
        timer.tic();
        eigs_sym(eigval,eigvec,S,nev);
        timer.toc();
        return
                std::make_pair(
                        Eigen::Map<Eigen::Array<double, Eigen::Dynamic, 1             >> (eigval.begin(), eigval.n_rows),
                        Eigen::Map<Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>> (eigvec.begin(), eigvec.n_rows, eigvec.n_cols)
                );
    }


    std::pair<Eigen::Array<double, Eigen::Dynamic, 1>,
              Eigen::Array<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>>
            eigs_Armadillo_real_nsym_sparse(Eigen::MatrixXcd matrix, int nev, class_tic_toc &timer)
    {
        using namespace arma;
        Eigen::MatrixXd matrixreal = matrix.real();
        long L = matrixreal.rows();
        mat D(matrixreal.data(), L,L,false, false);
        sp_mat S(D);
        vec eigval;
        mat eigvec;
        timer.tic();
        eigs_sym(eigval,eigvec,S,nev);
        timer.toc();
        return
        std::make_pair(
                Eigen::Map<Eigen::Array<double, Eigen::Dynamic, 1             >> (eigval.begin(), eigval.n_rows),
                Eigen::Map<Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>> (eigvec.begin(), eigvec.n_rows, eigvec.n_cols)
        );
    }

};


#endif //TRAINING_BENCHMARK_ARMADILLO_H
