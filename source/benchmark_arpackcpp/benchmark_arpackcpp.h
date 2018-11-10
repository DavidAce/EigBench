//
// Created by david on 2018-10-24.
//

#ifndef BENCHMARK_ARPACKPP_H
#define BENCHMARK_ARPACKPP_H
#include <Eigen/Core>
#include <complex>
#include <general/class_tic_toc.h>
#include <benchmark_arpackcpp/class_eigsolver.h>
#include <general/matrix_generator.h>
#include <general/matrix_container.h>

class benchmark_arpackcpp {
public:

    benchmark_arpackcpp() = default;


    auto run_eigs_real_symm_dense(matrix_container &matrices, std::vector<int> nev_list){
        using namespace arpack_extra::modes;
        class_tic_toc timer(true,5,"");
        std::vector<std::tuple<int,double,int,double,double,double,double>>    results;
        for (size_t l = 0; l <  matrices.L_list.size() ; l++){
            for (size_t s = 0; s <  matrices.sparcity_list.size() ; s++) {
                for(auto nev : nev_list) {
                    Eigen::ArrayXd timings(matrices.R);
                    Eigen::ArrayXd iterations(matrices.R);
                    for(size_t r = 0; r < matrices.R; r++){
                        auto &matrix = matrices.get_real_symm(r,s,l);
                        assert(matrices.L_list[l]== matrix.rows());
                        int ncv = std::min((int)matrix.rows()/2, nev*8);
                        timer.tic();
                        class_eigsolver solver;
                        solver.eigs_dense(matrix.data(), matrix.rows(), nev, -1, NAN, Form::SYMMETRIC, Ritz::LM, Side::R,true,false);
                        timer.toc();
                        timings(r) = timer.get_last_time_interval();
                        iterations(r) = solver.solution.meta.iter;
                    }

                    double tavg = timings.mean();
                    double tstd = std::sqrt((timings - tavg).square().sum()/(timings.size()-1));
                    double iavg = iterations.mean();
                    double istd = std::sqrt((iterations - iavg).square().sum()/(iterations.size()-1));
                    results.emplace_back(matrices.L_list[l],matrices.sparcity_list[s],nev,tavg,tstd,iavg, istd);

                }

            }
        }
        return results;
    }

    auto run_eigs_real_nsym_dense(matrix_container &matrices, std::vector<int> nev_list){
        using namespace arpack_extra::modes;
        class_tic_toc timer(true,5,"");
        std::vector<std::tuple<int,double,int,double,double,double,double>>    results;
        for (size_t l = 0; l <  matrices.L_list.size() ; l++){
            for (size_t s = 0; s <  matrices.sparcity_list.size() ; s++) {
                for(auto nev : nev_list) {
                    Eigen::ArrayXd timings(matrices.R);
                    Eigen::ArrayXd iterations(matrices.R);
                    for(size_t r = 0; r < matrices.R; r++){
                        auto &matrix = matrices.get_real_nsym(r,s,l);
                        assert(matrices.L_list[l]== matrix.rows());
//                        int ncv = std::min((int)matrix.rows()/2, nev*8);
                        timer.tic();
                        class_eigsolver solver;
                        solver.eigs_dense(matrix.data(), matrix.rows(), nev, -1, NAN, Form::NONSYMMETRIC, Ritz::LM, Side::R,true,false);
                        timer.toc();
                        timings(r) = timer.get_last_time_interval();
                        iterations(r) = solver.solution.meta.iter;
//                        std::cout << "iter: " <<  solver.solution.meta.iter << std::endl;

                    }

                    double tavg = timings.mean();
                    double tstd = std::sqrt((timings - tavg).square().sum()/(timings.size()-1));
                    double iavg = iterations.mean();
                    double istd = std::sqrt((iterations - iavg).square().sum()/(iterations.size()-1));
                    results.emplace_back(matrices.L_list[l],matrices.sparcity_list[s],nev,tavg,tstd,iavg, istd);

                }

            }
        }
        return results;
    }


    auto run_eigs_real_symm_sparse(matrix_container &matrices, std::vector<int> nev_list){
        using namespace arpack_extra::modes;
        class_tic_toc timer(true,5,"");
        std::vector<std::tuple<int,double,int,double,double,double,double>>    results;
        for (size_t l = 0; l <  matrices.L_list.size() ; l++){
            for (size_t s = 0; s <  matrices.sparcity_list.size() ; s++) {
                for(auto nev : nev_list) {
                    Eigen::ArrayXd timings(matrices.R);
                    Eigen::ArrayXd iterations(matrices.R);
                    for(size_t r = 0; r < matrices.R; r++){
                        auto &matrix = matrices.get_real_symm(r,s,l);
                        assert(matrices.L_list[l]== matrix.rows());
                        int ncv = std::min((int)matrix.rows()/2, nev*8);
                        timer.tic();
                        class_eigsolver solver;
                        solver.eigs_sparse(matrix.data(), matrix.rows(), nev, -1, NAN, Form::SYMMETRIC, Ritz::LM, Side::R,true,false);
                        timer.toc();
                        timings(r) = timer.get_last_time_interval();
                        iterations(r) = solver.solution.meta.iter;
//                        std::cout << "iter: " <<  solver.solution.meta.iter << std::endl;
                    }

                    double tavg = timings.mean();
                    double tstd = std::sqrt((timings - tavg).square().sum()/(timings.size()-1));
                    double iavg = iterations.mean();
                    double istd = std::sqrt((iterations - iavg).square().sum()/(iterations.size()-1));
                    results.emplace_back(matrices.L_list[l],matrices.sparcity_list[s],nev,tavg,tstd,iavg, istd);

                }

            }
        }
        return results;
    }




    auto run_eigs_real_nsym_sparse(matrix_container &matrices, std::vector<int> nev_list){
        using namespace arpack_extra::modes;
        class_tic_toc timer(true,5,"");
        std::vector<std::tuple<int,double,int,double,double,double,double>>    results;
        for (size_t l = 0; l <  matrices.L_list.size() ; l++){
            for (size_t s = 0; s <  matrices.sparcity_list.size() ; s++) {
                for(auto nev : nev_list) {
                    Eigen::ArrayXd timings(matrices.R);
                    Eigen::ArrayXd iterations(matrices.R);
                    for(size_t r = 0; r < matrices.R; r++){
                        auto &matrix = matrices.get_real_nsym(r,s,l);
                        assert(matrices.L_list[l]== matrix.rows());
//                        int ncv = std::min((int)matrix.rows()/2, nev*8);
                        timer.tic();
                        class_eigsolver solver;
                        solver.eigs_sparse(matrix.data(), matrix.rows(), nev, -1, NAN, Form::NONSYMMETRIC, Ritz::LM, Side::R,true,false);
                        timer.toc();
                        timings(r) = timer.get_last_time_interval();
                        iterations(r) = solver.solution.meta.iter;
//                        std::cout << "iter: " <<  solver.solution.meta.iter << std::endl;
                    }

                    double tavg = timings.mean();
                    double tstd = std::sqrt((timings - tavg).square().sum()/(timings.size()-1));
                    double iavg = iterations.mean();
                    double istd = std::sqrt((iterations - iavg).square().sum()/(iterations.size()-1));
                    results.emplace_back(matrices.L_list[l],matrices.sparcity_list[s],nev,tavg,tstd,iavg, istd);

                }

            }
        }
        return results;
    }


    auto run_eigs_cplx_symm_dense(matrix_container &matrices, std::vector<int> nev_list){
        using namespace arpack_extra::modes;
        class_tic_toc timer(true,5,"");
        std::vector<std::tuple<int,double,int,double,double,double,double>>    results;
        for (size_t l = 0; l <  matrices.L_list.size() ; l++){
            for (size_t s = 0; s <  matrices.sparcity_list.size() ; s++) {
                for(auto nev : nev_list) {
                    Eigen::ArrayXd timings(matrices.R);
                    Eigen::ArrayXd iterations(matrices.R);
                    for(size_t r = 0; r < matrices.R; r++){
                        auto &matrix = matrices.get_real_symm(r,s,l);
                        assert(matrices.L_list[l]== matrix.rows());
                        int ncv = std::min((int)matrix.rows()/2, nev*8);
                        timer.tic();
                        class_eigsolver solver;
                        solver.eigs_dense(matrix.data(), matrix.rows(), nev, -1, NAN, Form::SYMMETRIC, Ritz::LM, Side::R,true,false);
                        timer.toc();
                        timings(r) = timer.get_last_time_interval();
                        iterations(r) = solver.solution.meta.iter;
                    }

                    double tavg = timings.mean();
                    double tstd = std::sqrt((timings - tavg).square().sum()/(timings.size()-1));
                    double iavg = iterations.mean();
                    double istd = std::sqrt((iterations - iavg).square().sum()/(iterations.size()-1));
                    results.emplace_back(matrices.L_list[l],matrices.sparcity_list[s],nev,tavg,tstd,iavg, istd);

                }

            }
        }
        return results;
    }




    auto run_eigs_cplx_nsym_dense(matrix_container &matrices, std::vector<int> nev_list){
        using namespace arpack_extra::modes;
        class_tic_toc timer(true,5,"");
        std::vector<std::tuple<int,double,int,double,double,double,double>>    results;
        for (size_t l = 0; l <  matrices.L_list.size() ; l++){
            for (size_t s = 0; s <  matrices.sparcity_list.size() ; s++) {
                for(auto nev : nev_list) {
                    Eigen::ArrayXd timings(matrices.R);
                    Eigen::ArrayXd iterations(matrices.R);
                    for(size_t r = 0; r < matrices.R; r++){
                        auto &matrix = matrices.get_cplx_nsym(r,s,l);
                        assert(matrices.L_list[l]== matrix.rows());
//                        int ncv = std::min((int)matrix.rows()/2, nev*8);
                        timer.tic();
                        class_eigsolver solver;
                        solver.eigs_dense(matrix.data(), matrix.rows(), nev, -1, NAN, Form::NONSYMMETRIC, Ritz::LM, Side::R,true,false);
                        timer.toc();
                        timings(r) = timer.get_last_time_interval();
                        iterations(r) = solver.solution.meta.iter;
//                        std::cout << "iter: " <<  solver.solution.meta.iter << std::endl;

                    }

                    double tavg = timings.mean();
                    double tstd = std::sqrt((timings - tavg).square().sum()/(timings.size()-1));
                    double iavg = iterations.mean();
                    double istd = std::sqrt((iterations - iavg).square().sum()/(iterations.size()-1));
                    results.emplace_back(matrices.L_list[l],matrices.sparcity_list[s],nev,tavg,tstd,iavg, istd);

                }

            }
        }
        return results;
    }





    auto run_eigs_cplx_symm_sparse(matrix_container &matrices, std::vector<int> nev_list){
        using namespace arpack_extra::modes;
        class_tic_toc timer(true,5,"");
        std::vector<std::tuple<int,double,int,double,double,double,double>>    results;
        for (size_t l = 0; l <  matrices.L_list.size() ; l++){
            for (size_t s = 0; s <  matrices.sparcity_list.size() ; s++) {
                for(auto nev : nev_list) {
                    Eigen::ArrayXd timings(matrices.R);
                    Eigen::ArrayXd iterations(matrices.R);
                    for(size_t r = 0; r < matrices.R; r++){
                        auto &matrix = matrices.get_cplx_symm(r,s,l);
                        assert(matrices.L_list[l]== matrix.rows());
                        int ncv = std::min((int)matrix.rows()/2, nev*8);
                        timer.tic();
                        class_eigsolver solver;
                        solver.eigs_sparse(matrix.data(), matrix.rows(), nev, -1, NAN, Form::SYMMETRIC, Ritz::LM, Side::R,true,false);
                        timer.toc();
                        timings(r) = timer.get_last_time_interval();
                        iterations(r) = solver.solution.meta.iter;
//                        std::cout << "iter: " <<  solver.solution.meta.iter << std::endl;
                    }

                    double tavg = timings.mean();
                    double tstd = std::sqrt((timings - tavg).square().sum()/(timings.size()-1));
                    double iavg = iterations.mean();
                    double istd = std::sqrt((iterations - iavg).square().sum()/(iterations.size()-1));
                    results.emplace_back(matrices.L_list[l],matrices.sparcity_list[s],nev,tavg,tstd,iavg, istd);

                }

            }
        }
        return results;
    }



    auto run_eigs_cplx_nsym_sparse(matrix_container &matrices, std::vector<int> nev_list){
        using namespace arpack_extra::modes;
        class_tic_toc timer(true,5,"");
        std::vector<std::tuple<int,double,int,double,double,double,double>>    results;
        for (size_t l = 0; l <  matrices.L_list.size() ; l++){
            for (size_t s = 0; s <  matrices.sparcity_list.size() ; s++) {
                for(auto nev : nev_list) {
                    Eigen::ArrayXd timings(matrices.R);
                    Eigen::ArrayXd iterations(matrices.R);
                    for(size_t r = 0; r < matrices.R; r++){
                        auto &matrix = matrices.get_cplx_nsym(r,s,l);
                        assert(matrices.L_list[l]== matrix.rows());
//                        int ncv = std::min((int)matrix.rows()/2, nev*8);
                        timer.tic();
                        class_eigsolver solver;
                        solver.eigs_sparse(matrix.data(), matrix.rows(), nev, -1, NAN, Form::NONSYMMETRIC, Ritz::LM, Side::R,true,false);
                        timer.toc();
                        timings(r) = timer.get_last_time_interval();
                        iterations(r) = solver.solution.meta.iter;
//                        std::cout << "iter: " <<  solver.solution.meta.iter << std::endl;
                    }

                    double tavg = timings.mean();
                    double tstd = std::sqrt((timings - tavg).square().sum()/(timings.size()-1));
                    double iavg = iterations.mean();
                    double istd = std::sqrt((iterations - iavg).square().sum()/(iterations.size()-1));
                    results.emplace_back(matrices.L_list[l],matrices.sparcity_list[s],nev,tavg,tstd,iavg, istd);

                }

            }
        }
        return results;
    }






    template<typename Scalar>
    std::pair<Eigen::Array<std::complex<double>,Eigen::Dynamic,1             >,
              Eigen::Array<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic>>
    eigs_Arpackcpp_auto(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix, int nev, class_tic_toc &timer)
    {
        using namespace arpack_extra::modes;
        class_eigsolver autosolver;
        Scalar sigma = NAN;
        timer.tic();
        autosolver.eigs_auto(matrix.data(), matrix.rows(), nev, -1, sigma, Ritz::LM, Side::R,true,true);
        timer.toc();
        return
                std::make_pair(
                        Eigen::Map<const Eigen::Array<std::complex<double>,Eigen::Dynamic,1             >>(autosolver.solution.eigvals.data(),autosolver.solution.meta.nev_found),
                        Eigen::Map<const Eigen::Array<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic>>(autosolver.solution.eigvecs.data(),autosolver.solution.meta.rows, autosolver.solution.meta.cols)
                );
    }


};


#endif //BENCHMARK_ARPACKPP_H
