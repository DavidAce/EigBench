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
#include <general/results_container.h>

class benchmark_arpackcpp {
public:

    benchmark_arpackcpp() = default;


    template<
            eigutils::eigSetting::Type    type,
            eigutils::eigSetting::Form    form,
            eigutils::eigSetting::Storage storage,
            eigutils::eigSetting::Ritz    ritz = eigutils::eigSetting::Ritz::LM,
            eigutils::eigSetting::Side    side = eigutils::eigSetting::Side::R
            >
    auto run_eigs(matrix_container &matrices, std::vector<int> nev_list, std::vector<int> ncv_list){
        using namespace eigutils::eigSetting;
        class_tic_toc timer(true,5,"");
        results_container results(matrices,nev_list, ncv_list);
        for (size_t l = 0; l <  matrices.L_list.size() ; l++){
            for (size_t s = 0; s <  matrices.sparcity_list.size() ; s++) {
                for(size_t ne = 0; ne < nev_list.size(); ne++) {
                    for(size_t nc = 0; nc < ncv_list.size(); nc++) {
                        int nev = nev_list[ne];
                        int ncv = ncv_list[nc];
                        Eigen::ArrayXd timings(matrices.R);
                        Eigen::ArrayXd iterations(matrices.R);
                        Eigen::ArrayXd ncvparams(matrices.R);
                        Eigen::ArrayXd matvecops(matrices.R);
                        for (int r = 0; r < matrices.R; r++) {
                            auto &matrix = matrices.get_matrix<type, form>(r, s, l);
                            assert(matrices.L_list[l] == matrix.rows());
                            timer.tic();
                            class_eigsolver solver;
                            solver.eigs<storage>(matrix.data(), matrix.rows(), nev, ncv, NAN, form, ritz, side, true,
                                                 false);
                            timer.toc();
                            results.add_entry(r, ne,nc, s, l, matrix.rows(), matrices.sparcity_list[s], nev,
                                              solver.solution.meta.ncv_used, timer.get_last_time_interval(),
                                              solver.solution.meta.iter, solver.solution.meta.counter);
                        }
                    }
                }

            }
        }
        return results;
    }


    template<
            eigutils::eigSetting::Type    type,
            eigutils::eigSetting::Form    form,
            eigutils::eigSetting::Storage storage,
            eigutils::eigSetting::Ritz    ritz = eigutils::eigSetting::Ritz::LM,
            eigutils::eigSetting::Side    side = eigutils::eigSetting::Side::R
    >
    auto run_eigs_shift_invert(matrix_container &matrices, std::vector<int> nev_list, std::vector<int> ncv_list){
        using namespace eigutils::eigSetting;
        class_tic_toc timer(true,5,"");
        results_container results(matrices,nev_list, ncv_list);
        for (size_t l = 0; l <  matrices.L_list.size() ; l++){
            for (size_t s = 0; s <  matrices.sparcity_list.size() ; s++) {
                for(size_t ne = 0; ne < nev_list.size(); ne++) {
                    for(size_t nc = 0; nc < ncv_list.size(); nc++) {
                        int nev = nev_list[ne];
                        int ncv = ncv_list[nc];
                        Eigen::ArrayXd timings(matrices.R);
                        Eigen::ArrayXd iterations(matrices.R);
                        Eigen::ArrayXd ncvparams(matrices.R);
                        Eigen::ArrayXd matvecops(matrices.R);
                        for (int r = 0; r < matrices.R; r++) {
                            auto &matrix = matrices.get_matrix<type, form>(r, s, l);
                            assert(matrices.L_list[l] == matrix.rows());
                            timer.tic();
                            class_eigsolver solver;
                            solver.eigs<storage>(matrix.data(), matrix.rows(), nev, ncv, 1.0, form, ritz, side, true,
                                                 false);
                            timer.toc();
                            results.add_entry(r, ne,nc, s, l, matrix.rows(), matrices.sparcity_list[s], nev,
                                              solver.solution.meta.ncv_used, timer.get_last_time_interval(),
                                              solver.solution.meta.iter, solver.solution.meta.counter);
                        }
                    }
                }

            }
        }
        return results;
    }


};


#endif //BENCHMARK_ARPACKPP_H
