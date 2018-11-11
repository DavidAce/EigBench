//
// Created by david on 2018-11-11.
//

#ifndef EIGBENCH_RESULTS_CONTAINER_H
#define EIGBENCH_RESULTS_CONTAINER_H
#include <iostream>
#include <iomanip>
#include <Eigen/Core>
#include <vector>
#include "matrix_container.h"

class results_container {
private:
    int R;
    int L;
    int S;
    Eigen::ArrayXd timings;
    Eigen::ArrayXd iterations;
    Eigen::ArrayXd ncvparams;
    Eigen::ArrayXd matvecops;
    const matrix_container &matrices;
    const std::vector<int> &nev_list;

    std::vector<std::tuple<int,double,int,double,double,double,double,double,double,double>>    results;

public:


    results_container(const matrix_container &matrices_, const std::vector<int> &nev_list_)
                 : matrices(matrices_), nev_list(nev_list_)
    {
        R = matrices.R;
        L = matrices.L;
        S = matrices.S;
        timings.resize(R);
        iterations.resize(R);
        ncvparams.resize(R);
        matvecops.resize(R);
    }

    void add_entry(int r, int n, int sp, int l, int L, double Sparcity, int Nev, int Ncv, double Time, int Iter, int Ops){
        assert (r       <  R                          and "Too many reps");
        assert(Nev      == nev_list[n]                and "Nev mismatch");
        assert(Sparcity == matrices.sparcity_list[sp] and "Sparcity mismatch");
        assert(L        == matrices.L_list[l]         and "L mismatch");
        timings(r)      = Time;
        iterations(r)   = Iter;
        ncvparams(r)    = Ncv;
        matvecops(r)    = Ops;

        if (r == R - 1){
            double tavg = timings.mean();
            double tstd = std::sqrt((timings - tavg).square().sum()/(timings.size()-1));
            double iavg = iterations.mean();
            double istd = std::sqrt((iterations - iavg).square().sum()/(iterations.size()-1));
            double navg = ncvparams.mean();
//            double nstd = std::sqrt((ncvparams - navg).square().sum()/(ncvparams.size()-1));
            double oavg = ncvparams.mean();
            double ostd = std::sqrt((matvecops - oavg).square().sum()/(matvecops.size()-1));
            results.emplace_back(L,Sparcity,Nev,tavg,tstd,iavg, istd,navg,oavg,ostd);
        }

    }



    void print_results(std::string header){
        std::cout << header << std::endl;
        std::cout << std::setprecision(8)
                  << "      " <<std::setw(6)   << std::left << "L"
                  << "      " <<std::setw(12)  << std::left << "Sparcity"
                  << "      " <<std::setw(6)   << std::left << "Nev"
                  << "      " <<std::setw(18)  << std::left << "Average [seconds]"
                  << "      " <<std::setw(18)  << std::left << "std"
                  << "      " <<std::setw(18)  << std::left << "Average [iter]"
                  << "      " <<std::setw(18)  << std::left << "std"
                  << "      " <<std::setw(18)  << std::left << "Ncv"
                  << std::endl;
        for(auto &res : results){
            std::cout << std::setprecision(8)
                      << "      " << std::setw(6)  << std::fixed << std::left << std::get<0>(res)
                      << "      " << std::setw(12) << std::fixed << std::left << std::get<1>(res)
                      << "      " << std::setw(6)  << std::fixed << std::left << std::get<2>(res)
                      << "      " << std::setw(18) << std::fixed << std::left << std::get<3>(res)
                      << "      " << std::setw(18) << std::fixed << std::left << std::get<4>(res)
                      << "      " << std::setw(18) << std::fixed << std::left << std::get<5>(res)
                      << "      " << std::setw(18) << std::fixed << std::left << std::get<6>(res)
                      << "      " << std::setw(18) << std::fixed << std::left << std::get<7>(res)
                      << std::endl;
        }
        std::cout << std::endl << std::flush;
    }



};


#endif //EIGBENCH_RESULTS_CONTAINER_H
