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
    const std::vector<int> &ncv_list;

    std::vector<std::tuple<
            int,      //L
            double,   //Sparcity
            int,      //Nev
            int,      //Ncv
            double,   //tavg
            double,   //tmed
            double,   //tstd
            double,   //iavg
            double,   //imed
            double,   //istd
            double,   //oavg
            double,   //omed
            double    //ostd
            >>    results;

public:


    results_container(const matrix_container &matrices_, const std::vector<int> &nev_list_, const std::vector<int> &ncv_list_)
                 : matrices(matrices_), nev_list(nev_list_), ncv_list(ncv_list_)
    {
        R = matrices.R;
        L = matrices.L;
        S = matrices.S;
        timings.resize(R);
        iterations.resize(R);
        ncvparams.resize(R);
        matvecops.resize(R);
    }

    void add_entry(int r, int nev,int ncv, int sp, int l, int L, double Sparcity, int Nev, int Ncv, double Time, int Iter, int Ops){
        assert (r       <  R                          and "Too many reps");
        assert(Nev      == nev_list[nev]              and "Nev mismatch");
        assert(Ncv      == ncv_list[ncv]              and "Ncv mismatch");
        assert(Sparcity == matrices.sparcity_list[sp] and "Sparcity mismatch");
        assert(L        == matrices.L_list[l]         and "L mismatch");
        timings(r)      = Time;
        iterations(r)   = Iter;
        ncvparams(r)    = Ncv;
        matvecops(r)    = Ops;

        if (r == R - 1){
            std::sort(timings.data()   , timings.data()   +timings.size());
            std::sort(iterations.data(), iterations.data()+iterations.size());
            std::sort(ncvparams.data() , ncvparams.data() +ncvparams.size());
            std::sort(matvecops.data() , matvecops.data() +matvecops.size());


            double tavg = timings.mean();
            double tmed = timings((int)(R/2));
            double tstd = std::sqrt((timings - tavg).square().sum()/(timings.size()-1));
            double iavg = iterations.mean();
            double imed = iterations((int)(R/2));
            double istd = std::sqrt((iterations - iavg).square().sum()/(iterations.size()-1));
            int    navg = (int)ncvparams.mean();
//            double nstd = std::sqrt((ncvparams - navg).square().sum()/(ncvparams.size()-1));
            double oavg = matvecops.mean();
            double omed = matvecops((int)(R/2));
            double ostd = std::sqrt((matvecops - oavg).square().sum()/(matvecops.size()-1));
            results.emplace_back(L,Sparcity,Nev,navg,tavg,tmed,tstd,iavg,imed,istd,oavg,omed,ostd);
        }

    }



    void print_results(std::string header){
        std::cout << header << std::endl;
        std::cout
                << std::setw(6)   << std::right << "L"
                << std::setw(10)  << std::right << "sparcity"
                << std::setw(6)   << std::right << "nev"
                << std::setw(6)   << std::right << "ncv"
                << std::setw(14)  << std::right << "avg time[ms]"
                << std::setw(10)  << std::right << "med"
                << std::setw(10)  << std::right << "std"
                << std::setw(10)  << std::right << "avg iter"
                << std::setw(8)   << std::right << "med"
                << std::setw(8)   << std::right << "std"
                << std::setw(8)   << std::right << "avg ops"
                << std::setw(8)   << std::right << "med"
                << std::setw(8)   << std::right << "std"
                << std::endl;
        for(auto &res : results){
            std::cout << std::setprecision(2)
                      << std::setw(6)  << std::fixed << std::right << std::get<0>(res)   // L
                      << std::setw(10) << std::fixed << std::right << std::get<1>(res)   // sparcity
                      << std::setw(6)  << std::fixed << std::right << std::get<2>(res)   // nev
                      << std::setw(6)  << std::fixed << std::right << std::get<3>(res)   // ncv
                      << std::setprecision(3)
                      << std::setw(14) << std::fixed << std::right << std::get<4>(res)*1000 // avg time[ms]
                      << std::setw(10) << std::fixed << std::right << std::get<5>(res)*1000 // med
                      << std::setw(10) << std::fixed << std::right << std::get<6>(res)*1000 // std
                      << std::setprecision(1)
                      << std::setw(10) << std::fixed << std::right << std::get<7>(res)   // avg iter
                      << std::setw(8)  << std::fixed << std::right << std::get<8>(res)   // med
                      << std::setw(8)  << std::fixed << std::right << std::get<9>(res)   // std
                      << std::setw(8)  << std::fixed << std::right << std::get<10>(res)  // avg ops
                      << std::setw(8)  << std::fixed << std::right << std::get<11>(res)  // med
                      << std::setw(8)  << std::fixed << std::right << std::get<12>(res)  // std
                      << std::endl;
        }
        std::cout << std::endl << std::flush;
    }



};


#endif //EIGBENCH_RESULTS_CONTAINER_H
