//
// Created by david on 2018-11-08.
//

#ifndef EIGBENCH_MATRIX_CONTAINER_H
#define EIGBENCH_MATRIX_CONTAINER_H
#include <vector>
#include <complex>
#include <Eigen/Core>
#include "matrix_generator.h"

class matrix_container{
public:
    std::vector<int>     L_list        ;
    std::vector<double>  sparcity_list ;
    int R;
    int L;
    int S;
    std::vector<Eigen::MatrixXd >  matrix_real_symm_list;
    std::vector<Eigen::MatrixXd >  matrix_real_nsym_list;
    std::vector<Eigen::MatrixXcd>  matrix_cplx_symm_list;
    std::vector<Eigen::MatrixXcd>  matrix_cplx_nsym_list;
    matrix_container(std::vector<int>     L_list_,
                    std::vector<double>  sparcity_list_,
                    int reps_
    )
            : L_list(L_list_), sparcity_list(sparcity_list_), R(reps_)
    {
        matrix_generator matgen;
        L = (int)L_list.size();
        S = (int)sparcity_list.size();
        matrix_real_symm_list.resize(R *L * S);
        matrix_real_nsym_list.resize(R *L * S);
        matrix_cplx_symm_list.resize(R *L * S);
        matrix_cplx_nsym_list.resize(R *L * S);
        int i = 0;
        for(auto L : L_list){
            for (auto s: sparcity_list) {
                for (int r = 0 ; r < R; r++) {
                    matrix_real_symm_list[i] = matgen.make_symmetric_dense<double>(L, s);
                    matrix_real_nsym_list[i] = matgen.make_nonsymmetric_dense<double>(L, s);
                    matrix_cplx_symm_list[i] = matgen.make_symmetric_dense < std::complex < double >> (L, s);
                    matrix_cplx_nsym_list[i] = matgen.make_nonsymmetric_dense < std::complex < double >> (L, s);
                    i++;
                }
            }
        }
    }


    const auto & get_real_symm(size_t r, size_t s, size_t l){return matrix_real_symm_list[r + s*R + l * R*S];}
    const auto & get_real_nsym(size_t r, size_t s, size_t l){return matrix_real_nsym_list[r + s*R + l * R*S];}
    const auto & get_cplx_symm(size_t r, size_t s, size_t l){return matrix_cplx_symm_list[r + s*R + l * R*S];}
    const auto & get_cplx_nsym(size_t r, size_t s, size_t l){return matrix_cplx_nsym_list[r + s*R + l * R*S];}
//    const auto & get_rows     (size_t r, size_t s, size_t l){return matrix_real_nsym_list[r + s*R + l * R*S].rows();}
//    const auto & get_sparcity (size_t r, size_t s, size_t l){return matrix_real_symm_list[r + s*R + l * R*S];}



};

#endif //EIGBENCH_MATRIX_CONTAINER_H
