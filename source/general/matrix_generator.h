//
// Created by david on 2018-10-24.
//

#ifndef TRAINING_MATRIX_GENERATOR_H
#define TRAINING_MATRIX_GENERATOR_H
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <iomanip>
#include <general/nmspc_random_numbers.h>
class matrix_generator {
public:
    matrix_generator () = default;

    template <typename Scalar>
    Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> make_symmetric_dense(int L, double sparcity = 0.1, bool verbose = false){
        Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> mat_dense(L,L);
        mat_dense.setZero();
        int max_nonzero_elems = sparcity * L * L;
        int counter = 0;
        if constexpr (std::is_same<Scalar,std::complex<double>>::value){
            while(counter < max_nonzero_elems){
                int i = rn::uniform_integer(0,L-1);
                int j = rn::uniform_integer(0,L-1);
                if (mat_dense(i,j) != 0.0){ continue;}
                if (mat_dense(i,j) == 0.0 and i != j){
                    mat_dense(i,j) = rn::uniform_complex_1();
                    mat_dense(j,i) = mat_dense.conjugate()(i,j);
                    counter +=2;
                 }
                 else if (mat_dense(i,j) == 0.0 and i == j){
                    mat_dense(i,j) = rn::uniform_double_1();
                    counter += 1;
                }
            }
        }else{
            while(counter < max_nonzero_elems){
                int i = rn::uniform_integer(0,L-1);
                int j = rn::uniform_integer(0,L-1);
                if (mat_dense(i,j) != 0.0){ continue;}
                if (mat_dense(i,j) == 0.0 and i != j){
                    mat_dense(i,j) = rn::uniform_double_1();
                    mat_dense(j,i) = mat_dense(i,j);
                    counter +=2;
                }
                else if (mat_dense(i,j) == 0 and i == j){
                    mat_dense(i,j) = rn::uniform_double_1();
                    counter +=1;
                }
            }
        }
        if (verbose){
            // mat_dense = (mat_dense + mat_dense.adjoint()).eval();
            double final_sparcity = (double) (mat_dense.array().cwiseAbs() > 1e-15 )
                    .select(Eigen::MatrixXd::Ones(L,L),0).sum() / mat_dense.size();
            auto determinant = mat_dense.determinant();
            std::cout << "Matrix size = "      << std::setprecision(2)<< std::left << std::setw(6)  << mat_dense.rows()
                      << " Target sparcity = " << std::setprecision(6)<< std::left << std::setw(10) << sparcity
                      << " Final sparcity  = " << std::setprecision(6)<< std::left << std::setw(10) << final_sparcity
                      << " Determinant = "     << std::setprecision(6)<< std::left << std::setw(10) << determinant
                      << std::endl;
        }

        return mat_dense;
    }

    template <typename Scalar>
    Eigen::SparseMatrix<Scalar> make_hermitian_sparse(int L, double sparcity = 0.01){
        return make_symmetric_dense<Scalar>(L, sparcity).sparseView();
    }

    template <typename Scalar>
    Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> make_nonsymmetric_dense(int L, double sparcity = 0.01, bool verbose = false){
        Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> mat_dense(L,L);
        mat_dense.setZero();
        int max_nonzero_elems = sparcity * L * L;
        int counter = 0;
        if constexpr (std::is_same<Scalar,std::complex<double>>::value){
            while(counter < max_nonzero_elems){
                int i = rn::uniform_integer(0,L-1);
                int j = rn::uniform_integer(0,L-1);
                if(mat_dense(i,j) != 0.0){continue;}
                mat_dense(i,j) = rn::uniform_complex_1();
                counter += 1;
            }
        }else{
            while(counter < max_nonzero_elems){
                int i = rn::uniform_integer(0,L-1);
                int j = rn::uniform_integer(0,L-1);
                if(mat_dense(i,j) != 0.0){continue;}
                mat_dense(i,j) = rn::uniform_double_1();
                counter +=1;
            }
        }
        if (verbose){
            double final_sparcity = (double) (mat_dense.array().cwiseAbs() > 1e-15 )
                    .select(Eigen::MatrixXd::Ones(L,L),0).sum() / mat_dense.size();
            auto determinant = mat_dense.determinant();
            std::cout << "Matrix size = "      << std::setprecision(2)<< std::left << std::setw(6)  << mat_dense.rows()
                      << " Target sparcity = " << std::setprecision(6)<< std::left << std::setw(10) << sparcity
                      << " Final sparcity  = " << std::setprecision(6)<< std::left << std::setw(10) << final_sparcity
                      << " Determinant = "     << std::setprecision(6)<< std::left << std::setw(10) << determinant
                      << std::endl;
        }
        return mat_dense;


    }



};


#endif //TRAINING_MATRIX_GENERATOR_H
