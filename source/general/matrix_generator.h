//
// Created by david on 2018-10-24.
//

#ifndef TRAINING_MATRIX_GENERATOR_H
#define TRAINING_MATRIX_GENERATOR_H
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <general/nmspc_random_numbers.h>
class matrix_generator {
public:
    matrix_generator () = default;

    template <typename Scalar>
    Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> make_hermitian_dense(int L, double sparcity = 0.01){
        Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> mat_dense(L,L);
        mat_dense.setZero();
        int max_nonzero_elems = sparcity * L * L;
        int counter = 0;
        if constexpr (std::is_same<Scalar,std::complex<double>>::value){
            while(counter < max_nonzero_elems){
                int i = rn::uniform_integer(0,L-1);
                int j = rn::uniform_integer(0,L-1);
                if (mat_dense(i,j) != 0){ continue;}
                if (mat_dense(i,j) == 0 and i != j){
                    mat_dense(i,j) = rn::uniform_complex_1();
                    mat_dense(j,i) = mat_dense.conjugate()(i,j);
                    counter +=2;
                 }
                 else if (mat_dense(i,j) == 0 and i == j){
                    mat_dense(i,j) = rn::uniform_double_1();
                    counter += 1;
                }
            }
        }else{
            while(counter < max_nonzero_elems){
                int i = rn::uniform_integer(0,L-1);
                int j = rn::uniform_integer(0,L-1);
                if (mat_dense(i,j) != 0){ continue;}
                if (mat_dense(i,j) == 0 and i != j){
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

        return mat_dense;
    }

    template <typename Scalar>
    Eigen::SparseMatrix<Scalar> make_hermitian_sparse(int L, double sparcity = 0.01){
        return make_hermitian_dense<Scalar>(L,sparcity).sparseView();
    }

};


#endif //TRAINING_MATRIX_GENERATOR_H
