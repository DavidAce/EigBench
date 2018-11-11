//
// Created by david on 2018-10-25.
//

#ifndef EIGBENCH_CLASS_SPARSIFYER_H
#define EIGBENCH_CLASS_SPARSIFYER_H
#include <complex>
#include <vector>

template <typename Scalar>
class class_sparsifyer {
public:
    double              sparcity;
    int                 nnz;
    int                 nx;         // Leading dimension of the matrix.
    int                 n;          // Size of the matrix.
    std::vector<int>    irow;       // pointer to an array that stores the row indices of the nonzeros in A.
    std::vector<int>    pcol;       // pointer to an array of pointers to the beginning of each column of A in valA.
    std::vector<Scalar> valA;       // pointer to an array that stores the nonzero values.

    class_sparsifyer(const Scalar *data, int L,double sparcity_threshold = 1e-15, bool prune = true);
};


#endif //EIGBENCH_CLASS_SPARSIFYER_H
