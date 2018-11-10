//
// Created by david on 2018-10-30.
//
#include <assert.h>
#include <iterator>
#include "arpack_sparse_lu.h"
#include <arpackpp/arlsmat.h>
#include <arpackpp/arlsnsym.h>
#include <arpackpp/arlscomp.h>
#include <arpackpp/arlscomp.h>
#include <arpackpp/arlssym.h>
#include <arpackpp/arseig.h>
#include "lcmatrxa.h"
#include "lnmatrxb.h"


using namespace arpack_extra::modes;
/*! \brief Prints the content of a vector nicely */
template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
    if (!v.empty()) {
//        out << "[ ";
        std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, "  "));
        out << '\n';
    }
    return out;
}

template<typename Scalar>
void arpack_sparse_lu<Scalar>::eigs() {
    solution.meta.eigvecs_found = false;
    solution.meta.eigvals_found = false;
    solution.eigvals.clear();
    solution.eigvecs.clear();
    nev_internal = std::min(matrix.L/2,solverConf.eigMaxNev);
    ncv_internal = std::max(solverConf.eigMaxNcv, 2+solverConf.eigMaxNev);
    ncv_internal = std::min(ncv_internal, matrix.L);
    assert(ncv_internal >= solverConf.eigMaxNev + 2 and ncv_internal <= matrix.L);
    assert(nev_internal >= 1 and nev_internal <= matrix.L / 2);

    // Stupidly enough, the only difference between the github and "apt" versions of arpack++, is that the apt version only accepts char*
    solverConf.writeRitzChar();
    if constexpr (std::is_same<Scalar,std::complex<double>>::value){
        this->eigs_comp();
    }else{
        if(solverConf.form == Form::SYMMETRIC){this->eigs_sym();}else {this->eigs_nsym();}
    }

    if (solverConf.remove_phase) {
        this->subtract_phase();
    }

}



template<>
void arpack_sparse_lu<double>::eigs_sym() {
    ARluSymMatrix<double>
            matrix_arpack(matrix.L,
            matrix.nnz,
            matrix.vals.data(),
            matrix.irow.data(),
            matrix.pcol.data());
    ARluSymStdEig<double> solver(
            nev_internal,
            matrix_arpack,
            solverConf.ritz_char,
            ncv_internal,
            solverConf.eigThreshold,
            solverConf.eigMaxIter,
            nullptr,
            true);

    switch (solverConf.shift){
            case Shift::OFF :
                solver.SetRegularMode();
                break;
            case Shift::ON :
                solver.SetShiftInvertMode(std::real(solverConf.sigma));
                break;
        }

    find_solution(solver, nev_internal);
}


template<>
void arpack_sparse_lu<double>::eigs_nsym() {

    int     nx;
    int     n;          // Dimension of the problem.
    int     nnz;        // Number of nonzero elements in A.
    int*    irow;       // pointer to an array that stores the row
    // indices of the nonzeros in A.
    int*    pcol;       // pointer to an array of pointers to the
    // beginning of each column of A in vector A.
    double* A;          // pointer to an array that stores the
    // nonzero elements of A.

    // Creating a 100x100 matrix.

    nx = 10;
    BlockTridMatrix(nx, n, nnz, A, irow, pcol);
    ARluNonSymMatrix<double, double> matrixT(n, nnz, A, irow, pcol);

    // Defining what we need: the four eigenvectors of A with largest magnitude.

    ARluNonSymStdEig<double> dprob(4, matrixT);



    ARluNonSymMatrix<double,double>
            matrix_arpack(matrix.L,
                          matrix.nnz,
                          matrix.vals.data(),
                          matrix.irow.data(),
                          matrix.pcol.data());

    ARluNonSymStdEig<Scalar> solver(
            nev_internal,
            matrix_arpack,
            solverConf.ritz_char,
            ncv_internal,
            solverConf.eigThreshold,
            solverConf.eigMaxIter,
            nullptr,
            true);

    switch (solverConf.shift){
        case Shift::OFF :
            solver.SetRegularMode();
            break;
        case Shift::ON :
            solver.SetShiftInvertMode(std::real(solverConf.sigma));
            break;
    }

    find_solution(solver, nev_internal);
}





template<>
void arpack_sparse_lu<std::complex<double>>::eigs_comp() {
    std::cout << "nonzeros : " << matrix.nnz  << std::endl;
    std::cout << "L        : " << matrix.L  << std::endl;
    std::cout << "N        : " << matrix.N  << std::endl;
    std::cout << "vals size: " << matrix.vals.size()  << std::endl;
    std::cout << "irow size: " << matrix.irow.size()  << std::endl;
    std::cout << "pcol size: " << matrix.pcol.size()  << std::endl;
    std::cout << "irow \n"     << matrix.irow << std::endl;
    std::cout << "pcol \n"     << matrix.pcol << std::endl;
    std::cout << "vals \n"     << matrix.vals << std::endl;



    int                nx;
    int                n;     // Dimension of the problem.
    int                nnz;   // Number of nonzero elements in A.
    int*               irow;  // pointer to an array that stores the row
    // indices of the nonzeros in A.
    int*               pcol;  // pointer to an array of pointers to the
    // beginning of each column of A in valA.
    arcomplex<double>* valA;  // pointer to an array that stores the
    // nonzero elements of A.
    // Creating a complex matrix.

    nx = 16;
    n  = nx*nx;
    CompMatrixA(nx, nnz, valA, irow, pcol);
    ARluNonSymMatrix<arcomplex<double>, double> A(n, nnz, valA, irow, pcol);

    // Defining what we need: the four eigenvectors of A with largest magnitude.

    ARluCompStdEig<double> solver(4L, A);
//
//    ARluNonSymMatrix<std::complex<double>, double>
//            matrix_arpack(matrix.L,
//                          matrix.nnz,
//                          matrix.vals.data(),
//                          matrix.irow.data(),
//                          matrix.pcol.data());
//
////    std::cout << "Arpack matrix: " << matrix_arpack.FactorAsI() << std::endl;
//    ARluCompStdEig<double> solver(
//            nev_internal,
//            matrix_arpack);

    switch (solverConf.shift){
        case Shift::OFF :
            solver.SetRegularMode();
            break;
        case Shift::ON :
            solver.SetShiftInvertMode(solverConf.sigma);
            break;
    }

    find_solution(solver, nev_internal);
}


template class arpack_sparse_lu<double>;
template class arpack_sparse_lu<std::complex<double>>;
