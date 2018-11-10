//
// Created by david on 2018-10-30.
//
#include <assert.h>
#include "arpack_dense.h"
#include <arpackpp/ardsnsym.h>
#include <arpackpp/ardscomp.h>
#include <arpackpp/ardgcomp.h>
#include <arpackpp/ardssym.h>
#include <arpackpp/arseig.h>




using namespace arpack_extra::modes;


template<typename Scalar>
void arpack_dense<Scalar>::eigs() {
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
void arpack_dense<double>::eigs_sym() {
    ARdsSymMatrix<double>
            matrix_arpack(matrix.L,
                          matrix.vals.data());
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
void arpack_dense<double>::eigs_nsym() {
    if (nev_internal == 1){nev_internal ++;}
    ARdsNonSymMatrix<double,double>
            matrix_arpack(matrix.L,
                          matrix.vals.data());

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
void arpack_dense<std::complex<double>>::eigs_comp() {
    ARdsNonSymMatrix<std::complex<double>, double>
            matrix_arpack(matrix.L,
                          matrix.vals.data());
    ARluCompStdEig<double> solver(
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
            solver.SetShiftInvertMode(solverConf.sigma);
            break;
    }

    find_solution(solver, nev_internal);
}



// Explicit instantiations

template class arpack_dense<double>;
template class arpack_dense<std::complex<double>>;
