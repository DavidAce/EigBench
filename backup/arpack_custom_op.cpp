//
// Created by david on 2018-10-30.
//
#include <assert.h>
#include <iostream>
#include "matrix_product_dense.h"
#include "matrix_product_sparse.h"
#include "arpack_custom_op.h"

#include <arpackpp/arssym.h>
#include <arpackpp/arsnsym.h>
#include <arpackpp/arscomp.h>
#include <general/nmspc_type_check.h>
namespace tc = TypeCheck;
using namespace eigutils::eigSetting;



template<typename MatrixType>
void arpack_custom_op<MatrixType>::eigs() {
    solution.meta.eigvecs_found = false;
    solution.meta.eigvals_found = false;
    solution.eigvals.clear();
    solution.eigvecs.clear();
    nev_internal = std::min(matrix.rows()/2,solverConf.eigMaxNev);
    ncv_internal = std::max(solverConf.eigMaxNcv, 2+solverConf.eigMaxNev);
    ncv_internal = std::min(ncv_internal, matrix.rows());
    assert(ncv_internal >= solverConf.eigMaxNev + 2 and ncv_internal <= matrix.rows());
    assert(nev_internal >= 1 and nev_internal <= matrix.rows() / 2);

    // Stupidly enough, the only difference between the github and "apt" versions of arpack++,
    // is that the apt version only accepts char*, whereas the github one accepts string and char*.
    // For this reason we have to convert the ritz to a format both can take.
    solverConf.writeRitzChar();
    matrix.set_mode(solverConf.form);
    matrix.set_side(solverConf.side);
    // Calculate shift-inverse mat-vec mult operator by LU decomposition
    if(solverConf.shift == eigutils::eigSetting::Shift::ON){
        matrix.set_shift(solverConf.sigma);
        matrix.FactorOP();
    }

    assert(solverConf.confOK and "solverConf isn't ready!");
    // Dispatch to symmetric or nonsymmetric. If complex, there's only a nonsymmetric option available.
    if constexpr (std::is_same<Scalar,std::complex<double>>::value){
        this->eigs_comp();
    }else{
        if(solverConf.form == Form::SYMMETRIC){this->eigs_sym();}else {this->eigs_nsym();}
    }

    // The solution to  the eigenvalue equation Av = l*v is determined up to a constant phase factor.
    // By checking the first element in the v, one can compute the phase and remove it from all elements
    // of v.
    if (solverConf.remove_phase) {
        this->subtract_phase();
    }
}



template<typename MatrixType>
void arpack_custom_op<MatrixType>::eigs_sym() {
    if constexpr(std::is_same<Scalar,double>::value) {
        assert(solverConf.form       == Form::SYMMETRIC and "ERROR: solverConf not SYMMETRIC");
        assert(matrix.get_form()     == Form::SYMMETRIC and "ERROR: matrix not SYMMETRIC");
        ARSymStdEig<double, MatrixType> solver(
                matrix.rows(),
                nev_internal,
                &matrix,
                &MatrixType::MultAx,
                solverConf.ritz_char,
                ncv_internal,
                solverConf.eigThreshold,
                solverConf.eigMaxIter,
                nullptr);
        switch (solverConf.shift) {
            case Shift::OFF :
                break;
            case Shift::ON :
                solver.SetShiftInvertMode(std::real(solverConf.sigma), &matrix, &MatrixType::MultOPv);
                break;
        }


        this->find_solution(solver, nev_internal);
        this->copy_solution_symm(solver);
    }else{
        std::cerr << "ERROR: Called with wrong eigs_nsym() with wrong type: " << tc::type_name<MatrixType>() << '\n';
        exit(1);
    }
}


template<typename MatrixType>
void arpack_custom_op<MatrixType>::eigs_nsym() {
    if constexpr(std::is_same<Scalar, double>::value) {
        assert(solverConf.form == Form::NONSYMMETRIC and "ERROR: solverConf not NONSYMMETRIC");
        assert(matrix.get_form() == Form::NONSYMMETRIC and "ERROR: matrix not NONSYMMETRIC");
        if (nev_internal == 1) { nev_internal++; }
        ARNonSymStdEig<double, MatrixType> solver(
                matrix.rows(),
                nev_internal,
                &matrix,
                &MatrixType::MultAx,
                solverConf.ritz_char,
                ncv_internal,
                solverConf.eigThreshold,
                solverConf.eigMaxIter,
                nullptr);


        switch (solverConf.shift) {
            case Shift::OFF :
                break;
            case Shift::ON :
                solver.SetShiftInvertMode(std::real(solverConf.sigma), &matrix, &MatrixType::MultOPv);
                break;
        }

        this->find_solution(solver, nev_internal);
        this->copy_solution_nsym(solver);
    }else{
        std::cerr << "ERROR: Called with wrong eigs_nsym() with wrong type: " << tc::type_name<MatrixType>() << '\n';
        exit(1);
    }
}





template<typename MatrixType>
void arpack_custom_op<MatrixType>::eigs_comp() {
    if constexpr(std::is_same<Scalar, std::complex<double>>::value){
        ARCompStdEig<double, MatrixType> solver(
                matrix.rows(),
                nev_internal,
                &matrix,
                &MatrixType::MultAx,
                solverConf.ritz_char,
                ncv_internal,
                solverConf.eigThreshold,
                solverConf.eigMaxIter,
                nullptr);

        switch (solverConf.shift){
            case Shift::OFF :
                break;
            case Shift::ON :
                solver.SetShiftInvertMode(solverConf.sigma, &matrix, &MatrixType::MultOPv);
                break;
        }

        this->find_solution(solver, nev_internal);
        this->copy_solution_symm(solver);
    }else{
        std::cerr << "ERROR: Called with wrong eigs_nsym() with wrong type: " << tc::type_name<MatrixType>() << '\n';
        exit(1);
    }
}


// Explicit instantiations

template class arpack_custom_op<DenseMatrixProduct<double>>;
template class arpack_custom_op<DenseMatrixProduct<std::complex<double>>>;
template class arpack_custom_op<SparseMatrixProduct<double>>;
template class arpack_custom_op<SparseMatrixProduct<std::complex<double>>>;
