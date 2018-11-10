//
// Created by david on 2018-10-30.
//

#ifndef CLASS_EIGSOLVER_BASE_H
#define CLASS_EIGSOLVER_BASE_H
#include <complex>
#include <vector>
#include <map>
#include <memory>
#include <general/class_tic_toc.h>
#include "nmspc_arpack_extra.h"
#define profile_arpack 0


template<typename MatrixType>
class arpack_base {
public:

    using Scalar = typename MatrixType::Scalar;
    class_tic_toc t_sol;
    class_tic_toc t_get;
    class_tic_toc t_sub;
    class_tic_toc t_all;


//    void shift_invert_eigvals(Scalar sigma);
    void subtract_phase();


    template <typename Derived>
    void find_solution(Derived &solver, int nev) {
        if (solverConf.compute_eigvecs) {
            solver.FindEigenvectors();
            solution.meta.eigvals_found = solver.EigenvaluesFound();
            solution.meta.eigvecs_found = solver.EigenvectorsFound();
            solution.meta.iter = solver.GetIter();
            solution.meta.n = solver.GetN();
            solution.meta.nev_found = std::min(nev, solver.GetNev());
            solution.meta.rows = solver.GetN();
            solution.meta.cols = std::min(nev, solver.GetNev());
        }else{
            solver.FindEigenvalues();
            solution.meta.eigvals_found = solver.EigenvaluesFound();
            solution.meta.iter = solver.GetIter();
            solution.meta.n = solver.GetN();
            solution.meta.nev_found = std::min(nev, solver.GetNev());
            solution.meta.rows = solver.GetN();
            solution.meta.cols = std::min(nev, solver.GetNev());

        }
    }

    template <typename Derived>
    void copy_solution_symm(Derived &solver) {
        solution.eigvals.assign(solver.RawEigenvalues() , solver.RawEigenvalues() + solution.meta.nev_found);
        if (solverConf.compute_eigvecs) {
            solution.eigvecs.assign(solver.RawEigenvectors(), solver.RawEigenvectors() + solution.meta.n * solution.meta.nev_found);
        }

    }

    template <typename Derived>
    void copy_solution_nsym(Derived &solver) {
        for (int j = 0; j < solver.ConvergedEigenvalues(); j++) {
            solution.eigvals.emplace_back(std::complex<double>(solver.EigenvalueReal(j), solver.EigenvalueImag(j)));
        }
        if(solverConf.compute_eigvecs){
            for (int j = 0; j < solver.ConvergedEigenvalues(); j++){
                for (int i = 0; i < solver.GetN(); i++){
                    solution.eigvecs.emplace_back(std::complex<double>(solver.EigenvectorReal(j,i), solver.EigenvectorImag(j,i)));
                }
            }
        }
    }


    MatrixType               &matrix;
    arpack_extra::SolverConf &solverConf;
    arpack_extra::Solution   &solution;

    arpack_base(
            MatrixType               &matrix_,
            arpack_extra::SolverConf &solverConf_,
            arpack_extra::Solution   &solution_);
};


#endif //EIGBENCH_CLASS_EIGSOLVER_BASE_H
