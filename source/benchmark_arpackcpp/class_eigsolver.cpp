//
// Created by david on 2018-10-29.
//

#include "class_eigsolver.h"
//#include "arpack_extra/arpack_custom_op.h"
//#include "arpack_extra/arpack_dense.h"
//#include "arpack_extra/arpack_dense_matProd.h"
//#include "arpack_extra/arpack_sparse_lu.h"

//using namespace arpack_extra::modes;




void class_eigsolver::conf_init  (const int L,
                                  const int nev,
                                  const int ncv,
                                  const std::complex<double> sigma          ,
                                  const arpack_extra::modes::Type type      ,
                                  const arpack_extra::modes::Form form      ,
                                  const arpack_extra::modes::Ritz ritz      ,
                                  const arpack_extra::modes::Side side      ,
                                  const arpack_extra::modes::Storage storage,
                                  const bool compute_eigvecs_               ,
                                  const bool remove_phase_
)
{
    using namespace arpack_extra::modes;
    bool is_shifted             = sigma == sigma;
    solverConf.compute_eigvecs  = compute_eigvecs_;
    solverConf.remove_phase     = remove_phase_;
    solverConf.eigMaxNev        = nev;
    solverConf.eigMaxNcv        = ncv;
    solverConf.shift            = is_shifted ? Shift::ON : Shift::OFF;
    solverConf.sigma            = sigma;
    solverConf.type             = type;
    solverConf.form             = form;
    solverConf.ritz             = ritz;
    solverConf.side             = side;
    solverConf.storage          = storage;

    if (solverConf.form == Form::NONSYMMETRIC){
        if (solverConf.eigMaxNev == 1) {
            solverConf.eigMaxNev = 2;
        }
    }

    if (ncv <= 2*nev ){
        solverConf.eigMaxNcv = std::max(L/32, nev*16);
        solverConf.eigMaxNcv = std::min(L/2, solverConf.eigMaxNcv);
//            solverConf.eigMaxNcv = std::min(L/2, 64);
    }

    solverConf.confOK = true;
}


















