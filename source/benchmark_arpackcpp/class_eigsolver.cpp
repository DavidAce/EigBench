//
// Created by david on 2018-10-29.
//

#include "class_eigsolver.h"


//using namespace eigSetting;




void class_eigsolver::conf_init  (const int L,
                                  const int nev,
                                  const int ncv,
                                  const std::complex<double> sigma          ,
                                  const eigutils::eigSetting::Type type      ,
                                  const eigutils::eigSetting::Form form      ,
                                  const eigutils::eigSetting::Ritz ritz      ,
                                  const eigutils::eigSetting::Side side      ,
                                  const eigutils::eigSetting::Storage storage,
                                  const bool compute_eigvecs_               ,
                                  const bool remove_phase_
)
{
    using namespace eigutils::eigSetting;
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



    if (ncv <= nev or ncv >= L ){
        solverConf.eigMaxNcv = std::min(L, nev*3);
        solverConf.eigMaxNcv = std::max(8, solverConf.eigMaxNcv);
        solverConf.eigMaxNcv = std::min(L, solverConf.eigMaxNcv);
    }

    if (solverConf.form == Form::NONSYMMETRIC){
        if (solverConf.eigMaxNev == 1) {
            solverConf.eigMaxNev = 2;
        }
        solverConf.eigMaxNcv = std::min(L, solverConf.eigMaxNcv*2);

    }


    solverConf.confOK = true;
}


















