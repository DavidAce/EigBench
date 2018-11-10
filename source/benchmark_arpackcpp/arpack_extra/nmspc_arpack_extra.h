//
// Created by david on 2018-06-07.
//

#ifndef ARPACK_EXTRA_H
#define ARPACK_EXTRA_H

#include <iostream>
#include <array>
#include <map>
#include <complex>
#include <vector>


namespace arpack_extra{
    namespace modes{
        enum class Form{SYMMETRIC, NONSYMMETRIC};       // Real Symmetric, Real General or Complex General
        enum class Storage {DENSE,SPARSE};              // Dense or sparse matrix
        enum class Shift {ON,OFF};                      // Enable or disable shift invert
        enum class Ritz {LA,SA,LM,SM,LR,SR,LI,SI,BE};   // Choice of eigenvalue. LA is largest algebraic, and so on.
        enum class Side {L,R};                          // Left or right eigenvectors
        enum class Type {REAL,CPLX};                    // Real or complex, i.e. double or std::complex<double> matrix
  }




    template<typename Scalar_>
    class DenseType{
    public:
        using Scalar = Scalar_;
        int                       L;                // Leading dimension of the matrix.
        int                       N;                // Size of the matrix.
        std::vector<Scalar>       vals;              // pointer to an array that stores all the mat.

        void clear(){
            vals.clear();
        }
    };

    template<typename Scalar_>
    class SparseType{
    public:
        using Scalar = Scalar_;
        double              sparcity;
        int                 nnz;
        int                 L;            // Leading dimension of the matrix.
        int                 N;            // Size of the matrix.
        std::vector<int>    irow;         // pointer to an array that stores the row indices of the nonzeros in A.
        std::vector<int>    pcol;         // pointer to an array of pointers to the beginning of each column of A in valA.
        std::vector<Scalar> vals;         // pointer to an array that stores the nonzero values.
        void clear(){
            irow.clear();
            pcol.clear();
            vals.clear();
        }

    };








    class SolverConf{
    private:
//        void setup_strings(){
//
//        }
    public:
        bool confOK = false;
        using  MapType = std::map<arpack_extra::modes::Ritz, std::string>;
        MapType RitzToString;
        char ritz_char[3];


        arpack_extra::modes::Form           form        = arpack_extra::modes::Form::NONSYMMETRIC;
        arpack_extra::modes::Storage        storage     = arpack_extra::modes::Storage::DENSE;
        arpack_extra::modes::Shift          shift       = arpack_extra::modes::Shift::OFF;
        arpack_extra::modes::Side           side        = arpack_extra::modes::Side::R;
        arpack_extra::modes::Ritz           ritz        = arpack_extra::modes::Ritz::LM;
        arpack_extra::modes::Type           type        = arpack_extra::modes::Type::REAL;
        bool    compute_eigvecs                = false;
        bool    remove_phase                   = false;
        double  eigThreshold                   = 1e-12;
        int     eigMaxIter                     = 2000;
        int     eigMaxNev                      = 1;
        int     eigMaxNcv                      = 16;
        std::complex<double>  sigma            = std::numeric_limits<std::complex<double>>::quiet_NaN();     // Sigma value for shift-invert mode.

        SolverConf() {
            RitzToString = {
                    {arpack_extra::modes::Ritz::LA, "LA"},
                    {arpack_extra::modes::Ritz::SA, "SA"},
                    {arpack_extra::modes::Ritz::LM, "LM"},
                    {arpack_extra::modes::Ritz::SM, "SM"},
                    {arpack_extra::modes::Ritz::LR, "LR"},
                    {arpack_extra::modes::Ritz::SR, "SR"},
                    {arpack_extra::modes::Ritz::LI, "LI"},
                    {arpack_extra::modes::Ritz::SI, "SI"},
                    {arpack_extra::modes::Ritz::BE, "BE"}
            };
        }

        void writeRitzChar()
        // Writes ritz to string and checks that it is valid for the given problem.
        // The valid ritzes are stated in the arpack++ manual page 78.
        // Except, for REAL matrices it can't be LA or SA
        {
//            std::cout << "WRITING RITZCHAR! " << std::endl;
            using namespace arpack_extra::modes;
            if (type==Type::CPLX or form==Form::NONSYMMETRIC){
//                std::cout << "CPLX or NONSYMMETRIC" << std::endl;
                if (ritz==Ritz::LA or
                    ritz==Ritz::SA or
                    ritz==Ritz::BE
                    )
                {
                    std::cerr << "WARNING: Invalid ritz for nonsym problem: " << RitzToString.at(ritz) << std::endl;
                    if (ritz==Ritz::LA){ritz = Ritz::LR;}
                    if (ritz==Ritz::SA){ritz = Ritz::SR;}
                    if (ritz==Ritz::BE){ritz = Ritz::LM;}
                    std::cerr << "         Changed ritz to : " << RitzToString.at(ritz)<< std::endl;
                }
            }else if (type==Type::REAL and form==Form::SYMMETRIC) {
//                std::cout << "REAL AND SYMMETRIC" << std::endl;
                if (ritz==Ritz::LR or
                    ritz==Ritz::SR or
                    ritz==Ritz::LI or
                    ritz==Ritz::SI
                    )
                {
                    std::cerr << "WARNING: Invalid ritz for nonsym problem: " << RitzToString.at(ritz)<< std::endl;
                    if (ritz==Ritz::LR){ritz = Ritz::LA;}
                    if (ritz==Ritz::SR){ritz = Ritz::SA;}
                    if (ritz==Ritz::LI){ritz = Ritz::LM;}
                    if (ritz==Ritz::SI){ritz = Ritz::SM;}
                    std::cerr << "         Changed ritz to : " << RitzToString.at(ritz)<< std::endl;
                }
            }
//
//            switch (ritz)
//            {
//                case Ritz::LA:  ("LA").copy(ritz_char, 2);  break;
//                case Ritz::SA:  ("SA").copy(ritz_char, 2);  break;
//                case Ritz::LM:  ("LM").copy(ritz_char, 2);  break;
//                case Ritz::SM:  ("SM").copy(ritz_char, 2);  break;
//                case Ritz::LR:  ("LR").copy(ritz_char, 2);  break;
//                case Ritz::SR:  ("SR").copy(ritz_char, 2);  break;
//                case Ritz::LI:  ("LI").copy(ritz_char, 2);  break;
//                case Ritz::SI:  ("SI").copy(ritz_char, 2);  break;
//                case Ritz::BE:  ("BE").copy(ritz_char, 2);  break;
//
////                default:       ritz_char = 'LM';  break;
//            }


            RitzToString.at(ritz).copy(ritz_char, 3);
            confOK = true;
        }

    };


//    template<typename Scalar>
    class Solution{
    public:
        using Scalar = std::complex<double>;
        std::vector<Scalar> eigvecs;
        std::vector<Scalar> eigvals;

        struct Meta{
            int     rows            = 0;
            int     cols            = 0;
            int     iter            = 0;
            int     nev_found       = 0; // Found eigenvectors. aka cols.
            int     n               = 0; // Linear dimension of the input matrix to diagonalize, aka rows.
            int     counter         = 0;
            bool    eigvecs_found   = false;
            bool    eigvals_found   = false;
        } meta;


        const auto & ref_eig_vecs_vals() const{
            return std::make_pair(eigvals, eigvecs);
        }
        auto get_eig_vecs_vals() const{
            return std::make_pair(eigvals, eigvecs);
        }

//
//        const std::vector<Scalar> & ref_eigvecs() const{
//            if(meta.eigvecs_found) {
//                return eigvecs;
//            }else{
//                std::cerr << "Eigenvectors haven't been computed yet. Exiting " << std::endl;
//                exit(1);
//            }
//        }
//
//
//        const std::vector<Scalar> & ref_eigvals() const {
//            if(meta.eigvals_found) {
//                return eigvals;
//            }else{
//                std::cerr << "Eigenvalues haven't been computed yet. Exiting " << std::endl;
//                exit(1);
//            }
//        }
//
//
//        const std::vector<Scalar>  get_eigvecs() const{
//            if(meta.eigvecs_found) {
////        return std::vector<Scalar> (solution.RawEigenvectors(), solution.RawEigenvectors() + L * nev_found);
//
//                return eigvecs;
//            }else{
//                std::cerr << "Eigenvectors haven't been computed yet. Exiting " << std::endl;
//                exit(1);
//            }
//        }
//
//
//        const std::vector<Scalar> get_eigvals() const {
//            if(meta.eigvals_found) {
//                return eigvals;
//            }else{
//                std::cerr << "Eigenvalues haven't been computed yet. Exiting " << std::endl;
//                exit(1);
//            }
//        }


    };

}



#endif //ARPACK_EXTRA_H
