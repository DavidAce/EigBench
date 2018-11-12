//#include <IO/class_hdf5_file.h>

/*
 *  This example writes data to the HDF5 file.
 *  Number of processes is assumed to be 1 or multiples of 2 (up to 8)
 */




#include <iostream>
#include <iomanip>
#include <iterator>
#include <vector>
#include <list>
#include <complex>
#include <ctime>
#include <cstring>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/LU>

//#include <El.hpp>
#include <omp.h>
#ifdef OpenBLAS_AVAILABLE
#include <cblas.h>
#endif

#ifdef MKL_AVAILABLE
#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>
#include <mkl_service.h>
#include <mkl.h>
#endif

#include <general/nmspc_eigutils.h>
#include <general/class_tic_toc.h>
#include <benchmark_armadillo/benchmark_armadillo.h>
#include <benchmark_arpackcpp/benchmark_arpackcpp.h>
#include <benchmark_eigen3/benchmark_eigen3.h>
#include <benchmark_lapacke/benchmark_lapacke.h>
#include <benchmark_flens/benchmark_flens.h>
#include <general/matrix_generator.h>
#include <general/matrix_container.h>
#undef complex
using namespace std::complex_literals;





int main () {
    int num_threads = 1;
#ifdef OpenBLAS_AVAILABLE
    openblas_set_num_threads(num_threads);
    std::cout << "Using OpenBLAS with " << openblas_get_num_threads() << " thread(s)" << std::endl;
#endif

#ifdef MKL_AVAILABLE
    mkl_set_num_threads(num_threads);
        std::cout << "Using Intel MKL with " << num_threads << " thread(s)" << std::endl;
#endif

//    omp_set_num_threads(num_threads);
    Eigen::setNbThreads(num_threads);
    using namespace eigutils::eigSetting;
    using real = double;
    using cplx = std::complex<double>;
    rn::seed(3);

    int reps = 10;
    std::vector<int>     L_list           = {256,512,1024,2048};
    std::vector<double>  sparcity_list    = {0.1};
    std::vector<int>     nev_list         = {1,8,16,32};
    std::vector<int>     ncv_list         = {8,16,24,32,64,96,128,160,196};
    matrix_container matrices(L_list,sparcity_list,reps);

    benchmark_arpackcpp bench_arpackcpp;
    benchmark_armadillo bench_armadillo;
    std::cout << "Eigenvectors: true" << std::endl;
    matrices.generate_matrices<Type::REAL, Form::SYMMETRIC>();

    auto results_arpackcpp_real_symm_dense_si   = bench_arpackcpp.run_eigs_shift_invert<Type::REAL,Form::SYMMETRIC,Storage::DENSE>  (matrices,nev_list,ncv_list);
    results_arpackcpp_real_symm_dense_si.print_results("Arpack++  eigs real symmetric dense w shift invert");

    auto results_arpackcpp_real_symm_sparse_si   = bench_arpackcpp.run_eigs_shift_invert<Type::REAL,Form::SYMMETRIC,Storage::SPARSE>  (matrices,nev_list,ncv_list);
    results_arpackcpp_real_symm_sparse_si.print_results("Arpack++  eigs real symmetric sparse w shift invert");


    auto results_arpackcpp_real_symm_dense   = bench_arpackcpp.run_eigs<Type::REAL,Form::SYMMETRIC,Storage::DENSE>  (matrices,nev_list,ncv_list);
    results_arpackcpp_real_symm_dense.print_results("Arpack++  eigs real symmetric dense");
    auto results_arpackcpp_real_symm_sparse  = bench_arpackcpp.run_eigs<Type::REAL,Form::SYMMETRIC,Storage::SPARSE> (matrices,nev_list,ncv_list);
    results_arpackcpp_real_symm_sparse.print_results("Arpack++  eigs real symmetric sparse");
    matrices.generate_matrices<Type::REAL, Form::NONSYMMETRIC>();
    auto results_arpackcpp_real_nsym_dense   = bench_arpackcpp.run_eigs<Type::REAL,Form::NONSYMMETRIC,Storage::DENSE>  (matrices,nev_list,ncv_list);
    results_arpackcpp_real_nsym_dense.print_results("Arpack++  eigs real nonsymmetric dense" );
    auto results_arpackcpp_real_nsym_sparse  = bench_arpackcpp.run_eigs<Type::REAL,Form::NONSYMMETRIC,Storage::SPARSE> (matrices,nev_list,ncv_list);
    results_arpackcpp_real_nsym_sparse.print_results("Arpack++  eigs real nonsymmetric sparse");
    matrices.generate_matrices<Type::CPLX, Form::SYMMETRIC>();
    auto results_arpackcpp_cplx_symm_dense   = bench_arpackcpp.run_eigs<Type::CPLX,Form::SYMMETRIC,Storage::DENSE>  (matrices,nev_list,ncv_list);
    results_arpackcpp_cplx_symm_dense.print_results("Arpack++  eigs cplx symmetric dense" );
    auto results_arpackcpp_cplx_symm_sparse  = bench_arpackcpp.run_eigs<Type::CPLX,Form::SYMMETRIC,Storage::SPARSE> (matrices,nev_list,ncv_list);
    results_arpackcpp_cplx_symm_sparse.print_results("Arpack++  eigs cplx symmetric sparse");
    matrices.generate_matrices<Type::CPLX, Form::NONSYMMETRIC>();
    auto results_arpackcpp_cplx_nsym_dense   = bench_arpackcpp.run_eigs<Type::CPLX,Form::NONSYMMETRIC,Storage::DENSE>  (matrices,nev_list,ncv_list);
    results_arpackcpp_cplx_nsym_dense.print_results("Arpack++  eigs cplx nonsymmetric dense" );
    auto results_arpackcpp_cplx_nsym_sparse  = bench_arpackcpp.run_eigs<Type::CPLX,Form::NONSYMMETRIC,Storage::SPARSE> (matrices,nev_list,ncv_list);
    results_arpackcpp_cplx_nsym_sparse.print_results("Arpack++  eigs cplx nonsymmetric sparse");






//    auto results_armadillo_cplx_symm_sparse  = bench_armadillo.run_eigs_real_symm_sparse(reps,nev_list,L_list,sparcity_list);
//    print_results("Armadillo eigs real symmetric dense"    , results_armadillo_cplx_symm_sparse);


    return 0;




//
//
//
//
//
//    // Now do full diagonalization
//
//
//    std::cout << std::setprecision(14);
//    std::vector<std::tuple<int, double,double,double,double,double,double>> results_eig;
//    for (auto &matrix : matrix_list ){
//        int L = (int)matrix.rows();  //Linear size of the matrix
//
//
//        std::pair<Eigen::ArrayXcd,Eigen::ArrayXXcd> results_Eigen_sparse;
//        std::pair<Eigen::ArrayXcd,Eigen::ArrayXXcd> results_Eigen_dense ;
//        std::pair<Eigen::ArrayXcd,Eigen::ArrayXXcd> results_Armad_dense ;
//        std::pair<Eigen::ArrayXcd,Eigen::ArrayXXcd> results_zheev_dense ;
//        std::pair<Eigen::ArrayXcd,Eigen::ArrayXXcd> results_Flens_dense ;
//
//        double sparcity = (double) (matrix.array().cwiseAbs() > 1e-15 )
//                .select(Eigen::MatrixXd::Ones(L,L),0).sum() /matrix.size();
//        auto   matrix_sparse            = denseToSparse(matrix);
//        Eigen::MatrixXcd matrix_dense   = matrix_sparse;
//
//
//        results_Eigen_sparse            = bench_eigen3.eig_sym_Eigen(matrix_sparse , t_eigen_cx_sparse);
//        results_Eigen_dense             = bench_eigen3.eig_sym_Eigen(matrix_dense  , t_eigen_cx_dense);
//        results_Armad_dense             = bench_armadi.eig_sym_Armadillo_dense(matrix_dense, t_armadillo_cx_dense);
//        results_zheev_dense             = bench_lapack.eig_sym_Lapacke_zheev(matrix_dense, t_zheev);
////        results_Flens_dense             = bench_flens .eig_sym_Flens(matrix_dense, t_flens);
//
//        results_eig.emplace_back(
//                L,
//                sparcity,
//                t_eigen_cx_sparse.get_last_time_interval(),
//                t_eigen_cx_dense .get_last_time_interval(),
//                t_armadillo_cx_dense .get_last_time_interval(),
//                t_zheev   .get_last_time_interval(),
//                t_flens   .get_last_time_interval()
//        );
//    }


    std::cout << "Times for full diagonalization [seconds]" << std::endl;
//    std::cout << std::setprecision(8)
//            << "      " <<std::setw(6)  << std::left << "L"
//            << "      " <<std::setw(12) << std::left << "target_sparcity"
//            << "      " <<std::setw(42) << std::left << "eigen3 sparse SelfAdjointEigenSolver(cplx)"
//            << "      " <<std::setw(42) << std::left << "eigen3 dense  SelfAdjointEigenSolver(cplx)"
//            << "      " <<std::setw(42) << std::left << "armadillo dense eigs_sym(cplx)"
//            << "      " <<std::setw(42) << std::left << "dense zheev(cplx)"
////            << "      " <<std::setw(42) << std::left << "dense flens(cplx)"
//            << std::endl;
//    for(auto &res : results_eig){
//        std::cout << std::setprecision(8)
//                << "      " << std::setw(6)  << std::left << std::get<0>(res)
//                << "      " << std::setw(12) << std::left << std::get<1>(res)
//                << "      " << std::setw(42) << std::left << std::get<2>(res)
//                << "      " << std::setw(42) << std::left << std::get<3>(res)
//                << "      " << std::setw(42) << std::left << std::get<4>(res)
//                << "      " << std::setw(42) << std::left << std::get<5>(res)
//                << "      " << std::setw(42) << std::left << std::get<6>(res)
//                << std::endl;
//    }
    std::cout << std::endl << std::flush;

    return 0;
}



//Eigen::IOFormat niceFormat(Eigen::StreamPrecision,0," ", "\L", "","","","");
//Eigen::ArrayXcd eig_Elemental(Eigen::MatrixXcd matrix, class_tic_toc &timer){
//
////    Set a dummy communicator for libelemental, which works with mpi. This is used in xDMRG.
////    Note that this program is single threaded anyway.
//    El::Environment env;
//    [[maybe_unused]] El::mpi::Comm comm = El::mpi::COMM_WORLD;
//
//    using Scalar        = El::Complex<double>;
//    int L = matrix.rows();
//    El::HermitianEigCtrl<Scalar> ctrl;
////    void eig(int L, std::complex<double> *hermitian_matrix,  double * eigvals, std::complex<double>  * eigvecs){
//    El::Matrix<Scalar> elmatrix(L,L, static_cast<Scalar*>(matrix.data()),L);
//    El::Matrix<Scalar> vecs;
//    El::Matrix<double> vals;
//    timer.tic();
//    El::HermitianEig( El::LOWER, elmatrix, vals, vecs,ctrl );
//    timer.toc();
//    El::Print(vals);
//    Eigen::ArrayXcd eigvals(L);
////    std::move(vecs.Buffer(), vecs.Buffer() + vecs.Height()*vecs.Width(), eigvecs.data());
//    std::move(vals.Buffer(), vals.Buffer() + vals.Height()*vals.Width(), eigvals.data());
//    return eigvals;
//
//}