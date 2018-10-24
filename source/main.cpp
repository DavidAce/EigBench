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
#include <general/class_tic_toc.h>
//#include <El.hpp>
//#include <flens/flens.cxx>
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


#include <unsupported/Eigen/KroneckerProduct>
#include <arpackpp/class_eigsolver_arpack.h>
#include <eigen3/benchmark_eigen3.h>
#include <lapacke/benchmark_lapacke.h>
#include <general/matrix_generator.h>
#undef complex

using namespace std::complex_literals;


/*! \brief Prints the content of a vector nicely */
template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
    if (!v.empty()) {
//        out << "[ ";
        std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, "\n"));
        out << '\n';
    }
    return out;
}




inline std::vector<Eigen::MatrixXcd> gen_manybody_spin(const Eigen::MatrixXcd &s, int sites) {
    std::vector<Eigen::MatrixXcd> S;
    Eigen::MatrixXcd I = Eigen::MatrixXcd::Identity(s.rows(),s.cols());
    Eigen::MatrixXcd tmp;
    for (int i = 0; i < sites; i++) {
        tmp = i == 0 ? s : I;
        for (int j = 1; j < sites; j++) {
            tmp = kroneckerProduct(tmp, i == j ? s : I).eval();
        }
        S.emplace_back(tmp);
    }
    return S;
}

template<typename T1, typename T2>
inline auto mod(const T1 x, const T2 y) {
//        return x >= 0 ? x % y : x % y + y;
    return (x % y + y) % y;

}

Eigen::MatrixXcd defineMatrix(int sites){
    double J = 0.9;
    double g = 1.1;
    Eigen::Matrix2cd sx = (Eigen::Matrix2cd() << 0.0, 1.0, 1.0, 0.0).finished();
    Eigen::Matrix2cd sy = (Eigen::Matrix2cd() << 0.0, -1.0i, 1.0i, 0.0).finished();
    Eigen::Matrix2cd sz = (Eigen::Matrix2cd() << 1.0, 0.0, 0.0, -1.0).finished();
    Eigen::Matrix2cd I  = (Eigen::Matrix2cd() << 1.0, 0.0, 0.0, 1.0).finished();
    auto SX = gen_manybody_spin(sx, sites);
    auto SY = gen_manybody_spin(sy, sites);
    auto SZ = gen_manybody_spin(sz, sites);
    Eigen::MatrixXcd H((int)std::pow(2,sites),(int)std::pow(2,sites));
    H.setZero();
    for (int pos=0; pos < sites; pos++){
        int i = mod(pos    , sites);
        int j = mod(pos + 1, sites);
        H += -J * SZ[i] * SZ[j] - 0.5 * g * (SX[i] + SX[j]);
    }
//    std::cout <<  "Matrix: \n" << std::setprecision(8) << H+H.adjoint() << std::endl;
    return H + H.adjoint();
}


template<typename Derived>
Eigen::SparseMatrix<typename Derived::Scalar> denseToSparse (Eigen::MatrixBase<Derived> &matrix, double prune_threshold = 0, bool prune = false){
    if (prune){
        return (matrix.array().cwiseAbs() < 1-prune_threshold).select(matrix,0).sparseView();

    }else{
        return matrix.sparseView();
    }
}








//
//Eigen::ArrayXcd eig_Flens(Eigen::MatrixXd matrix, class_tic_toc &timer){
//    // Read the example at
//    // http://apfel.mathematik.uni-ulm.de/~lehn/flens/flens/examples/tut04-page04.html
//    using namespace std;
//    using namespace flens;
//    int L = matrix.rows();
//    assert(L == matrix.cols());
//    typedef GeMatrix<FullStorage<double, ColMajor> >   Matrix;
//    typedef DenseVector<Array<double> >                Vector;
//
//    Matrix   VL(L, L), VR(L, L);// Left and right eigenvectors
//    Vector   wr(L), wi(L);      // The real and imaginary parts of the eigenvalues
//
//    typedef FullStorageView<double, ColMajor>  FSView;
//    typedef GeMatrix<FSView>                   GeMatrixView;
//    GeMatrixView  A    = FSView(L, L, matrix.data(), L);
//    timer.tic();
//    flens::lapack::ev(true, true, A, wr, wi, VL, VR);
//    timer.toc();
//
//    Eigen::ArrayXcd eigvals(L);
//    eigvals.real() = Eigen::Map<Eigen::ArrayXd>(wr.data(), L);
//    eigvals.imag() = Eigen::Map<Eigen::ArrayXd>(wi.data(), L);
//    return eigvals;
//}






// Lib elemental only does hermitian eigenvalue decomposition

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


int main () {
    int num_threads = 8;
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

    class_tic_toc t_eigen_rl_sparse       (true,6,"Eigen real sparse        ");
    class_tic_toc t_eigen_rl_dense        (true,6,"Eigen real dense         ");
    class_tic_toc t_eigen_cx_sparse       (true,6,"Eigen cplx sparse        ");
    class_tic_toc t_eigen_cx_dense        (true,6,"Eigen cplx dense         ");
    class_tic_toc t_armadillo_rl_sparse   (true,6,"Armadillo real sparse    ");
    class_tic_toc t_armadillo_rl_dense    (true,6,"Armadillo real dense     ");
    class_tic_toc t_armadillo_cx_sparse   (true,6,"Armadillo cplx sparse    ");
    class_tic_toc t_armadillo_cx_dense    (true,6,"Armadillo cplx dense     ");
    class_tic_toc t_arpackcpp_rl_sparse   (true,6,"Arpack++ real  sparse    ");
    class_tic_toc t_arpackcpp_rl_dense    (true,6,"Arpack++ real  dense     ");
    class_tic_toc t_arpackcpp_cx_sparse   (true,6,"Arpack++ cplx  sparse    ");
    class_tic_toc t_arpackcpp_cx_dense    (true,6,"Arpack++ cplx  dense     ");
//    class_tic_toc t_eleme(true,6,"Elemental ");
//    class_tic_toc t_flens(true,6,"Flens     ");
    class_tic_toc t_zheev(true,6,"zheev     ");
    class_tic_toc t_dgeev(true,6,"dgeev     ");
    class_tic_toc t_zheevd(true,6,"zheevd    ");
    class_tic_toc t_lapak4(true,6,"zheevd(nv)");


    matrix_generator matgen;
    std::vector<int> L_list     = {128,256,512};
    std::vector<int> nev_list   = {1,8};
    std::vector<Eigen::MatrixXcd> matrix_list;
    //Generate matrices
    for(auto L : L_list){
        matrix_list.emplace_back(matgen.make_hermitian_dense<std::complex<double>>(L));
    }


    std::vector<int> site_list  = {7,8,9};
    std::vector<int> nev_list   = {1,8};
    std::vector<Eigen::MatrixXcd> matrix_list;
    //Generate matrices
    for(auto sites : site_list){
        matrix_list.emplace_back(defineMatrix(sites));
    }

    using real = double;
    using cplx = std::complex<double>;

    // Start benchmarking eigs, i.e. few eigenvalues
    std::cout << std::setprecision(14);
    std::vector<std::tuple<int, int,double,double,double,double>> results_eigs;
    for (auto &matrix : matrix_list ){
        for (auto &nev : nev_list) {
            int L = matrix.rows();  //Linear size of the matrix
//            Eigen::ArrayXXcd eigvals_Flens(L,nev);
//            Eigen::ArrayXXcd eigvecs_dgeev(L,nev);
            Eigen::ArrayXXcd eigvecs_Armadillo_sparse(L,nev);
            Eigen::ArrayXXcd eigvecs_Arpackcpp_dense(L,nev);
            double sparcity = (double) (matrix.array().cwiseAbs() > 1e-12 )
                    .select(Eigen::MatrixXd::Ones(L,L),0).sum() / matrix.size();
            auto   matrix_sparse            = denseToSparse(matrix);
            Eigen::MatrixXcd matrix_dense   = matrix_sparse;
            eigvecs_Armadillo_sparse            = eigs_Armadillo_sparse      (matrix_dense  ,nev,t_armadillo_rl_sparse);
            eigvecs_Arpackcpp_dense             = eigs_Arpackcpp_dense<real> (matrix_dense  ,nev,t_arpackcpp_rl_dense);
            eigvecs_Arpackcpp_dense             = eigs_Arpackcpp_dense<cplx> (matrix_dense  ,nev,t_arpackcpp_cx_dense);

            results_eigs.emplace_back(
                    L,
                    nev,
                    sparcity,
                    t_armadillo_rl_sparse.get_last_time_interval(),
                    t_arpackcpp_rl_dense.get_last_time_interval(),
                    t_arpackcpp_cx_dense.get_last_time_interval()
            );
        }

    }


    std::cout << "Times for partial eigenvalue decomposition [seconds]" << std::endl;
    std::cout << std::setprecision(8)
            << "      " <<std::setw(6)  << std::left << "L"
            << "      " <<std::setw(6)  << std::left << "nev"
            << "      " <<std::setw(12) << std::left << "sparcity"
            << "      " <<std::setw(32) << std::left << "armadillo eigs_sym(real)"
            << "      " <<std::setw(32) << std::left << "arpackcpp eigs_sym(real)"
            << "      " <<std::setw(32) << std::left << "arpackcpp eigs_sym(cplx)"
              << std::endl;
    for(auto &res : results_eigs){
        std::cout << std::setprecision(8)
                << "      " << std::setw(6)  << std::left << std::get<0>(res)
                << "      " << std::setw(6)  << std::left << std::get<1>(res)
                << "      " << std::setw(12) << std::left << std::get<2>(res)
                << "      " << std::setw(32) << std::left << std::get<3>(res)
                << "      " << std::setw(32) << std::left << std::get<4>(res)
                << "      " << std::setw(32) << std::left << std::get<5>(res)
                << std::endl;
    }
    std::cout << std::endl << std::flush;


    std::cout << std::setprecision(14);
    std::vector<std::tuple<int, double,double,double,double>> results_eig;
    for (auto &matrix : matrix_list ){
        int L = matrix.rows();  //Linear size of the matrix
//            Eigen::ArrayXXcd eigvals_Flens(L,nev);
//            Eigen::ArrayXXcd eigvecs_dgeev(L,nev);
        Eigen::ArrayXXcd eigvecs_Eigen_sparse(L,L);
        Eigen::ArrayXXcd eigvecs_Eigen_dense (L,L);
        Eigen::ArrayXXcd eigvecs_Armad_dense (L,L);

        double sparcity = (double) (matrix.array().cwiseAbs() > 1e-12 )
                .select(Eigen::MatrixXd::Ones(L,L),0).sum()
                          / (double)matrix.size();
        auto   matrix_sparse            = denseToSparse(matrix);
        Eigen::MatrixXcd matrix_dense   = matrix_sparse;
        eigvecs_Eigen_sparse            = eig_Eigen(matrix_sparse , t_eigen_cx_sparse);
        eigvecs_Eigen_dense             = eig_Eigen(matrix_dense  , t_eigen_cx_dense);
        eigvecs_Armad_dense             = eig_Armadillo_dense(matrix_dense, t_armadillo_cx_dense);

        results_eig.emplace_back(
                L,
                sparcity,
                t_eigen_cx_sparse.get_last_time_interval(),
                t_eigen_cx_dense .get_last_time_interval(),
                t_armadillo_cx_dense .get_last_time_interval()
        );


    }


    std::cout << "Times for full diagonalization [seconds]" << std::endl;
    std::cout << std::setprecision(8)
            << "      " <<std::setw(6)  << std::left << "L"
            << "      " <<std::setw(12) << std::left << "sparcity"
            << "      " <<std::setw(48) << std::left << "eigen sparse SelfAdjointEigenSolver(cplx)"
            << "      " <<std::setw(48) << std::left << "eigen dense  SelfAdjointEigenSolver(cplx)"
            << "      " <<std::setw(48) << std::left << "armadillo eig_sym(cplx)"
            << std::endl;
    for(auto &res : results_eig){
        std::cout << std::setprecision(8)
                << "      " << std::setw(6)  << std::left << std::get<0>(res)
                << "      " << std::setw(12) << std::left << std::get<1>(res)
                << "      " << std::setw(48) << std::left << std::get<2>(res)
                << "      " << std::setw(48) << std::left << std::get<3>(res)
                << "      " << std::setw(48) << std::left << std::get<4>(res)
                << std::endl;
    }
    std::cout << std::endl << std::flush;





    // Now benchmark full diag
//
//    std::cout << std::setprecision(14);
//    std::vector<std::tuple<int,double,double,double,double>> results_eig;
//    for (auto sites : site_list){
//        Eigen::MatrixXcd matrix = defineMatrix(sites);
//        int L = matrix.rows();  //Linear size of the matrix
//        assert(L == std::pow(2,sites));
//
////        Eigen::ArrayXXcd eigvals_Flens(2*sites,gamVals);
//        std::cout << "sites = " << sites << std::endl;
//        Eigen::ArrayXXcd eig_Eigen_sparse(L,prunevals);
//        Eigen::ArrayXXcd eig_Eigen_dense (L,prunevals);
//        Eigen::ArrayXXcd eig_Armad_sparse(L,prunevals);
//        Eigen::ArrayXXcd eig_Armad_dense (L,prunevals);
//        Eigen::ArrayXXcd eigvals_dgeev   (L,prunevals);
//        double sparcity = (double) (matrix.array().cwiseAbs() > 1e-12 )
//                .select(Eigen::MatrixXd::Ones(matrix.rows(),matrix.cols()),0).sum()
//                          / (double)matrix.size();
//        std::cout << "sparcity      : " << std::setprecision(16) << sparcity << std::endl;
//        auto   matrix_sparse            = denseToSparse(matrix);
//        Eigen::MatrixXcd matrix_dense   = matrix_sparse;
//
//        double sparcity_dense = (double) (matrix_dense.array().cwiseAbs() > 1e-12 )
//                .select(Eigen::MatrixXd::Ones(matrix_dense.rows(),matrix_dense.cols()),0).sum()
//                                / (double)matrix_dense.size();
//        std::cout << "sparcity dense: " << std::setprecision(16) << sparcity_dense << std::endl;
//
//
////            std::cout << "Matrix : \n" << matrix_dense << std::endl;
//        assert(matrix_sparse.rows() == matrix_sparse.cols());
//        eig_Eigen_sparse.col(0)  = eig_Eigen_sparse            (matrix_sparse, t_eigen_sparse);
//        eig_Eigen_dense .col(0)  = eig_Eigen_dense             (matrix_dense , t_eigen_dense);
//        eig_Armad_dense.col(0)   = eigs_Armadillo_sparse(matrix_dense, t_armad_sparse);
//        eig_Armad_sparse .col(0) = eig_Armadillo_dense(matrix_dense, t_armad_dense);
//        results.emplace_back(
//                L,
//                t_eigen_sparse.get_last_time_interval(),
//                t_eigen_dense .get_last_time_interval(),
//                t_armad_sparse.get_last_time_interval(),
//                t_armad_dense .get_last_time_interval()
//                );
//    }
//
//
//
//    std::cout << std::setprecision(8)
//              << "      " <<std::setw(6)  << std::left << "L"
//              << "      " <<std::setw(12) << std::left << "eigen_sparse"
//              << "      " <<std::setw(12) << std::left << "eigen_dense "
//              << "      " <<std::setw(12) << std::left << "armad_sparse"
//              << "      " <<std::setw(12) << std::left << "armad_dense "
//              << std::endl;
//    for(auto &res : results){
//        std::cout << std::setprecision(8)
//                  << "      " << std::setw(6)  << std::left << std::get<0>(res)
//                  << "      " << std::setw(12) << std::left << std::get<1>(res)
//                  << "      " << std::setw(12) << std::left << std::get<2>(res)
//                  << "      " << std::setw(12) << std::left << std::get<3>(res)
//                  << "      " << std::setw(12) << std::left << std::get<4>(res)
//                  << std::endl;
//    }
//    std::cout << std::endl << std::flush;


    return 0;
}



//Eigen::IOFormat niceFormat(Eigen::StreamPrecision,0," ", "\n", "","","","");


//DISREGARD THIS BELOW

// numthreads = 1, using openblas, L=1000

//arma_config::blas             = true
//arma_config::mp_threads       = 10
//arma_config::openmp           = false
//arma_config::mp_threshold     = 320
//Eigen     : 205.952579 s
//        Armadillo : 27.623776 s
//        Elemental : 41.100007 s


//        Eigen     : 11.588257 s
//        Armadillo : 0.986238 s
//        Elemental : 1.237153 s
//        Flens     : 24.690452 s
//        zheev     : 11.857386 s
//        zgeev     : 8.149970 s
//        zheevd    : 3.084791 s
//        zheevd(nv): 0.677886 s

// numthreads 1
//        Eigen     : 13.381641 s
//        Armadillo : 1.802595 s
//        Elemental : 2.568189 s
//        Flens     : 25.920710 s
//        zheev     : 11.950945 s
//        zgeev     : 11.980822 s
//        zheevd    : 4.466797 s
//        zheevd(nv): 1.378681 s

//        Eigen     : 8.589689 s
//        Armadillo : 1.559550 s
//        Elemental : 2.235573 s
//        Flens     : 24.555642 s
//        zheev     : 11.822214 s
//        zgeev     : 11.620216 s
//        zheevd    : 4.213716 s
//        zheevd(nv): 1.345663 s