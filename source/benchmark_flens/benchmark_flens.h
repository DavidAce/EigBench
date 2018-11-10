//
// Created by david on 2018-10-24.
//

#ifndef EIGBENCH_BENCHMARK_FLENS_H
#define EIGBENCH_BENCHMARK_FLENS_H
#include <complex>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <flens/flens.cxx>
#include <general/class_tic_toc.h>
class benchmark_flens {
public:
    benchmark_flens() = default;

    auto eig_sym_Flens(Eigen::MatrixXcd matrix, class_tic_toc &timer){
        // Read the example at
        // http://apfel.mathematik.uni-ulm.de/~lehn/flens/flens/examples/tut04-page04.html
        using namespace std;
        using namespace flens;
        using ScalarType = std::complex<double>;
        int L = matrix.rows();
        assert(L == matrix.cols());

        typedef GeMatrix<FullStorage<ScalarType, ColMajor>>             Matrix;
        typedef DenseVector<Array<ScalarType>>                          ZVector;
        typedef DenseVector<Array<double> >                             DVector;
        typedef FullStorageView<ScalarType, ColMajor>                   FSView;
        typedef GeMatrix<FSView>                                        GeMatrixView;
        Matrix   VL(L, L), VR(L, L);// Left and right eigenvectors
        ZVector  w(L), wr(L), wi(L);
        ZVector    work;
        DVector    rWork;

        GeMatrixView  A    = FSView(L, L, matrix.data(), L);
        timer.tic();
//        flens::lapack::ev(true, true, A, wr, wi, VL, VR);
//        flens::lapack::ev(true,true, A.upper().hermitian() ,w, VL,VR);


//        template <typename MA, typename VW, typename MVL, typename MVR>
//        typename RestrictTo<IsComplexGeMatrix<MA>::value
//                            && IsComplexDenseVector<VW>::value
//                            && IsComplexGeMatrix<MVL>::value
//                            && IsComplexGeMatrix<MVR>::value,
//                typename RemoveRef<MA>::Type::IndexType>::Type
//        ev(bool     computeVL,
//           bool     computeVR,
//           MA       &&A,
//           VW       &&w,
//           MVL      &&VL,
//           MVR      &&VR);






//        lapack::ev(true,true, A, wr, wi, VL, VR, work);


        timer.toc();
//        Eigen::ArrayXXcd eigvecs(L);
//        eigvecs.real() = Eigen::Map<Eigen::ArrayXXd>(VL.data(), L,L);
//        eigvecs.imag() = Eigen::Map<Eigen::ArrayXXd>(VR.data(), L,L);
        return std::make_pair(Eigen::Map<Eigen::ArrayXcd>(w.data(), L),
                              Eigen::Map<Eigen::ArrayXXcd> (VR.data(), L,L));

//        return Eigen::Map<Eigen::ArrayXd>(w.data(), L);;





//
//        typedef complex<double>                             ZDouble;
//        typedef GeMatrix<FullStorage<ZDouble, ColMajor> >   ZGeMatrix;
//        typedef DenseVector<Array<ZDouble> >                ZDenseVector;
//        typedef DenseVector<Array<double> >                 DDenseVector;
//
//        const int L = 4;
//
//        ZGeMatrix      A(L, L);
//        DDenseVector   w(L);
//
//
//        A.upper() = ZDouble(1,0), ZDouble(1,1), ZDouble(2,1), ZDouble(5,2),
//                ZDouble(2,0), ZDouble(1,4), ZDouble(2,7),
//                ZDouble(3,0), ZDouble(1,4),
//                ZDouble(4,0);
//
//        cerr << "A.upper().hermitian() = "
//             << A.upper().hermitian() << endl;
//
//        ZDenseVector    work;
//        DDenseVector    rWork;
//        lapack::ev(true, A.upper().hermitian(), w, work, rWork);
//
//        cerr << "A = " << A << endl;
//        cerr << "w = " << w << endl;






    }
};


#endif //EIGBENCH_BENCHMARK_FLENS_H
