/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre 
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.inria.fr)

The OpenMEEG software is a C++ package for solving the forward/inverse
problems of electroencephalography and magnetoencephalography.

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's authors,  the holders of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/

#include <cstdlib>
#include <sstream>
#include <limits>
#include <algorithm>

#include <matrix.h>
#include <symmatrix.h>
#include <sparse_matrix.h>
#include <vector.h>

namespace OpenMEEG {

    const Matrix& Matrix::set(const double d) {
        const size_t sz = size();
        for (size_t i=0; i<sz; i++)
            data()[i] = d;
        return *this;
    }

    Matrix::Matrix(const SymMatrix& A): LinOp(A.nlin(),A.ncol(),FULL,2),value(size()) {
        for (Index j=0; j<ncol(); ++j)
            for (Index i=0; i<nlin(); ++i)
                (*this)(i,j) = A(i,j);
    }

    Matrix::Matrix(const SparseMatrix& A): LinOp(A.nlin(),A.ncol(),FULL,2),value(size()) {
        set(0.0);
        for (SparseMatrix::const_iterator it=A.begin(); it!=A.end(); ++it)
            (*this)(it->first.first,it->first.second) = it->second;
    }

    Matrix::Matrix(const Vector& v,const Dimension M,const Dimension N): LinOp(M,N,FULL,2) {
        om_assert(M*N==v.size());
        value = v.value;
    }

    /// pseudo inverse
    Matrix Matrix::pinverse(const double tolrel) const {
        if (ncol()>nlin())
            return transpose().pinverse().transpose();

        Matrix U,V;
        SparseMatrix S;
        svd(U,S,V,false);
        // Following LAPACK The singular values of A, sorted so that S(i) >= S(i+1).
        const double atol = S(0,0)*((tolrel==0.0) ? std::numeric_limits<double>::epsilon() : tolrel);
        const double tol = std::max(nlin(),ncol())* atol;

        // Since here nlin()>=ncol(), S minimum dimension is ncol().

        Dimension rank = 0;
        for (Index i=0; i<ncol(); ++i)
            if (S(i,i)>tol)
                ++rank;

        if (rank==0) {
            Matrix result(ncol(),nlin());
            result.set(0.0);
            return result;
        }

        Matrix D(rank,rank);
        D.set(0.0);
        for (Index i=0; i<rank; ++i)
            D(i,i) = 1.0/S(i,i);

        // Keep only the first rank columns of U and V.

        const Matrix Vr(V.transpose(),rank);
        const Matrix Ur(U,rank);

        return Vr*D*Ur.transpose();
    }

    Matrix Matrix::transpose() const {
        Matrix result(ncol(),nlin());
        for (Index i=0; i<nlin(); ++i)
            for (Index j=0; j<ncol(); ++j)
                result(j,i) = (*this)(i,j);
        return result;
    }

    void Matrix::svd(Matrix& U,SparseMatrix& S,Matrix& V,const bool complete) const {
        // XXX output is now V.transpose() not V (now as in Lapack and numpy)
    #ifdef HAVE_LAPACK
        Matrix cpy(*this,DEEP_COPY);
        const Index mini = std::min(nlin(),ncol());
        const Index maxi = std::max(nlin(),ncol());
        U = Matrix(nlin(),nlin());
        U.set(0.0);
        S = SparseMatrix(nlin(),ncol());
        V = Matrix(ncol(),ncol());
        V.set(0.0);
        double* s = new double[mini];
        // int lwork = 4 *mini*mini + maxi + 9*mini; 
        // http://www.netlib.no/netlib/lapack/double/dgesdd.f :
        BLAS_INT* iwork = new BLAS_INT[8*mini];
        BLAS_INT lwork = 4*mini*mini+std::max(maxi,4*mini*mini+4*mini);
        BLAS_INT Info = 0;
        double *work = new double[lwork];
        if (complete) { // complete SVD
            DGESDD('A',sizet_to_int(nlin()),sizet_to_int(ncol()),cpy.data(),sizet_to_int(nlin()),s,U.data(),sizet_to_int(U.nlin()),V.data(),sizet_to_int(V.nlin()),work,lwork,iwork,Info);
        } else { // only first min(m,n)
            DGESDD('S',sizet_to_int(nlin()),sizet_to_int(ncol()),cpy.data(),sizet_to_int(nlin()),s,U.data(),sizet_to_int(U.nlin()),V.data(),sizet_to_int(V.nlin()),work,lwork,iwork,Info);
        }
        for (Index i=0; i<mini; ++i)
            S(i,i) = s[i];
        delete[] s;
        delete[] work;
        delete[] iwork;
        if (Info<0) {
            std::cout << "in svd: the "<< -Info << "-th argument had an illegal value." << std::endl;
        } else if (Info>0) {
            std::cout << "in svd: DBDSDC did not converge, updating process failed." << std::endl;
        }
    #else
        std::cerr<<"svd not implemented without blas/lapack"<<std::endl;
    #endif
    }

    Matrix Matrix::operator*(const SparseMatrix& mat) const {
        om_assert(ncol()==mat.nlin());
        Matrix out(nlin(),mat.ncol());
        out.set(0.0);

        for (SparseMatrix::const_iterator it=mat.begin(); it!=mat.end(); ++it) {
            const Index i = it->first.first;
            const Index j = it->first.second;
            const double val = it->second;
            for (Index k=0; k<nlin(); ++k)
                out(k,j) += (*this)(k,i)*val;
        }
        return out;
    }

    Matrix Matrix::operator*(const double x) const {
        Matrix C(nlin(),ncol());
        const size_t sz = size();
        for (size_t k=0; k<sz; ++k)
            C.data()[k] = data()[k]*x;
        return C;
    }

    Matrix Matrix::operator/(const double x) const {
        return (*this)*(1.0/x);
    }

    void Matrix::operator*=(const double x) {
        const size_t sz = size();
        for (size_t k=0; k<sz; ++k)
            data()[k] *= x;
    }

    void Matrix::operator/=(const double x) {
        *this *= 1.0/x;
    }

    void Matrix::info() const {
        if (nlin()==0 && ncol()==0) {
            std::cout << "Empty matrix" << std::endl;
            return;
        }

        std::cout << "Dimensions : " << nlin() << " x " << ncol() << std::endl;

        double minv = (*this)(0,0);
        double maxv = (*this)(0,0);
        Index mini = 0;
        Index maxi = 0;
        Index minj = 0;
        Index maxj = 0;

        //#pragma omp parallel for reduction()
        for (Index i=0; i<nlin(); ++i)
            for (Index j=0; j<ncol(); ++j) {
                const double value = (*this)(i,j);
                if (minv>value) {
                    minv = value;
                    mini = i;
                    minj = j;
                } else if (maxv<value) {
                    maxv = value;
                    maxi = i;
                    maxj = j;
                }
            }

        std::cout << "Min Value : " << minv << " (" << mini << "," << minj << ")" << std::endl;
        std::cout << "Max Value : " << maxv << " (" << maxi << "," << maxj << ")" << std::endl;
        std::cout << "First Values" << std::endl;

        for (Index i=0; i<std::min(nlin(),5U); ++i) {
            for (Index j=0; j<std::min(ncol(),5U); ++j)
                std::cout << (*this)(i,j) << " ";
            std::cout << std::endl;
        }
    }

    // =======
    // = IOs =
    // =======

    void Matrix::load(const char* filename) {
        maths::ifstream ifs(filename);
        try {
            ifs >> maths::format(filename,maths::format::FromSuffix) >> *this;
        } catch (maths::Exception& e) {
            ifs >> *this;
        }
    }

    void Matrix::save(const char* filename) const {
        maths::ofstream ofs(filename);
        try {
            ofs << maths::format(filename,maths::format::FromSuffix) << *this;
        } catch (maths::Exception& e) {
            ofs << *this;
        }
    }
}
