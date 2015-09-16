// keycpp.h -- Common Matlab functions implemented in C++
/** @file */

#ifndef KEYCPP_H_
#define KEYCPP_H_

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <string>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <boost/range.hpp>
#include <utility>
#include <algorithm>
#include <limits>
#include <ctime>

#ifdef __linux__ 
    #include <sys/time.h>
#elif _WIN32
    #include <time.h>
#else
    static_assert(false,"Unsupported Operating System");
#endif

#include <stdarg.h>
#include "vector_k.h"
#include "Matrix.h"
#include "_kiss_fft_guts.h"
#include "Spline.h"
#include "Figure.h"



/** \brief The keycpp namespace prevents KeyCpp functions and classes from interfering with
 *         other C++ libraries, for instance the std library.
 */
namespace keycpp
{
	static constexpr double pi = 3.1415926535897932384626433832795;
	static constexpr double eps = std::numeric_limits<double>::epsilon();
	static constexpr double Inf = std::numeric_limits<double>::infinity();
	static constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
}
	
namespace keycpp
{
	class KeyCppException : public std::runtime_error
	{
		public:
			KeyCppException(const std::string &msg) : std::runtime_error(msg){}
	};
	
	extern "C"{	
	/** \brief This provides a C interface to BLAS's double dot product function. */
	double ddot_(const int *N, const double *a, const int *inca, const double *b, const int *incb);
	
	/** \brief This provides a C interface to BLAS's complex double dot product function. */
	void zdotu_(std::complex<double> *result, const int *N, const std::complex<double> *a, const int *inca, const std::complex<double> *b, const int *incb);
	
	/** \brief This provides a C interface to BLAS's double vector addition function. */
	void daxpy_(const int *N, const double *alpha, const double *x, const int *incx, double *y, const int *incy);
	
	/** \brief This provides a C interface to BLAS's complex double vector addition function. */
	void zaxpy_(const int *N, const std::complex<double> *alpha, const std::complex<double> *x, const int *incx, std::complex<double> *y, const int *incy);
	
	/** \brief This provides a C interface to BLAS's double scalar-vector multiplication function. */
	void dscal_(const int *N, const double *alpha, double *x, const int *incx);
	
	/** \brief This provides a C interface to BLAS's complex double scalar-vector multiplication function. */
	void zscal_(const int *N, const std::complex<double> *alpha, std::complex<double> *x, const int *incx);
	
	/** \brief This provides a C interface to LAPACK's complex generalized eigenvalue solver. */
	void zggev_(const char *jobvl, const char *jobvr, const int *n, std::complex<double> *a,
		        const int *lda, std::complex<double> *b, const int *ldb, std::complex<double> *alpha,
		        std::complex<double> *beta, std::complex<double> *vl,
		        const int *ldvl, std::complex<double> *vr, const int *ldvr,
		        std::complex<double> *work, const int *lwork, double *rwork, int *info);
		        
	/** \brief This provides a C interface to LAPACK's double precision eigenvalue solver for a general matrix. */
	void dgeev_(const char *jobvl, const char *jobvr, const int *n, double *A, const int *lda,
	            double *wr, double *wi, double *vl, const int *ldvl, double *vr, const int *ldvr,
	            double *work, const int *lwork, int *info);
	            
	/** \brief This provides a C interface to LAPACK's complex eigenvalue solver for a general matrix. */
	void zgeev_(const char *jobvl, const char *, const int *n, std::complex<double> *A, const int *lda,
	            std::complex<double> *w, std::complex<double> *VL, const int *ldvl, std::complex<double> *VR, const int *ldvr,
	            std::complex<double> *work, const int *lwork, double *rwork, int *info);
		    
    /** \brief This provides a C interface to LAPACK's double precision reciprocal condition number estimator. */
    void dgecon_(const char *norm, const int *n, double *a,
                 const int *lda, const double *anorm, double *rcond,
                 double *work, int *iwork, int *info);
    
    /** \brief This provides a C interface to LAPACK's double precision LU decomposition function. */
    void dgetrf_(const int *m, const int *n, double *a, const int *lda,
                 int *lpiv, int *info);
    
    /** \brief This provides a C interface to LAPACK's double precision LU solver. */
    void dgetrs_(const char *trans, int *n, int *nrhs, double *a, const int *lda,
                 int *ipiv, double *b, int *ldb, int *info);
    
    /** \brief This provides a C interface to LAPACK's double precision norm function. */
    double dlange_(const char *norm, const int *m, const int *n,
                   const double *a, const int *lda, double *work);
		    
    /** \brief This provides a C interface to LAPACK's complex-valued reciprocal condition number estimator. */
    void zgecon_(const char *norm, const int *n, std::complex<double> *a,
                 const int *lda, const double *anorm, double *rcond,
                 std::complex<double> *work, double *rwork, int *info);
    
    /** \brief This provides a C interface to LAPACK's complex LU decomposition function. */
    void zgetrf_(const int *m, const int *n, std::complex<double> *a, const int *lda,
                 int *lpiv, int *info);
    
    /** \brief This provides a C interface to LAPACK's complex LU solver. */
    void zgetrs_(const char *trans, int *n, int *nrhs, std::complex<double> *a, const int *lda,
                 int *ipiv, std::complex<double> *b, int *ldb, int *info);
    
    /** \brief This provides a C interface to LAPACK's complex norm function. */
    double zlange_(const char *norm, const int *m, const int *n,
                   const std::complex<double> *a, const int *lda, double *work);
                   
    /** \brief This provides a C interface to LAPACK's double precision linear system solver. */
    void dgesv_(const int *n, const int *nrhs, double *a, const int *lda, int *ipiv,
                double *b, const int *ldb, const int *info);
                
    /** \brief This provides a C interface to LAPACK's double precision matrix inverse function. */
    void dgetri_(const int* n, double* A, const int* lda, const int* ipiv, double* work, const int* lwork, int* info);
    
    /** \brief This provides a C interface to LAPACK's complex matrix inverse function. */
    void zgetri_(const int* n, std::complex<double>* A, const int* lda, const int* ipiv,
                 std::complex<double>* work, const int* lwork, int* info);
                 
    /** \brief This provides a C interface to LAPACK's double precision SVD function. */
    void dgesvd_(const char *jobu, const char *jobvt, const int* m, const int* n, double* A,
                 const int* lda, double* S, double* U, const int* ldu, double* VT, const int* ldvt,
                 double* work, const int* lwork, int* info); 
                 
    /** \brief This provides a C interface to LAPACK's complex SVD function. */
    void zgesvd_(const char *jobu, const char *jobvt, const int* m, const int* n, std::complex<double>* A,
                 const int* lda, double* S, std::complex<double>* U, const int* ldu, std::complex<double>* VT, const int* ldvt,
                 std::complex<double>* work, const int* lwork, double *rwork, int* info);
	}

	matrix<std::complex<double> > eig(const matrix<std::complex<double> > &A,
	                                       const matrix<std::complex<double> > &B,
	                                       matrix<std::complex<double> > *vr_return = NULL,
	                                       matrix<std::complex<double> > *vl_return = NULL);

	matrix<std::complex<double> > eig(const matrix<std::complex<double> > &A,
	                                       matrix<std::complex<double> > *vr_return = NULL,
	                                       matrix<std::complex<double> > *vl_return = NULL);
	                                       
	matrix<std::complex<double>> eig(const matrix<double> &A,
                                          matrix<std::complex<double>> *vr_return = NULL,
                                          matrix<std::complex<double>> *vl_return = NULL);
	
	double rcond(const matrix<double> &A);
	double rcond(const matrix<std::complex<double>> &A);
	matrix<std::complex<double>> linsolve(const matrix<std::complex<double>>& A_in,
                                          const matrix<std::complex<double>>& b_in);
	matrix<double> linsolve(const matrix<double>& A_in,
                            const matrix<double>& b_in);
    matrix<double> inv(const matrix<double>& A_in);
    matrix<std::complex<double>> inv(const matrix<std::complex<double>>& A_in);
	
	template<class T, size_t dim> matrix<size_t,2> size(const matrix<T,dim> &A);
	template<class T> bool isnan(const T &a);
	template<> bool isnan<>(const std::complex<double> &a);
	template<class T,size_t dim> matrix<bool,dim> isnan(const matrix<T,dim> &A);
	
	template<class T,size_t dim>
	matrix<T,dim> eop(const matrix<T,dim> &A, T (*f)(const T&))
	{
	    matrix<T,dim> B;
	    B.resize(size(A));
	    for(size_t ii = 0; ii < A.numel(); ii++)
	    {
	        B(ii) = (*f)(A(ii));
	    }
	    return B;
	}
	
	template<class T,size_t dim>
	matrix<T,dim> eop(const matrix<std::complex<T>,dim> &A, T (*f)(const std::complex<T>&))
	{
	    matrix<T,dim> B;
	    B.resize(size(A));
	    for(size_t ii = 0; ii < A.numel(); ii++)
	    {
	        B(ii) = (*f)(A(ii));
	    }
	    return B;
	}
	
	template<class T, size_t dim>
	matrix<T,dim> eop(const matrix<T,dim> &A, T (*f)(T))
	{
	    matrix<T,dim> B;
	    B.resize(size(A));
	    for(size_t ii = 0; ii < A.numel(); ii++)
	    {
	        B(ii) = (*f)(A(ii));
	    }
	    return B;
	}
	
	template<class T,size_t dim>
	matrix<T,dim> eop(const matrix<std::complex<T>,dim> &A, T (*f)(std::complex<T>))
	{
	    matrix<T,dim> B;
	    B.resize(size(A));
	    for(size_t ii = 0; ii < A.numel(); ii++)
	    {
	        B(ii) = (*f)(A(ii));
	    }
	    return B;
	}
	
	template<class T, class U> struct observe
	{
		matrix<T,2>& y;
		matrix<U,2>& x_ode;
		size_t ii;

		observe(matrix<T,2> &p_y, matrix<U,2> &p_x_ode) : y(p_y), x_ode(p_x_ode), ii(0) { };

		void operator()(const matrix<T,2> &y_temp, U x_temp)
		{
		    if(ii >= y.size(1))
		    {
		        matrix<size_t,2> temp_size = size(y);
		        temp_size(0)++;
		        y.resize(temp_size);
		        x_ode.resize(temp_size(0));
		    }
		    x_ode(ii) = x_temp;
		    for(size_t jj = 0; jj < y_temp.length(); jj++)
		    {
			    y(ii,jj) = y_temp(jj);
			}
			ii++;
		}
	};
	
	template<class T>
	std::ostream& operator<<(std::ostream &out, const vector_k<T> &v1)
	{
        for(size_t ii = 0; ii < v1.size(); ii++)
        {
            out << v1[ii];
            out << std::endl;
        }
        return out;
    }

    /** \brief Returns a vector containing the product of all the elements in each
     *         column of the matrix A.
     */
	template<class T> matrix<T> prod(const matrix<T> &A)
	{
	    if(A.empty())
	    {
	        matrix<T> x(1);
	        x(0) = 1.0;
	        return x;
	    }
	    if(A.isVec())
	    {
		    matrix<T> B(1);
		    B(0) = A(0);
		    for(size_t jj = 1; jj < A.numel(); jj++)
		    {
			        B(0) *= A(jj);
		    }
		    return B;
	    }
	    else
	    {
		    matrix<T> B = A.row(0);
		    for(size_t jj = 0; jj < A.size(2); jj++)
		    {
		        for(size_t ii = 1; ii < A.size(1); ii++)
		        {
			        B(jj) *= A(ii,jj);
			    }
		    }
		    return B;
		}
	}
    
    template<class T, class U>
    decltype(std::declval<T>() % std::declval<U>()) mod(const T& a, const U& b)
    {
        decltype(std::declval<T>() % std::declval<U>()) r = a % b;
        return r < 0 ? r + b : r;
    }
	
	/** \brief Returns a matrix of row differences between adjacent rows.
	 *
	 *   TODO: Add recursive functionality and make sure it picks first non-singleton dimension.
	 *         Also, accept dimension as argument. See MATLAB docs.
     */
	template<class T> matrix<T> diff(const matrix<T> &A)
	{
	    if(A.size(1) <= 0 || A.size(2) <= 0)
	    {
	        throw KeyCppException("Cannot compute diff() on empty matrix!");
	    }
		
		if(!A.isVec())
		{
		    matrix<T> B(A.size(1)-1,A.size(2));
		    for(size_t ii = 0; ii < B.size(1); ii++)
		    {
		        for(size_t jj = 0; jj < B.size(2); jj++)
		        {
			        B(ii,jj) = A(ii+1,jj) - A(ii,jj);
			    }
		    }
		
		    return B;
        }
        else
        {
		    matrix<T> B(A.numel()-1);
		    for(size_t ii = 0; ii < B.numel(); ii++)
		    {
			    B(ii) = A(ii+1) - A(ii);
		    }
		
		    return B;
        }
	}

	template<class T,size_t dim> matrix<std::complex<T>,dim> conj(const matrix<std::complex<T>,dim> &A)
	{
		return eop(A,static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::conj));
	}

	template<class T,size_t dim> matrix<T,dim> real(const matrix<std::complex<T>,dim> &A)
	{
		return eop(A,static_cast<T (*)(const std::complex<T> &)>(&std::real));
	}

	template<class T, size_t dim> matrix<T,dim> imag(const matrix<std::complex<T>,dim> &A)
	{
		return eop(A,static_cast<T (*)(const std::complex<T> &)>(&std::imag));
	}

	template<class T,size_t dim> matrix<T,dim> abs(const matrix<T,dim> &A)
	{
		return eop(A,static_cast<T (*)(T)>(&std::abs));
	}

	template<class T,size_t dim> matrix<T,dim> abs(const matrix<std::complex<T>,dim> &A)
	{
		return eop(A,static_cast<T (*)(const std::complex<T> &)>(&std::abs));
	}
	
	inline std::complex<double> csqrt(const double& a)
	{
	    return std::sqrt((std::complex<double>)a);
	}
	
	inline std::complex<double> csqrt(const std::complex<double>& a)
	{
	    return std::sqrt(a);
	}
	
	/** \brief Return a vector containing the sine of each element of A.
	 */
	template<class T, size_t dim>
    matrix<std::complex<T>,dim> sin(const matrix<std::complex<T>,dim> &A)
    {
        return eop(A, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::sin<T>));
    }
	
	/** \brief Return a vector containing the sine of each element of A.
	 */
	template<class T, size_t dim> matrix<T,dim> sin(const matrix<T,dim> &A)
	{
		return eop(A,static_cast<T (*)(T)>(&std::sin));
	}
	
	/** \brief Return a vector containing the cosine of each element of A.
	 */
	template<class T, size_t dim>
    matrix<std::complex<T>,dim> cos(const matrix<std::complex<T>, dim> &A)
    {
        return eop(A, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::cos<T>));
    }
	
	/** \brief Return a vector containing the cos of each element of A.
	 */
	template<class T, size_t dim> matrix<T, dim> cos(const matrix<T,dim> &A)
	{
		return eop(A,static_cast<T (*)(T)>(&std::cos));
	}
	
	/** \brief Return a vector containing the tangent of each element of A.
	 */
	template<class T, size_t dim>
    matrix<std::complex<T>,dim> tan(const matrix<std::complex<T>,dim> &A)
    {
        return eop(A, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::tan<T>));
    }
	
	/** \brief Return a vector containing the tangent of each element of A.
	 */
	template<class T, size_t dim> matrix<T,dim> tan(const matrix<T,dim> &A)
	{
		return eop(A,static_cast<T (*)(T)>(&std::tan));
	}
	
	/** \brief Return a vector containing the arc cosine of each element of A.
	 */
	template<class T, size_t dim>
    matrix<std::complex<T>,dim> acos(const matrix<std::complex<T>,dim> &A)
    {
        return eop(A, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::acos<T>));
    }
	
	/** \brief Return a vector containing the arc cosine of each element of A.
	 */
	template<class T, size_t dim> matrix<T,dim> acos(const matrix<T,dim> &A)
	{
		return eop(A,static_cast<T (*)(T)>(&std::acos));
	}
	
	/** \brief Return a vector containing the arc sine of each element of A.
	 */
	template<class T, size_t dim>
    matrix<std::complex<T>,dim> asin(const matrix<std::complex<T>,dim> &A)
    {
        return eop(A, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::asin<T>));
    }
	
	/** \brief Return a vector containing the arc sine of each element of A.
	 */
	template<class T, size_t dim> matrix<T,dim> asin(const matrix<T,dim> &A)
	{
		return eop(A,static_cast<T (*)(T)>(&std::asin));
	}
	
	/** \brief Return a vector containing the exponential of each element of A.
	 */
	template<class T, size_t dim>
    matrix<std::complex<T>,dim> exp(const matrix<std::complex<T>,dim> &A)
    {
        return eop(A, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::exp<T>));
    }
	
	/** \brief Return a vector containing the exponential of each element of A.
	 */
	template<class T,size_t dim> matrix<T,dim> exp(const matrix<T,dim> &A)
	{
		return eop(A,static_cast<T (*)(T)>(&std::exp));
	}
	
	/** \brief Return a vector containing the natural logarithm of each element of A.
	 */
	template<class T, size_t dim>
    matrix<std::complex<T>,dim> log(const matrix<std::complex<T>,dim> &A)
    {
        return eop(A, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::log<T>));
    }
	
	/** \brief Return a vector containing the natural logarithm of each element of A.
	 */
	template<class T, size_t dim> matrix<T,dim> log(const matrix<T,dim> &A)
	{
		return eop(A,static_cast<T (*)(T)>(&std::log));
	}
	
	/** \brief Return a vector containing the base 10 logarithm of each element of A.
	 */
	template<class T, size_t dim>
    matrix<std::complex<T>,dim> log10(const matrix<std::complex<T>,dim> &A)
    {
        return eop(A, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::log10<T>));
    }
	
	/** \brief Return a vector containing the base 10 logarithm of each element of A.
	 */
	template<class T, size_t dim> matrix<T,dim> log10(const matrix<T,dim> &A)
	{
		return eop(A,static_cast<T (*)(T)>(&std::log10));
	}
	
	/** \brief Return a vector containing the sqrt of each element of A.
	 */
	template<class T, size_t dim>
    matrix<std::complex<T>,dim> sqrt(const matrix<std::complex<T>,dim> &A)
    {
        return eop(A, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::sqrt<T>));
    }
	
	/** \brief Return a vector containing the sqrt of each element of A.
	 */
	template<class T, size_t dim> matrix<T,dim> sqrt(const matrix<T,dim> &A)
	{
		return eop(A,static_cast<T (*)(T)>(&std::sqrt));
	}
	
	/** \brief Return a vector containing the csqrt of each element of A.
	 */
	template<class T, size_t dim>
    matrix<std::complex<T>,dim> csqrt(const matrix<std::complex<T>,dim> &A)
    {
        return eop(A, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&csqrt));
    }
	
	/** \brief Return a vector containing the csqrt of each element of A.
	 */
	template<class T, size_t dim> matrix<T,dim> csqrt(const matrix<T,dim> &A)
	{
		return eop(A,static_cast<T (*)(T)>(&csqrt));
	}
	
	template<class T> matrix<T> eye(const int &N)
	{
		matrix<T> A(N,N);
		for(int ii = 0; ii < N; ii++)
		{
			A(ii,ii) = 1.0;
		}
	
		return A;
	}
	
    /**  \brief Returns the number of elements along dimension dim.
     *   @details
     *    Returns the number of elements along dimension dim.
     *   @param[in] A matrix for which you want to know the size.
     *   @param[in] dim Dimension along which you want the number of elements.
     *              1 = number of rows, 2 = number of columns
     *   @return An integer value of the number of elements along the desired dimension. 
     */
	template<class T, size_t dim> int size(const matrix<T,dim> &A, const int &n)
	{
		return A.size(n);
	}
	
    /**  \brief Returns the size of matrix A.
     *   @details
     *    Returns the size of matrix A in a matrix_size_type variable.
     *   @usage auto msize = size(A); int num_rows = msize.rows; int num_cols = msize.cols;
     *   @param[in] A matrix for which you want to know the size.
     *   @return A matrix_size_type variable containing the number of rows and cols.
     */
	template<class T, size_t dim> matrix<size_t,2> size(const matrix<T,dim> &A)
	{
	    matrix<size_t,2> msize(dim);
	    for(size_t ii = 0; ii < dim; ii++)
	    {
	        msize(ii) = A.size(ii+1);
	    }
		return msize;
	}
	
    /**  \brief Returns a matrix of size M x N containing all zeros.
     *   @details
     *    Returns a matrix of size M x N containing all zeros.
     *   @param[in] M Number of rows.
     *   @param[in] N Number of columns.
     *   @return An M x N matrix containing zeros for each element. 
     */
	template<class T> matrix<T> zeros(const int &M, const int &N)
	{
		matrix<T> A(M,N);
		return A;
	}
	
    /**  \brief Returns a matrix of size N x N containing all zeros.
     *   @details
     *    Returns a matrix of size N x N containing all zeros.
     *   @param[in] N Number of rows and columns.
     *   @return An N x N matrix containing zeros for each element. 
     */
	template<class T> matrix<T,2> zeros(const int &N)
	{
		matrix<T,2> A(N,N);
		return A;
	}
	
    /**  \brief Returns a matrix of size M x N containing all ones.
     *   @details
     *    Returns a matrix of size M x N containing all ones.
     *   @param[in] M Number of rows.
     *   @param[in] N Number of columns.
     *   @return An M x N matrix containing ones for each element. 
     */
	template<class T> matrix<T> ones(const int &M, const int &N)
	{
		matrix<T> A(M,N);
		for(int ii = 0; ii < M; ii++)
		{
		    for(int jj = 0; jj < N; jj++)
		    {
		        A(ii,jj) = 1.0;
		    }
		}
		return A;
	}
	
    /**  \brief Returns a matrix of size N x N containing all ones.
     *   @details
     *    Returns a matrix of size N x N containing all ones.
     *   @param[in] N Number of rows and columns.
     *   @return An N x N matrix containing ones for each element. 
     */
	template<class T> matrix<T,2> ones(const int &N)
	{
		matrix<T> A(N,N);
		for(int ii = 0; ii < N; ii++)
		{
		    for(int jj = 0; jj < N; jj++)
		    {
		        A(ii,jj) = 1.0;
		    }
		}
		return A;
	}
	
	template<class T> matrix<T> diag(const std::initializer_list<T>& lst, const int &d = 0)
	{
		matrix<T> A(std::abs(d)+lst.size(),std::abs(d)+lst.size());
		int ii = 0;
		if(d != 0)
		{
		    for(const auto& l : lst)
		    {
			    if(d < 0)
	            {
		            A(ii+std::abs(d),ii) = l;
		        }
		        else
	            {
		            A(ii,ii+std::abs(d)) = l;
		        }
			    ii++;
		    }
		}
		else
		{
		    for(const auto& l : lst)
		    {
		        A(ii,ii) = l;
			    ii++;
		    }
		}
		return A;
	}
	
	template<class T> matrix<T,2> diag(const matrix<T,2> &A, const int &d = 0)
	{
	    if(A.empty())
	    {
	        throw KeyCppException("Cannot compute diagonal of empty matrix!");
	    }
	    if(A.isVec())
	    {
	        matrix<T,2> v1(std::abs(d)+A.numel(),std::abs(d)+A.numel());
		    if(d != 0)
		    {
		        for(int ii = 0; ii < A.numel(); ii++)
		        {
		            if(d < 0)
		            {
			            v1(ii+std::abs(d),ii) = A(ii);
			        }
			        else
		            {
			            v1(ii,ii+std::abs(d)) = A(ii);
			        }
		        }
		    }
		    else
		    {
		        for(int ii = 0; ii < A.numel(); ii++)
		        {
			        v1(ii,ii) = A(ii);
		        }
		    }
		    return v1;
	    }
	    else
	    {
	        int min_dim;
		    matrix<T> v1;
		    if(d == 0)
		    {
	            if(A.size(1) < A.size(2))
	            {
	                min_dim = A.size(1);
	            }
	            else
	            {
	                min_dim = A.size(2);
	            }
		        v1 = matrix<T>(min_dim);
		        for(int ii = 0; ii < min_dim; ii++)
		        {
			        v1(ii) = A(ii,ii);
		        }
		    }
		    else
		    {
		        if(d > 0)
		        {
		            min_dim = A.size(2) - std::abs(d);
		            v1 = matrix<T>(min_dim);
		            for(int ii = 0; ii < min_dim; ii++)
		            {
			            v1(ii) = A(ii,ii+std::abs(d));
		            }
		        }
		        else
		        {
		            min_dim = A.size(1) - std::abs(d);
		            v1 = matrix<T>(min_dim);
		            for(int ii = 0; ii < min_dim; ii++)
		            {
			            v1(ii) = A(ii+std::abs(d),ii);
		            }
		        }
		    }
		    return v1;
		}
	}
	
	template<class T> matrix<T,2> repmat(const matrix<T,2> &A, const int &m, const int &n)
	{
		matrix<T,2> B(m*A.size(1), n*A.size(2));
		for(int ii = 0; ii < m; ii++)
		{
			for(int jj = 0; jj < n; jj++)
			{
				for(int kk = 0; kk < A.size(1); kk++)
				{
					for(int mm = 0; mm < A.size(2); mm++)
					{
						B(ii*A.size(1) + kk, jj*A.size(2) + mm) = A(kk,mm);
					}
				}
			}
		}
		return B;
	}
	
	/** \brief Performs array multiplication on matrices A and B.
	 *
	 *  Each element of A is multiplied by each element of B. The matrix that is
	 *  returned is the same size as A and B.
	 */
	template<class T, class U, size_t dim> matrix<decltype(std::declval<T>()*std::declval<U>()),dim> times(const matrix<T,dim>& A, const matrix<U,dim>& B)
	{
	    if(A.empty() || B.empty())
	    {
	        throw KeyCppException("Cannot multiply an empty matrix!");
	    }
	    if(size(A) != size(B))
	    {
	        throw KeyCppException("Matrix dimensions must agree in times().");
	    }
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> C;
		C.resize(size(A));
		for(int ii = 0; ii < A.numel(); ii++)
		{
		    C(ii) = A(ii)*B(ii);
		}
		
		return C;
	}
	
	template<class T, class U, class V, size_t dim> matrix<decltype(std::declval<T>()*std::declval<U>()*std::declval<V>()),dim> times(const matrix<T,dim>& A, const matrix<U,dim>& B, const matrix<V,dim>& C)
	{
	    if(A.empty() || B.empty() || C.empty())
	    {
	        throw KeyCppException("Cannot multiply an empty matrix!");
	    }
	    if(size(A) != size(B))
	    {
	        throw KeyCppException("Matrix dimensions must agree in times().");
	    }
	    if(size(A) != size(C))
	    {
	        throw KeyCppException("Matrix dimensions must agree in times().");
	    }
		matrix<decltype(std::declval<T>()*std::declval<U>()*std::declval<V>()),dim> D;
		D.resize(size(A));
		for(int ii = 0; ii < A.numel(); ii++)
		{
		    D(ii) = A(ii)*B(ii)*C(ii);
		}
		
		return D;
	}
	
	template<class T, class U, class V, class W, size_t dim> matrix<decltype(std::declval<T>()*std::declval<U>()*std::declval<V>()*std::declval<W>()),dim> times(const matrix<T,dim>& A, const matrix<U,dim>& B, const matrix<V,dim>& C, const matrix<W,dim>& D)
	{
	    if(A.empty() || B.empty() || C.empty() || D.empty())
	    {
	        throw KeyCppException("Cannot multiply an empty matrix!");
	    }
	    if(size(A) != size(B))
	    {
	        throw KeyCppException("Matrix dimensions must agree in times().");
	    }
	    if(size(A) != size(C))
	    {
	        throw KeyCppException("Matrix dimensions must agree in times().");
	    }
	    if(size(A) != size(D))
	    {
	        throw KeyCppException("Matrix dimensions must agree in times().");
	    }
		matrix<decltype(std::declval<T>()*std::declval<U>()*std::declval<V>()*std::declval<W>()),dim> E;
		E.resize(size(A));
		for(int ii = 0; ii < A.numel(); ii++)
		{
		    E(ii) = A(ii)*B(ii)*C(ii)*D(ii);
		}
		
		return E;
	}
	
	template<class T, class U, class V, class W, class X, size_t dim> matrix<decltype(std::declval<T>()*std::declval<U>()*std::declval<V>()*std::declval<W>()*std::declval<X>()),dim> times(const matrix<T,dim>& A, const matrix<U,dim>& B, const matrix<V,dim>& C, const matrix<W,dim>& D, const matrix<X,dim>& E)
	{
	    if(A.empty() || B.empty() || C.empty() || D.empty() || E.empty())
	    {
	        throw KeyCppException("Cannot multiply an empty matrix!");
	    }
	    if(size(A) != size(B))
	    {
	        throw KeyCppException("Matrix dimensions must agree in times().");
	    }
	    if(size(A) != size(C))
	    {
	        throw KeyCppException("Matrix dimensions must agree in times().");
	    }
	    if(size(A) != size(D))
	    {
	        throw KeyCppException("Matrix dimensions must agree in times().");
	    }
	    if(size(A) != size(E))
	    {
	        throw KeyCppException("Matrix dimensions must agree in times().");
	    }
		matrix<decltype(std::declval<T>()*std::declval<U>()*std::declval<V>()*std::declval<W>()*std::declval<X>()),dim> F;
		F.resize(size(A));
		for(int ii = 0; ii < A.numel(); ii++)
		{
		    F(ii) = A(ii)*B(ii)*C(ii)*D(ii)*E(ii);
		}
		
		return F;
	}
	
	/** \brief Performs right array division on matrices A and B.
	 *
	 *  Each element of A is divided by each element of B. The matrix that is
	 *  returned is the same size as A and B. Equivalent to A./B in MATLAB.
	 */
	template<class T, class U, size_t dim> matrix<decltype(std::declval<T>()*std::declval<U>()),dim> rdivide(const matrix<T,dim>& A, const matrix<U,dim>& B)
	{
	    if(A.empty() || B.empty())
	    {
	        throw KeyCppException("Cannot rdivide an empty matrix!");
	    }
	    if(size(A) != size(B))
	    {
	        throw KeyCppException("Matrix dimensions must agree in rdivide().");
	    }
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> C;
		C.resize(size(A));
		for(int ii = 0; ii < A.numel(); ii++)
		{
		    C(ii) = A(ii)/B(ii);
		}
		return C;
	}
	
	/** \brief Performs left array division on matrices B and A.
	 *
	 *  Each element of A is divided by each element of B. The matrix that is
	 *  returned is the same size as B and A. Equivalent to B.\A in MATLAB.
	 */
	template<class T, class U> matrix<decltype(std::declval<T>()*std::declval<U>())> ldivide(const matrix<T>& B, const matrix<U>& A)
	{
	    if(A.size(1) <= 0 || A.size(2) <= 0 || B.size(1) <= 0 || B.size(2) <= 0)
	    {
	        throw KeyCppException("Cannot divide an empty matrix!");
	    }
	    if(A.size(1) != B.size(1) || A.size(2) != B.size(2))
	    {
	        throw KeyCppException("Matrix dimensions must agree in ldivide().");
	    }
		matrix<decltype(std::declval<T>()*std::declval<U>())> C(A.size(1),B.size(2));
		for(int ii = 0; ii < A.size(1); ii++)
		{
			for(int jj = 0; jj < A.size(2); jj++)
			{
				C(ii,jj) = A(ii,jj)/B(ii,jj);
			}
		}
		return C;
	}

	template<class T> double sign(const T &val)
	{
	    return (T(0) < val) - (val < T(0));
	}
	
	template<class T> matrix<double> sign(const matrix<T> &A)
	{
	    matrix<double> B(A.size(1),A.size(2));
	    for(size_t ii = 0; ii < B.numel(); ii++)
	    {
	        B(ii) = (T(0) < A(ii)) - (A(ii) < T(0));
	    }
	    return B;
	}
	
	template<class T> T angle(const std::complex<T> &x)
	{
		return arg(x);
	}
	
	template<class T,size_t dim>
    matrix<T,dim> angle(const matrix<std::complex<T>,dim> &A)
    {
        return eop(A, static_cast<T (*)(const std::complex<T> &)>(&std::arg<T>));
    }
	
	template<class T> matrix<T,2> max(const matrix<T,2> &A)
	{
	    if(!A.isVec())
	    {
	        matrix<T,2> v(A.size(2));
	        
	        for(size_t ii = 0; ii < A.size(2); ii++)
	        {
	            T temp = T();
	            bool flag = true;
	            for(size_t jj = 0; jj < A.size(1); jj++)
	            {
	                if(!isnan(A(jj,ii)) && (A(jj,ii) > temp || flag))
	                {
	                    temp = A(jj,ii);
	                    flag = false;
	                }
	            }
	            v(ii) = temp;
	        }
	        return v;
	    }
	    else
	    {
	        T a = T();
	        bool flag = true;
		    size_t index = 0;
		    for(size_t ii = 0; ii < A.numel(); ii++)
		    {
			    if(!isnan(A(ii)) && (A(ii) > a || flag))
			    {
				    a = A(ii);
				    index = ii;
				    flag = false;
			    }
		    }
	        matrix<T,2> v(1);
	        v(0) = A(index);
		    return v;
	    }
	}
	
	inline matrix<std::complex<double>,2> max(const matrix<std::complex<double>,2> &A)
	{
	    if(!A.isVec())
	    {
	        matrix<std::complex<double>,2> v(A.size(2));
	        
	        for(size_t ii = 0; ii < A.size(2); ii++)
	        {
                double a = nan("");
                double b = nan("");
                size_t index = 0;
	            for(size_t jj = 0; jj < A.size(1); jj++)
	            {
	                if(!keycpp::isnan(A(jj,ii)) && (std::abs(A(jj,ii)) > a || ((std::abs(A(jj,ii)) - std::abs(a)) < eps && angle(A(jj,ii)) > b) || (std::isnan(a) || std::isnan(b))))
	                {
                        a = abs(A(jj,ii));
                        b = angle(A(jj,ii));
                        index = jj;
	                }
	            }
	            v(ii) = A(index,ii);
	        }
	        return v;
	    }
	    else
	    {
            double a = nan("");
            double b = nan("");
            size_t index = 0;
		    for(size_t ii = 0; ii < A.numel(); ii++)
		    {
                if(!keycpp::isnan(A(ii)) && (std::abs(A(ii)) > a || ((std::abs(A(ii)) - std::abs(a)) < eps && angle(A(ii)) > b) || (std::isnan(a) || std::isnan(b))))
                {
                    a = abs(A(ii));
                    b = angle(A(ii));
                    index = ii;
                }
		    }
	        matrix<std::complex<double>,2> v(1);
	        v(0) = A(index);
		    return v;
	    }
	}
	
	template<class T> matrix<T,2> min(const matrix<T,2> &A)
	{
	    if(!A.isVec())
	    {
	        matrix<T,2> v(A.size(2));
	        
	        for(size_t ii = 0; ii < A.size(2); ii++)
	        {
	            T temp = T();
	            bool flag = true;
	            for(size_t jj = 0; jj < A.size(1); jj++)
	            {
	                if(!isnan(A(jj,ii)) && (A(jj,ii) < temp || flag))
	                {
	                    temp = A(jj,ii);
	                    flag = false;
	                }
	            }
	            v(ii) = temp;
	        }
	        return v;
	    }
	    else
	    {
	        T a = T();
	        bool flag = true;
		    size_t index = 0;
		    for(size_t ii = 0; ii < A.numel(); ii++)
		    {
			    if(!isnan(A(ii)) && (A(ii) < a || flag))
			    {
				    a = A(ii);
				    index = ii;
				    flag = false;
			    }
		    }
	        matrix<T,2> v(1);
	        v(0) = A(index);
		    return v;
	    }
	}
	
	inline matrix<std::complex<double>,2> min(const matrix<std::complex<double>,2> &A)
	{
	    if(!A.isVec())
	    {
	        matrix<std::complex<double>,2> v(A.size(2));
	        
	        for(size_t ii = 0; ii < A.size(2); ii++)
	        {
                double a = nan("");
                double b = nan("");
                size_t index = 0;
	            for(size_t jj = 0; jj < A.size(1); jj++)
	            {
	                if(!keycpp::isnan(A(jj,ii)) && (std::abs(A(jj,ii)) < a || ((std::abs(A(jj,ii)) - std::abs(a)) < eps && angle(A(jj,ii)) < b) || (std::isnan(a) || std::isnan(b))))
	                {
                        a = abs(A(jj,ii));
                        b = angle(A(jj,ii));
                        index = jj;
	                }
	            }
	            v(ii) = A(index,ii);
	        }
	        return v;
	    }
	    else
	    {
            double a = nan("");
            double b = nan("");
            size_t index = 0;
		    for(size_t ii = 0; ii < A.numel(); ii++)
		    {
                if(!keycpp::isnan(A(ii)) && (std::abs(A(ii)) < a || ((std::abs(A(ii)) - std::abs(a)) < eps && angle(A(ii)) < b) || (std::isnan(a) || std::isnan(b))))
                {
                    a = abs(A(ii));
                    b = angle(A(ii));
                    index = ii;
                }
		    }
	        matrix<std::complex<double>,2> v(1);
	        v(0) = A(index);
		    return v;
	    }
	}

	/** \brief Returns the transpose of matrix A.
	 */
	template<class T> matrix<T> transpose(const matrix<T> &A)
	{
		matrix<T> B(A.size(2),A.size(1));
		for(size_t ii = 0; ii < A.size(1); ii++)
		{
			for(size_t jj = 0; jj < A.size(2); jj++)
			{
				B(jj,ii) = A(ii,jj);
			}
		}
		return B;
	}

	/** \brief Returns the complex-conjugate transpose of matrix A.
	 */
	template<class T> matrix<T> ctranspose(const matrix<T> &A)
	{
		matrix<T> B(A.size(2),A.size(1));
		for(size_t ii = 0; ii < A.size(1); ii++)
		{
			for(size_t jj = 0; jj < A.size(2); jj++)
			{
				B(jj,ii) = conj(A(ii,jj));
			}
		}
		return B;
	}

	/** \brief Returns the complex-conjugate transpose of matrix A.
	 */
	template<> inline matrix<double> ctranspose<double>(const matrix<double> &A)
	{
		return transpose(A);
	}
	
	/** \brief Computes the sum of each column of A.
	 */
	template<class T> matrix<T,2> sum(const matrix<T,2> &A)
	{
	    if(A.empty())
	    {
	        throw KeyCppException("Cannot compute sum of empty matrix!");
	    }
	    if(!A.isVec())
	    {
		    matrix<T,2> v1(A.size(2));
		    for(size_t ii = 0; ii < v1.numel(); ii++)
		    {
			    v1(ii) = sum(A.col(ii));
		    }
		    return v1;
	    }
	    else
	    {
	        matrix<T,2> v1(1);
	        for(size_t ii = 0; ii < A.numel(); ii++)
		    {
			    v1(0) += A(ii);
		    }
		    return v1;
	    }
	}
	
	template<class T> matrix<T,2> linspace(const T &x1, const T &x2, const size_t &N)
	{
		matrix<T,2> x(N);
		if(N == 1)
		{
			x(0) = x2;
			return x;
		}

		T delta_x = (x2-x1)/(N-1);

		for(size_t ii = 0; ii < N; ii++)
		{
			x(ii) = x1 + ii*delta_x;
		}

		x(N-1) = x2;

		return x;
	}

	
    /**  \brief Produces a vector containing N values logarithmically spaced between
     *          10^(x1) and 10^(x2), inclusively.
     *   @details
     *    Produces a vector containing N values logarithmically spaced between
     *    10^(x1) and 10^(x2), inclusively.
     *   @param[in] x1 The base 10 logarithm of the minimum value.
     *   @param[in] x2 The base 10 logarithm of the maximum value.
     *   @param[in] N The number of values between 10^(x1) and 10^(x2). 
     *   @return A vector containing N logarithmically spaced values between
     *          10^(x1) and 10^(x2), inclusively. 
     */
	template<class T> matrix<T,2> logspace(const T &x1, const T &x2, const int &N)
	{
		matrix<T,2> x(N);
		if(N == 1)
		{
			x(0) = x2;
			return x;
		}

		T delta_x = (x2-x1)/(N-1);

		for(size_t ii = 0; ii < N; ii++)
		{
			x(ii) = pow(10.0,(x1 + delta_x*ii));
		}

		return x;
	}
	
	template<class T> matrix<T,2> unwrap(const matrix<T,2>& v1, const T &tol = pi)
	{
	    if(v1.size(1) > 1)
	    {
	        throw KeyCppException("This function requires the number of rows to be equal to 1.");
	    }
		matrix<T,2> v2(v1.size(2));
		v2(0) = v1(0);
		int correction = 0;
		for(size_t ii = 1; ii < v1.size(2); ii++)
		{
			if((v1(ii) - v1(ii-1)) > tol)
			{
				correction -= 1;
			}
			else if((v1(ii) - v1(ii-1)) < -tol)
			{
				correction += 1;
			}
			v2(ii) = v1(ii) + correction*2*pi;
		}
		return v2;
	}
	
	/** \brief Computes the mean of vector v1.
	 */
	template<class T> matrix<T,2> mean(const matrix<T,2>& v1)
	{
	    if(v1.isVec())
	    {
		    matrix<T,2> m(1);
		    double tot = 0.0;
		    for(size_t ii = 0; ii < v1.numel(); ii++)
		    {
			    m(0) += v1(ii);
			    tot += 1.0;
		    }
		    return m/tot;
		}
		else
		{
		    matrix<T,2> m(v1.size(2));
		    for(size_t ii = 0; ii < v1.size(2); ii++)
		    {
		        m(ii) = mean(v1.col(ii));
		    }
		    return m;
		}
	}
	
	template<class T, class U> T interp1_vec(const matrix<U,2> &x, const matrix<T,2> &y, const U &x_interp, std::string method = "linear", Extrap extrap = Extrap())
	{
		if(x.empty() || y.empty())
		{
			throw KeyCppException("Error in interp1! Cannot interpolate on an empty vector!");
		}
		if(x.length() != y.length())
		{
			throw KeyCppException("Error in interp1! Variables `x` and `y` have incompatible sizes!");
		}
		if(!x.isVec() || !y.isVec())
		{
			throw KeyCppException("Error in interp1! Variables `x` and `y` must be vectors!");
		}
		std::transform(method.begin(), method.end(), method.begin(), ::tolower);
		size_t N = x.length();
		T y2;
		if(method.compare("spline") == 0)
		{
			Spline<U,T> spline(N,x,y,extrap);
			spline.compute_spline();
			y2 = spline.J(x_interp);
		}
		else if(method.compare("linear") == 0)
		{
			for(size_t ii = 0; ii < (N-1); ii++)
			{
				if(x_interp == x(ii))
				{
					return y(ii);
				}
				else if((x_interp > x(ii) && x_interp < x(ii+1)) || (x_interp < x(ii) && x_interp > x(ii+1)))
				{
					return (y(ii) + (x_interp - x(ii))*(y(ii+1) - y(ii))/(x(ii+1) - x(ii)));
				}
			}

			if(x_interp == x(N-1))
			{
				return y(N-1);
			}
			
			if(extrap.isString)
			{
				if(extrap.extrap_string.compare("extrap") == 0)
				{
					if(x_interp < x(0))
					{
						return (y(0) + (x_interp - x(0))*(y(1) - y(0))/(x(1) - x(0)));
					}
					else
					{
						return (y(N-2) + (x_interp - x(N-2))*(y(N-1) - y(N-2))/(x(N-1) - x(N-2)));
					}
				}
				else
				{
					throw KeyCppException("ERROR! Could not interpolate!! Unknown std::string in interp1!");
					return nan("");
				}
			}
			else if(extrap.isDouble)
			{
				return extrap.extrap_val;
			}
			else
			{
				throw KeyCppException("ERROR! Could not interpolate!!");
				return nan("");
			}
		}
		else if(method.compare("nearest") == 0)
		{
			U min_val = 1e15;
			int index = -1;
			for(size_t ii = 0; ii < N; ii++)
			{
				if(std::abs(x(ii) - x_interp) < std::abs(min_val))
				{
					min_val = x(ii) - x_interp;
					index = ii;
				}
				else if(std::abs(x(ii) - x_interp) == std::abs(min_val) && (x(ii) - x_interp) > min_val)
				{
					min_val = x(ii) - x_interp;
					index = ii;
				}
			}
			if(index >= 0)
			{
				return y(index);
			}
			
			if(extrap.isString)
			{
				if(extrap.extrap_string.compare("extrap") == 0)
				{
					throw KeyCppException("ERROR in interp1! Cannot extrapolate using method `nearest`!");
				}
				else
				{
					throw KeyCppException("ERROR! Could not interpolate!! Unknown std::string in interp1!");
					return nan("");
				}
			}
			else if(extrap.isDouble)
			{
				return extrap.extrap_val;
			}
			else
			{
				throw KeyCppException("ERROR! Could not interpolate!!");
				return nan("");
			}
		}
		else
		{
			throw KeyCppException("Error in interp1! Unrecognized interpolation method!");
		}

		return y2;
	}

	template<class T, class U> matrix<T,2> interp1(const matrix<U,2> &x, const matrix<T,2> &y, const matrix<U,2> &x_interp, std::string method = "linear", Extrap extrap = Extrap())
	{
		if(x.empty() || y.empty() || x_interp.empty())
		{
			throw KeyCppException("Error in interp1! Cannot interpolate on an empty vector!");
		}
		if((x.length() != y.size(1) && !y.isVec()) || (x.length() != y.length() && y.isVec()))
		{
			throw KeyCppException("Error in interp1! Variables `x` and `y` have incompatible sizes!");
		}
		std::transform(method.begin(), method.end(), method.begin(), ::tolower);
		if(y.isVec())
		{
		    if(x_interp.isVec())
		    {
		        size_t N = x.numel();
		        size_t N_int = x_interp.numel();
		        matrix<T,2> y2(N_int);

		        if(method.compare("spline") == 0)
		        {	
			        Spline<U,T> spline(N,x,y,extrap);
			        spline.compute_spline();
			        for(size_t ii = 0; ii < N_int; ii++)
			        {
				        y2(ii) = spline.J(x_interp(ii));
			        }
		        }
		        else if(method.compare("linear") == 0 || method.compare("nearest") == 0)
		        {
			        for(size_t ii = 0; ii < N_int; ii++)
			        {
				        y2(ii) = interp1_vec(x, y, x_interp(ii),method,extrap);
			        }
		        }
		        else
		        {
			        throw KeyCppException("Error in interp1! Unrecognized interpolation method!");
		        }

		        return y2;
	        }
	        else
	        {
		        matrix<T,2> y2(x_interp.size(1),x_interp.size(2));

		        for(size_t kk = 0; kk < x_interp.size(2); kk++)
		        {
			        y2.col(kk) = interp1(x,y,x_interp.col(kk),method, extrap);
		        }

		        return y2;
	        }
	    }
	    else
	    {
		    if(x_interp.isVec())
		    {
		        matrix<T,2> y2(x_interp.numel(),y.size(2));

		        for(size_t kk = 0; kk < y.size(2); kk++)
		        {
			        y2.col(kk) = interp1(x,y.col(kk),x_interp,method, extrap);
		        }

		        return y2;
		    }
		    else
		    {
		        throw KeyCppException("This mode of interp1 is currently not supported!");
		    }
	    }
	}

	template<class T, class U> matrix<T,2> interp1(const matrix<U,2> &x, const matrix<T,2> &y, const U &x_interp, std::string method = "linear", Extrap extrap = Extrap())
	{
		if(x.empty() || y.empty())
		{
			throw KeyCppException("Error in interp1! Cannot interpolate on an empty vector!");
		}
		if(!x.isVec())
		{
		    throw KeyCppException("First input to interp1 must be vector!");
		}
		if(!y.isVec())
		{
		    if(x.length() != y.size(1))
		    {
			    throw KeyCppException("Error in interp1! Variables `x` and `y` have incompatible sizes!");
		    }
		
            matrix<T,2> y2(y.size(2));

            for(size_t kk = 0; kk < y.size(2); kk++)
            {
	            y2(kk) = interp1_vec(x,y.col(kk),x_interp,method, extrap);
            }

            return y2;
        }
        else
        {
		    if(x.length() != y.length())
		    {
			    throw KeyCppException("Error in interp1! Variables `x` and `y` have incompatible sizes!");
		    }
            matrix<T,2> y2(1);
            y2(0) = interp1_vec(x,y,x_interp,method, extrap);
            return y2;
        }
	}

	template<class U, class T> matrix<T,2> trapz(const matrix<U,2> &eta, const matrix<T,2> &integrand)
	{
		if(eta.empty() || integrand.empty())
		{
			throw KeyCppException("Error in trapz()! Empty vector/matrix supplied!");
		}
		if(!eta.isVec())
		{
		    throw KeyCppException("Error in trapz()! First input must be a vector!");
		}
		if(eta.numel() != integrand.size(1) && integrand.size(1) > 1)
		{
			throw KeyCppException("Error in trapz()! Vector and matrix sizes are not compatible!");
		}
		if(eta.numel() != integrand.size(2) && integrand.size(1) <= 1)
		{
			throw KeyCppException("Error in trapz()! Vector and matrix sizes are not compatible!");
		}
		
		size_t N;
		matrix<T,2> z;
		if(eta.numel() == integrand.size(1))
		{
		    N = integrand.size(2);
		    z = matrix<T,2>(N);
		    for(size_t ii = 0; ii < N; ii++)
		    {
			    size_t M = eta.numel();
                T sum = 0.0;

                for(size_t jj = 0; jj < M-1; jj++)
                {
                    sum += (eta(jj+1) - eta(jj))*(integrand(jj+1,ii) + integrand(jj,ii));
                }
                z(ii) = 0.5*sum;
		    }
		}
		else
		{
		    N = 1;
		    z = matrix<T,2>(N);
			
		    size_t M = eta.numel();
            T sum = 0.0;

            for(size_t jj = 0; jj < M-1; jj++)
            {
                sum += (eta(jj+1) - eta(jj))*(integrand(jj+1) + integrand(jj));
            }
            z(0) = 0.5*sum;
		}
		return z;
	}
	
    /// Method by: Fornberg B. Calculation of weights in finite difference formulas. SIAM Review 1998; 40(3):685–691.
    inline matrix<double> calc_weights(double x_loc, matrix<double> x, int q, int m)
    {
        matrix<double> c(q+1,m+1);

        double c1 = 1.0;
        double c4 = x(0) - x_loc;

        c(0,0) = 1.0;

        for(int ii = 1; ii < q+1; ii++)
        {
            int mn;
            if(ii < m)
            {
                mn = ii;
            }
            else
            {
                mn = m;
            }
            double c2 = 1.0;
            double c5 = c4;
            c4 = x(ii) - x_loc;
            for(int jj = 0; jj < ii; jj++)
            {
                double c3 = x(ii) - x(jj);
                c2 = c2*c3;
                if(jj == ii-1)
                {
                    for(int kk = mn; kk > 0; kk--)
                    {
                        c(ii,kk) = c1*(((double)kk)*c(ii-1,kk-1) - c5*c(ii-1,kk))/c2;
                    }
                    c(ii,0) = -c1*c5*c(ii-1,0)/c2;
                }
                for(int kk = mn; kk > 0; kk--)
                {
                    c(jj,kk) = (c4*c(jj,kk) - ((double)kk)*c(jj,kk-1))/c3;
                }
                c(jj,0) = c4*c(jj,0)/c3;
            }
            c1 = c2;
        }
        
        return c;
    }
    
    inline matrix<size_t> calc_s(int N, int q, bool periodic)
    {
        matrix<size_t> s(N+1);
        if(!periodic)
        {
            if(mod(q,2) == 0)
            {
                for(int ii = 0; ii < q/2; ii++)
                {
                    s(ii) = 0;
                }
                for(int ii = q/2; ii < N-q/2+1; ii++)
                {
                    s(ii) = ii - q/2;
                }
                for(int ii = N-q/2+1; ii < N+1; ii++)
                {
                    s(ii) = N-q;
                }
            }
            else
            {
                for(int ii = 0; ii < (q-1)/2; ii++)
                {
                    s(ii) = 0;
                }
                for(int ii = (q-1)/2; ii < N-(q+1)/2+1; ii++)
                {
                    s(ii) = ii - (q-1)/2;
                }
                for(int ii = N-(q+1)/2+1; ii < N+1; ii++)
                {
                    s(ii) = N-q;
                }
            }
        }
        else
        {
            if(mod(q,2) == 0)
            {
                for(int ii = 0; ii < N+1; ii++)
                {
                    s(ii) = mod(ii - q/2,N+1);
                }
            }
            else
            {
                for(int ii = 0; ii < N+1; ii++)
                {
                    s(ii) = mod(ii - (q-1)/2,N+1);
                }
            }
        }
        return s;
    }

    inline matrix<double> calc_weights_mat(matrix<double> x, size_t q, size_t n)
    {
        matrix<double> C(x.numel(),x.numel());
        auto s = calc_s(x.numel()-1,q,false);
        for(size_t ii = 0; ii < x.numel(); ii++)
        {
            matrix<double> temp1(q+1);
            for(size_t jj = 0; jj < q+1; jj++)
            {
                temp1(jj) = x(s(ii) + jj);
            }
            auto temp2 = calc_weights(x(ii),temp1,q,n);
            
            for(size_t jj = 0; jj < q+1; jj++)
            {
                C(ii,s(ii) + jj) = temp2(jj,n);
            }
        }
        return C;
    }


	template<class T, class U> matrix<T,2> diffxy(const matrix<U,2> &eta, const matrix<T,2> &u, const int &index = 2)
	{
		if(eta.empty() || u.empty())
		{
			throw KeyCppException("Error in diffxy()! Empty matrix supplied!");
		}
		
		if(!eta.isVec())
		{
		    if(eta.size(1) != u.size(1) && eta.size(2) != u.size(2))
		    {
			    throw KeyCppException("Error in diffxy()! Matrix sizes are not compatible!");
		    }
            
		    size_t N = u.size(1);
		    size_t P = u.size(2);

		    matrix<T,2> du(N,P);
		    if(index == 1)
		    {
			    for(size_t p = 0; p < P; p++)
			    {
			        int q;
			        if(N > 5)
			        {
			            q = 4;
			        }
			        else
			        {
			            q = N-1;
			        }
			    
                    auto C = calc_weights_mat(eta.col(p),q,1);
                    
                    du.col(p) = C*u.col(p);
			    }
		    }
		    else
		    {
			    for(size_t n = 0; n < N; n++)
			    {
			        int q;
			        if(P > 5)
			        {
			            q = 4;
			        }
			        else
			        {
			            q = P-1;
			        }
			        
                    auto C = calc_weights_mat(eta.row(n),q,1);
                    
                    du.row(n) = transpose(C*transpose(u.row(n)));
			    }
		    }

		    return du;
	    }
	    else
	    {
	        if(eta.numel() != u.size(1) && eta.numel() != u.size(2))
            {
                throw KeyCppException("Error in diffxy()! Vector and matrix sizes are not compatible!");
            }
            
	        int q;
	        if(numel(eta) > 5)
	        {
	            q = 4;
	        }
	        else
	        {
	            q = numel(eta)-1;
	        }
            
            auto C = calc_weights_mat(eta,q,1);
            
            size_t N = u.size(1);
            size_t P = u.size(2);

            matrix<T,2> du(N,P);
            if(N == eta.numel())
            {
                for(size_t p = 0; p < P; p++)
                {
                    du.col(p) = C*u.col(p);
                }
            }
            else
            {
                for(size_t ii = 0; ii < N; ii++)
                {
                    du.row(ii) = transpose(C*transpose(u.row(ii)));
                }
            }

            return du;
	    }
	}

	template<class T> matrix<std::complex<double>,2> fft(const matrix<T,2> &u, int N = -1)
	{
		if(u.empty())
		{
			throw KeyCppException("Error in fft()! Empty vector supplied!");
		}
		if(!u.isVec())
		{
			throw KeyCppException("Error in fft()! u must have a singleton dimension!");
		}
		
		if(N < 0)
		{
			N = u.length();
		}
		
		kiss_fft_cpx *cx_in = new kiss_fft_cpx[N];
		kiss_fft_cpx *cx_out = new kiss_fft_cpx[N];

		matrix<std::complex<double>,2> u_hat(N);

		for(int ii = 0; ii < N; ii++)
		{
			cx_in[ii].r = real((std::complex<double>)u(ii));
			cx_in[ii].i = imag((std::complex<double>)u(ii));
		}

		kiss_fft_cfg cfg = kiss_fft_alloc(N,false,NULL,NULL);
		kiss_fft(cfg,cx_in,cx_out);

		for(int ii = 0; ii < N; ii++)
		{
			u_hat(ii) = std::complex<double>((double)cx_out[ii].r,(double)cx_out[ii].i);
		}

		free(cfg);
		delete [] cx_in;
		delete [] cx_out;

		return u_hat;
	}

	template<class T> matrix<std::complex<double>,2> ifft(const matrix<T,2> &u, int N = -1)
	{
		if(u.empty())
		{
			throw KeyCppException("Error in ifft()! Empty vector supplied!");
		}
		if(!u.isVec())
		{
			throw KeyCppException("Error in ifft()! u must have a singleton dimension!");
		}
		
		if(N < 0)
		{
			N = u.length();
		}
		
		kiss_fft_cpx *cx_in = new kiss_fft_cpx[N];
		kiss_fft_cpx *cx_out = new kiss_fft_cpx[N];

		matrix<std::complex<double>,2> u_hat(N);

		for(int ii = 0; ii < N; ii++)
		{
			cx_in[ii].r = real((std::complex<double>)u(ii));
			cx_in[ii].i = imag((std::complex<double>)u(ii));
		}

		kiss_fft_cfg cfg = kiss_fft_alloc(N,true,NULL,NULL);
		kiss_fft(cfg,cx_in,cx_out);

		for(int ii = 0; ii < N; ii++)
		{
			u_hat(ii) = std::complex<double>((double)cx_out[ii].r,(double)cx_out[ii].i)/((double)N);
		}

		free(cfg);
		delete [] cx_in;
		delete [] cx_out;

		return u_hat;
	}
	
	template<class T>
	matrix<T> fftshift(matrix<T> A, size_t dim = 0)
	{
	    if((dim == 0 && !A.isVec()) || dim > 2)
	    {
	        throw KeyCppException("Error in fftshift!");
	    }
	    if(dim == 0)
	    {
	        size_t N = A.numel();
	        matrix<T> B(N);
	        size_t m = std::ceil(N/2);
	        for(size_t ii = 0; ii < N; ii++)
	        {
	            if(ii >= N-m)
	            {
	                B(ii-(N-m)) = A(ii);
	            }
	            else
	            {
	                B(ii+m) = A(ii);
	            }
	        }
	        return B;
	    }
	    else
	    {
	        matrix<T> B(A.size(1),A.size(2));
	        if(dim == 1)
	        {
	            for(size_t ii = 0; ii < A.size(2); ii++)
	            {
	                B.col(ii) = fftshift(A.col(ii));
	            }
	        }
	        else
	        {
	            for(size_t ii = 0; ii < A.size(1); ii++)
	            {
	                B.row(ii) = fftshift(A.row(ii));
	            }
	        }
	        return B;
	    }
	}
	
	template<class T>
	matrix<T> ifftshift(matrix<T> A, size_t dim = 0)
	{
	    if((dim == 0 && !A.isVec()) || dim > 2)
	    {
	        throw KeyCppException("Error in ifftshift!");
	    }
	    if(dim == 0)
	    {
	        size_t N = A.numel();
	        matrix<T> B(N);
	        size_t m = std::floor(N/2);
	        for(size_t ii = 0; ii < N; ii++)
	        {
	            if(ii >= m)
	            {
	                B(ii-m) = A(ii);
	            }
	            else
	            {
	                B(ii+N-m) = A(ii);
	            }
	        }
	        return B;
	    }
	    else
	    {
	        matrix<T> B(A.size(1),A.size(2));
	        if(dim == 1)
	        {
	            for(size_t ii = 0; ii < A.size(2); ii++)
	            {
	                B.col(ii) = ifftshift(A.col(ii));
	            }
	        }
	        else
	        {
	            for(size_t ii = 0; ii < A.size(1); ii++)
	            {
	                B.row(ii) = ifftshift(A.row(ii));
	            }
	        }
	        return B;
	    }
	}
}
	
namespace boost
{
    namespace numeric
    {
        namespace odeint
        {
            template<class T>
            struct is_resizeable<keycpp::matrix<T,2> >
            {
                typedef boost::true_type type;
                const static bool value = type::value;
            };
        }
    }
}

namespace keycpp
{
	template<class T,class Y>
	struct ODE_type
	{
		matrix<T,2> t;
		matrix<Y,2> y;
	};

	template<class T, class U, class F>
	ODE_type<U,T> ode45(F odeClass, const std::initializer_list<U> &x_span, matrix<T,2> ICs, double abs_tol = 1.0e-10, double rel_tol = 1.0e-6)
	{
		if(x_span.size() <= 0)
		{
			throw KeyCppException("Error in ode45()! x_span cannot be empty!");
		}
		if(ICs.empty())
		{
			throw KeyCppException("Error in ode45()! Must provide initial conditions!");
		}
		if(!ICs.isVec())
		{
			throw KeyCppException("Error in ode45()! Initial conditions must have one singleton dimension!");
		}
		if(x_span.size() != 2)
		{
			throw KeyCppException("Error in ode45()! Invalid vector x_span!");
		}

		U x0 = *(x_span.begin());
		U xf = *(x_span.end()-1);
		U delta_x0 = (xf-x0)/1000.0;
		matrix<T,2> y_temp(2,ICs.length());
		matrix<U,2> x_temp(2);

		{
			using namespace boost::numeric::odeint;
			integrate_adaptive(make_controlled<runge_kutta_dopri5<matrix<T,2> > >(abs_tol, rel_tol), odeClass, ICs, x0, xf, delta_x0, observe<T,U>(y_temp,x_temp));
		}

        ODE_type<U,T> ans;
        ans.t = x_temp;
        ans.y = y_temp;
		
		return ans;
	}
	
	template<class T, class U, class F>
	matrix<T> ode45(F odeClass, matrix<U,2> x_ode, matrix<T,2> ICs, double abs_tol = 1.0e-10, double rel_tol = 1.0e-6)
	{
		if(x_ode.empty())
		{
			throw KeyCppException("Error in ode45()! Vector x_ode cannot be empty!");
		}
		if(ICs.empty())
		{
			throw KeyCppException("Error in ode45()! Must provide initial conditions!");
		}
		if(!ICs.isVec())
		{
			throw KeyCppException("Error in ode45()! Initial conditions must be a vector!");
		}
		if(x_ode.length() < 2 || !x_ode.isVec())
		{
			throw KeyCppException("Error in ode45()! Invalid vector x_ode!");
		}

		U delta_x0 = x_ode(1) - x_ode(0);
		matrix<T,2> y(x_ode.length(),ICs.length());
		matrix<U,2> x(x_ode.length());

		{
			using namespace boost::numeric::odeint;
			integrate_times(make_dense_output<runge_kutta_dopri5<matrix<T,2> > >(abs_tol, rel_tol), odeClass, ICs, x_ode.begin(), x_ode.end(), delta_x0, observe<T,U>(y,x));
		}
		
		return y;
	}
	
	inline void set(Figure &h, std::string property, double val)
	{
		h.set(property,val);
	}
	
	inline void set(Figure &h, std::string property, std::string val)
	{
		h.set(property,val);
	}
	
	inline void set(Figure &h, std::string property, std::initializer_list<size_t> list)
	{
		h.set(property,list);
	}
	
	inline void print(Figure &h, std::string pterm, std::string pfilename)
	{
		h.print(pterm,pfilename);
	}
	
	template<class T>
	struct Sort_Matrix
	{
		matrix<T> B;
		matrix<size_t> index;
	};
	
	template<class T> Sort_Matrix<T> sort(const matrix<T> &AA, const size_t &dim = 2, std::string method = "ascend")
	{
		std::transform(method.begin(), method.end(), method.begin(), ::tolower);
		if(method.compare("ascend") != 0 && method.compare("descend") != 0)
		{
			throw KeyCppException("Invalid sort method!");
		}
		if(AA.empty())
		{
			throw KeyCppException("Tried to sort empty matrix!");
		}
		matrix<T> A(AA);
		
		if(!A.isVec())
		{
		    bool swapped = true;
		    T temp;
		    size_t temp_i;
		    matrix<T> B(A.size(1), A.size(2));
		    matrix<size_t> index(A.size(1),A.size(2));
		    if(dim == 2)
		    {
			    for(size_t ii = 0; ii < A.size(1); ii++)
			    {
				    for(size_t jj = 0; jj < A.size(2); jj++)
				    {
					    index(ii,jj) = ii;
				    }
			    }
			    for(size_t jj = 0; jj < A.size(2); jj++)
			    {
				    swapped = true;
				    while(swapped)
				    {
					    swapped = false;
					    for(size_t ii = 1; ii < A.size(1); ii++)
					    {
						    if(((A(ii-1,jj)) > (A(ii,jj)) && method.compare("ascend") == 0) || ((A(ii-1,jj)) < (A(ii,jj)) && method.compare("descend") == 0))
						    {
							    temp = A(ii-1,jj);
							    A(ii-1,jj) = A(ii,jj);
							    A(ii,jj) = temp;
							    temp_i = index(ii-1,jj);
							    index(ii-1,jj) = index(ii,jj);
							    index(ii,jj) = temp_i;
							    swapped = true;
						    }
					    }
				    }
			    }
			    for(size_t ii = 0; ii < A.size(2); ii++)
			    {
				    for(size_t jj = 0; jj < A.size(2); jj++)
				    {
					    B(ii,jj) = A(index(ii,jj),jj);
				    }
			    }
		    }
		    else if(dim == 1)
		    {
			    for(size_t ii = 0; ii < A.size(1); ii++)
			    {
				    for(size_t jj = 0; jj < A.size(2); jj++)
				    {
					    index(ii,jj) = jj;
				    }
			    }
			    for(size_t ii = 0; ii < A.size(1); ii++)
			    {
				    swapped = true;
				    while(swapped)
				    {     
					    swapped = false;
					    for(size_t jj = 1; jj < A.size(2); jj++)
					    {
						    if(((A(ii,jj-1)) > (A(ii,jj)) && method.compare("ascend") == 0) || ((A(ii,jj-1)) < (A(ii,jj)) && method.compare("descend") == 0))
						    {
							    temp = A(ii,jj-1);
							    A(ii,jj-1) = A(ii,jj);
							    A(ii,jj) = temp;
							    temp_i = index(ii,jj-1);
							    index(ii,jj-1) = index(ii,jj);
							    index(ii,jj) = temp_i;
							    swapped = true;
						    }
					    }
				    }
			    }
			    for(size_t ii = 0; ii < A.size(1); ii++)
			    {
				    for(size_t jj = 0; jj < A.size(2); jj++)
				    {
					    B(ii,jj) = A(ii,index(ii,jj));
				    }
			    }
		    }
		    else
		    {
			    throw KeyCppException("Invalid dimension in sort().");
		    }
		
		    Sort_Matrix<T> sort_matrix;
		    sort_matrix.B = B;
		    sort_matrix.index = index;

		    return sort_matrix;
		}
		else
		{
		    bool swapped = true;
		    T temp;
		    size_t temp_i;
		    matrix<size_t,2> index(A.numel());
		    for(size_t ii = 0; ii < A.numel(); ii++)
		    {
			    index(ii) = ii;
		    }
		    matrix<T,2> B(A.numel());
		    for(size_t ii = 0; ii < B.numel(); ii++)
		    {
			    B(ii) = A(ii);
		    }
		    swapped = true;
		    while(swapped)
		    {     
			    swapped = false;
			    for(size_t ii = 1; ii < A.numel(); ii++)
			    {
				    if(((B(ii-1)) > (B(ii)) && method.compare("ascend") == 0) || ((B(ii-1)) < (B(ii)) && method.compare("descend") == 0))
				    {
					    temp = B(ii-1);
					    B(ii-1) = B(ii);
					    B(ii) = temp;
					    temp_i = index(ii-1);
					    index(ii-1) = index(ii);
					    index(ii) = temp_i;
					    swapped = true;
				    }
			    }
		    }
		
		    Sort_Matrix<T> sort_matrix;
		    sort_matrix.B = B;
		    sort_matrix.index = index;

		    return sort_matrix;
		}
	}
	
    /** \brief Displays on standard output any parameter passed to it provided 
     *         the operator << is defined for its type.
     */
	template<class T>
	void disp(const T &x, std::ostream& outStream = std::cout)
	{
	    outStream << x << std::endl;
	    return;
	}
	
	/** \brief Prints the prompt to the screen and then waits for user input.
	 *         Currently the option must be supplied as "s" because C++ is a 
	 *         statically typed language.
	 */
	inline std::string input(const std::string &prompt, std::string option)
	{
	    if(option.empty())
	    {
	        throw KeyCppException("Evaluating input expressions is currently unsupported. Use option \"s\" instead.");
	    }
		std::transform(option.begin(), option.end(), option.begin(), ::tolower);
	    if(option.compare("s") != 0)
	    {
	        throw KeyCppException("Unknown option provided to input()!");
	    }
	    std::cout << prompt;
	    std::string in;
	    std::cin >> in;
	    
	    return in;
	}
	
	/** \brief Converts a std::string to a double. Currently only works on single numbers.
	 *         In the future this should be expanded to work on vectors and matrices. (see MATLAB docs)
	 */
	inline double str2num(const std::string &in)
	{
	    return atof(in.c_str());
	}
	
	/** \brief Returns the length of the largest dimension of A.
	 */
	template<class T>
	size_t length(const matrix<T> &A)
	{
	    return A.length();
	}
	
	template<class T>
	size_t numel(const matrix<T> &A)
	{
	    return A.numel();
	}
	
	/** \brief Finds and returns the indices of non-zero elements of v1.
	 */
	template<class T,size_t dim>
	matrix<size_t,2> find(const matrix<T,dim> &v1, const size_t &k = -1, std::string start = "")
	{
		std::transform(start.begin(), start.end(), start.begin(), ::tolower);
	    matrix<size_t,2> v2;
	    if(v1.empty())
	    {
	        return v2;
	    }
	    if(k < 0 || k > v1.numel() || start.empty() || (start.compare("first") != 0 &&
	       start.compare("last") != 0))
	    {
	        v2.reserve(v1.numel());
	        for(size_t ii = 0; ii < v1.numel(); ii++)
	        {
	            if(std::abs(v1(ii)) > eps)
	            {
	                auto temp = numel(v2);
	                v2.resize(temp + 1);
	                v2(temp) = ii;
	            }
	        }
	    }
	    else if(start.compare("first") == 0 && k > 0 && k < v1.numel())
	    {
	        v2.reserve(k);
	        size_t count = 0;
	        for(size_t ii = 0; ii < v1.numel(); ii++)
	        {
	            if(std::abs(v1(ii)) > eps)
	            {
	                auto temp = numel(v2);
	                v2.resize(temp + 1);
	                v2(temp) = ii;
	                count++;
	                if(count >= k)
	                {
	                    break;
	                }
	            }
	        }
	    }
	    else if(start.compare("last") == 0 && k > 0 && k < v1.numel())
	    {
	        v2.reserve(k);
	        size_t count = 0;
	        for(size_t ii = v1.numel()-1; ii >= 0; ii--)
	        {
	            if(std::abs(v1(ii)) > eps)
	            {
	                auto temp = numel(v2);
	                v2.resize(temp + 1);
	                v2(temp) = ii;
	                count++;
	                if(count >= k)
	                {
	                    break;
	                }
	            }
	        }
	    }
	    else
	    {
	        throw KeyCppException("Unknown arguments in find()!");
	    }
	    
	    return v2;
	}
	
	template<class T>
	matrix<T> reshape(const matrix<T> &A, const size_t &m, const size_t &n)
	{
	    if(A.empty())
	    {
	        throw KeyCppException("Cannot reshape empty matrix!");
	    }
	    if(numel(A) != m*n)
	    {
	        throw KeyCppException("To reshape the number of elements must not change.");
	    }
	    
	    matrix<T> B(m,n);
	    size_t count = 0;
	    size_t m0 = A.size(1);
	    size_t n0 = A.size(2);
	    for(size_t ii = 0; ii < m; ii++)
	    {
	        for(size_t jj = 0; jj < n; jj++)
	        {
	            B(ii,jj) = A(count % m0, count/m0);
	            count++;
	        }
	    }
	    return B;
	}
	
	/** \brief Computes the dot product between the first non-singleton dimension of A and B.
	 */
	template<class T, class U>
	matrix<decltype(std::declval<T>()*std::declval<U>())> dot(const matrix<T> &A, const matrix<U> &B, const size_t &dim = -1)
	{
	    if(A.empty() || B.empty())
	    {
	        throw KeyCppException("Cannot dot multiply empty matrices!");
	    }
	    if(!A.isVec())
	    {
	        if(A.size(1) != B.size(1) && A.size(2) != B.size(2))
	        {
	            throw KeyCppException("Matrices must be same size in dot()!");
	        }
	        matrix<decltype(std::declval<T>()*std::declval<U>())> result;
	        if((A.size(1) > 1 || dim == 1) && dim != 2)
	        {
	            result = matrix<decltype(std::declval<T>()*std::declval<U>())>(A.size(2));
	            for(size_t ii = 0; ii < result.numel(); ii++)
	            {
	                result(ii) = dot(A.col(ii),B.col(ii));
	            }
	        }
	        else
	        {
	            result = matrix<decltype(std::declval<T>()*std::declval<U>())>(A.size(1));
	            for(size_t ii = 0; ii < result.numel(); ii++)
	            {
	                result(ii) = dot(A.row(ii),B.row(ii));
	            }
	        }
	        return result;
	    }
	    else
	    {
	        if(!A.isVec() || !B.isVec())
	        {
	            throw KeyCppException("Inputs must be vectors in dot()!");
	        }
	        if(A.numel() != B.numel())
	        {
	            throw KeyCppException("Vectors must be same size in dot()!");
	        }
	        matrix<decltype(std::declval<T>()*std::declval<U>())> result(1);
	        for(size_t ii = 0; ii < A.numel(); ii++)
	        {
	            result(0) += A(ii)*B(ii);
	        }
	        return result;
	    }
	}
	
	/** \brief Computes the cross product between vectors v1 and v2. Both vectors
	 *         must have exactly 3 elements.
	 */
	template<class T, class U>
	matrix<decltype(std::declval<T>()*std::declval<U>())> cross(const matrix<T> &v1, const matrix<U> &v2)
	{
	    if(v1.empty() || v2.empty())
	    {
	        throw KeyCppException("Cannot cross multiply an empty vector!");
	    }
        if(!v1.isVec() || !v2.isVec())
        {
            throw KeyCppException("Inputs must be vectors in cross()!");
        }
	    if(v1.numel() != 3 || v2.numel() != 3)
	    {
	        throw KeyCppException("Vectors must be have length of 3 in cross()!");
	    }
	    matrix<decltype(std::declval<T>()*std::declval<U>())> result(3);
	    result(0) = v1(1)*v2(2) - v1(2)*v2(1);
	    result(1) = v1(2)*v2(0) - v1(0)*v2(2);
	    result(2) = v1(0)*v2(1) - v1(1)*v2(0);
	    
	    return result;
	}
	
	
	template<class T,class X>
	class SVD_type
	{
	    public:
		    matrix<X> S;
		    matrix<T> U;
		    matrix<T> V;
		SVD_type() : S(), U(), V() {};
	};
	
	double norm(const matrix<double,2,0> &A_in, std::string method = "2");
	SVD_type<double,double> svd(const matrix<double> &A_in, std::string method = "");
	double norm(const matrix<std::complex<double>,2,0> &A_in, std::string method = "2");
	SVD_type<std::complex<double>,double> svd(const matrix<std::complex<double>> &A_in, std::string method = "");
	
	/** \brief Estimates the rank of a matrix by counting the singular values
	 *         whose absolute value is greater than epsilon.
	 */
	template<class T>
	int rank(const matrix<T> &A)
	{
	    auto output = svd(A);
	    return length(find(diag(output.S)));
	}
	
	/** \brief Computes the nullspace of matrix A.
	 */
	template<class T>
	matrix<T> null(const matrix<T> &A)
	{
	    auto output = svd(A);
	    vector_k<int> index;
	    auto S = diag(output.S);
	    index.reserve(S.size());
	    size_t max_dim = (A.size(1) > A.size(2))?A.size(1):A.size(2);
	    for(int ii = 0; ii < S.size(); ii++)
	    {
	        if(std::abs(S[ii]) < eps*max_dim)
	        {
	            index.push_back(ii);
	        }
	    }
	    matrix<T> B(output.V.size(1),index.size());
	    for(size_t ii = 0; ii < B.size(1); ii++)
	    {
	        for(size_t jj = 0; jj < B.size(2); jj++)
	        {
	            B(ii,jj) = output.V(ii,index[jj]);
	        }
	    }
	    
	    return B;
	}
	
	/** \brief Returns true if a is nonzero.
	 */
	template<class T>
	bool any(const T &a)
	{
        if(std::abs(a) < eps)
        {
            return false;
        }
	    return true;
	}
	
	/** \brief Returns true if any elements of A are nonzero.
	 */
	template<class T, size_t dim>
	bool any(const matrix<T,dim> &A)
	{
	    for(size_t ii = 0; ii < A.numel(); ii++)
	    {
            if(std::abs(A(ii)) > eps)
            {
                return true;
            }
	    }
	    return false;
	}
	
	/** \brief Returns true if a is nonzero.
	 */
	template<class T>
	bool all(const T &a)
	{
        if(std::abs(a) < eps)
        {
            return false;
        }
	    return true;
	}
	
	/** \brief Returns true if all elements of A are nonzero.
	 */
	template<class T,size_t dim>
	bool all(const matrix<T,dim> &A)
	{
	    for(size_t ii = 0; ii < A.numel(); ii++)
	    {
            if(!all(A(ii)))
            {
                return false;
            }
	    }
	    return true;
	}
	
	/** \brief Returns boolean value that is true if
	 *         a is finite.
	 */
	template<class T>
	bool finite(const T &a)
	{
	    bool out;
        if(std::isfinite(a))
        {
            out = true;
        }
        else
        {
            out = false;
        }
	    return out;
	}
	
	/** \brief Returns matrix containing boolean values that are true if
	 *         corresponding elements of A are finite.
	 */
	template<class T,size_t dim>
	matrix<bool,dim> finite(const matrix<T,dim> &A)
	{
	    matrix<bool,dim> out;
	    out.resize(size(A));
	    for(size_t ii = 0; ii < A.numel(); ii++)
	    {
	        out(ii) = finite(A(ii));
	    }
	    return out;
	}
	
	/** \brief Returns boolean value that is true if
	 *         a is infinite.
	 */
	template<class T>
	bool isinf(const T &a)
	{
	    return std::isinf(a);
	}
	
	/** \brief Returns matrix containing boolean values that are true if
	 *         corresponding elements of A are infinite.
	 */
	template<class T,size_t dim>
	matrix<bool,dim> isinf(const matrix<T,dim> &A)
	{
	    matrix<bool,dim> out;
	    out.resize(size(A));
	    for(size_t ii = 0; ii < A.numel(); ii++)
	    {
	        out(ii) = isinf(A(ii));
	    }
	    return out;
	}
	
	
	/** \brief Returns boolean value that is true if
	 *         a is NaN.
	 */
	template<>
	inline bool isnan<>(const std::complex<double> &a)
	{
	    return (std::isnan(real(a)) || std::isnan(imag(a)));
	}
	
	/** \brief Returns boolean value that is true if
	 *         a is NaN.
	 */
	template<class T>
	bool isnan(const T &a)
	{
	    return std::isnan(a);
	}
	
	/** \brief Returns matrix containing boolean values that are true if
	 *         corresponding elements of A are NaN.
	 */
	template<class T,size_t dim>
	matrix<bool,dim> isnan(const matrix<T,dim> &A)
	{
	    matrix<bool,dim> out;
	    out.resize(size(A));
	    for(size_t ii = 0; ii < A.numel(); ii++)
	    {
	        out(ii) = keycpp::isnan(A(ii));
	    }
	    return out;
	}
	
	/** \brief Returns true if matrix is empty.
	 */
	template<class T>
	matrix<bool> isempty(const matrix<T> &A)
	{
	    return A.empty();
	}
	
	/** \brief Returns true if a is real.
	 */
	template<class T>
	bool isreal(const T &a)
	{
        if(std::abs(imag(a)) < eps)
        {
            return true;
        }
	    return false;
	}
	
	/** \brief Returns true if all elements of A are real.
	 */
	template<class T,size_t dim>
	bool isreal(const matrix<T,dim> &A)
	{
	    for(int ii = 0; ii < A.numel(); ii++)
	    {
            if(!isreal(A(ii)))
            {
                return false;
            }
	    }
	    return true;
	}
	
	/** \brief Rounds the real and imaginary parts of std::complex<double> a towards
	 *         positive infinity seperately.
	 */
	inline std::complex<double> ceil(const std::complex<double> &a)
	{
	    std::complex<double> b;
	    b = ceil(real(a)) + std::complex<double>(0.0,1.0)*ceil(imag(a));
	    return b;
	}
	
	/** \brief Rounds the elements of A towards positive infinity.
	 */
	template<class T, size_t dim>
	matrix<T,dim> ceil(const matrix<T,dim> &A)
	{
	    return eop(A,static_cast<T (*)(T)>(&std::ceil));
	}
	
	/** \brief Rounds the elements of A towards the closest integer.
	 */
	template<class T, size_t dim>
	matrix<T,dim> round(const matrix<T,dim> &A)
	{
	    return eop(A,static_cast<T (*)(T)>(&std::round));
	}
	
	/** \brief Rounds the real and imaginary parts of std::complex<double> a towards
	 *         negative infinity seperately.
	 */
	inline std::complex<double> floor(const std::complex<double> &a)
	{
	    std::complex<double> b;
	    b = floor(real(a)) + std::complex<double>(0.0,1.0)*floor(imag(a));
	    return b;
	}
	
	/** \brief Rounds the elements of A towards negative infinity.
	 */
	template<class T, size_t dim>
	matrix<T,dim> floor(const matrix<T,dim> &A)
	{
	    return eop(A,static_cast<T (*)(T)>(&std::floor));
	}
	
	template<class T, class U>
	decltype(std::declval<T>()*std::declval<U>()) polyval(const matrix<T,2> &p, const U &x)
	{
	    decltype(std::declval<T>()*std::declval<U>()) val = 0.0;
	    
        for(size_t ii = 0; ii < p.numel()-1; ii++)
        {
            val += p(ii)*pow(x,p.numel()-ii-1);
        }
	    val += p(p.numel()-1);
	    return val;
	}
	
	/** \brief Computes all roots of polynomial p by solving for the eigenvalues
	 *         of the companion matrix.
	 */
	template<class T>
	matrix<std::complex<double>,2> roots(const matrix<T,2> &p)
	{
	    size_t n = p.numel()-1;
	    matrix<T> A = diag(diag(ones<T>(n-1)),-1);
	    for(size_t ii = 0; ii < A.size(2); ii++)
	    {
            A(0,ii) = -p(ii+1)/p(0);
        }
        matrix<std::complex<double>,2> v = transpose(diag(eig(A)));
        return v;
	}
	
	/** \brief Data type for using the tic() and toc(tictoc_type Timer) commands.
	 */
	struct tictoc_type
	{
	    timeval start, stop, elapsed;
	};
	
	/** \brief Start the timer.
	 */
	inline tictoc_type tic()
	{
	    tictoc_type Timer;
	    gettimeofday(&Timer.start,NULL);
	    return Timer;
	}
	
	/** \brief Stop the timer. The number of elapsed seconds is returned.
	 */
	inline double toc(tictoc_type &Timer)
	{
	    gettimeofday(&Timer.stop,NULL);
	    timersub(&Timer.stop,&Timer.start,&Timer.elapsed);
	    return ((double)Timer.elapsed.tv_sec + 1e-6*((double)Timer.elapsed.tv_usec));
	}
	
	/** \brief Overload of the C++ function sprintf(). This overload provides a more MATLAB-like
	 *         interface. Specifically, the output is returned instead of passed by reference.
	 */
	inline std::string sprintf(const std::string &fmt, ...)
	{
        int size = 100;
        va_list ap;
        char* s;
        while(1)
        {
            s = new char[size + 1];
            va_start(ap, fmt);
            int n = vsnprintf(s, size, fmt.c_str(), ap);
            va_end(ap);
            if(n > -1 && n < size)
            {
                std::string str(s);
                delete [] s;
                str.resize(n);
                return str;
            }
            if (n > -1)
            {
                size = n + 1;
            }
            else
            {
                size *= 2;
            }
            delete [] s;
        }
        std::string str(s);
        return str;
    }
    
    /** \brief Returns a vector of integers containing the current: year, month, day,
     *         hour, minute, and second. This is based on the system clock. The number of hours
     *         is based on the 24-hour clock.
     */
    inline matrix<size_t,2> clock()
    {
        time_t t = time(0);
        struct tm * now = localtime(&t);
        matrix<size_t,2> dt(6);
        dt(0) = (now->tm_year + 1900); // year
        dt(1) = (now->tm_mon + 1); // month
        dt(2) = (now->tm_mday); // day
        dt(3) = (now->tm_hour); // hour
        dt(4) = (now->tm_min); // minute
        dt(5) = (now->tm_sec); // seconds
        
        return dt;
    }
    
    /** \brief Returns the Moore-Penrose Pseudoinverse of matrix A. Currently only the SVD method is implemented.
     *         This restricts matrix A to be only square matrices. This is currently slower than inv(), use with care.
     */
    template<class T>
    matrix<T> pinv(const matrix<T> &A)
    {
        if(A.empty())
        {
            throw KeyCppException("Cannot compute pseudoinverse of empty matrix!");
        }
        if(A.size(1) != A.size(2))
        {
            throw KeyCppException("A must be a square matrix when computing pseudoinverse using SVD.");
        }
        auto svd_out = svd(A);
        
        matrix<T> s_inv = diag(svd_out.S);
        for(size_t ii = 0; ii < s_inv.numel(); ii++)
        {
            if(any(s_inv(ii)))
            {
                s_inv(ii) = 1.0/s_inv(ii);
            }
        }
        matrix<T> Ap = svd_out.V*diag(s_inv)*ctranspose(svd_out.U);
        return Ap;
    }
    
    inline std::string removeWhiteSpace(std::string in)
    {
        bool stop = false;
        while(in.length() > 0 && !stop)
        {
            if(isspace(in.at(0)))
            {
                in = in.substr(1, in.length() - 1);
            }
            else if(isspace(in.at(in.length() - 1)))
            {
                in = in.substr(0, in.length() - 1);
            }
            else
            {
                stop = true;
            }
        }
        return in;
    }

    inline std::string removeQuotes(std::string in)
    {
	    bool stop = false;
	    while(in.length() > 0 && !stop)
	    {
		    if(in.at(0) == '\"')
		    {
			    in = in.substr(1,in.length()-1);
		    }
		    else if(in.at(in.length()-1) == '\"')
		    {
			    in = in.substr(0,in.length()-1);
		    }
		    else
		    {
			    stop = true;
		    }
	    }

	    return in;
    }
    
    /** \brief Returns a matrix containing the data read from a text file. Values must be white space separated. */
    inline matrix<double> importdata(std::string filename)
    {
        std::ifstream in(filename.c_str());
        
        if(!in.is_open())
        {
            throw KeyCppException("ERROR!! Could not open file in importdata!");
        }
        
        matrix<double> A;
        A.reserve(1000);
        
        std::stringstream ss;
	    std::string dummy = "";
	    double temp;
	
        while(std::getline(in, dummy))
        {
            vector_k<double> b;
            b.reserve(100);
	        ss.str("");
	        ss.clear();
	        ss << dummy;
	        std::getline(ss,dummy,'/');
	        ss.str("");
	        ss.clear();
	        dummy = removeWhiteSpace(dummy);
            if(!dummy.empty())
            {
	            ss << dummy;
	            while(ss >> temp)
	            {
                    b.push_back(temp);
                }
                if(A.empty())
                {
                    A = matrix<double>(1,b.size());
                    A.row(0) = b;
                }
                else
                {
                    auto temp_size = size(A);
                    temp_size(0)++;
                    A.resize(temp_size);
                    A.row(temp_size(0)-1) = b;
                }
            }
        }
        in.close();
        return A;
    }
    
    /** \brief Returns a matrix containing the data read from a text file. Values must be white space separated. */
    inline matrix<std::complex<double>> importdata_complex(std::string filename)
    {
        std::ifstream in(filename.c_str());
        
        if(!in.is_open())
        {
            throw KeyCppException("ERROR!! Could not open file in importdata!");
        }
        
        matrix<std::complex<double>> A;
        A.reserve(1000);
        
        std::stringstream ss;
	    std::string dummy = "";
	    std::complex<double> temp;
	
        while(std::getline(in, dummy))
        {
            vector_k<std::complex<double>> b;
            b.reserve(100);
	        ss.str("");
	        ss.clear();
	        ss << dummy;
	        std::getline(ss,dummy,'/');
	        ss.str("");
	        ss.clear();
	        dummy = removeWhiteSpace(dummy);
            if(!dummy.empty())
            {
	            ss << dummy;
	            while(ss >> temp)
	            {
                    b.push_back(temp);
                }
                if(A.empty())
                {
                    A = matrix<std::complex<double>>(1,b.size());
                    A.row(0) = b;
                }
                else
                {
                    auto temp_size = size(A);
                    temp_size(0)++;
                    A.resize(temp_size);
                    A.row(temp_size(0)-1) = b;
                }
            }
        }
        in.close();
        return A;
    }
    
    /** \brief Returns the standard deviation of inputed vector. */
    template<class T>
    matrix<T,2> stdev(matrix<T,2> v1)
    {
        if(v1.isVec())
        {
            T v_bar = mean(v1);
            matrix<T,2> temp(1);
            for(size_t ii = 0; ii < v1.numel(); ii++)
            {
                temp(0) += std::abs(v1(ii) - v_bar)*std::abs(v1(ii) - v_bar);
            }
            temp(0) = std::sqrt(temp(0)/((double)v1.numel()-1.0));
            
            return temp;
        }
        else
        {
            matrix<T,2> temp(v1.size(2));
            for(size_t ii = 0; ii < v1.size(2); ii++)
            {
                temp(ii) = stdev(v1.col(ii));
            }
            return temp;
        }
    }
    
    /** \brief Returns the standard deviation of inputed vector. */
    inline matrix<double,2> stdev(matrix<std::complex<double>,2> v1)
    {
        if(v1.isVec())
        {
            std::complex<double> v_bar = mean(v1);
            matrix<double,2> temp(1);
            for(size_t ii = 0; ii < v1.numel(); ii++)
            {
                temp(0) += std::abs((v1(ii) - v_bar)*conj(v1(ii) - v_bar));
            }
            temp(0) = std::sqrt(temp(0)/((double)v1.numel()-1.0));
            
            return temp;
        }
        else
        {
            matrix<double,2> temp(v1.size(2));
            for(size_t ii = 0; ii < v1.size(2); ii++)
            {
                temp(ii) = stdev(v1.col(ii));
            }
            return temp;
        }
    }
    
    /** \brief Returns the variance (square of standard deviation) for inputed vector. */
    template<class T>
    T var(matrix<T,2> v1)
    {
        return pow(stdev(v1),2);
    }
    
    /** \brief Returns the variance (square of standard deviation) for inputed vector. */
    inline double var(matrix<std::complex<double>,2> v1)
    {
        return pow(stdev(v1),2);
    }
    
    template<class T> matrix<double> rms(const matrix<T,2> &A, size_t dimension)
    {
        if(dimension > 2 || dimension < 1)
        {
            throw KeyCppException("Invalid dimension in rms()!");
        }
        size_t N;
        if(dimension == 2)
        {
            N = A.size(1);
            matrix<double> B(N);
            for(size_t jj = 0; jj < N; jj++)
            {
                for(size_t ii = 0; ii < A.size(dimension); ii++)
                {
                    B(jj) += std::abs(A(jj,ii))*std::abs(A(jj,ii));
                }
                B(jj) = std::sqrt(B(jj)/A.size(dimension));
            }
            return B;
        }
        else
        {
            N = A.size(2);
            matrix<double> B(N);
            for(size_t jj = 0; jj < N; jj++)
            {
                for(size_t ii = 0; ii < A.size(dimension); ii++)
                {
                    B(jj) += std::abs(A(ii,jj))*std::abs(A(ii,jj));
                }
                B(jj) = std::sqrt(B(jj)/A.size(dimension));
            }
            return B;
        }
    }
    
    template<class T> matrix<double> rms(const matrix<T,2> &A)
    {
        size_t dimension = 0;
        for(size_t ii = 0; ii < 2; ii++)
        {
            dimension = ii+1;
            if(A.size(ii+1) > 1)
            {
                break;
            }
        }
        return rms(A,dimension);
    }
    
    namespace rng_ns
    {
        inline int& getChoice()
        {
            static int choice = 1;
            return choice;
        }
        
        inline std::mt19937& get_mt_rng()
        {
            // Choice = 1;
            static std::mt19937 mt_rng(std::random_device{}());
            return mt_rng;
        }
        
        inline std::ranlux24_base& get_lagfib_rng()
        {
            // Choice = 2;
            static std::ranlux24_base lagfib_rng(std::random_device{}());
            return lagfib_rng;
        }
    }
    
    inline void rng(size_t seed = 0, std::string generator = "twister")
    {
		std::transform(generator.begin(), generator.end(), generator.begin(), ::tolower);
		if(generator.compare("twister") == 0)
		{
		    rng_ns::getChoice() = 1;
		    rng_ns::get_mt_rng() = std::mt19937(seed);
		}
		else if(generator.compare("multFibonacci") == 0)
		{
		    rng_ns::getChoice() = 2;
		    rng_ns::get_lagfib_rng() = std::ranlux24_base(seed);
		}
    }
    
    inline void rng(std::string shuffle = "", std::string generator = "")
    {
		std::transform(generator.begin(), generator.end(), generator.begin(), ::tolower);
		std::transform(shuffle.begin(), shuffle.end(), shuffle.begin(), ::tolower);
		if(generator.empty())
		{
		    if(shuffle.compare("shuffle") == 0)
		    {
		        if(rng_ns::getChoice() == 2)
		        {
		            rng(std::random_device{}(),"multFibonacci");
		        }
		        else
                {
		            rng(std::random_device{}(),"twister");
                }
		    }
		    else if(shuffle.compare("default") == 0)
		    {
		        rng_ns::getChoice() = 1;
		        rng_ns::get_mt_rng() = std::mt19937(0);
		    }
		}
		else
		{
		    if(generator.compare("twister") == 0)
		    {
		        rng_ns::getChoice() = 1;
		        rng_ns::get_mt_rng() = std::mt19937(std::random_device{}());
		    }
		    else if(generator.compare("multFibonacci") == 0)
		    {
		        rng_ns::getChoice() = 2;
		        rng_ns::get_lagfib_rng() = std::ranlux24_base(std::random_device{}());
		    }
		}
    }
    
    /** \brief Returns a random double between 0 and 1.0
	 */
	inline double rand()
	{
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        double randomNumber;
        if(rng_ns::getChoice() == 1)
        {
            randomNumber = distribution(rng_ns::get_mt_rng());
        }
        else if(rng_ns::getChoice() == 2)
        {
            randomNumber = distribution(rng_ns::get_lagfib_rng());
        }
        else
        {
            throw KeyCppException("Unknown random number generator!");
        }
		return randomNumber;
	}
	
	/** \brief Returns an N x N matrix of random doubles between 0 and 1.0
	 */
	inline matrix<double> rand(const unsigned int &N)
	{
	    matrix<double> A(N,N);
	    for(unsigned int ii = 0; ii < N; ii++)
	    {
	        for(unsigned int jj = 0; jj < N; jj++)
	        {
	            A(ii,jj) = keycpp::rand();
	        }
        }
	    
		return A;
	}
	
	/** \brief Returns an M x N matrix of random doubles between 0 and 1.0
	 */
	inline matrix<double> rand(const unsigned int &M, const unsigned int &N)
	{
	    matrix<double> A(M,N);
	    for(unsigned int ii = 0; ii < M; ii++)
	    {
	        for(unsigned int jj = 0; jj < N; jj++)
	        {
	            A(ii,jj) = keycpp::rand();
	        }
        }
	    
		return A;
	}
	
	inline double randn()
	{
        std::normal_distribution<double> distribution(0.0, 1.0);
        double randomNumber;
        if(rng_ns::getChoice() == 1)
        {
            randomNumber = distribution(rng_ns::get_mt_rng());
        }
        else if(rng_ns::getChoice() == 2)
        {
            randomNumber = distribution(rng_ns::get_lagfib_rng());
        }
        else
        {
            throw KeyCppException("Unknown random number generator!");
        }
		return randomNumber;
	}

	inline matrix<double> randn(const unsigned int &N)
	{
	    matrix<double> A(N,N);
	    for(unsigned int ii = 0; ii < N; ii++)
	    {
	        for(unsigned int jj = 0; jj < N; jj++)
	        {
	            A(ii,jj) = keycpp::randn();
	        }
        }
	    
		return A;
	}

	inline matrix<double> randn(const unsigned int &M, const unsigned int &N)
	{
	    matrix<double> A(M,N);
	    for(unsigned int ii = 0; ii < M; ii++)
	    {
	        for(unsigned int jj = 0; jj < N; jj++)
	        {
	            A(ii,jj) = keycpp::randn();
	        }
        }
	    
		return A;
	}

	/** \brief Generalized complex-valued eigenvalue solver using LAPACK function call. 
	 *  
	 *  This function returns the eigenvalues(lambda) of the complex-valued generalized
	 *  eigenvalue problem: Ax_r = lambda*Bx_r or x_l^T*A = lambda*x_l^T*B. The eigenvalues
	 *  are returned by default. To return the right or left eigenvectors, supply the
	 *  function with a std::complex<double> matrix object in the 3rd or 4th parameters, respectively.
	 */
    inline matrix<std::complex<double> > eig(const matrix<std::complex<double> > &A, const matrix<std::complex<double> > &B, matrix<std::complex<double> > *vr_return, matrix<std::complex<double> > *vl_return)
	{
		unsigned int n;
		int nn, lda, ldb, ldvl, ldvr, lwork, info;
		n = (unsigned)A.size(1);
		lda = ldb = (int)A.size(1);
		nn = n;
		lwork = n*n + 64;
		char jobvl, jobvr;

		if(vl_return == NULL)
		{
			jobvl = 'N';
			ldvl = 1;
		}
		else
		{
			jobvl = 'V';
			ldvl = n;
		}

		if(vr_return == NULL)
		{
			jobvr = 'N';
			ldvr = 1;
		}
		else
		{
			jobvr = 'V';
			ldvr = n;
		}

		std::complex<double> *a = new std::complex<double>[n*n];
		std::complex<double> *b = new std::complex<double>[n*n];
		std::complex<double> *vr = new std::complex<double>[n*n];
		std::complex<double> *vl = new std::complex<double>[n*n];
		std::complex<double> *alpha = new std::complex<double>[n];
		std::complex<double> *beta = new std::complex<double>[n];
		std::complex<double> *work = new std::complex<double>[lwork];
		double *rwork = new double[8*n];
		for(unsigned int ii = 0; ii < n; ii++)
		{
			for(unsigned int jj = 0; jj < n; jj++)
			{
				a[ii*n + jj] = A(jj,ii);
				b[ii*n + jj] = B(jj,ii);
			}
		}

		zggev_(&jobvl, &jobvr, &nn, a, &lda, b, &ldb, alpha, beta, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);

		matrix<std::complex<double> > lambda(n,n);
		for(unsigned int ii = 0; ii < n; ii++)
		{
			lambda(ii,ii) = alpha[ii]/beta[ii];
		}
		if(jobvr == 'V')
		{
			(*vr_return) = matrix<std::complex<double> >(n,n);
			for(unsigned int ii = 0; ii < n; ii++)
			{
				for(unsigned int jj = 0; jj < n; jj++)
				{
					(*vr_return)(jj,ii) = vr[ii*n + jj];
				}
			}
		}
		if(jobvl == 'V')
		{
			(*vl_return) = matrix<std::complex<double> >(n,n);
			for(unsigned int ii = 0; ii < n; ii++)
			{
				for(unsigned int jj = 0; jj < n; jj++)
				{
					(*vl_return)(jj,ii) = vl[ii*n + jj];
				}
			}
		}

		delete [] a;
		delete [] b;
		delete [] vr;
		delete [] vl;
		delete [] alpha;
		delete [] beta;
		delete [] work;
		delete [] rwork;

		return lambda;
	}
	
	/** \brief Complex-valued eigenvalue solver using LAPACK function call. 
	 *  
	 *  This function returns the eigenvalues(lambda) of the complex-valued
	 *  eigenvalue problem: Ax_r = lambda*x_r or x_l^T*A = lambda*x_l^T. The eigenvalues
	 *  are returned by default. To return the right or left eigenvectors, supply the
	 *  function with a std::complex<double> matrix object in the 2nd or 3rd parameters, respectively.
	 */
    inline matrix<std::complex<double>> eig(const matrix<std::complex<double> > &A, matrix<std::complex<double> > *vr_return, matrix<std::complex<double> > *vl_return)
	{
		unsigned int n;
		int nn, lda, ldb, ldvl, ldvr, lwork, info;
		n = (unsigned)A.size(1);
		lda = ldb = (int)A.size(1);
		nn = n;
		lwork = 2*n;
		char jobvl, jobvr;

		if(vl_return == NULL)
		{
			jobvl = 'N';
			ldvl = 1;
		}
		else
		{
			jobvl = 'V';
			ldvl = n;
		}

		if(vr_return == NULL)
		{
			jobvr = 'N';
			ldvr = 1;
		}
		else
		{
			jobvr = 'V';
			ldvr = n;
		}

		std::complex<double> *a = new std::complex<double>[n*n];
		std::complex<double> *vr = new std::complex<double>[n*n];
		std::complex<double> *vl = new std::complex<double>[n*n];
		std::complex<double> *w = new std::complex<double>[n];
		std::complex<double> *work = new std::complex<double>[lwork];
		double *rwork = new double[lwork];
		for(unsigned int ii = 0; ii < n; ii++)
		{
			for(unsigned int jj = 0; jj < n; jj++)
			{
				a[ii*n + jj] = A(jj,ii);
			}
		}

	    zgeev_(&jobvl, &jobvr, &nn, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);

		matrix<std::complex<double> > lambda(n,n);
		for(unsigned int ii = 0; ii < n; ii++)
		{
			lambda(ii,ii) = w[ii];
		}
		if(jobvr == 'V')
		{
			(*vr_return) = matrix<std::complex<double> >(n,n);
			for(unsigned int ii = 0; ii < n; ii++)
			{
				for(unsigned int jj = 0; jj < n; jj++)
				{
					(*vr_return)(jj,ii) = vr[ii*n + jj];
				}
			}
		}
		if(jobvl == 'V')
		{
			(*vl_return) = matrix<std::complex<double> >(n,n);
			for(unsigned int ii = 0; ii < n; ii++)
			{
				for(unsigned int jj = 0; jj < n; jj++)
				{
					(*vl_return)(jj,ii) = vl[ii*n + jj];
				}
			}
		}

		delete [] a;
		delete [] vr;
		delete [] vl;
		delete [] w;
		delete [] work;
		delete [] rwork;

		return lambda;
	}
	
	
	/** \brief Double precision eigenvalue solver using LAPACK function call. 
	 *  
	 *  This function returns the eigenvalues(lambda) of the
	 *  eigenvalue problem: Ax_r = lambda*x_r or x_l^T*A = lambda*x_l^T. The eigenvalues
	 *  are returned by default. To return the right or left eigenvectors, supply the
	 *  function with a std::complex<double> matrix object in the 2nd or 3rd parameters, respectively.
	 */
    inline matrix<std::complex<double> > eig(const matrix<double> &A, matrix<std::complex<double> > *vr_return, matrix<std::complex<double> > *vl_return)
	{
	    matrix<std::complex<double>> B(A.size(1),A.size(2));
	    for(unsigned int ii = 0; ii < B.size(1); ii++)
	    {
	        for(unsigned int jj = 0; jj < B.size(2); jj++)
	        {
	            B(ii,jj) = A(ii,jj);
	        }
	    }
	    return eig(B, vr_return, vl_return);
	}
	
	inline double rcond(const matrix<double> &A)
	{
	    if(A.size(1) != A.size(2))
	    {
	        throw KeyCppException("In rcond(), input must be a square matrix.");
	    }
	    
	    int info, n, lda;
        double anorm, rcond;
        
        int *iw = new int[A.size(1)];
        double *w = new double[A.size(1)*A.size(2) + 64];
        double *x = new double[A.size(1)*A.size(2)];
        for(unsigned int ii = 0; ii < A.size(2); ii++)
        {
            for(unsigned int jj = 0; jj < A.size(1); jj++)
            {
                x[ii*A.size(1) + jj] = A(jj,ii);
            }
        }
        n = (int)A.size(1);
        lda = n;

        /* Computes the norm of x */
        anorm = dlange_("1", &n, &n, x, &lda, w);

        /* Modifies x in place with a LU decomposition */
        dgetrf_(&n, &n, x, &lda, iw, &info);
        if(info != 0)
        {
            throw KeyCppException("Unknown error in rcond()!");
        }

        /* Computes the reciprocal norm */
        dgecon_("1", &n, x, &lda, &anorm, &rcond, w, iw, &info);
        if(info != 0)
        {
            throw KeyCppException("Unknown error in rcond()!");
        }
        
        delete [] iw;
        delete [] w;
        delete [] x;

        return rcond;
	}
	
	inline double rcond(const matrix<std::complex<double>> &A)
	{
	    if(A.size(1) != A.size(2))
	    {
	        throw KeyCppException("In rcond(), input must be a square matrix.");
	    }
	    
	    int info, n, lda;
        double anorm, rcond;
        
        int *iw = new int[A.size(1)];
        double *w1 = new double[A.size(1)*A.size(2) + 64];
        std::complex<double> *w2 = new std::complex<double>[A.size(1)*A.size(2) + 64];
        std::complex<double> *x = new std::complex<double>[A.size(1)*A.size(2)];
        for(unsigned int ii = 0; ii < A.size(2); ii++)
        {
            for(unsigned int jj = 0; jj < A.size(1); jj++)
            {
                x[ii*A.size(1) + jj] = A(jj,ii);
            }
        }
        n = (int)A.size(1);
        lda = n;

        /* Computes the norm of x */
        anorm = zlange_("1", &n, &n, x, &lda, w1);

        /* Modifies x in place with a LU decomposition */
        zgetrf_(&n, &n, x, &lda, iw, &info);
        if(info != 0)
        {
            throw KeyCppException("Unknown error in rcond()!");
        }

        /* Computes the reciprocal norm */
        zgecon_("1", &n, x, &lda, &anorm, &rcond, w2, w1, &info);
        if(info != 0)
        {
            throw KeyCppException("Unknown error in rcond()!");
        }
        
        delete [] iw;
        delete [] w1;
        delete [] w2;
        delete [] x;

        return rcond;
	}
	
	inline matrix<std::complex<double>> linsolve(const matrix<std::complex<double>>& A_in,
	                                      const matrix<std::complex<double>>& b_in)
	{
		if(b_in.empty() || A_in.empty())
		{
			throw KeyCppException("Error in linsolve()! Empty matrix or vector supplied!");
		}
		if(A_in.size(2) != b_in.size(1) || b_in.size(2) != 1)
		{
			throw KeyCppException("Error in linsolve()! Matrix and vector sizes are not compatible!");
		}
		
		unsigned int n = (unsigned)b_in.size(1);
		int nn = n;
		int m = (int)A_in.size(2);
		int nrhs = 1;

        int info, lda;
        double anorm, rcond;
        
        int *iw = new int[A_in.size(1)];
        std::complex<double> *w1 = new std::complex<double>[A_in.size(1)*A_in.size(2) + 64];
        double *w2 = new double[A_in.size(1)*A_in.size(2) + 64];
        std::complex<double> *A = new std::complex<double>[A_in.size(1)*A_in.size(2)];
        for(unsigned int ii = 0; ii < A_in.size(2); ii++)
        {
            for(unsigned int jj = 0; jj < A_in.size(1); jj++)
            {
                A[ii*A_in.size(1) + jj] = A_in(jj,ii);
            }
        }
        lda = n;

        /* Computes the norm of A */
        anorm = zlange_("1", &nn, &nn, A, &lda, w2);

        /* Modifies A in place with a LU decomposition */
        zgetrf_(&nn, &nn, A, &lda, iw, &info);
        if(info != 0)
        {
            if(info > 0)
            {
                std::cerr << "Warning: Matrix is singular. Results may be inaccurate.\n";
            }
            else
            {
                throw KeyCppException("Unknown error in linsolve()!");
            }
        }

        /* Computes the reciprocal norm */
        zgecon_("1", &nn, A, &lda, &anorm, &rcond, w1, w2, &info);
        if(info != 0)
        {
            throw KeyCppException("Unknown error in linsolve()!");
        }
        
        if(rcond < 1e-15)
        {
            std::cerr << "Warning: Matrix is close to singular or badly scaled. Results may be inaccurate.\nrcond = " << rcond << std::endl;
        }
        
        matrix<std::complex<double>> x_out(b_in.size(1),1);
        for(size_t ii = 0; ii < b_in.size(1); ii++)
        {
            x_out(ii) = b_in(ii,0);
        }
        zgetrs_("N", &m, &nrhs, A, &lda, iw, &x_out(0), &nn, &info);
        if(info != 0)
        {
            throw KeyCppException("Unknown error in linsolve()!");
        }
        
        delete [] iw;
        delete [] w1;
        delete [] w2;
        delete [] A;
        
		return x_out;
	}
	
	inline matrix<double> linsolve(const matrix<double>& A_in,
	                                      const matrix<double>& b_in)
	{
		if(b_in.empty() || A_in.size(1) <= 0 || A_in.size(2) <= 0)
		{
			throw KeyCppException("Error in linsolve()! Empty matrix or vector supplied!");
		}
		if(A_in.size(2) != b_in.size(1) || b_in.size(2) != 1)
		{
			throw KeyCppException("Error in linsolve()! Matrix and vector sizes are not compatible!");
		}
		
		unsigned int n = (unsigned)b_in.size(1);
		unsigned int m = (unsigned)A_in.size(2);
		int nrhs = 1, nn = n, mm = m;

        int info = 0, lda;
        double anorm, rcond;
        
        int *iw = new int[A_in.size(1)];
        double *w1 = new double[A_in.size(1)*A_in.size(2) + 64];
        int *w2 = new int[A_in.size(1)*A_in.size(2) + 64];
        double *A = new double[A_in.size(1)*A_in.size(2)];
        for(unsigned int ii = 0; ii < A_in.size(2); ii++)
        {
            for(unsigned int jj = 0; jj < A_in.size(1); jj++)
            {
                A[ii*A_in.size(1) + jj] = A_in(jj,ii);
            }
        }
        lda = n;

        /* Computes the norm of A */
        anorm = dlange_("1", &nn, &nn, A, &lda, w1);

        /* Modifies A in place with a LU decomposition */
        dgetrf_(&nn, &nn, A, &lda, iw, &info);
        if(info != 0)
        {
            if(info > 0)
            {
                std::cerr << "Warning: Matrix is singular. Results may be inaccurate.\n";
            }
            else
            {
                throw KeyCppException("Unknown error in linsolve()!");
            }
        }

        /* Computes the reciprocal norm */
        dgecon_("1", &nn, A, &lda, &anorm, &rcond, w1, w2, &info);
        if(info != 0)
        {
            throw KeyCppException("Unknown error in linsolve()!");
        }
        
        if(rcond < 1e-15)
        {
            std::cerr << "Warning: Matrix is close to singular or badly scaled. Results may be inaccurate.\nrcond = " << rcond << std::endl;
        }
        
        matrix<double> x_out(b_in.size(1),1);
        for(size_t ii = 0; ii < b_in.size(1); ii++)
        {
            x_out(ii) = b_in(ii,0);
        }
        dgetrs_("N", &mm, &nrhs, A, &lda, iw, &x_out(0), &nn, &info);
        if(info != 0)
        {
            throw KeyCppException("Unknown error in linsolve()!");
        }
        
        delete [] iw;
        delete [] w1;
        delete [] w2;
        delete [] A;
        
		return x_out;
	}
	
	inline matrix<double> inv(const matrix<double>& A_in)
	{
	    if(A_in.size(1) <= 0 || A_in.size(2) <= 0)
		{
			throw KeyCppException("Error in inv()! Empty matrix supplied!");
		}
		if(A_in.size(1) != A_in.size(2))
		{
		    throw KeyCppException("Error in inv()! Matrix must be square!");
		}
		
		unsigned int n = (unsigned)A_in.size(1);
		int nn = (int)n;

        int info, lda;
        double anorm, rcond;
        
        int *iw = new int[A_in.size(1)];
        int lwork = (int)(A_in.size(1)*A_in.size(2)) + 64;
        double *w1 = new double[lwork];
        int *w2 = new int[lwork];
        double *A = new double[A_in.size(1)*A_in.size(2)];
        for(unsigned int ii = 0; ii < A_in.size(2); ii++)
        {
            for(unsigned int jj = 0; jj < A_in.size(1); jj++)
            {
                A[ii*A_in.size(1) + jj] = A_in(jj,ii);
            }
        }
        lda = n;

        /* Computes the norm of A */
        anorm = dlange_("1", &nn, &nn, A, &lda, w1);

        /* Modifies A in place with a LU decomposition */
        dgetrf_(&nn, &nn, A, &lda, iw, &info);
        if(info != 0)
        {
            throw KeyCppException("Unknown error in inv()!");
        }

        /* Computes the reciprocal norm */
        dgecon_("1", &nn, A, &lda, &anorm, &rcond, w1, w2, &info);
        if(info != 0)
        {
            throw KeyCppException("Unknown error in inv()!");
        }
        
        if(rcond < 1e-15)
        {
            std::cerr << "Warning: Matrix is close to singular or badly scaled. Results may be inaccurate.\nrcond = " << rcond << std::endl;
        }

        dgetri_(&nn,A,&nn,iw,w1,&lwork,&info);

        matrix<double> A_out(A_in.size(1),A_in.size(2));
        for(unsigned int ii = 0; ii < n; ii++)
        {
            for(unsigned int jj = 0; jj < n; jj++)
            {
                A_out(jj,ii) = A[ii*n + jj];
            }
        }

        delete [] iw;
        delete [] w1;
        delete [] w2;
        delete [] A;
        
        return A_out;
    }
    
    inline matrix<std::complex<double>> inv(const matrix<std::complex<double>>& A_in)
	{
		if(A_in.size(1) <= 0 || A_in.size(2) <= 0)
		{
			throw KeyCppException("Error in inv()! Empty matrix supplied!");
		}
		if(A_in.size(1) != A_in.size(2))
		{
			throw KeyCppException("Error in inv()! Matrix must be square!");
		}
		
		unsigned int n = (unsigned)A_in.size(1);
		int nn = n;

        int info, lda;
        double anorm, rcond;
        
        int *iw = new int[A_in.size(1)];
        int lwork = (int)(A_in.size(1)*A_in.size(2)) + 64;
        std::complex<double> *w1 = new std::complex<double>[lwork];
        double *w2 = new double[lwork];
        std::complex<double> *A = new std::complex<double>[A_in.size(1)*A_in.size(2)];
        for(unsigned int ii = 0; ii < A_in.size(2); ii++)
        {
            for(unsigned int jj = 0; jj < A_in.size(1); jj++)
            {
                A[ii*A_in.size(1) + jj] = A_in(jj,ii);
            }
        }
        lda = n;

        /* Computes the norm of A */
        anorm = zlange_("1", &nn, &nn, A, &lda, w2);

        /* Modifies A in place with a LU decomposition */
        zgetrf_(&nn, &nn, A, &lda, iw, &info);
        if(info != 0)
        {
            throw KeyCppException("Unknown error in inv()!");
        }

        /* Computes the reciprocal norm */
        zgecon_("1", &nn, A, &lda, &anorm, &rcond, w1, w2, &info);
        if(info != 0)
        {
            throw KeyCppException("Unknown error in inv()!");
        }
        
        if(rcond < 1e-15)
        {
            std::cerr << "Warning: Matrix is close to singular or badly scaled. Results may be inaccurate.\nrcond = " << rcond << std::endl;
        }
        
        zgetri_(&nn,A,&nn,iw,w1,&lwork,&info);

        matrix<std::complex<double>> A_out(A_in.size(1),A_in.size(2));
        for(unsigned int ii = 0; ii < n; ii++)
        {
            for(unsigned int jj = 0; jj < n; jj++)
            {
                A_out(jj,ii) = A[ii*n + jj];
            }
        }
        
        delete [] iw;
        delete [] w1;
        delete [] w2;
        delete [] A;
        
		return A_out;
	}
	
	inline double norm(const matrix<double,2,0> &A_in, std::string method)
	{
	    if(A_in.empty())
		{
			throw KeyCppException("Error in norm()! Empty matrix supplied!");
		}
		std::transform(method.begin(), method.end(), method.begin(), ::tolower);
		if(method.compare("1") != 0 && method.compare("2") != 0 && method.compare("inf") != 0 && method.compare("fro") != 0)
		{
		    throw KeyCppException("Unknown norm type!");
		}
		if(method.compare("fro") == 0)
		{
		    method = "f";
		}
		else if(method.compare("inf") == 0)
		{
		    method = "i";
		}
		
        double anorm;
		if(method.compare("2") == 0)
		{
		    if(A_in.isVec())
		    {
		        anorm = sqrt(sum(times(abs(A_in),abs(A_in))));
		    }
		    else
		    {
		        auto svd_out = svd(A_in);
		        anorm = max(max(svd_out.S));
		    }
		}
		else
		{
	        int m = (int)A_in.size(1);
	        int n = (int)A_in.size(2);

            int lda;
            int lwork = (int)(A_in.size(1)*A_in.size(2)) + 64;
            double *w1 = new double[lwork];
            lda = m;
            anorm = dlange_(method.c_str(), &m, &n, &A_in.mData[0], &lda, w1);
            
            delete [] w1;
        }
        
        return anorm;
	}
	
	/** \brief Computes the singular value decomposition of matrix A_in.
	 */
	inline SVD_type<double,double> svd(const matrix<double> &A_in, std::string method)
	{
	    if(A_in.empty())
		{
			throw KeyCppException("Error in svd()! Empty matrix supplied!");
		}
		std::transform(method.begin(), method.end(), method.begin(), ::tolower);
		
	    unsigned int m = (unsigned)A_in.size(1);
	    unsigned int n = (unsigned)A_in.size(2);
	    int nn = n, mm = m;
        
        unsigned int ldvt;
        int info, lda, ldu, ldvt2;
        
        unsigned int col_u;
        
        SVD_type<double,double> out;
        int lwork = (int)(A_in.size(1)*A_in.size(2)) + 64;
        double *work = new double[lwork];
        double *A = new double[A_in.size(1)*A_in.size(2)];
        double *U;
        double *S;
        double *VT;
        unsigned int min_dim;
        if(m > n)
        {
            min_dim = n;
        }
        else
        {
            min_dim = m;
        }
        S = new double[min_dim];
        for(unsigned int ii = 0; ii < A_in.size(2); ii++)
        {
            for(unsigned int jj = 0; jj < A_in.size(1); jj++)
            {
                A[ii*A_in.size(1) + jj] = A_in(jj,ii);
            }
        }
        lda = m;
        ldu = m;
        std::string jobu, jobvt;
        if(method.compare("0") == 0)
        {
            jobu = "S";
            jobvt = "A";
            U = new double[m*min_dim];
            col_u = min_dim;
            VT = new double[n*n];
            ldvt = n;
            if(m > n)
            {
                out.S = matrix<double>(n,n);
            }
            else
            {
                out.S = matrix<double>(m,n);
            }
        }
        else if(method.compare("econ") == 0)
        {
            jobu = "S";
            jobvt = "S";
            U = new double[m*min_dim];
            col_u = min_dim;
            VT = new double[n*min_dim];
            if(m > n)
            {
                ldvt = n;
            }
            else
            {
                ldvt = m;
            }
            if(m >= n)
            {
                out.S = matrix<double>(n,n);
            }
            else if(m < n)
            {
                out.S = matrix<double>(m,m);
            }
            else
            {
                out.S = matrix<double>(m,n);
            }
        }
        else if(method.empty())
        {
            jobu = "A";
            jobvt = "A";
            U = new double[m*m];
            col_u = m;
            VT = new double[n*n];
            ldvt = n;
            out.S = matrix<double>(m,n);
        }
        else
        {
            throw KeyCppException("Unknown argument to svd()!");
        }
		ldvt2 = ldvt;
		dgesvd_(jobu.c_str(), jobvt.c_str(), &mm, &nn, A, &lda, S, U, &ldu, VT, &ldvt2, work, &lwork, &info);
                 
        if(info != 0)
        {
            throw KeyCppException("Unknown error in svd()!");
        }
        
        out.U = matrix<double>(m,col_u);
        out.V = matrix<double>(n,ldvt);
        
        for(unsigned int ii = 0; ii < min_dim; ii++)
        {
            out.S(ii,ii) = S[ii];
        }
        for(unsigned int ii = 0; ii < m; ii++)
        {
            for(unsigned int jj = 0; jj < col_u; jj++)
            {
                out.U(ii,jj) = U[jj*m + ii];
            }
        }
        
        for(unsigned int ii = 0; ii < n; ii++)
        {
            for(unsigned int jj = 0; jj < ldvt; jj++)
            {
                out.V(ii,jj) = VT[ii*ldvt + jj];
            }
        }
                 
        delete [] work;
        delete [] A;
        delete [] U;
        delete [] S;
        delete [] VT;
        
        return out;
	}
	
	inline double norm(const matrix<std::complex<double>,2,0> &A_in, std::string method)
	{
	    if(A_in.size(1) <= 0 || A_in.size(2) <= 0)
		{
			throw KeyCppException("Error in norm()! Empty matrix supplied!");
		}
		std::transform(method.begin(), method.end(), method.begin(), ::tolower);
		if(method.compare("1") != 0 && method.compare("2") != 0 && method.compare("inf") != 0 && method.compare("fro") != 0)
		{
		    throw KeyCppException("Unknown norm type!");
		}
		if(method.compare("fro") == 0)
		{
		    method = "f";
		}
		else if(method.compare("inf") == 0)
		{
		    method = "i";
		}
		
        double anorm;
		if(method.compare("2") == 0)
		{
		    if(A_in.isVec())
		    {
		        anorm = sqrt(sum(times(abs(A_in),abs(A_in))));
		    }
		    else
		    {
		        auto svd_out = svd(A_in);
		        anorm = max(max(svd_out.S));
		    }
		}
		else
		{
	        int m = (int)A_in.size(1);
	        int n = (int)A_in.size(2);

            int lda;
            int lwork = (int)(A_in.size(1)*A_in.size(2)) + 64;
            double *w1 = new double[lwork];
            lda = m;
            anorm = zlange_(method.c_str(), &m, &n, &A_in.mData[0], &lda, w1);
            
            delete [] w1;
        }
        
        return anorm;
	}
	
	/** \brief Computes the singular value decomposition of matrix A_in.
	 */
	inline SVD_type<std::complex<double>, double> svd(const matrix<std::complex<double>> &A_in, std::string method)
	{
	    if(A_in.size(1) <= 0 || A_in.size(2) <= 0)
		{
			throw KeyCppException("Error in svd()! Empty matrix supplied!");
		}
		std::transform(method.begin(), method.end(), method.begin(), ::tolower);
		
	    unsigned int m = (unsigned)A_in.size(1);
	    unsigned int n = (unsigned)A_in.size(2);
	    int mm = (int)m, nn = (int)n;

        unsigned int ldvt;
        int info, lda, ldu, ldvt2;
        
        unsigned int col_u;
        
        SVD_type<std::complex<double>,double> out;
        int lwork = (int)(A_in.size(1)*A_in.size(2)) + 64;
        double *rwork = new double[lwork];
        std::complex<double> *work = new std::complex<double>[lwork];
        std::complex<double> *A = new std::complex<double>[A_in.size(1)*A_in.size(2)];
        std::complex<double> *U;
        double *S;
        std::complex<double> *VT;
        unsigned int min_dim;
        if(m > n)
        {
            min_dim = n;
        }
        else
        {
            min_dim = m;
        }
        S = new double[min_dim];
        for(unsigned int ii = 0; ii < A_in.size(2); ii++)
        {
            for(unsigned int jj = 0; jj < A_in.size(1); jj++)
            {
                A[ii*A_in.size(1) + jj] = A_in(jj,ii);
            }
        }
        lda = m;
        ldu = m;
        std::string jobu, jobvt;
        if(method.compare("0") == 0)
        {
            jobu = "S";
            jobvt = "A";
            U = new std::complex<double>[m*min_dim];
            col_u = min_dim;
            VT = new std::complex<double>[n*n];
            ldvt = n;
            if(m > n)
            {
                out.S = matrix<double>(n,n);
            }
            else
            {
                out.S = matrix<double>(m,n);
            }
        }
        else if(method.compare("econ") == 0)
        {
            jobu = "S";
            jobvt = "S";
            U = new std::complex<double>[m*min_dim];
            col_u = min_dim;
            VT = new std::complex<double>[n*min_dim];
            if(m > n)
            {
                ldvt = n;
            }
            else
            {
                ldvt = m;
            }
            if(m >= n)
            {
                out.S = matrix<double>(n,n);
            }
            else if(m < n)
            {
                out.S = matrix<double>(m,m);
            }
            else
            {
                out.S = matrix<double>(m,n);
            }
        }
        else if(method.empty())
        {
            jobu = "A";
            jobvt = "A";
            U = new std::complex<double>[m*m];
            col_u = m;
            VT = new std::complex<double>[n*n];
            ldvt = n;
            out.S = matrix<double>(m,n);
        }
        else
        {
            throw KeyCppException("Unknown argument to svd()!");
        }
		ldvt2 = ldvt;
		zgesvd_(jobu.c_str(), jobvt.c_str(), &mm, &nn, A, &lda, S, U, &ldu, VT, &ldvt2, work, &lwork, rwork, &info);
                 
        if(info != 0)
        {
            throw KeyCppException("Unknown error in SVD()!");
        }
        
        out.U = matrix<std::complex<double>>(m,col_u);
        out.V = matrix<std::complex<double>>(n,ldvt);
        
        for(unsigned int ii = 0; ii < min_dim; ii++)
        {
            out.S(ii,ii) = S[ii];
        }
        for(unsigned int ii = 0; ii < m; ii++)
        {
            for(unsigned int jj = 0; jj < col_u; jj++)
            {
                out.U(ii,jj) = U[jj*m + ii];
            }
        }
        
        for(unsigned int ii = 0; ii < n; ii++)
        {
            for(unsigned int jj = 0; jj < ldvt; jj++)
            {
                out.V(ii,jj) = VT[ii*ldvt + jj];
            }
        }
                 
        delete [] rwork;
        delete [] work;
        delete [] A;
        delete [] U;
        delete [] S;
        delete [] VT;
        
        return out;
	}
	
	inline matrix<std::complex<double>> lu(const matrix<std::complex<double>>& A_in)
	{
		if(A_in.empty())
		{
			throw KeyCppException("Error in lu()! Empty matrix supplied!");
		}
		
		int n = (int)A_in.size(2);
		int m = (int)A_in.size(1);

        int info, lda;
        
        int *iw = new int[A_in.size(1)];
        std::complex<double> *A = new std::complex<double>[A_in.size(1)*A_in.size(2)];
        for(unsigned int ii = 0; ii < A_in.size(2); ii++)
        {
            for(unsigned int jj = 0; jj < A_in.size(1); jj++)
            {
                A[ii*A_in.size(1) + jj] = A_in(jj,ii);
            }
        }
        lda = m;
        
        /* Modifies A in place with a LU decomposition */
        zgetrf_(&m, &n, A, &lda, iw, &info);
        if(info != 0)
        {
            if(info > 0)
            {
                std::cerr << "Warning: Matrix is singular. Results may be inaccurate.\n";
            }
            else
            {
                throw KeyCppException("Unknown error in lu()!");
            }
        }
        
        matrix<std::complex<double>> A_out(A_in.size(1),A_in.size(2));
        for(unsigned int ii = 0; ii < A_out.size(2); ii++)
        {
            for(unsigned int jj = 0; jj < A_out.size(1); jj++)
            {
                 A_out(jj,ii) = A[ii*A_out.size(1) + jj];
            }
        }
        
        delete [] A;
        delete [] iw;
        
        return A_out;
    }
	
	inline matrix<double> lu(const matrix<double>& A_in)
	{
		if(A_in.empty())
		{
			throw KeyCppException("Error in lu()! Empty matrix supplied!");
		}
		
		int n = (int)A_in.size(2);
		int m = (int)A_in.size(1);

        int info, lda;
        
        int *iw = new int[A_in.size(1)];
        double *A = new double[A_in.size(1)*A_in.size(2)];
        for(size_t ii = 0; ii < A_in.size(2); ii++)
        {
            for(size_t jj = 0; jj < A_in.size(1); jj++)
            {
                A[ii*A_in.size(1) + jj] = A_in(jj,ii);
            }
        }
        lda = m;
        
        /* Modifies A in place with a LU decomposition */
        dgetrf_(&m, &n, A, &lda, iw, &info);
        if(info != 0)
        {
            if(info > 0)
            {
                std::cerr << "Warning: Matrix is singular. Results may be inaccurate.\n";
            }
            else
            {
                throw KeyCppException("Unknown error in lu()!");
            }
        }
        
        matrix<double> A_out(A_in.size(1),A_in.size(2));
        for(size_t ii = 0; ii < A_out.size(2); ii++)
        {
            for(size_t jj = 0; jj < A_out.size(1); jj++)
            {
                 A_out(jj,ii) = A[ii*A_out.size(1) + jj];
            }
        }
        
        delete [] A;
        delete [] iw;
        
        return A_out;
    }
	
	inline matrix<std::complex<double>> lu(const matrix<std::complex<double>>& A_in, int *iw)
	{
		if(A_in.empty())
		{
			throw KeyCppException("Error in lu()! Empty matrix supplied!");
		}
		
		int n = (int)A_in.size(2);
		int m = (int)A_in.size(1);

        int info, lda;
        
        std::complex<double> *A = new std::complex<double>[A_in.size(1)*A_in.size(2)];
        for(unsigned int ii = 0; ii < A_in.size(2); ii++)
        {
            for(unsigned int jj = 0; jj < A_in.size(1); jj++)
            {
                A[ii*A_in.size(1) + jj] = A_in(jj,ii);
            }
        }
        lda = m;
        
        /* Modifies A in place with a LU decomposition */
        zgetrf_(&m, &n, A, &lda, iw, &info);
        if(info != 0)
        {
            if(info > 0)
            {
                std::cerr << "Warning: Matrix is singular. Results may be inaccurate.\n";
            }
            else
            {
                throw KeyCppException("Unknown error in lu()!");
            }
        }
        
        matrix<std::complex<double>> A_out(A_in.size(1),A_in.size(2));
        for(unsigned int ii = 0; ii < A_out.size(2); ii++)
        {
            for(unsigned int jj = 0; jj < A_out.size(1); jj++)
            {
                 A_out(jj,ii) = A[ii*A_out.size(1) + jj];
            }
        }
        
        delete [] A;
        
        return A_out;
    }
	
	inline matrix<double> lu(const matrix<double>& A_in, int *iw)
	{
		if(A_in.empty())
		{
			throw KeyCppException("Error in lu()! Empty matrix supplied!");
		}
		
		int n = (int)A_in.size(2);
		int m = (int)A_in.size(1);

        int info, lda;
        
        double *A = new double[A_in.size(1)*A_in.size(2)];
        for(size_t ii = 0; ii < A_in.size(2); ii++)
        {
            for(size_t jj = 0; jj < A_in.size(1); jj++)
            {
                A[ii*A_in.size(1) + jj] = A_in(jj,ii);
            }
        }
        lda = m;
        
        /* Modifies A in place with a LU decomposition */
        dgetrf_(&m, &n, A, &lda, iw, &info);
        if(info != 0)
        {
            if(info > 0)
            {
                std::cerr << "Warning: Matrix is singular. Results may be inaccurate.\n";
            }
            else
            {
                throw KeyCppException("Unknown error in lu()!");
            }
        }
        
        matrix<double> A_out(A_in.size(1),A_in.size(2));
        for(size_t ii = 0; ii < A_out.size(2); ii++)
        {
            for(size_t jj = 0; jj < A_out.size(1); jj++)
            {
                 A_out(jj,ii) = A[ii*A_out.size(1) + jj];
            }
        }
        
        delete [] A;
        
        return A_out;
    }
    
    template<class T, class U>
    struct meshgrid_type
    {
        matrix<T> X;
        matrix<U> Y;
    };
    
    template<class T, class U>
    meshgrid_type<T,U> meshgrid(const matrix<T> &x, const matrix<U> &y)
    {
        meshgrid_type<T,U> return_struct;
        if(!x.isVec() || !y.isVec())
        {
            throw KeyCppException("Error in meshgrid! Inputs must be vectors!");
        }
        
        if(x.size(1) == 1)
        {
            return_struct.X = repmat(x,numel(y),1);
        }
        else
        {
            return_struct.X = repmat(transpose(x),numel(y),1);
        }
        return_struct.Y = repmat(y,1,numel(x));
        if(y.size(2) == 1)
        {
            return_struct.Y = repmat(y,1,numel(x));
        }
        else
        {
            return_struct.Y = repmat(transpose(y),1,numel(x));
        }
        return return_struct;
    }
    
    inline matrix<double> gallery(std::string type, int n, int k)
    {
		std::transform(type.begin(), type.end(), type.begin(), ::tolower);
		
		if(type.compare("orthog") == 0)
		{
		    if(k == -1)
		    {
		        matrix<double> Q(n,n);
		        for(int ii = 0; ii < n; ii++)
		        {
		            for(int jj = 0; jj < n; jj++)
		            {
		                double iii = ii, jjj = jj;
		                Q(ii,jj) = std::cos((iii)*(jjj)*pi/(n-1.0));
		            }
		        }
		        return Q;
		    }
		    else
		    {
		        throw KeyCppException("Unrecognized option in gallery!");
		    }
		}
    }
    
    /*inline std::complex<double> str2complex(std::string str)
    {
        std::string dummy = removeWhiteSpace(str);
        std::stringstream ss;
        ss << dummy;
        std::getline(ss,dummy,'+');
        dummy = removeWhiteSpace(dummy);
        
        // Check if that + was in front of the first number
        if(dummy.empty())
        {
            std::stringstream ss2;
            ss2 << ss.str();
            
            std::getline(ss2,dummy,'+');
            if(ss2.rdbuf()->in_avail() == 0)
            {
		        ss2.str("");
		        ss2.clear();
		        ss2 << ss.str();
                
                std::getline(ss2,dummy,'-');
                if(ss2.rdbuf()->in_avail() == 0)
                {
                    std::cerr << "Invalid string in str2complex.\n";
                    return std::complex<double>();
                }
                // dummy should now contain the first number which is positive:
                
            }
        }
    }*/
    
    template<class T>
    matrix<T> pow(const matrix<T> &A, double n)
    {
        matrix<T> B(A.size(1), A.size(2));
        for(size_t ii = 0; ii < A.numel(); ii++)
        {
            B(ii) = std::pow(A(ii),n);
        }
        
        return B;
    }
}

#include "znaupd.h"

#endif
