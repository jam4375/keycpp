// Matlab.h -- Common Matlab functions implemented in C++
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

	vector_k<std::complex<double> > eig(const matrix<std::complex<double> > &A,
	                                       const matrix<std::complex<double> > &B,
	                                       matrix<std::complex<double> > *vr_return = NULL,
	                                       matrix<std::complex<double> > *vl_return = NULL);

	vector_k<std::complex<double> > eig(const matrix<std::complex<double> > &A,
	                                       matrix<std::complex<double> > *vr_return = NULL,
	                                       matrix<std::complex<double> > *vl_return = NULL);
	                                       
	vector_k<std::complex<double>> eig(const matrix<double> &A,
                                          matrix<std::complex<double>> *vr_return = NULL,
                                          matrix<std::complex<double>> *vl_return = NULL);
	
	double rcond(const matrix<double> &A);
	double rcond(const matrix<std::complex<double>> &A);
	vector_k<std::complex<double>> linsolve(const matrix<std::complex<double>>& A_in,
                                               const vector_k<std::complex<double>>& b_in);
	vector_k<double> linsolve(const matrix<double>& A_in,
                                 const vector_k<double>& b_in);
    matrix<double> inv(const matrix<double>& A_in);
    matrix<std::complex<double>> inv(const matrix<std::complex<double>>& A_in);
	
	template<class T, size_t dim> matrix<size_t,1> size(const matrix<T,dim> &A);
	
	template<class T,size_t dim>
	matrix<T,dim> eop(const matrix<T,dim> &A, T (*f)(const T&))
	{
	    matrix<T,dim> B;
	    B.resize(size(A));
	    for(size_t ii = 0; ii < A.mData.size(); ii++)
	    {
	        B.mData[ii] = (*f)(A.mData[ii]);
	    }
	    return B;
	}
	
	template<class T,size_t dim>
	matrix<T,dim> eop(const matrix<std::complex<T>,dim> &A, T (*f)(const std::complex<T>&))
	{
	    matrix<T,dim> B;
	    B.resize(size(A));
	    for(size_t ii = 0; ii < A.mData.size(); ii++)
	    {
	        B.mData[ii] = (*f)(A.mData[ii]);
	    }
	    return B;
	}
	
	template<class T>
	vector_k<T> eop(const vector_k<T> &v1, T (*f)(const T&))
	{
	    vector_k<T> v2(v1.size());
	    for(size_t ii = 0; ii < v1.size(); ii++)
	    {
	        v2[ii] = (*f)(v1[ii]);
	    }
	    return v2;
	}
	
	template<class T, size_t dim>
	matrix<T,dim> eop(const matrix<T,dim> &A, T (*f)(T))
	{
	    matrix<T,dim> B;
	    B.resize(size(A));
	    for(size_t ii = 0; ii < A.mData.size(); ii++)
	    {
	        B.mData[ii] = (*f)(A.mData[ii]);
	    }
	    return B;
	}
	
	template<class T>
	vector_k<T> eop(const vector_k<T> &v1, T (*f)(T))
	{
	    vector_k<T> v2(v1.size());
	    for(size_t ii = 0; ii < v1.size(); ii++)
	    {
	        v2[ii] = (*f)(v1[ii]);
	    }
	    return v2;
	}
	
	template<class T>
	vector_k<T> eop(const vector_k<std::complex<T>> &v1, T (*f)(const std::complex<T>&))
	{
	    vector_k<T> v2(v1.size());
	    for(size_t ii = 0; ii < v1.size(); ii++)
	    {
	        v2[ii] = (*f)(v1[ii]);
	    }
	    return v2;
	}
	
	template<class T,size_t dim>
	matrix<T,dim> eop(const matrix<std::complex<T>,dim> &A, T (*f)(std::complex<T>))
	{
	    matrix<T,dim> B;
	    B.resize(size(A));
	    for(size_t ii = 0; ii < A.mData.size(); ii++)
	    {
	        B.mData[ii] = (*f)(A.mData[ii]);
	    }
	    return B;
	}
	
	template<class T>
	vector_k<T> eop(const vector_k<std::complex<T>> &v1, T (*f)(std::complex<T>))
	{
	    vector_k<T> v2(v1.size());
	    for(size_t ii = 0; ii < v1.size(); ii++)
	    {
	        v2[ii] = (*f)(v1[ii]);
	    }
	    return v2;
	}
	
	
	template<class T, class U> struct observe
	{
		matrix<T,2>& y;
		matrix<U,1>& x_ode;
		size_t ii;

		observe(matrix<T,2> &p_y, matrix<U,1> &p_x_ode) : y(p_y), x_ode(p_x_ode), ii(0) { };

		void operator()(const matrix<T,1> &y_temp, U x_temp)
		{
		    if(ii >= y.size(1))
		    {
		        matrix<size_t,1> temp_size = size(y);
		        temp_size(0)++;
		        y.resize(temp_size);
		        x_ode.resize(temp_size(0));
		    }
		    x_ode.mData[ii] = x_temp;
		    for(size_t jj = 0; jj < y_temp.size(1); jj++)
		    {
			    y(ii,jj) = y_temp(jj);
			}
			ii++;
		}
	};
	
	template<class T, size_t dim>
	std::ostream& operator<<(std::ostream &out, const matrix<T,dim> &A)
	{
	    static_assert((dim == 2 || dim == 1),"This function is only available for matrices of dimension 1 or 2.");
	    
	    if(dim == 2)
	    {
            for(size_t ii = 0; ii < A.size(1); ii++)
            {
                for(size_t jj = 0; jj < A.size(2); jj++)
                {
                    out << A(ii,jj) << " ";
                }
                out << std::endl;
            }
        }
        else if(dim == 1)
        {
            for(size_t ii = 0; ii < A.size(1); ii++)
            {
                out << A.mData[ii] << std::endl;
            }
        }
        
        return out;
    }
	
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

    /** \brief Returns the product of all the elements of the vector x.
     */
	template<class T> T prod(const vector_k<T> &x)
	{
	    if(x.empty())
	    {
	        return 1.0;
	    }
		T x_out;
		x_out = x[0];
		for(size_t ii = 1; ii < x.size(); ii++)
		{
			x_out *= x[ii];
		}
		return x_out;
	}

    /** \brief Returns a vector containing the product of all the elements in each
     *         column of the matrix A.
     */
	template<class T> vector_k<T> prod(const matrix<T> &A)
	{
	    if(A.size(1) <= 0 || A.size(2) <= 0)
	    {
	        vector_k<T> x(1);
	        x[0] = 1.0;
	        return x;
	    }
		vector_k<T> B = A.row(0);
		for(size_t jj = 0; jj < A.size(2); jj++)
		{
		    for(size_t ii = 1; ii < A.size(1); ii++)
		    {
			    B[jj] *= A(ii,jj);
			}
		}
		return B;
	}
	
	/** \brief Returns a vector of differences between adjacent elements.
     */
	template<class T> vector_k<T> diff(const vector_k<T> &v1)
	{
	    if(v1.empty())
	    {
	        throw KeyCppException("Cannot compute diff() on empty vector!");
	    }
		vector_k<T> v2(v1.size()-1);
		for(size_t ii = 0; ii < v2.size(); ii++)
		{
		    v2[ii] = v1[ii+1] - v1[ii];
		}
		return v2;
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

	template<class T> vector_k<std::complex<T>> conj(const vector_k<std::complex<T>> &v1)
	{
		return eop(v1,static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::conj));
	}

	template<class T,size_t dim> matrix<std::complex<T>,dim> conj(const matrix<std::complex<T>,dim> &A)
	{
		return eop(A,static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::conj));
	}

	template<class T> vector_k<T> real(const vector_k<std::complex<T>> &v1)
	{
		return eop(v1,static_cast<T (*)(const std::complex<T> &)>(&std::real));
	}

	template<class T,size_t dim> matrix<T,dim> real(const matrix<std::complex<T>,dim> &A)
	{
		return eop(A,static_cast<T (*)(const std::complex<T> &)>(&std::real));
	}

	template<class T> vector_k<T> imag(const vector_k<std::complex<T>> &v1)
	{
		return eop(v1,static_cast<T (*)(const std::complex<T> &)>(&std::imag));
	}

	template<class T, size_t dim> matrix<T,dim> imag(const matrix<std::complex<T>,dim> &A)
	{
		return eop(A,static_cast<T (*)(const std::complex<T> &)>(&std::imag));
	}

	template<class T> vector_k<T> abs(const vector_k<T> &v1)
	{
		return eop(v1,static_cast<T (*)(T)>(&std::abs));
	}

	template<class T> vector_k<T> abs(const vector_k<std::complex<T>> &v1)
	{
		return eop(v1,static_cast<T (*)(const std::complex<T> &)>(&std::abs));
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

	template<class T, class U> vector_k<decltype(std::declval<T>()*std::declval<U>())> operator+(const vector_k<T>& v1, const vector_k<U>& v2)
	{
	    if(v1.empty() || v2.empty())
	    {
	        throw KeyCppException("Cannot add empty vector!");
	    }
	    if(v1.size() != v2.size())
	    {
	        throw KeyCppException("Cannot add vectors of different sizes!");
	    }
		vector_k<decltype(std::declval<T>()*std::declval<U>())> result(v1.size());
		for(int ii = 0; ii < result.size(); ii++)
		{
			result[ii] = v1[ii]+v2[ii];
		}
		return result;
	}
	
	template<>
	inline vector_k<double> operator+(const vector_k<double>& v1, const vector_k<double>& v2)
	{
	    if(v1.empty() || v2.empty())
	    {
	        throw KeyCppException("Cannot add empty vector!");
	    }
	    if(v1.size() != v2.size())
	    {
	        throw KeyCppException("Cannot add vectors of different sizes!");
	    }
	    int N = (int)v1.size(), incx = (int)v1.get_inc(), incy = 1;
	    double alpha = 1.0;
		vector_k<double> result(v2);
	    daxpy_(&N, &alpha, &v1[0], &incx, &result[0], &incy);
		return result;
	}
	
	template<>
	inline vector_k<std::complex<double>> operator+(const vector_k<std::complex<double>>& v1, const vector_k<std::complex<double>>& v2)
	{
	    if(v1.empty() || v2.empty())
	    {
	        throw KeyCppException("Cannot add empty vector!");
	    }
	    if(v1.size() != v2.size())
	    {
	        throw KeyCppException("Cannot add vectors of different sizes!");
	    }
	    int N = (int)v1.size(), incx = (int)v1.get_inc(), incy = 1;
	    std::complex<double> alpha = 1.0;
		vector_k<std::complex<double>> result(v2);
	    zaxpy_(&N, &alpha, &v1[0], &incx, &result[0], &incy);
		return result;
	}

	template<class T, class U> vector_k<decltype(std::declval<T>()*std::declval<U>())> operator+(const vector_k<T>& v1, const U& a)
	{
		vector_k<decltype(std::declval<T>()*std::declval<U>())> result(v1.size());
		for(size_t ii = 0; ii < result.size(); ii++)
		{
			result[ii] = v1[ii]+a;
		}
		return result;
	}

	template<class T, class U> vector_k<decltype(std::declval<T>()*std::declval<U>())> operator+(const U& a, const vector_k<T>& v2)
	{
		vector_k<decltype(std::declval<T>()*std::declval<U>())> result(v2.size());
		for(size_t ii = 0; ii < result.size(); ii++)
		{
			result[ii] = v2[ii]+a;
		}
		return result;
	}

	template<class T, class U, size_t dim> matrix<decltype(std::declval<T>()*std::declval<U>()),dim> operator+(const matrix<T,dim>& A, const U& a)
	{
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> B;
		matrix<size_t,1> temp = size(A);
		B.resize(temp);
		for(size_t ii = 0; ii < B.mData.size(); ii++)
		{
			B.mData[ii] = a+A.mData[ii];
		}
		return B;
	}

	template<class T, class U> matrix<decltype(std::declval<T>()*std::declval<U>())> operator+(const U& a, const matrix<T>& A)
	{
		matrix<decltype(std::declval<T>()*std::declval<U>())> result(A.size(1),A.size(2));
		for(size_t ii = 0; ii < result.size(1); ii++)
		{
		    for(size_t jj = 0; jj < result.size(2); jj++)
		    {
			    result(ii,jj) = A(ii,jj)+a;
			}
		}
		return result;
	}
	
	

	template<class T, class U> vector_k<decltype(std::declval<T>()*std::declval<U>())> operator-(const vector_k<T>& v1, const U& a)
	{
		vector_k<decltype(std::declval<T>()*std::declval<U>())> result(v1.size());
		for(size_t ii = 0; ii < result.size(); ii++)
		{
			result[ii] = v1[ii]-a;
		}
		return result;
	}

	template<class T, class U> vector_k<decltype(std::declval<T>()*std::declval<U>())> operator-(const U& a, const vector_k<T>& v2)
	{
		vector_k<decltype(std::declval<T>()*std::declval<U>())> result(v2.size());
		for(size_t ii = 0; ii < result.size(); ii++)
		{
			result[ii] = a - v2[ii];
		}
		return result;
	}

	template<class T, class U,size_t dim> matrix<decltype(std::declval<T>()*std::declval<U>()),dim> operator-(const matrix<T,dim>& A, const U& a)
	{
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> result;
		result.resize(size(A));
		for(size_t ii = 0; ii < result.mData.size(); ii++)
		{
			result.mData[ii] = A.mData[ii]-a;
		}
		return result;
	}

	template<class T, class U,size_t dim> matrix<decltype(std::declval<T>()*std::declval<U>()),dim> operator-(const U& a, const matrix<T,dim>& A)
	{
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> result;
		result.resize(size(A));
		for(size_t ii = 0; ii < result.mData.size(); ii++)
		{
			result.mData[ii] = a-A.mData[ii];
		}
		return result;
	}

	template<class T, class U> vector_k<decltype(std::declval<T>()*std::declval<U>())> operator-(const vector_k<T>& v1, const vector_k<U>& v2)
	{
	    if(v1.empty() || v2.empty())
	    {
	        throw KeyCppException("Cannot subtract empty vector!");
	    }
	    if(v1.size() != v2.size())
	    {
	        throw KeyCppException("Cannot subtract vectors of different sizes!");
	    }
		vector_k<decltype(std::declval<T>()*std::declval<U>())> result(v1.size());
		for(size_t ii = 0; ii < result.size(); ii++)
		{
			result[ii] = v1[ii]-v2[ii];
		}
		return result;
	}
	
	template<>
	inline vector_k<double> operator-(const vector_k<double>& v1, const vector_k<double>& v2)
	{
	    if(v1.empty() || v2.empty())
	    {
	        throw KeyCppException("Cannot subtract empty vector!");
	    }
	    if(v1.size() != v2.size())
	    {
	        throw KeyCppException("Cannot subtract vectors of different sizes!");
	    }
	    int N = (int)v1.size(), incx = (int)v2.get_inc(), incy = 1;
	    double alpha = -1.0;
		vector_k<double> result(v1);
	    daxpy_(&N, &alpha, &v2[0], &incx, &result[0], &incy);
		return result;
	}
	
	template<>
	inline vector_k<std::complex<double>> operator-(const vector_k<std::complex<double>>& v1, const vector_k<std::complex<double>>& v2)
	{
	    if(v1.empty() || v2.empty())
	    {
	        throw KeyCppException("Cannot subtract empty vector!");
	    }
	    if(v1.size() != v2.size())
	    {
	        throw KeyCppException("Cannot subtract vectors of different sizes!");
	    }
	    int N = (int)v1.size(), incx = (int)v2.get_inc(), incy = 1;
	    std::complex<double> alpha = -1.0;
		vector_k<std::complex<double>> result(v1);
	    zaxpy_(&N, &alpha, &v2[0], &incx, &result[0], &incy);
		return result;
	}

	template<class T, class U,size_t dim> matrix<decltype(std::declval<T>()*std::declval<U>()),dim> operator*(const T& a, const matrix<U,dim>& A)
	{
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> B;
		matrix<size_t,1> temp = size(A);
		B.resize(temp);
		for(size_t ii = 0; ii < B.mData.size(); ii++)
		{
			B.mData[ii] = a*A.mData[ii];
		}
		return B;
	}
	
	template<class T, class U,size_t dim> matrix<decltype(std::declval<T>()*std::declval<U>()),dim> operator*(const matrix<U,dim>& A, const T& a)
	{
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> B;
		matrix<size_t,1> temp = size(A);
		B.resize(temp);
		for(size_t ii = 0; ii < B.mData.size(); ii++)
		{
			B.mData[ii] = a*A.mData[ii];
		}
		return B;
	}
	
	template<class T, class U> matrix<decltype(std::declval<T>()*std::declval<U>())> operator*(const vector_k<T>& v1, const matrix<U>& A)
	{
		if(A.size(1) == 1)
		{
			matrix<decltype(std::declval<T>()*std::declval<U>())> B(v1.size(),A.size(2));
			for(size_t ii = 0; ii < B.size(1); ii++)
			{
				for(size_t jj = 0; jj < B.size(2); jj++)
				{
					B(ii,jj) = v1[ii]*A(0,jj);
				}
			}
			return B;	
		}
		else if(A.size(1) == v1.size())
		{
			matrix<decltype(std::declval<T>()*std::declval<U>())> B(1,A.size(2));
			for(size_t ii = 0; ii < A.size(2); ii++)
			{
				for(size_t jj = 0; jj < A.size(1); jj++)
				{
					B(0,ii) = v1[jj]*A(jj,ii);
				}
			}
			return B;
		}
		else
		{
			return matrix<decltype(std::declval<T>()*std::declval<U>())>();
		}
	}

	template<class T, class U> vector_k<decltype(std::declval<T>()*std::declval<U>())> operator*(const T& a, const vector_k<U>& v1)
	{
		vector_k<decltype(std::declval<T>()*std::declval<U>())> v2(v1.size());
		for(size_t ii = 0; ii < v2.size(); ii++)
		{
			v2[ii] = a*v1[ii];
		}
		return v2;
	}

	template<>
	inline vector_k<double> operator*(const double& a, const vector_k<double>& v1)
	{
		vector_k<double> v2(v1);
		int N = (int)v1.size(), incx = (int)v1.get_inc();
		dscal_(&N, &a, &v2[0], &incx);
		return v2;
	}

	template<>
	inline vector_k<std::complex<double>> operator*(const std::complex<double>& a, const vector_k<std::complex<double>>& v1)
	{
		vector_k<std::complex<double>> v2(v1);
		int N = (int)v1.size(), incx = (int)v2.get_inc();
		zscal_(&N, &a, &v2[0], &incx);
		return v2;
	}
	
	template<class T, class U> vector_k<decltype(std::declval<T>()*std::declval<U>())> operator*(const vector_k<T>& v1, const U& a)
	{
		vector_k<decltype(std::declval<T>()*std::declval<U>())> v2(v1.size());
		for(size_t ii = 0; ii < v2.size(); ii++)
		{
			v2[ii] = a*v1[ii];
		}
		return v2;
	}

	template<>
	inline vector_k<double> operator*(const vector_k<double>& v1, const double& a)
	{
		vector_k<double> v2(v1);
		int N = (int)v1.size(), incx = (int)v2.get_inc();
		dscal_(&N, &a, &v2[0], &incx);
		return v2;
	}

	template<>
	inline vector_k<std::complex<double>> operator*(const vector_k<std::complex<double>>& v1, const std::complex<double>& a)
	{
		vector_k<std::complex<double>> v2(v1);
		int N = (int)v1.size(), incx = (int)v2.get_inc();
		zscal_(&N, &a, &v2[0], &incx);
		return v2;
	}
	
	template<class T> vector_k<T> operator-(const vector_k<T>& v1)
	{
		vector_k<T> v2(v1.size());
		for(size_t ii = 0; ii < v2.size(); ii++)
		{
			v2[ii] = -v1[ii];
		}
		return v2;
	}
	
	template<>
	inline vector_k<double> operator-(const vector_k<double>& v1)
	{
		vector_k<double> v2(v1);
		double a = -1.0;
		int N = (int)v1.size(), incx = (int)v2.get_inc();
		dscal_(&N, &a, &v2[0], &incx);
		return v2;
	}
	
	template<>
	inline vector_k<std::complex<double>> operator-(const vector_k<std::complex<double>>& v1)
	{
		vector_k<std::complex<double>> v2(v1);
		std::complex<double> a = -1.0;
		int N = (int)v1.size(), incx = (int)v2.get_inc();
		zscal_(&N, &a, &v2[0], &incx);
		return v2;
	}
	
	template<class T,size_t dim> matrix<T,dim> operator-(const matrix<T,dim>& A)
	{
		matrix<T,dim> B;
		B.resize(size(A));
		for(size_t ii = 0; ii < B.mData.size(); ii++)
		{
			B.mData[ii] = -A.mData[ii];
		}
		return B;
	}
	
	template<class T> vector_k<T> operator+(const vector_k<T>& v1)
	{
		vector_k<T> v2(v1);
		return v2;
	}
	
	template<class T> matrix<T> operator+(const matrix<T>& A)
	{
		matrix<T> B(A.size(1),A.size(2));
		for(size_t ii = 0; ii < B.size(1); ii++)
		{
			for(size_t jj = 0; jj < B.size(2); jj++)
			{
				B(ii,jj) = A(ii,jj);
			}
		}
		return B;
	}
	
	template<class T, class U, size_t dim> matrix<decltype(std::declval<T>()*std::declval<U>()),dim> operator/(const matrix<T,dim>& A, const U& a)
	{
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> B;
		matrix<size_t,1> temp = size(A);
		B.resize(temp);
		
		for(size_t ii = 0; ii < A.mData.size(); ii++)
		{
			B.mData[ii] = A.mData[ii]/a;
		}
		return B;
	}
	
	template<class T, class U> vector_k<decltype(std::declval<T>()*std::declval<U>())> operator/(const vector_k<T>& v1, const U& a)
	{
		vector_k<decltype(std::declval<T>()*std::declval<U>())> v2(v1.size());
		for(size_t ii = 0; ii < v2.size(); ii++)
		{
			v2[ii] = v1[ii]/a;
		}
		return v2;
	}
	
	template<>
	inline vector_k<double> operator/(const vector_k<double>& v1, const double& a)
	{
		vector_k<double> v2(v1);
		double aa = 1.0/a;
		int N = (int)v1.size(), incx = (int)v2.get_inc();
		dscal_(&N, &aa, &v2[0], &incx);
		return v2;
	}
	
	template<>
	inline vector_k<std::complex<double>> operator/(const vector_k<std::complex<double>>& v1, const std::complex<double>& a)
	{
		vector_k<std::complex<double>> v2(v1);
		std::complex<double> aa = 1.0/a;
		int N = (int)v1.size(), incx = (int)v2.get_inc();
		zscal_(&N, &aa, &v2[0], &incx);
		return v2;
	}
	
	template<class T, class U,size_t dim> matrix<decltype(std::declval<T>()*std::declval<U>()),dim> operator/(const U& a, const matrix<T,dim>& A)
	{
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> B;
		matrix<size_t,1> temp = size(A);
		B.resize(temp);
		
		for(size_t ii = 0; ii < A.mData.size(); ii++)
		{
			B.mData[ii] = a/A.mData[ii];
		}
		return B;
	}
	
	template<class T, class U> vector_k<decltype(std::declval<T>()*std::declval<U>())> operator/(const U& a, const vector_k<T>& v1)
	{
		vector_k<decltype(std::declval<T>()*std::declval<U>())> v2(v1.size());
		for(size_t ii = 0; ii < v2.size(); ii++)
		{
			v2[ii] = a/v1[ii];
		}
		return v2;
	}
	
	/** \brief Return a vector containing the sine of each element of v1.
	 */
	template<class T>
    vector_k<std::complex<T>> sin(const vector_k<std::complex<T>> & v1)
    {
        return eop(v1, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::sin<T>));
    }
	
	/** \brief Return a vector containing the sine of each element of v1.
	 */
	template<class T> vector_k<T> sin(const vector_k<T> &v1)
	{
		return eop(v1,static_cast<T (*)(T)>(&std::sin));
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
	
	/** \brief Return a vector containing the cosine of each element of v1.
	 */
	template<class T>
    vector_k<std::complex<T>> cos(const vector_k<std::complex<T>> & v1)
    {
        return eop(v1, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::cos<T>));
    }
	
	/** \brief Return a vector containing the cosine of each element of v1.
	 */
	template<class T> vector_k<T> cos(const vector_k<T> &v1)
	{
		return eop(v1,static_cast<T (*)(T)>(&std::cos));
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
	
	/** \brief Return a vector containing the tangent of each element of v1.
	 */
	template<class T>
    vector_k<std::complex<T>> tan(const vector_k<std::complex<T>> & v1)
    {
        return eop(v1, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::tan<T>));
    }
	
	/** \brief Return a vector containing the tangent of each element of v1.
	 */
	template<class T> vector_k<T> tan(const vector_k<T> &v1)
	{
		return eop(v1,static_cast<T (*)(T)>(&std::tan));
	}
	
	/** \brief Return a vector containing the tangent of each element of A.
	 */
	template<class T>
    vector_k<std::complex<T>> tan(const matrix<std::complex<T>> &A)
    {
        return eop(A, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::tan<T>));
    }
	
	/** \brief Return a vector containing the tangent of each element of A.
	 */
	template<class T> vector_k<T> tan(const matrix<T> &A)
	{
		return eop(A,static_cast<T (*)(T)>(&std::tan));
	}
	
	/** \brief Return a vector containing the arc cosine of each element of v1.
	 */
	template<class T>
    vector_k<std::complex<T>> acos(const vector_k<std::complex<T>> & v1)
    {
        return eop(v1, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::acos<T>));
    }
	
	/** \brief Return a vector containing the arc cosine of each element of v1.
	 */
	template<class T> vector_k<T> acos(const vector_k<T> &v1)
	{
		return eop(v1,static_cast<T (*)(T)>(&std::acos));
	}
	
	/** \brief Return a vector containing the arc cosine of each element of A.
	 */
	template<class T>
    vector_k<std::complex<T>> acos(const matrix<std::complex<T>> &A)
    {
        return eop(A, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::acos<T>));
    }
	
	/** \brief Return a vector containing the arc cosine of each element of A.
	 */
	template<class T> vector_k<T> acos(const matrix<T> &A)
	{
		return eop(A,static_cast<T (*)(T)>(&std::acos));
	}
	
	/** \brief Return a vector containing the arc sine of each element of v1.
	 */
	template<class T>
    vector_k<std::complex<T>> asin(const vector_k<std::complex<T>> & v1)
    {
        return eop(v1, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::asin<T>));
    }
	
	/** \brief Return a vector containing the arc sine of each element of v1.
	 */
	template<class T> vector_k<T> asin(const vector_k<T> &v1)
	{
		return eop(v1,static_cast<T (*)(T)>(&std::asin));
	}
	
	/** \brief Return a vector containing the arc sine of each element of A.
	 */
	template<class T>
    vector_k<std::complex<T>> asin(const matrix<std::complex<T>> &A)
    {
        return eop(A, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::asin<T>));
    }
	
	/** \brief Return a vector containing the arc sine of each element of A.
	 */
	template<class T> vector_k<T> asin(const matrix<T> &A)
	{
		return eop(A,static_cast<T (*)(T)>(&std::asin));
	}
	
	/** \brief Return a vector containing the exponential of each element of v1.
	 */
	template<class T>
    vector_k<std::complex<T>> exp(const vector_k<std::complex<T>> & v1)
    {
        return eop(v1, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::exp<T>));
    }
	
	/** \brief Return a vector containing the exponential of each element of v1.
	 */
	template<class T> vector_k<T> exp(const vector_k<T> &v1)
	{
		return eop(v1,static_cast<T (*)(T)>(&std::exp));
	}
	
	/** \brief Return a vector containing the exponential of each element of A.
	 */
	template<class T,size_t dim>
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
	
	/** \brief Return a vector containing the natural logarithm of each element of v1.
	 */
	template<class T>
    vector_k<std::complex<T>> log(const vector_k<std::complex<T>> & v1)
    {
        return eop(v1, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::log<T>));
    }
	
	/** \brief Return a vector containing the natural logarithm of each element of v1.
	 */
	template<class T> vector_k<T> log(const vector_k<T> &v1)
	{
		return eop(v1,static_cast<T (*)(T)>(&std::log));
	}
	
	/** \brief Return a vector containing the natural logarithm of each element of A.
	 */
	template<class T>
    vector_k<std::complex<T>> log(const matrix<std::complex<T>> &A)
    {
        return eop(A, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::log<T>));
    }
	
	/** \brief Return a vector containing the natural logarithm of each element of A.
	 */
	template<class T> vector_k<T> log(const matrix<T> &A)
	{
		return eop(A,static_cast<T (*)(T)>(&std::log));
	}
	
	/** \brief Return a vector containing the base 10 logarithm of each element of v1.
	 */
	template<class T>
    vector_k<std::complex<T>> log10(const vector_k<std::complex<T>> & v1)
    {
        return eop(v1, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::log10<T>));
    }
	
	/** \brief Return a vector containing the base 10 logarithm of each element of v1.
	 */
	template<class T> vector_k<T> log10(const vector_k<T> &v1)
	{
		return eop(v1,static_cast<T (*)(T)>(&std::log10));
	}
	
	/** \brief Return a vector containing the base 10 logarithm of each element of A.
	 */
	template<class T>
    vector_k<std::complex<T>> log10(const matrix<std::complex<T>> &A)
    {
        return eop(A, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::log10<T>));
    }
	
	/** \brief Return a vector containing the base 10 logarithm of each element of A.
	 */
	template<class T> vector_k<T> log10(const matrix<T> &A)
	{
		return eop(A,static_cast<T (*)(T)>(&std::log10));
	}
	
	/** \brief Return a vector containing the sqrt of each element of v1.
	 */
	template<class T>
    vector_k<std::complex<T>> sqrt(const vector_k<std::complex<T>> & v1)
    {
        return eop(v1, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::sqrt<T>));
    }
	
	/** \brief Return a vector containing the sqrt of each element of v1.
	 */
	template<class T> vector_k<T> sqrt(const vector_k<T> &v1)
	{
		return eop(v1,static_cast<T (*)(T)>(&std::sqrt));
	}
	
	/** \brief Return a vector containing the sqrt of each element of A.
	 */
	template<class T>
    vector_k<std::complex<T>> sqrt(const matrix<std::complex<T>> &A)
    {
        return eop(A, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&std::sqrt<T>));
    }
	
	/** \brief Return a vector containing the sqrt of each element of A.
	 */
	template<class T> vector_k<T> sqrt(const matrix<T> &A)
	{
		return eop(A,static_cast<T (*)(T)>(&std::sqrt));
	}
	
	/** \brief Return a vector containing the csqrt of each element of v1.
	 */
	template<class T>
    vector_k<std::complex<T>> csqrt(const vector_k<std::complex<T>> & v1)
    {
        return eop(v1, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&csqrt<T>));
    }
	
	/** \brief Return a vector containing the csqrt of each element of v1.
	 */
	template<class T> vector_k<T> csqrt(const vector_k<T> &v1)
	{
		return eop(v1,static_cast<T (*)(T)>(&csqrt));
	}
	
	/** \brief Return a vector containing the csqrt of each element of A.
	 */
	template<class T>
    vector_k<std::complex<T>> csqrt(const matrix<std::complex<T>> &A)
    {
        return eop(A, static_cast<std::complex<T> (*)(const std::complex<T> &)>(&csqrt<T>));
    }
	
	/** \brief Return a vector containing the csqrt of each element of A.
	 */
	template<class T> vector_k<T> csqrt(const matrix<T> &A)
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
	template<class T> int size(const matrix<T> &A, const int &dim)
	{
		return A.size(dim);
	}
	
    /**  \brief Returns the size of matrix A.
     *   @details
     *    Returns the size of matrix A in a matrix_size_type variable.
     *   @usage auto msize = size(A); int num_rows = msize.rows; int num_cols = msize.cols;
     *   @param[in] A matrix for which you want to know the size.
     *   @return A matrix_size_type variable containing the number of rows and cols.
     */
	template<class T, size_t dim> matrix<size_t,1> size(const matrix<T,dim> &A)
	{
	    matrix<size_t,1> msize(dim);
	    for(size_t ii = 0; ii < dim; ii++)
	    {
	        msize.mData[ii] = A.size(ii+1);
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
	
    /**  \brief Returns a vector of length N containing all zeros.
     *   @details
     *    Returns a vector of length N containing all zeros.
     *   @param[in] N Number of elements.
     *   @return A vector of length N containing zeros for each element. 
     */
	template<class T> vector_k<T> zeros(const int &N)
	{
		vector_k<T> v1(N);
		return v1;
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
	
    /**  \brief Returns a vector of length N containing all ones.
     *   @details
     *    Returns a vector of length N containing all ones.
     *   @param[in] N Number of elements.
     *   @return A vector of length N containing ones for each element. 
     */
	template<class T> vector_k<T> ones(const int &N)
	{
		vector_k<T> v1(N);
		for(int ii = 0; ii < N; ii++)
		{
		    v1[ii] = 1.0;
		}
		return v1;
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
	
	template<class T> matrix<T,2> diag(const matrix<T,1> &v1, const int &d = 0)
	{
		matrix<T,2> A(std::abs(d)+v1.size(1),std::abs(d)+v1.size(1));
		if(d != 0)
		{
		    for(int ii = 0; ii < v1.size(1); ii++)
		    {
		        if(d < 0)
		        {
			        A(ii+std::abs(d),ii) = v1.mData[ii];
			    }
			    else
		        {
			        A(ii,ii+std::abs(d)) = v1.mData[ii];
			    }
		    }
		}
		else
		{
		    for(int ii = 0; ii < v1.size(1); ii++)
		    {
			    A(ii,ii) = v1.mData[ii];
		    }
		}
		return A;
	}
	
	template<class T> vector_k<T> diag(const matrix<T> &A, const int &d = 0)
	{
	    if(A.empty())
	    {
	        throw KeyCppException("Cannot compute diagonal of empty matrix!");
	    }
	    int min_dim;
		vector_k<T> v1;
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
		    v1 = vector_k<T>(min_dim);
		    for(int ii = 0; ii < min_dim; ii++)
		    {
			    v1[ii] = A(ii,ii);
		    }
		}
		else
		{
		    if(d > 0)
		    {
		        min_dim = A.size(2) - std::abs(d);
		        v1 = vector_k<T>(min_dim);
		        for(int ii = 0; ii < min_dim; ii++)
		        {
			        v1[ii] = A(ii,ii+std::abs(d));
		        }
		    }
		    else
		    {
		        min_dim = A.size(1) - std::abs(d);
		        v1 = vector_k<T>(min_dim);
		        for(int ii = 0; ii < min_dim; ii++)
		        {
			        v1[ii] = A(ii+std::abs(d),ii);
		        }
		    }
		}
		return v1;
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
	
	template<class T> matrix<T,2> repmat(const vector_k<T> v1, const int &m, const int &n)
	{
		matrix<T,2> B(m, n*v1.size());
		for(int ii = 0; ii < m; ii++)
		{
			for(int jj = 0; jj < n; jj++)
			{
				for(int kk = 0; kk < v1.size(); kk++)
				{
					B(ii*v1.size() + kk, jj) = v1[kk];
				}
			}
		}
		return B;
	}
	
	template<class T> matrix<T,2> repmat(const matrix<T,1> v1, const int &m, const int &n)
	{
		matrix<T,2> B(m, n*v1.size(1));
		for(int ii = 0; ii < m; ii++)
		{
			for(int jj = 0; jj < n; jj++)
			{
				for(int kk = 0; kk < v1.size(1); kk++)
				{
					B(ii,jj*v1.size(1) + kk) = v1(kk);
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
		for(int ii = 0; ii < A.mData.size(); ii++)
		{
		    C.mData[ii] = A.mData[ii]*B.mData[ii];
		}
		
		return C;
	}
	
	/** \brief Performs array multiplication on vectors v1 and v2.
	 *
	 *  Each element of v1 is multiplied by each element of v2. The vector that is
	 *  returned is the same size as v1 and v2.
	 */
	template<class T, class U> vector_k<decltype(std::declval<T>()*std::declval<U>())> times(const vector_k<T>& v1, const vector_k<U>& v2)
	{
	    if(v1.empty() || v2.empty())
	    {
	        throw KeyCppException("Cannot multiply an empty vector!");
	    }
	    if(v1.size() != v2.size())
	    {
	        throw KeyCppException("Vector dimensions must agree in times().");
	    }
		vector_k<decltype(std::declval<T>()*std::declval<U>())> v3(v1.size());
		for(int ii = 0; ii < v1.size(); ii++)
		{
			v3[ii] = v1[ii]*v2[ii];
		}
		return v3;
	}
	
	/** \brief Performs right array division on matrices A and B.
	 *
	 *  Each element of A is divided by each element of B. The matrix that is
	 *  returned is the same size as A and B. Equivalent to A./B in MATLAB.
	 */
	template<class T, class U> matrix<decltype(std::declval<T>()*std::declval<U>())> rdivide(const matrix<T>& A, const matrix<U>& B)
	{
	    if(A.size(1) <= 0 || A.size(2) <= 0 || B.size(1) <= 0 || B.size(2) <= 0)
	    {
	        throw KeyCppException("Cannot rdivide an empty matrix!");
	    }
	    if(A.size(1) != B.size(1) || A.size(2) != B.size(2))
	    {
	        throw KeyCppException("Matrix dimensions must agree in rdivide().");
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
	
	/** \brief Performs right array division on vectors v1 and v2.
	 *
	 *  Each element of v1 is divided by each element of v2. The vector that is
	 *  returned is the same size as v1 and v2. Equivalent to v1./v2 in MATLAB.
	 */
	template<class T, class U> vector_k<decltype(std::declval<T>()*std::declval<U>())> rdivide(const vector_k<T>& v1, const vector_k<U>& v2)
	{
	    if(v1.empty() || v2.empty())
	    {
	        throw KeyCppException("Cannot divide an empty vector!");
	    }
	    if(v1.size() != v2.size())
	    {
	        throw KeyCppException("Vector dimensions must agree in rdivide().");
	    }
		vector_k<decltype(std::declval<T>()*std::declval<U>())> v3(v1.size());
		for(int ii = 0; ii < v1.size(); ii++)
		{
			v3[ii] = v1[ii]/v2[ii];
		}
		return v3;
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
	
	/** \brief Performs left array division on vectors v2 and v1.
	 *
	 *  Each element of v1 is divided by each element of v2. The vector that is
	 *  returned is the same size as v2 and v1. Equivalent to v2.\v1 in MATLAB.
	 */
	template<class T, class U> vector_k<decltype(std::declval<T>()*std::declval<U>())> ldivide(const vector_k<T>& v2, const vector_k<U>& v1)
	{
	    if(v1.empty() || v2.empty())
	    {
	        throw KeyCppException("Cannot divide an empty vector!");
	    }
	    if(v1.size() != v2.size())
	    {
	        throw KeyCppException("Vector dimensions must agree in ldivide().");
	    }
		vector_k<decltype(std::declval<T>()*std::declval<U>())> v3(v1.size());
		for(int ii = 0; ii < v1.size(); ii++)
		{
			v3[ii] = v1[ii]/v2[ii];
		}
		return v3;
	}

	template<class T> int sign(const T &val)
	{
	    return (T(0) < val) - (val < T(0));
	}
	
	template<class T> T angle(const std::complex<T> &x)
	{
		return arg(x);
	}
	
	template<class T>
    vector_k<T> angle(const vector_k<std::complex<T>> &v1)
    {
        return eop(v1, static_cast<T (*)(const std::complex<T> &)>(&std::arg<T>));
    }
	
	template<class T,size_t dim>
    matrix<T,dim> angle(const matrix<std::complex<T>,dim> &A)
    {
        return eop(A, static_cast<T (*)(const std::complex<T> &)>(&std::arg<T>));
    }

	template<class T> T max(const vector_k<T> &x)
	{
		double a = nan("");
		size_t index = 0;
		for(size_t ii = 0; ii < x.size(); ii++)
		{
			if(!std::isnan(x[ii]) && (x[ii] > a || std::isnan(a)))
			{
				a = x[ii];
				index = ii;
			}
		}
		return x[index];
	}

	inline std::complex<double> max(const vector_k<std::complex<double> > &x)
	{
		double a = nan("");
		double b = nan("");
		size_t index = 0;
		for(size_t ii = 0; ii < x.size(); ii++)
		{
			if(!std::isnan(real(x[ii])) && !std::isnan(imag(x[ii])) && ((abs(x[ii]) > a && angle(x[ii]) > b) || (std::isnan(a) || std::isnan(b))))
			{
				a = abs(x[ii]);
				b = angle(x[ii]);
				index = ii;
			}
		}
		return x[index];
	}
	
	template<class T,size_t dim> matrix<T,dim-1> max(const matrix<T,dim> &A)
	{
	    static_assert(dim == 2,"This function is only available for matrices of dimension 2.");
	    matrix<T,dim-1> v(A.size(2));
	    
	    for(size_t jj = 0; jj < A.size(2); jj++)
	    {
	        vector_k<T> v2(A.size(1));
	        for(size_t ii = 0; ii < A.size(1); ii++)
	        {
	            v2[ii] = A(ii,jj);
	        }
	        v(jj) = max(v2);
	    }
	    return v;
	}
	
	template<class T> T max(const matrix<T,1> &A)
	{
	    T v;
	    
	    vector_k<T> v2(A.size(1));
        for(size_t ii = 0; ii < A.size(1); ii++)
        {
            v2[ii] = A.mData[ii];
        }
        v = max(v2);
	    return v;
	}

	template<class T> T min(const matrix<T,1> &x)
	{
		double a = nan("");
		int index = 0;
		for(int ii = 0; ii < x.size(1); ii++)
		{
			if(!std::isnan(x(ii)) && (x(ii) < a || std::isnan(a)))
			{
				a = x(ii);
				index = ii;
			}
		}
		return x(index);
	}

	inline std::complex<double> min(const vector_k<std::complex<double> > &x)
	{
		double a = nan("");
		double b = nan("");
		size_t index = 0;
		for(size_t ii = 0; ii < x.size(); ii++)
		{
			if(!std::isnan(real(x[ii])) && !std::isnan(imag(x[ii])) && ((abs(x[ii]) < a && angle(x[ii]) < b) || (std::isnan(a) || std::isnan(b))))
			{
				a = abs(x[ii]);
				b = angle(x[ii]);
				index = ii;
			}
		}
		return x[index];
	}
	
	
	template<class T> matrix<T,1> min(const matrix<T,2> &A)
	{
	    matrix<T,1> v(A.size(2));
	    for(size_t ii = 0; ii < v.size(1); ii++)
	    {
	        v(ii) = min(A.col(ii));
	    }
	    return v;
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
	
	template<class T> matrix<T,2> transpose(const matrix<T,1> &v1)
	{
		matrix<T,2> B(v1.size(1),1);
		for(size_t ii = 0; ii < v1.size(1); ii++)
		{
			B(ii,0) = v1.mData[ii];
		}
		return B;
	}
	
	/** \brief Returns the transpose of vector v1.
	 */
	template<class T> matrix<T> transpose(const vector_k<T> &v1)
	{
		matrix<T> B(1,v1.size(1));
		for(size_t ii = 0; ii < v1.size(1); ii++)
		{
			B(0,ii) = v1[ii];
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
	
	template<class T> matrix<T,2> ctranspose(const matrix<T,1> &v1)
	{
		matrix<T,2> B(v1.size(1),1);
		for(size_t ii = 0; ii < v1.size(1); ii++)
		{
			B(ii,0) = conj(v1.mData[ii]);
		}
		return B;
	}
	
	template<> inline matrix<double,2> ctranspose(const matrix<double,1> &v1)
	{
		return transpose(v1);
	}
	
	/** \brief Returns the complex-conjugate transpose of vector v1.
	 */
	template<class T> matrix<T> ctranspose(const vector_k<T> &v1)
	{
		matrix<T> B(1,v1.size(1));
		for(size_t ii = 0; ii < v1.size(1); ii++)
		{
			B(0,ii) = conj(v1[ii]);
		}
		return B;
	}
	
	/** \brief Computes the sum of vector v1.
	 */
	template<class T> T sum(const vector_k<T> &v1)
	{
		T a = 0.0;
		for(size_t ii = 0; ii < v1.size(); ii++)
		{
			a += v1[ii];
		}
		return a;
	}

	template<class T> T sum(const matrix<T,1> &v1)
	{
		T a = 0.0;
		for(size_t ii = 0; ii < v1.size(1); ii++)
		{
			a += v1.mData[ii];
		}
		return a;
	}
	
	/** \brief Computes the sum of each column of A.
	 */
	template<class T> matrix<T,1> sum(const matrix<T,2> &A)
	{
		matrix<T,1> v1(A.size(2));
		for(size_t ii = 0; ii < v1.size(1); ii++)
		{
			v1.mData[ii] = sum(A.col(ii));
		}
		return v1;
	}
	
	/** \brief Converts matrix A to a column vector.
	 */
	template<class T> vector_k<T> mat2vec(const matrix<T> &A)
	{
	    if(A.empty())
	    {
	        throw KeyCppException("Cannot convert empty matrix to a vector!");
	    }
		vector_k<T> v1(A.size(1));
		for(size_t ii = 0; ii < v1.size(); ii++)
		{
			v1[ii] = A(ii,0);
		}
		return v1;
	}
	
	/** \brief Converts a column vector to a 1 x length(v1) matrix.
	 */
	template<class T> matrix<T> vec2mat(const vector_k<T> &v1)
	{
	    if(v1.empty())
	    {
	        throw KeyCppException("Cannot convert empty vector to a matrix!");
	    }
		matrix<T> A(v1.size(),1);
		for(size_t ii = 0; ii < v1.size(); ii++)
		{
			A(ii,0) = v1[ii];
		}
		return A;
	}
	
	template<class T> matrix<T,1> linspace(const T &x1, const T &x2, const size_t &N)
	{
		matrix<T,1> x(N);
		if(N == 1)
		{
			x.mData[0] = x2;
			return x;
		}

		T delta_x = (x2-x1)/(N-1);

		for(size_t ii = 0; ii < N; ii++)
		{
			x.mData[ii] = x1 + ii*delta_x;
		}

		x.mData[N-1] = x2;

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
	template<class T> vector_k<T> logspace(const T &x1, const T &x2, const int &N)
	{
		vector_k<T> x(N);
		if(N == 1)
		{
			x[0] = x2;
			return x;
		}

		T delta_x = (x2-x1)/(N-1);

		for(size_t ii = 0; ii < N; ii++)
		{
			x[ii] = pow(10.0,(x1 + delta_x*ii));
		}

		return x;
	}
	
	template<class T> matrix<T,1> unwrap(const matrix<T,1>& v1, const T &tol = pi)
	{
		matrix<T,1> v2(v1.size(1));
		v2(0) = v1(0);
		int correction = 0;
		for(size_t ii = 1; ii < v1.size(1); ii++)
		{
			if((v1.mData[ii] - v1.mData[ii-1]) > tol)
			{
				correction -= 1;
			}
			else if((v1.mData[ii] - v1.mData[ii-1]) < -tol)
			{
				correction += 1;
			}
			v2.mData[ii] = v1.mData[ii] + correction*2*pi;
		}
		return v2;
	}
	
	/** \brief Computes the mean of vector v1.
	 */
	template<class T> T mean(const vector_k<T>& v1)
	{
		T m = T(0);
		double tot = 0.0;
		for(size_t ii = 0; ii < v1.size(); ii++)
		{
			m += v1[ii];
			tot += 1.0;
		}
		return m/tot;
	}
	
	template<class T, class U> T interp1(const matrix<U,1> &x, const matrix<T,1> &y, const U &x_interp, std::string method = "linear", Extrap extrap = Extrap())
	{
		if(x.empty() || y.empty())
		{
			throw KeyCppException("Error in interp1! Cannot interpolate on an empty vector!");
		}
		if(x.size(1) != y.size(1))
		{
			throw KeyCppException("Error in interp1! Variables `x` and `y` have incompatible sizes!");
		}
		std::transform(method.begin(), method.end(), method.begin(), ::tolower);
		size_t N = x.size(1);
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
				if(x_interp == x.mData[ii])
				{
					return y.mData[ii];
				}
				else if((x_interp > x.mData[ii] && x_interp < x.mData[ii+1]) || (x_interp < x.mData[ii] && x_interp > x.mData[ii+1]))
				{
					return (y.mData[ii] + (x_interp - x.mData[ii])*(y.mData[ii+1] - y.mData[ii])/(x.mData[ii+1] - x.mData[ii]));
				}
			}

			if(x_interp == x.mData[N-1])
			{
				return y.mData[N-1];
			}
			
			if(extrap.isString)
			{
				if(extrap.extrap_string.compare("extrap") == 0)
				{
					if(x_interp < x.mData[0])
					{
						return (y.mData[0] + (x_interp - x.mData[0])*(y.mData[1] - y.mData[0])/(x.mData[1] - x.mData[0]));
					}
					else
					{
						return (y.mData[N-2] + (x_interp - x.mData[N-2])*(y.mData[N-1] - y.mData[N-2])/(x.mData[N-1] - x.mData[N-2]));
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
				if(std::abs(x.mData[ii] - x_interp) < std::abs(min_val))
				{
					min_val = x.mData[ii] - x_interp;
					index = ii;
				}
				else if(std::abs(x.mData[ii] - x_interp) == std::abs(min_val) && (x.mData[ii] - x_interp) > min_val)
				{
					min_val = x.mData[ii] - x_interp;
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


	template<class T, class U> matrix<T,1> interp1(const matrix<U,1> &x, const matrix<T,1> &y, const matrix<U,1> &x_interp, std::string method = "linear", Extrap extrap = Extrap())
	{
		if(x.empty() || y.empty() || x_interp.empty())
		{
			throw KeyCppException("Error in interp1! Cannot interpolate on an empty vector!");
		}
		if(x.size(1) != y.size(1))
		{
			throw KeyCppException("Error in interp1! Variables `x` and `y` have incompatible sizes!");
		}
		std::transform(method.begin(), method.end(), method.begin(), ::tolower);
		size_t N = x.size(1);
		size_t N_int = x_interp.size(1);
		matrix<T,1> y2(N_int);

		if(method.compare("spline") == 0)
		{	
			Spline<U,T> spline(N,x,y,extrap);
			spline.compute_spline();
			for(size_t ii = 0; ii < N_int; ii++)
			{
				y2.mData[ii] = spline.J(x_interp.mData[ii]);
			}
		}
		else if(method.compare("linear") == 0 || method.compare("nearest") == 0)
		{
			for(size_t ii = 0; ii < N_int; ii++)
			{
				y2.mData[ii] = interp1(x, y, x_interp.mData[ii],method,extrap);
			}
		}
		else
		{
			throw KeyCppException("Error in interp1! Unrecognized interpolation method!");
		}

		return y2;
	}

	template<class T, class U> matrix<T,2> interp1(const matrix<U,1> &x, const matrix<T,2> &y, const matrix<U,1> &x_interp, std::string method = "linear", Extrap extrap = Extrap())
	{
		if(x.empty() || y.size(1) <= 0 || y.size(2) <= 0 || x_interp.empty())
		{
			throw KeyCppException("Error in interp1! Cannot interpolate on an empty vector or matrix!");
		}
		if(x.size(1) != y.size(1))
		{
			throw KeyCppException("Error in interp1! Variables `x` and `y` have incompatible sizes!");
		}
		matrix<T,2> y2(x_interp.size(1),y.size(2));

		for(size_t kk = 0; kk < y.size(2); kk++)
		{
			y2.col(kk) = interp1(x,y.col(kk),x_interp,method, extrap);
		}

		return y2;
	}
	
	template<class T, class U> matrix<T,2> interp1(const matrix<U,1> &x, const matrix<T,1> &y, const matrix<U,2> &x_interp, std::string method = "linear", Extrap extrap = Extrap())
	{
		if(x.empty() || y.empty() || x_interp.size(1) <= 0 || x_interp.size(2) <= 0)
		{
			throw KeyCppException("Error in interp1! Cannot interpolate on an empty vector or matrix!");
		}
		if(x.size(1) != y.size(1))
		{
			throw KeyCppException("Error in interp1! Variables `x` and `y` have incompatible sizes!");
		}
		matrix<T,2> y2(x_interp.size(1),x_interp.size(2));

		for(size_t kk = 0; kk < x_interp.size(2); kk++)
		{
			y2.col(kk) = interp1(x,y,x_interp.col(kk),method, extrap);
		}

		return y2;
	}	

	template<class U, class T> T trapz(const matrix<U,1> &eta, const matrix<T,1> &integrand)
	{
		if(eta.empty() || integrand.empty())
		{
			throw KeyCppException("Error in trapz()! Empty vector supplied!");
		}
		if(eta.size(1) != integrand.size(1))
		{
			throw KeyCppException("Error in trapz()! Vector sizes are not compatible!");
		}
		size_t N = eta.size(1);
		T sum = 0.0;

		for(size_t ii = 0; ii < N-1; ii++)
		{
			sum += (eta.mData[ii+1] - eta.mData[ii])*(integrand.mData[ii+1] + integrand.mData[ii]);
		}
		return 0.5*sum;
	}

	template<class U, class T> matrix<T,1> trapz(const matrix<U,1> &eta, const matrix<T,2> &integrand)
	{
		if(eta.empty() || integrand.empty())
		{
			throw KeyCppException("Error in trapz()! Empty vector/matrix supplied!");
		}
		if(eta.size(1) != integrand.size(1) && integrand.size(1) > 1)
		{
			throw KeyCppException("Error in trapz()! Vector and matrix sizes are not compatible!");
		}
		if(eta.size(1) != integrand.size(2) && integrand.size(1) <= 1)
		{
			throw KeyCppException("Error in trapz()! Vector and matrix sizes are not compatible!");
		}
		
		size_t N;
		matrix<T,1> z;
		if(eta.size(1) == integrand.size(1))
		{
		    N = integrand.size(2);
		    z = matrix<T,1>(N);
		    for(size_t ii = 0; ii < N; ii++)
		    {
			    z.mData[ii] = trapz(eta,integrand.col(ii));
		    }
		}
		else
		{
		    N = 1;
		    z = matrix<T,1>(N);
			z(0) = trapz(eta,integrand.row(0));
		}
		return z;
	}


	template<class T, class U> matrix<T,2> diffxy(const matrix<U,2> &eta, const matrix<T,2> &u, const int &index = 2)
	{
		if(eta.size(1) <= 0 || eta.size(2) <= 0 || u.size(1) <= 0 || u.size(2) <= 0)
		{
			throw KeyCppException("Error in diffxy()! Empty matrix supplied!");
		}
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
				for(size_t n = 0; n < N-1; n++)
				{
					du(n,p) = (u(n+1,p) - u(n,p))/(eta(n+1,p) - eta(n,p));
				}
				du(N-1,p) = (u(N-1,p) - u(N-2,p))/(eta(N-1,p) - eta(N-2,p));
			}
		}
		else
		{
			for(size_t n = 0; n < N; n++)
			{
				for(size_t p = 0; p < P-1; p++)
				{
					du(n,p) = (u(n,p+1) - u(n,p))/(eta(n,p+1) - eta(n,p));
				}
				du(n,P-1) = (u(n,P-1) - u(n,P-2))/(eta(n,P-1) - eta(n,P-2));
			}
		}

		return du;
	}

	template<class T, class U> matrix<T,2> diffxy(const matrix<U,1> &eta, const matrix<T,2> &u)
	{
		if(eta.empty() || u.size(1) <= 0 || u.size(2) <= 0)
		{
			throw KeyCppException("Error in diffxy()! Empty vector or matrix supplied!");
		}
		if(eta.size(1) != u.size(1) && eta.size(1) != u.size(2))
		{
			throw KeyCppException("Error in diffxy()! Vector and matrix sizes are not compatible!");
		}
		size_t N = u.size(1);
		size_t P = u.size(2);

		matrix<T,2> du(N,P);
		if(N == eta.size(1))
		{
			for(size_t p = 0; p < P; p++)
			{
				for(size_t ii = 0; ii < N-1; ii++)
				{
					du(ii,p) = (u(ii+1,p) - u(ii,p))/(eta.mData[ii+1] - eta.mData[ii]);
				}
				du(N-1,p) = (u(N-1,p) - u(N-2,p))/(eta.mData[N-1] - eta.mData[N-2]);
			}
		}
		else
		{
			for(size_t ii = 0; ii < N; ii++)
			{
				for(size_t p = 0; p < P-1; p++)
				{
					du(ii,p) = (u(ii,p+1) - u(ii,p))/(eta(p+1) - eta(p));
				}
				du(ii,P-1) = (u(ii,P-1) - u(ii,P-2))/(eta(P-1) - eta(P-2));
			}
		}

		return du;
	}

	template<class T, class U> matrix<T,1> diffxy(const matrix<U,1> &eta, const matrix<T,1> &u)
	{
		if(eta.empty() || u.empty())
		{
			throw KeyCppException("Error in diffxy()! Empty vector supplied!");
		}
		if(eta.size(1) != u.size(1))
		{
			throw KeyCppException("Error in diffxy()! Vector sizes are not compatible!");
		}
		size_t N = u.size(1);

		matrix<T,1> du(N);
		for(size_t ii = 0; ii < N-1; ii++)
		{
			du.mData[ii] = (u.mData[ii+1] - u.mData[ii])/(eta.mData[ii+1] - eta.mData[ii]);
		}
		du.mData[N-1] = (u.mData[N-1] - u.mData[N-2])/(eta.mData[N-1] - eta.mData[N-2]);

		return du;
	}

	template<class T> matrix<std::complex<double>,1> fft(const matrix<T,1> &u, int N = -1)
	{
		if(u.empty())
		{
			throw KeyCppException("Error in fft()! Empty vector supplied!");
		}
		
		if(N < 0)
		{
			N = u.size(1);
		}
		
		kiss_fft_cpx *cx_in = new kiss_fft_cpx[N];
		kiss_fft_cpx *cx_out = new kiss_fft_cpx[N];

		matrix<std::complex<double>,1> u_hat(N);

		for(int ii = 0; ii < N; ii++)
		{
			cx_in[ii].r = real((std::complex<double>)u.mData[ii]);
			cx_in[ii].i = imag((std::complex<double>)u.mData[ii]);
		}

		kiss_fft_cfg cfg = kiss_fft_alloc(N,false,NULL,NULL);
		kiss_fft(cfg,cx_in,cx_out);

		for(int ii = 0; ii < N; ii++)
		{
			u_hat.mData[ii] = std::complex<T>((T)cx_out[ii].r,(T)cx_out[ii].i);
		}

		free(cfg);
		delete [] cx_in;
		delete [] cx_out;

		return u_hat;
	}
}
	
namespace boost
{
    namespace numeric
    {
        namespace odeint
        {
            template<class T>
            struct is_resizeable<keycpp::matrix<T,1> >
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
		matrix<T,1> t;
		matrix<Y,2> y;
	};

	template<class T, class U, class F>
	ODE_type<U,T> ode45(F odeClass, const std::initializer_list<U> &x_span, matrix<T,1> ICs, double abs_tol = 1.0e-10, double rel_tol = 1.0e-6)
	{
		if(x_span.size() <= 0)
		{
			throw KeyCppException("Error in ode45()! x_span cannot be empty!");
		}
		if(ICs.empty())
		{
			throw KeyCppException("Error in ode45()! Must provide initial conditions!");
		}
		if(x_span.size() != 2)
		{
			throw KeyCppException("Error in ode45()! Invalid vector x_span!");
		}

		U x0 = *(x_span.begin());
		U xf = *(x_span.end()-1);
		U delta_x0 = (xf-x0)/1000.0;
		matrix<T,2> y_temp(2,ICs.size(1));
		matrix<U,1> x_temp(2);

		{
			using namespace boost::numeric::odeint;
			integrate_adaptive(make_controlled<runge_kutta_dopri5<matrix<T,1> > >(abs_tol, rel_tol), odeClass, ICs, x0, xf, delta_x0, observe<T,U>(y_temp,x_temp));
		}

        ODE_type<U,T> ans;
        ans.t = x_temp;
        ans.y = y_temp;
		
		return ans;
	}
	
	template<class T, class U, class F>
	matrix<T> ode45(F odeClass, matrix<U,1> x_ode, matrix<T,1> ICs, double abs_tol = 1.0e-10, double rel_tol = 1.0e-6)
	{
		if(x_ode.empty())
		{
			throw KeyCppException("Error in ode45()! Vector x_ode cannot be empty!");
		}
		if(ICs.empty())
		{
			throw KeyCppException("Error in ode45()! Must provide initial conditions!");
		}
		if(x_ode.size(1) < 2)
		{
			throw KeyCppException("Error in ode45()! Invalid vector x_ode!");
		}

		U delta_x0 = x_ode(1) - x_ode(0);
		matrix<T,2> y_temp(x_ode.size(1),ICs.size(1));
		matrix<U,1> x_temp(x_ode.size(1));

		{
			using namespace boost::numeric::odeint;
			integrate_times(make_dense_output<runge_kutta_dopri5<matrix<T,1> > >(abs_tol, rel_tol), odeClass, ICs, x_ode.begin(), x_ode.end(), delta_x0, observe<T,U>(y_temp,x_temp));
		}

		matrix<T> y(y_temp);
		
		return y;
	}
	
	inline void set(Figure &h, std::string property, double val)
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
	
	template<class T> Sort_Matrix<T> sort(const matrix<T> &A, const size_t &dim = 2, std::string method = "ascend")
	{
		std::transform(method.begin(), method.end(), method.begin(), ::tolower);
		if(method.compare("ascend") != 0 && method.compare("descend") != 0)
		{
			throw KeyCppException("Invalid sort method!");
		}
		if(A.size(1) <= 0 || A.size(2) <= 0)
		{
			throw KeyCppException("Tried to sort empty matrix!");
		}
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
	
	template<class T>
	struct Sort_Vector
	{
		vector_k<T> v;
		vector_k<size_t> index;
	};
	
	template<class T> Sort_Vector<T> sort(const vector_k<T> &v1, std::string method = "ascend")
	{
		std::transform(method.begin(), method.end(), method.begin(), ::tolower);
		if(method.compare("ascend") != 0 && method.compare("descend") != 0)
		{
			throw KeyCppException("Invalid sort method!");
		}
		if(v1.empty())
		{
			throw KeyCppException("Tried to sort empty vector!");
		}
		bool swapped = true;
		T temp;
		size_t temp_i;
		vector_k<size_t> index(v1.size());
		for(size_t ii = 0; ii < v1.size(); ii++)
		{
			index[ii] = ii;
		}
		vector_k<T> v2(v1.size());
		for(size_t ii = 0; ii < v2.size(); ii++)
		{
			v2[ii] = v1[ii];
		}
		swapped = true;
		while(swapped)
		{     
			swapped = false;
			for(size_t ii = 1; ii < v1.size(); ii++)
			{
				if(((v2[ii-1]) > (v2[ii]) && method.compare("ascend") == 0) || ((v2[ii-1]) < (v2[ii]) && method.compare("descend") == 0))
				{
					temp = v2[ii-1];
					v2[ii-1] = v2[ii];
					v2[ii] = temp;
					temp_i = index[ii-1];
					index[ii-1] = index[ii];
					index[ii] = temp_i;
					swapped = true;
				}
			}
		}
		
		Sort_Vector<T> sort_vector;
		sort_vector.v = v2;
		sort_vector.index = index;

		return sort_vector;
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
	
	/** \brief Returns the number of elements in a vector.
	 */
	template<class T>
	size_t length(const vector_k<T> &v1)
	{
	    return v1.size();
	}
	
	/** \brief Returns the length of the largest dimension of A.
	 */
	template<class T>
	size_t length(const matrix<T> &A)
	{
	    size_t m = A.size(1);
	    size_t n = A.size(2);
	    return ((m > n)?(m):(n));
	}
	
	template<class T>
	size_t numel(const matrix<T> &A)
	{
	    return (A.size(1)*A.size(2));
	}
	
	/** \brief Finds and returns the indices of non-zero elements of v1.
	 */
	template<class T>
	vector_k<size_t> find(const vector_k<T> &v1, const size_t &k = -1, std::string start = "")
	{
		std::transform(start.begin(), start.end(), start.begin(), ::tolower);
	    vector_k<size_t> v2;
	    if(v1.empty())
	    {
	        return v2;
	    }
	    if(k < 0 || k > v1.size() || start.empty() || (start.compare("first") != 0 &&
	       start.compare("last") != 0))
	    {
	        v2.reserve(v1.size());
	        for(size_t ii = 0; ii < v1.size(); ii++)
	        {
	            if(std::abs(v1[ii]) > eps)
	            {
	                v2.push_back(ii);
	            }
	        }
	    }
	    else if(start.compare("first") == 0 && k > 0 && k < v1.size())
	    {
	        v2.reserve(k);
	        size_t count = 0;
	        for(size_t ii = 0; ii < v1.size(); ii++)
	        {
	            if(std::abs(v1[ii]) > eps)
	            {
	                v2.push_back(ii);
	                count++;
	                if(count >= k)
	                {
	                    break;
	                }
	            }
	        }
	    }
	    else if(start.compare("last") == 0 && k > 0 && k < v1.size())
	    {
	        v2.reserve(k);
	        size_t count = 0;
	        for(size_t ii = v1.size()-1; ii >= 0; ii--)
	        {
	            if(std::abs(v1[ii]) > eps)
	            {
	                v2.push_back(ii);
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
	struct matrix_find_type
	{
		vector_k<size_t> rows;
		vector_k<size_t> cols;
		vector_k<T> v;
	};
	
	/** \brief Finds and returns the row and column indices and values of non-zero elements of A.
	 */
	template<class T>
	matrix_find_type<T> find(const matrix<T> &A)
	{
	    matrix_find_type<T> out;
	    if(A.empty())
	    {
	        return out;
	    }
	  
        out.rows.reserve(numel(A));
        out.cols.reserve(numel(A));
        out.v.reserve(numel(A));
        for(size_t ii = 0; ii < A.size(1); ii++)
        {
            for(size_t jj = 0; jj < A.size(2); jj++)
            {
                if(std::abs(A(ii,jj)) > eps)
                {
                    out.rows.push_back(ii);
                    out.cols.push_back(jj);
                    out.v.push_back(A(ii,jj));
                }
            }
        }
	    
	    return out;
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
	
	template<class T>
	matrix<T> reshape(const vector_k<T> &v1, const size_t &m, const size_t &n)
	{
	    return reshape(vec2mat(v1),m,n);
	}
	
	/** \brief Computes the dot product between vectors v1 and v2.
	 */
	template<class T, class U>
	decltype(std::declval<T>()*std::declval<U>()) dot(const vector_k<T> &v1, const vector_k<U> &v2)
	{
	    if(v1.empty() || v2.empty())
	    {
	        throw KeyCppException("Cannot dot multiply an empty vector!");
	    }
	    if(v1.size() != v2.size())
	    {
	        throw KeyCppException("Vectors must be same size in dot()!");
	    }
	    decltype(std::declval<T>()*std::declval<U>()) result = 0.0;
	    for(size_t ii = 0; ii < v1.size(); ii++)
	    {
	        result += v1[ii]*v2[ii];
	    }
	    return result;
	}
	
	/** \brief Computes the dot product between vectors v1 and v2.
	 */
	template<>
	inline double dot(const vector_k<double> &v1, const vector_k<double> &v2)
	{
	    if(v1.empty() || v2.empty())
	    {
	        throw KeyCppException("Cannot dot multiply an empty vector!");
	    }
	    if(v1.size() != v2.size())
	    {
	        throw KeyCppException("Vectors must be same size in dot()!");
	    }
	    int inca = 1, incb = 1, N = (int)v1.size();
	    return ddot_(&N, &v1[0], &inca, &v2[0], &incb);
	}
	
	/** \brief Computes the dot product between vectors v1 and v2.
	 */
	template<>
	inline std::complex<double> dot(const vector_k<std::complex<double>> &v1, const vector_k<std::complex<double>> &v2)
	{
	    if(v1.empty() || v2.empty())
	    {
	        throw KeyCppException("Cannot dot multiply an empty vector!");
	    }
	    if(v1.size() != v2.size())
	    {
	        throw KeyCppException("Vectors must be same size in dot()!");
	    }
	    int inca = 1, incb = 1, N = (int)v1.size();
	    std::complex<double> result;
	    zdotu_(&result, &N, &v1[0], &inca, &v2[0], &incb);
	    return result;
	}
	
	/** \brief Computes the dot product between the first non-singleton dimension of A and B.
	 */
	template<class T, class U>
	vector_k<decltype(std::declval<T>()*std::declval<U>())> dot(const matrix<T> &A, const matrix<U> &B, const size_t &dim = -1)
	{
	    if(A.empty() || B.empty())
	    {
	        throw KeyCppException("Cannot dot multiply empty matrices!");
	    }
	    if(A.size(1) != B.size(1) && A.size(2) != B.size(2))
	    {
	        throw KeyCppException("Matrices must be same size in dot()!");
	    }
	    vector_k<decltype(std::declval<T>()*std::declval<U>())> result;
	    if((A.size(1) > 1 || dim == 1) && dim != 2)
	    {
	        result = vector_k<decltype(std::declval<T>()*std::declval<U>())>(A.size(2));
	        for(size_t ii = 0; ii < result.size(); ii++)
	        {
	            result[ii] = dot(A.col(ii),B.col(ii));
	        }
	    }
	    else
	    {
	        result = vector_k<decltype(std::declval<T>()*std::declval<U>())>(A.size(1));
	        for(size_t ii = 0; ii < result.size(); ii++)
	        {
	            result[ii] = dot(A.row(ii),B.row(ii));
	        }
	    }
	    return result;
	}
	
	/** \brief Computes the cross product between vectors v1 and v2. Both vectors
	 *         must have exactly 3 elements.
	 */
	template<class T, class U>
	vector_k<decltype(std::declval<T>()*std::declval<U>())> cross(const vector_k<T> &v1, const vector_k<U> &v2)
	{
	    if(v1.empty() || v2.empty())
	    {
	        throw KeyCppException("Cannot cross multiply an empty vector!");
	    }
	    if(v1.size() != 3 || v2.size() != 3)
	    {
	        throw KeyCppException("Vectors must be have length of 3 in cross()!");
	    }
	    vector_k<decltype(std::declval<T>()*std::declval<U>())> result(3);
	    result[0] = v1[1]*v2[2] - v1[2]*v2[1];
	    result[1] = v1[2]*v2[0] - v1[0]*v2[2];
	    result[2] = v1[0]*v2[1] - v1[1]*v2[0];
	    
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
	
	template<class T>
	double norm(const vector_k<T> &v1, std::string method = "2")
	{
	    if(v1.empty())
	    {
	        throw KeyCppException("Cannot compute norm of empty vector!");
	    }
		std::transform(method.begin(), method.end(), method.begin(), ::tolower);
	    double anorm = 0.0;
	    if(!method.empty() && method.find_first_not_of("-+0123456789") == std::string::npos)
	    {
	        int p = atoi(method.c_str());
	        for(size_t ii = 0; ii < v1.size(); ii++)
	        {
	            anorm += pow(std::abs(v1[ii]),p);
	        }
	        anorm = pow(anorm,1.0/p);
	    }
	    else if(method.compare("inf") == 0)
	    {
	        anorm = 0.0;
	        for(size_t ii = 0; ii < v1.size(); ii++)
	        {
	            if(std::abs(v1[ii]) > anorm)
	            {
	                anorm = std::abs(v1[ii]);
	            }
	        }
	    }
	    else if(method.compare("-inf") == 0)
	    {
	        anorm = -1.0;
	        for(size_t ii = 0; ii < v1.size(); ii++)
	        {
	            if(std::abs(v1[ii]) < anorm || anorm < 0.0)
	            {
	                anorm = std::abs(v1[ii]);
	            }
	        }
	    }
	    else
	    {
	        throw KeyCppException("Error in norm! Unknown norm type!");
	    }
	    
	    return anorm;
	}
	
	
	template<class T>
	double norm(const matrix<T,1> &v1, std::string method = "2")
	{
	    if(v1.empty())
	    {
	        throw KeyCppException("Cannot compute norm of empty vector!");
	    }
		std::transform(method.begin(), method.end(), method.begin(), ::tolower);
	    double anorm = 0.0;
	    if(!method.empty() && method.find_first_not_of("-+0123456789") == std::string::npos)
	    {
	        int p = atoi(method.c_str());
	        for(size_t ii = 0; ii < v1.size(1); ii++)
	        {
	            anorm += pow(std::abs(v1(ii)),p);
	        }
	        anorm = pow(anorm,1.0/p);
	    }
	    else if(method.compare("inf") == 0)
	    {
	        anorm = 0.0;
	        for(size_t ii = 0; ii < v1.size(1); ii++)
	        {
	            if(std::abs(v1(ii)) > anorm)
	            {
	                anorm = std::abs(v1(ii));
	            }
	        }
	    }
	    else if(method.compare("-inf") == 0)
	    {
	        anorm = -1.0;
	        for(size_t ii = 0; ii < v1.size(1); ii++)
	        {
	            if(std::abs(v1(ii)) < anorm || anorm < 0.0)
	            {
	                anorm = std::abs(v1(ii));
	            }
	        }
	    }
	    else
	    {
	        throw KeyCppException("Error in norm! Unknown norm type!");
	    }
	    
	    return anorm;
	}
	
	double norm(const matrix<double> &A_in, std::string method = "2");
	SVD_type<double,double> svd(const matrix<double> &A_in, std::string method = "");
	double norm(const matrix<std::complex<double>> &A_in, std::string method = "2");
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
	    matrix<T> B(output.V.size(1),length(index));
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
	template<class T>
	bool any(const matrix<T> &A)
	{
	    for(size_t ii = 0; ii < A.size(1); ii++)
	    {
	        for(size_t jj = 0; jj < A.size(2); jj++)
	        {
	            if(std::abs(A(ii,jj)) > eps)
	            {
	                return true;
	            }
	        }
	    }
	    return false;
	}
	
	/** \brief Returns true if any elements of v1 are nonzero.
	 */
	template<class T>
	bool any(const vector_k<T> &v1)
	{
	    for(size_t ii = 0; ii < v1.size(); ii++)
	    {
            if(std::abs(v1[ii]) > eps)
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
	template<class T>
	bool all(const matrix<T> &A)
	{
	    for(size_t ii = 0; ii < A.size(1); ii++)
	    {
	        for(size_t jj = 0; jj < A.size(2); jj++)
	        {
	            if(std::abs(A(ii,jj)) < eps)
	            {
	                return false;
	            }
	        }
	    }
	    return true;
	}
	
	/** \brief Returns true if all elements of v1 are nonzero.
	 */
	template<class T>
	bool all(const vector_k<T> &v1)
	{
	    for(size_t ii = 0; ii < v1.size(); ii++)
	    {
            if(std::abs(v1[ii]) < eps)
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
	vector_k<bool> finite(const T &a)
	{
	    bool out;
        if(isfinite(a))
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
	template<class T>
	matrix<bool> finite(const matrix<T> &A)
	{
	    matrix<bool> out(A.size(1),A.size(2));
	    for(size_t ii = 0; ii < A.size(1); ii++)
	    {
	        for(size_t jj = 0; jj < A.size(2); jj++)
	        {
	            out(ii,jj) = finite(A(ii,jj));
            }
	    }
	    return out;
	}
	
	/** \brief Returns vector containing boolean values that are true if
	 *         corresponding elements of v1 are finite.
	 */
	template<class T>
	vector_k<bool> finite(const vector_k<T> &v1)
	{
	    vector_k<bool> out(v1.size());
	    for(size_t ii = 0; ii < v1.size(); ii++)
	    {
            out[ii] = finite(v1[ii]);
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
	template<class T>
	matrix<bool> isinf(const matrix<T> &A)
	{
	    matrix<bool> out(A.size(1),A.size(2));
	    for(size_t ii = 0; ii < A.size(1); ii++)
	    {
	        for(size_t jj = 0; jj < A.size(2); jj++)
	        {
	            out(ii,jj) = std::isinf(A(ii,jj));
            }
	    }
	    return out;
	}
	
	/** \brief Returns vector containing boolean values that are true if
	 *         corresponding elements of v1 are infinite.
	 */
	template<class T>
	vector_k<bool> isinf(const vector_k<T> &v1)
	{
	    vector_k<bool> out(v1.size());
	    for(size_t ii = 0; ii < v1.size(); ii++)
	    {
	        out[ii] = std::isinf(v1[ii]);
	    }
	    return out;
	}
	
	/** \brief Returns matrix containing boolean values that are true if
	 *         corresponding elements of A are NaN.
	 */
	template<class T>
	matrix<bool> isnan(const matrix<T> &A)
	{
	    matrix<bool> out(A.size(1),A.size(2));
	    for(size_t ii = 0; ii < A.size(1); ii++)
	    {
	        for(size_t jj = 0; jj < A.size(2); jj++)
	        {
	            out(ii,jj) = std::isnan(A(ii,jj));
            }
	    }
	    return out;
	}
	
	/** \brief Returns boolean value that is true if
	 *         a is NaN.
	 */
	template<class T>
	bool isnan(const T &a)
	{
	    return std::isnan(a);
	}
	
	/** \brief Returns vector containing boolean values that are true if
	 *         corresponding elements of v1 are NaN.
	 */
	template<class T>
	vector_k<bool> isnan(const vector_k<T> &v1)
	{
	    vector_k<bool> out(v1.size());
	    for(size_t ii = 0; ii < v1.size(); ii++)
	    {
	        out[ii] = std::isnan(v1[ii]);
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
	
	/** \brief Returns true if vector is empty.
	 */
	template<class T>
	vector_k<bool> isempty(const vector_k<T> &v1)
	{
	    return v1.empty();
	}
	
	/** \brief Returns true if a is real.
	 */
	template<class T>
	bool isreal(const T &a)
	{
        if(abs(imag(a)) < eps)
        {
            return true;
        }
	    return false;
	}
	
	/** \brief Returns true if all elements of A are real.
	 */
	template<class T>
	bool isreal(const matrix<T> &A)
	{
	    for(int ii = 0; ii < A.size(1); ii++)
	    {
	        for(int jj = 0; jj < A.size(2); jj++)
	        {
	            if(abs(imag(A(ii,jj))) > eps)
	            {
	                return false;
	            }
	        }
	    }
	    return true;
	}
	
	/** \brief Returns true if all elements of v1 are real.
	 */
	template<class T>
	bool isreal(const vector_k<T> &v1)
	{
	    for(size_t ii = 0; ii < v1.size(); ii++)
	    {
            if(abs(imag(v1[ii])) > eps)
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
	
	/** \brief Rounds the elements of v1 towards positive infinity.
	 */
	template<class T>
	vector_k<T> ceil(const vector_k<T> &v1)
	{
	    vector_k<T> v2(v1.size());
	    for(int ii = 0; ii < v2.size(); ii++)
	    {
	        v2[ii] = ceil(v1[ii]);
	    }
	    return v2;
	}
	
	/** \brief Rounds the elements of A towards positive infinity.
	 */
	template<class T, size_t dim>
	matrix<T,dim> ceil(const matrix<T,dim> &A)
	{
	    return eop(A,static_cast<T (*)(T)>(&std::ceil));
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
	
	/** \brief Rounds the elements of v1 towards negative infinity.
	 */
	template<class T>
	vector_k<T> floor(const vector_k<T> &v1)
	{
	    vector_k<T> v2(v1.size());
	    for(int ii = 0; ii < v2.size(); ii++)
	    {
	        v2[ii] = floor(v1[ii]);
	    }
	    return v2;
	}
	
	/** \brief Rounds the elements of A towards negative infinity.
	 */
	template<class T, size_t dim>
	matrix<T,dim> floor(const matrix<T,dim> &A)
	{
	    return eop(A,static_cast<T (*)(T)>(&std::floor));
	}
	
	template<class T, class U>
	decltype(std::declval<T>()*std::declval<U>()) polyval(const vector_k<T> &p, const U &x)
	{
	    decltype(std::declval<T>()*std::declval<U>()) val = 0.0;
	    
        for(size_t ii = 0; ii < p.size()-1; ii++)
        {
            val += p[ii]*pow(x,p.size()-ii-1);
        }
	    val += p[p.size()-1];
	    return val;
	}
	
	/** \brief Computes all roots of polynomial p by solving for the eigenvalues
	 *         of the companion matrix.
	 */
	template<class T>
	vector_k<T> roots(const vector_k<T> &p)
	{
	    size_t n = p.size()-1;
	    matrix<T> A = diag(ones<T>(n-1),-1);
	    for(size_t ii = 0; ii < A.size(2); ii++)
	    {
            A(0,ii) = -p[ii+1]/p[0];
        }
        vector_k<T> v = eig(A);
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
    inline vector_k<size_t> clock()
    {
        time_t t = time(0);
        struct tm * now = localtime(&t);
        vector_k<size_t> dt(6);
        dt[0] = (now->tm_year + 1900); // year
        dt[1] = (now->tm_mon + 1); // month
        dt[2] = (now->tm_mday); // day
        dt[3] = (now->tm_hour); // hour
        dt[4] = (now->tm_min); // minute
        dt[5] = (now->tm_sec); // seconds
        
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
        
        vector_k<T> s_inv = diag(svd_out.S);
        for(size_t ii = 0; ii < s_inv.size(); ii++)
        {
            if(any(s_inv[ii]))
            {
                s_inv[ii] = 1.0/s_inv[ii];
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
    
    /** \brief Returns the standard deviation of inputed vector. */
    template<class T>
    T stdev(vector_k<T> v1)
    {
        T v_bar = mean(v1);
        T temp = 0.0;
        for(size_t ii = 0; ii < v1.size(); ii++)
        {
            temp += pow((v1[ii] - v_bar),2.0);
        }
        temp = std::sqrt(temp/(v1.size()-1.0));
        
        return temp;
    }
    
    /** \brief Returns the standard deviation of inputed vector. */
    inline double stdev(vector_k<std::complex<double>> v1)
    {
        std::complex<double> v_bar = mean(v1);
        double temp = 0.0;
        for(size_t ii = 0; ii < v1.size(); ii++)
        {
            temp += std::abs((v1[ii] - v_bar)*conj(v1[ii] - v_bar));
        }
        temp = std::sqrt(temp/((double)v1.size()-1.0));
        
        return temp;
    }
    
    /** \brief Returns the variance (square of standard deviation) for inputed vector. */
    template<class T>
    T var(vector_k<T> v1)
    {
        return pow(stdev(v1),2.0);
    }
    
    /** \brief Returns the variance (square of standard deviation) for inputed vector. */
    inline double var(vector_k<std::complex<double>> v1)
    {
        return pow(stdev(v1),2.0);
    }
    
    namespace rng_ns
    {
        int choice = 1;
        std::default_random_engine default_rng(std::random_device{}()); // Choice = 0;
        std::mt19937 mt_rng(std::random_device{}()); // Choice = 1;
        std::ranlux24_base lagfib_rng(std::random_device{}()); // Choice = 2;
    }
    
    inline void rng(size_t seed = 0, std::string generator = "twister")
    {
		std::transform(generator.begin(), generator.end(), generator.begin(), ::tolower);
		if(generator.compare("twister") == 0)
		{
		    rng_ns::choice = 1;
		    rng_ns::mt_rng = std::mt19937(seed);
		}
		else if(generator.compare("multFibonacci") == 0)
		{
		    rng_ns::choice = 2;
		    rng_ns::lagfib_rng = std::ranlux24_base(seed);
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
		        if(rng_ns::choice == 2)
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
		        rng_ns::choice = 1;
		        rng_ns::mt_rng = std::mt19937(0);
		    }
		}
		else
		{
		    if(generator.compare("twister") == 0)
		    {
		        rng_ns::choice = 1;
		        rng_ns::mt_rng = std::mt19937(std::random_device{}());
		    }
		    else if(generator.compare("multFibonacci") == 0)
		    {
		        rng_ns::choice = 2;
		        rng_ns::lagfib_rng = std::ranlux24_base(std::random_device{}());
		    }
		}
    }
    
    /** \brief Returns a random double between 0 and 1.0
	 */
	inline double rand()
	{
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        double randomNumber;
        if(rng_ns::choice == 1)
        {
            randomNumber = distribution(rng_ns::mt_rng);
        }
        else if(rng_ns::choice == 2)
        {
            randomNumber = distribution(rng_ns::lagfib_rng);
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
        if(rng_ns::choice == 1)
        {
            randomNumber = distribution(rng_ns::mt_rng);
        }
        else if(rng_ns::choice == 2)
        {
            randomNumber = distribution(rng_ns::lagfib_rng);
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
    inline vector_k<std::complex<double> > eig(const matrix<std::complex<double> > &A, const matrix<std::complex<double> > &B, matrix<std::complex<double> > *vr_return, matrix<std::complex<double> > *vl_return)
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

		vector_k<std::complex<double> > lambda(n);
		for(unsigned int ii = 0; ii < n; ii++)
		{
			lambda[ii] = alpha[ii]/beta[ii];
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
    inline vector_k<std::complex<double> > eig(const matrix<std::complex<double> > &A, matrix<std::complex<double> > *vr_return, matrix<std::complex<double> > *vl_return)
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

		vector_k<std::complex<double> > lambda(n);
		for(unsigned int ii = 0; ii < n; ii++)
		{
			lambda[ii] = w[ii];
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
    inline vector_k<std::complex<double> > eig(const matrix<double> &A, matrix<std::complex<double> > *vr_return, matrix<std::complex<double> > *vl_return)
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
	
	inline vector_k<std::complex<double>> linsolve(const matrix<std::complex<double>>& A_in,
	                                      const vector_k<std::complex<double>>& b_in)
	{
		if(b_in.empty() || A_in.size(1) <= 0 || A_in.size(2) <= 0)
		{
			throw KeyCppException("Error in linsolve()! Empty matrix or vector supplied!");
		}
		if(A_in.size(2) != b_in.size())
		{
			throw KeyCppException("Error in linsolve()! Matrix and vector sizes are not compatible!");
		}
		
		unsigned int n = (unsigned)b_in.size();
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
        
        vector_k<std::complex<double>> x_out(b_in);
        zgetrs_("N", &m, &nrhs, A, &lda, iw, &x_out[0], &nn, &info);
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
	
	
	inline vector_k<double> linsolve(const matrix<double>& A_in,
	                                      const vector_k<double>& b_in)
	{
		if(b_in.empty() || A_in.size(1) <= 0 || A_in.size(2) <= 0)
		{
			throw KeyCppException("Error in linsolve()! Empty matrix or vector supplied!");
		}
		if(A_in.size(2) != b_in.size())
		{
			throw KeyCppException("Error in linsolve()! Matrix and vector sizes are not compatible!");
		}
		
		unsigned int n = (unsigned)b_in.size();
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
        
        vector_k<double> x_out(b_in);
        dgetrs_("N", &mm, &nrhs, A, &lda, iw, &x_out[0], &nn, &info);
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
	
	inline double norm(const matrix<double> &A_in, std::string method)
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
		    auto svd_out = svd(A_in);
		    anorm = max(max(svd_out.S));
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
	    if(A_in.size(1) <= 0 || A_in.size(2) <= 0)
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
	
	inline double norm(const matrix<std::complex<double>> &A_in, std::string method)
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
		    auto svd_out = svd(A_in);
		    anorm = max(max(svd_out.S));
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
}

#endif
