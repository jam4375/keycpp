// Matlab.h -- Common Matlab functions implemented in C++

#ifndef KEYCPP_H_
#define KEYCPP_H_

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <string>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <utility>
#include <algorithm>
#include "Matrix.h"
#include "kiss_fft.h"
#include "Spline.h"
#include "Figure.h"

namespace keycpp
{
	#define pi 3.1415926535897932384626433832795
	
	class KeyCppException : public std::runtime_error
	{
		public:
			KeyCppException(const std::string &msg) : std::runtime_error(msg){}
	};
	
	/** \brief This provides a C interface to the fortran LAPACK function "zggev_"
	 */
	extern "C"{
	void zggev_(char *jobvl, char *jobvr, int *n, std::complex<double> *a,
		    int *lda, std::complex<double> *b, int *ldb, std::complex<double> *alpha,
		    std::complex<double> *beta, std::complex<double> *vl,
		    int *ldvl, std::complex<double> *vr, int *ldvr,
		    std::complex<double> *work, int *lwork, double *rwork, int *info);
	}

	std::vector<std::complex<double> > eig(matrix<std::complex<double> > A, matrix<std::complex<double> > B, matrix<std::complex<double> > *vr_return = NULL, matrix<std::complex<double> > *vl_return = NULL);

	std::complex<double> interp1(double *x, std::complex<double> **y, double x_interp, int index, int N);
	double rand_limits(double a, double b);
	
	template<class T, class U> struct observe
	{
		std::vector<T>& y;
		std::vector<U>& x_ode;

		observe(std::vector<T> &p_y, std::vector<U> &p_x_ode) : y(p_y), x_ode(p_x_ode) { };

		void operator()(const T &y_temp, U x_temp)
		{
			y.push_back(y_temp);
			x_ode.push_back(x_temp);
		}
	};

	template<class T> std::vector<std::complex<double> > conj(const std::vector<T> x)
	{
		std::vector<T> x_out(x.size());
		for(int ii = 0; ii < x.size(); ii++)
		{
			x_out[ii] = conj(x[ii]);
		}
		return x_out;
	}

	template<class T> matrix<std::complex<double> > conj(const matrix<T> x)
	{
		matrix<T> x_out(x.size(1),x.size(2));
		for(int ii = 0; ii < x.size(1); ii++)
		{
			for(int jj = 0; jj < x.size(2); jj++)
			{
				x_out(ii,jj) = conj(x(ii,jj));
			}
		}
		return x_out;
	}

	template<class T> std::vector<double> real(std::vector<T> x)
	{
		std::vector<double> y(x.size());
		for(int ii = 0; ii < x.size(); ii++)
		{
			y[ii] = real(x[ii]);
		}
		return y;
	}

	template<class T> std::vector<double> imag(std::vector<T> x)
	{
		std::vector<double> y(x.size());
		for(int ii = 0; ii < x.size(); ii++)
		{
			y[ii] = imag(x[ii]);
		}
		return y;
	}

	template<class T> std::vector<double> abs(std::vector<T> x)
	{
		std::vector<double> y(x.size());
		for(int ii = 0; ii < x.size(); ii++)
		{
			y[ii] = abs(x[ii]);
		}
		return y;
	}

	template<class T> std::vector<double> arg(std::vector<T> x)
	{
		std::vector<double> y(x.size());
		for(int ii = 0; ii < x.size(); ii++)
		{
			y[ii] = arg(x[ii]);
		}
		return y;
	}

	template<class T, class U> std::vector<decltype(std::declval<T>()*std::declval<U>())> operator+(const std::vector<T>& v1, const std::vector<U>& v2)
	{
		std::vector<decltype(std::declval<T>()*std::declval<U>())> result(v1.size());
		for(int ii = 0; ii < result.size(); ii++)
		{
			result[ii] = v1[ii]+v2[ii];
		}
		return result;
	}

	template<class T, class U> matrix<decltype(std::declval<T>()*std::declval<U>())> operator*(const T& a, const matrix<U>& A)
	{
		matrix<decltype(std::declval<T>()*std::declval<U>())> B(A.size(1),A.size(2));
		for(int ii = 0; ii < B.size(1); ii++)
		{
			for(int jj = 0; jj < B.size(2); jj++)
			{
				B(ii,jj) = a*A(ii,jj);
			}
		}
		return B;
	}
	
	template<class T, class U> matrix<decltype(std::declval<T>()*std::declval<U>())> operator*(const matrix<U>& A, const T& a)
	{
		matrix<decltype(std::declval<T>()*std::declval<U>())> B(A.size(1),A.size(2));
		for(int ii = 0; ii < B.size(1); ii++)
		{
			for(int jj = 0; jj < B.size(2); jj++)
			{
				B(ii,jj) = a*A(ii,jj);
			}
		}
		return B;
	}
	
	template<class T, class U> matrix<decltype(std::declval<T>()*std::declval<U>())> operator*(const std::vector<T>& v1, const matrix<U>& A)
	{
		if(A.size(1) == 1)
		{
			matrix<decltype(std::declval<T>()*std::declval<U>())> B(v1.size(),A.size(2));
			for(int ii = 0; ii < B.size(1); ii++)
			{
				for(int jj = 0; jj < B.size(2); jj++)
				{
					B(ii,jj) = v1[ii]*A(0,jj);
				}
			}
			return B;	
		}
		else if(A.size(1) == v1.size())
		{
			matrix<decltype(std::declval<T>()*std::declval<U>())> B(1,A.size(2));
			for(int ii = 0; ii < A.size(2); ii++)
			{
				for(int jj = 0; jj < A.size(1); jj++)
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

	template<class T, class U> std::vector<decltype(std::declval<T>()*std::declval<U>())> operator*(const T& a, const std::vector<U>& v1)
	{
		std::vector<decltype(std::declval<T>()*std::declval<U>())> v2(v1.size());
		for(int ii = 0; ii < v2.size(); ii++)
		{
			v2[ii] = a*v1[ii];
		}
		return v2;
	}
	
	template<class T, class U> std::vector<decltype(std::declval<T>()*std::declval<U>())> operator*(const std::vector<T>& v1, const U& a)
	{
		std::vector<decltype(std::declval<T>()*std::declval<U>())> v2(v1.size());
		for(int ii = 0; ii < v2.size(); ii++)
		{
			v2[ii] = a*v1[ii];
		}
		return v2;
	}
	
	template<class T> std::vector<T> exp(const std::vector<T> v1)
	{
		std::vector<T> v2(v1.size());
		for(int ii = 0; ii < v1.size(); ii++)
		{
			v2[ii] = exp(v1[ii]);
		}
		return v2;
	}
	
	template<class T> matrix<T> exp(const matrix<T> A)
	{
		matrix<T> B(A.size(1),A.size(2));
		for(int ii = 0; ii < A.size(1); ii++)
		{
			for(int jj = 0; jj < A.size(2); jj++)
			{
				B(ii,jj) = exp(A(ii,jj));
			}
		}
		return B;
	}

	template<class T> matrix<T> eye(int N)
	{
		matrix<T> A(N,N);
		for(int ii = 0; ii < N; ii++)
		{
			A(ii,ii) = 1.0;
		}
	
		return A;
	}
	
	template<class T> matrix<T> repmat(const matrix<T> A, int m, int n)
	{
		matrix<T> B(m*A.size(1), n*A.size(2));
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
	
	template<class T> matrix<T> repmat(const std::vector<T> v1, int m, int n)
	{
		matrix<T> B(m*v1.size(), n);
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
	
	template<class T, class U> matrix<decltype(std::declval<T>()*std::declval<U>())> etimes(const matrix<T>& A, const matrix<U>& B)
	{
		matrix<decltype(std::declval<T>()*std::declval<U>())> C(A.size(1),B.size(2));
		for(int ii = 0; ii < A.size(1); ii++)
		{
			for(int jj = 0; jj < A.size(2); jj++)
			{
				C(ii,jj) = A(ii,jj)*B(ii,jj);
			}
		}
		return C;
	}
	
	template<class T, class U> std::vector<decltype(std::declval<T>()*std::declval<U>())> etimes(const std::vector<T>& v1, const std::vector<U>& v2)
	{
		std::vector<decltype(std::declval<T>()*std::declval<U>())> v3(v1.size());
		for(int ii = 0; ii < v1.size(); ii++)
		{
			v3[ii] = v1[ii]*v2[ii];
		}
		return v3;
	}

	template<class T> int sign(T val)
	{
	    return (T(0) < val) - (val < T(0));
	}

	template<class T> T max(std::vector<T> x)
	{
		double a = 0;
		int index = 0;
		for(int ii = 0; ii < x.size(); ii++)
		{
			if(abs(x[ii]) > a)
			{
				a = abs(x[ii]);
				index = ii;
			}
		}
		return x[index];
	}

	template<class T> matrix<T> transpose(matrix<T> A)
	{
		matrix<T> B(A.size(2),A.size(1));
		for(int ii = 0; ii < A.size(1); ii++)
		{
			for(int jj = 0; jj < A.size(2); jj++)
			{
				B(jj,ii) = A(ii,jj);
			}
		}
		return B;
	}
	
	template<class T> matrix<T> transpose(std::vector<T> v1)
	{
		matrix<T> B(1,v1.size());
		for(int ii = 0; ii < v1.size(); ii++)
		{
			B(0,ii) = v1[ii];
		}
		return B;
	}
	
	template<class T> std::vector<T> sum(matrix<T> A)
	{
		std::vector<T> v1(A.size(2));
		for(int ii = 0; ii < v1.size(); ii++)
		{
			v1[ii] = T(0);
			for(int jj = 0; jj < A.size(1); jj++)
			{
				v1[ii] += A(jj,ii);
			}
		}
		return v1;
	}
	
	template<class T> T sum(std::vector<T> v1)
	{
		T a;
		for(int ii = 0; ii < v1.size(); ii++)
		{
			a += v1[ii];
		}
		return a;
	}
	
	template<class T> std::vector<T> mat2vec(matrix<T> A)
	{
		std::vector<T> v1(A.size(1));
		for(int ii = 0; ii < v1.size(); ii++)
		{
			v1[ii] = A(ii,0);
		}
		return v1;
	}
	
	template<class T> matrix<T> vec2mat(std::vector<T> v1)
	{
		matrix<T> A(v1.size(),1);
		for(int ii = 0; ii < v1.size(); ii++)
		{
			A(ii,0) = v1[ii];
		}
		return A;
	}

	template<class T> std::vector<T> linspace(T x1, T x2, int N)
	{
		std::vector<T> x(N);
		if(N == 1)
		{
			x[0] = x2;
			return x;
		}

		T delta_x = (x2-x1)/(N-1);

		for(int ii = 0; ii < N; ii++)
		{
			x[ii] = x1 + ii*delta_x;
		}

		x[N-1] = x2;

		return x;
	}

	template<class T> std::vector<T> logspace(T x1, T x2, int N)
	{
		std::vector<T> x(N);
		if(N == 1)
		{
			x[0] = x2;
			return x;
		}

		T delta_x = (x2-x1)/(N-1);

		for(int ii = 0; ii < N; ii++)
		{
			x[ii] = pow(10.0,(x1 + delta_x*ii));
		}

		return x;
	}
	
	template<class T> std::vector<T> unwrap(const std::vector<T>& v1, T tol = pi)
	{
		std::vector<T> v2(v1.size());
		v2[0] = v1[0];
		int correction = 0;
		for(int ii = 1; ii < v1.size(); ii++)
		{
			if((v1[ii] - v1[ii-1]) > tol)
			{
				correction -= 1;
			}
			else if((v1[ii] - v1[ii-1]) < -tol)
			{
				correction += 1;
			}
			v2[ii] = v1[ii] + correction*2*pi;
		}
		return v2;
	}
	
	template<class T> T mean(const std::vector<T>& v1)
	{
		T m = T(0);
		double tot = 0.0;
		for(int ii = 0; ii < v1.size(); ii++)
		{
			m += v1[ii];
			tot += 1.0;
		}
		return m/tot;
	}
	
	template<class T, class U> T interp1(std::vector<U> x, std::vector<T> y, U x_interp, std::string method = "linear", Extrap extrap = Extrap())
	{
		if(x.empty() || y.empty())
		{
			throw KeyCppException("Error in interp1! Cannot interpolate on an empty vector!");
		}
		if(x.size() != y.size())
		{
			throw KeyCppException("Error in interp1! Variables `x` and `y` have incompatible sizes!");
		}
		std::transform(method.begin(), method.end(), method.begin(), ::tolower);
		int N = x.size();
		T y2;
		if(method.compare("spline") == 0)
		{
			Spline<U,T> spline(N,x,y,extrap);
			spline.compute_spline();
			y2 = spline.J(x_interp);
		}
		else if(method.compare("linear") == 0)
		{
			for(int ii = 0; ii < (N-1); ii++)
			{
				if(x_interp == x[ii])
				{
					return y[ii];
				}
				else if((x_interp > x[ii] && x_interp < x[ii+1]) || (x_interp < x[ii] && x_interp > x[ii+1]))
				{
					return (y[ii] + (x_interp - x[ii])*(y[ii+1] - y[ii])/(x[ii+1] - x[ii]));
				}
			}

			if(x_interp == x[N-1])
			{
				return y[N-1];
			}
			
			if(extrap.isString)
			{
				if(extrap.extrap_string.compare("extrap") == 0)
				{
					if(x_interp < x[0])
					{
						return (y[0] + (x_interp - x[0])*(y[1] - y[0])/(x[1] - x[0]));
					}
					else
					{
						return (y[N-2] + (x_interp - x[N-2])*(y[N-1] - y[N-2])/(x[N-1] - x[N-2]));
					}
				}
				else
				{
					throw KeyCppException("ERROR! Could not interpolate!! Unknown string in interp1!");
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
			for(int ii = 0; ii < N; ii++)
			{
				if(std::abs(x[ii] - x_interp) < std::abs(min_val))
				{
					min_val = x[ii] - x_interp;
					index = ii;
				}
				else if(std::abs(x[ii] - x_interp) == std::abs(min_val) && (x[ii] - x_interp) > min_val)
				{
					min_val = x[ii] - x_interp;
					index = ii;
				}
			}
			if(index >= 0)
			{
				return y[index];
			}
			
			if(extrap.isString)
			{
				if(extrap.extrap_string.compare("extrap") == 0)
				{
					throw KeyCppException("ERROR in interp1! Cannot extrapolate using method `nearest`!");
				}
				else
				{
					throw KeyCppException("ERROR! Could not interpolate!! Unknown string in interp1!");
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


	template<class T, class U> std::vector<T> interp1(std::vector<U> x, std::vector<T> y, std::vector<U> x_interp, std::string method = "linear", Extrap extrap = Extrap())
	{
		if(x.empty() || y.empty() || x_interp.empty())
		{
			throw KeyCppException("Error in interp1! Cannot interpolate on an empty vector!");
		}
		if(x.size() != y.size())
		{
			throw KeyCppException("Error in interp1! Variables `x` and `y` have incompatible sizes!");
		}
		std::transform(method.begin(), method.end(), method.begin(), ::tolower);
		int N = x.size();
		int N_int = x_interp.size();
		std::vector<T> y2(N_int);

		if(method.compare("spline") == 0)
		{	
			Spline<U,T> spline(N,x,y,extrap);
			spline.compute_spline();
			for(int ii = 0; ii < N_int; ii++)
			{
				y2[ii] = spline.J(x_interp[ii]);
			}
		}
		else if(method.compare("linear") == 0 || method.compare("nearest") == 0)
		{
			for(int ii = 0; ii < N_int; ii++)
			{
				y2[ii] = interp1(x, y, x_interp[ii],method,extrap);
			}
		}
		else
		{
			throw KeyCppException("Error in interp1! Unrecognized interpolation method!");
		}

		return y2;
	}

	template<class T, class U> matrix<T> interp1(std::vector<U> x, matrix<T> y, std::vector<U> x_interp, std::string method = "linear", Extrap extrap = Extrap())
	{
		if(x.empty() || y.size(1) <= 0 || y.size(2) <= 0 || x_interp.empty())
		{
			throw KeyCppException("Error in interp1! Cannot interpolate on an empty vector or matrix!");
		}
		if(x.size() != y.size(1))
		{
			throw KeyCppException("Error in interp1! Variables `x` and `y` have incompatible sizes!");
		}
		matrix<T> y2(x_interp.size(),y.size(2));

		for(int kk = 0; kk < y.size(2); kk++)
		{
			y2.setCol(interp1(x,y.getCol(kk),x_interp,method, extrap),kk);
		}

		return y2;
	}


	template<class T, class U> matrix<T> interp1(std::vector<U> x, std::vector<std::vector<T> > y, std::vector<U> x_interp, std::string method = "linear", Extrap extrap = Extrap())
	{
		if(x.empty() || y.empty() || y[0].empty() || x_interp.empty())
		{
			throw KeyCppException("Error in interp1! Cannot interpolate on an empty vector or matrix!");
		}
		if(x.size() != y.size())
		{
			throw KeyCppException("Error in interp1! Variables `x` and `y` have incompatible sizes!");
		}
		
		matrix<T> y2(x_interp.size(),y[0].size());
		for(int kk = 0; kk < y[0].size(); kk++)
		{
			std::vector<T> y_temp(x.size());
			for(int ii = 0; ii < x.size(); ii++)
			{
				y_temp[ii] = y[ii][kk];
			}
			y2.setCol(interp1(x,y_temp,x_interp,method, extrap),kk);
		}

		return y2;
	}
	
	template<class T, class U> matrix<T> interp1(std::vector<U> x, std::vector<T> y, matrix<U> x_interp, std::string method = "linear", Extrap extrap = Extrap())
	{
		if(x.empty() || y.empty() || x_interp.size(1) <= 0 || x_interp.size(2) <= 0)
		{
			throw KeyCppException("Error in interp1! Cannot interpolate on an empty vector or matrix!");
		}
		if(x.size() != y.size())
		{
			throw KeyCppException("Error in interp1! Variables `x` and `y` have incompatible sizes!");
		}
		matrix<T> y2(x_interp.size(1),x_interp.size(2));

		for(int kk = 0; kk < x_interp.size(2); kk++)
		{
			y2.setCol(interp1(x,y,x_interp.getCol(kk),method, extrap),kk);
		}

		return y2;
	}	
	
	template<class T, class U> matrix<T> interp1(std::vector<U> x, std::vector<T> y, std::vector<std::vector<U> > x_interp, std::string method = "linear", Extrap extrap = Extrap())
	{
		if(x.empty() || y.empty() || x_interp.empty() || x_interp[0].empty())
		{
			throw KeyCppException("Error in interp1! Cannot interpolate on an empty vector or matrix!");
		}
		if(x.size() != y.size())
		{
			throw KeyCppException("Error in interp1! Variables `x` and `y` have incompatible sizes!");
		}
		matrix<T> y2(x_interp.size(),x_interp[0].size());

		for(int kk = 0; kk < x_interp.size(); kk++)
		{
			y2.setRow(interp1(x,y,x_interp[kk],method, extrap),kk);
		}

		return y2;
	}

	template<class T> std::vector<T> linsolve(const matrix<T>& A_in, const std::vector<T>& b_in)
	{
		if(b_in.empty() | A_in.size(1) <= 0 || A_in.size(2) <= 0)
		{
			throw KeyCppException("Error in linsolve()! Empty matrix or vector supplied!");
		}
		if(A_in.size(2) != b_in.size())
		{
			throw KeyCppException("Error in linsolve()! Matrix and vector sizes are not compatible!");
		}
		int N = b_in.size();
		matrix<T> A = A_in;
		std::vector<T> b = b_in;
		std::vector<T> x_out(N);

		T ratio, temp;

		// Gaussian elimination
		for(int ii = 0; ii < (N-1); ii++)
		{
			if(std::abs(A(ii,ii)) == 0.0)
			{
				throw KeyCppException("ERROR! Divide by zero in linsolve!!");
			}
			for(int jj = (ii+1); jj < N; jj++)
			{
				ratio = A(jj,ii)/A(ii,ii);
				for(int count=ii; count<N; count++)
				{
					A(jj,count) -= (ratio*A(ii,count));
				}
				b[jj] -= (ratio*b[ii]);
			}
		}

		// Back substitution
		if(std::abs(A(N-1,N-1)) == 0.0)
		{
			throw KeyCppException("ERROR! Divide by zero in linsolve!!");
		}
		else
		{
			x_out[N-1] = b[N-1]/A(N-1,N-1);
		}
		for(int ii = (N-2); ii >= 0; ii--)
		{
			temp = b[ii];
			for(int jj=(ii+1); jj<N; jj++)
			{
				temp -= (A(ii,jj)*x_out[jj]);
			}
			x_out[ii] = temp/A(ii,ii);
		}

		return x_out;
	}

	template<class U, class T> T trapz(std::vector<U> eta, std::vector<T> integrand)
	{
		if(eta.empty() || integrand.empty())
		{
			throw KeyCppException("Error in trapz()! Empty vector supplied!");
		}
		if(eta.size() != integrand.size())
		{
			throw KeyCppException("Error in trapz()! Vector sizes are not compatible!");
		}
		int N = eta.size();
		T sum = 0.0;

		for(int ii = 0; ii < N-1; ii++)
		{
			sum += (eta[ii+1] - eta[ii])*(integrand[ii+1] + integrand[ii]);
		}
		return 0.5*sum;
	}


	template<class T, class U> matrix<T> diffxy(matrix<U> eta, matrix<T> u, int index = 2)
	{
		if(eta.size(1) <= 0 || eta.size(2) <= 0 || u.size(1) <= 0 || u.size(2) <= 0)
		{
			throw KeyCppException("Error in diffxy()! Empty matrix supplied!");
		}
		if(eta.size(1) != u.size(1) && eta.size(2) != u.size(2))
		{
			throw KeyCppException("Error in diffxy()! Matrix sizes are not compatible!");
		}
		int N = u.size(1);
		int P = u.size(2);

		matrix<T> du(N,P);
		if(index == 1)
		{
			for(int p = 0; p < P; p++)
			{
				for(int n = 0; n < N-1; n++)
				{
					du(n,p) = (u(n+1,p) - u(n,p))/(eta(n+1,p) - eta(n,p));
				}
				du(N-1,p) = (u(N-1,p) - u(N-2,p))/(eta(N-1,p) - eta(N-2,p));
			}
		}
		else
		{
			for(int n = 0; n < N; n++)
			{
				for(int p = 0; p < P-1; p++)
				{
					du(n,p) = (u(n,p+1) - u(n,p))/(eta(n,p+1) - eta(n,p));
				}
				du(n,P-1) = (u(n,P-1) - u(n,P-2))/(eta(n,P-1) - eta(n,P-2));
			}
		}

		return du;
	}

	template<class T, class U> matrix<T> diffxy(std::vector<U> eta, matrix<T> u)
	{
		if(eta.empty() || u.size(1) <= 0 || u.size(2) <= 0)
		{
			throw KeyCppException("Error in diffxy()! Empty vector or matrix supplied!");
		}
		if(eta.size() != u.size(1) && eta.size() != u.size(2))
		{
			throw KeyCppException("Error in diffxy()! Vector and matrix sizes are not compatible!");
		}
		int N = u.size(1);
		int P = u.size(2);

		matrix<T> du(N,P);
		if(N == eta.size())
		{
			for(int p = 0; p < P; p++)
			{
				for(int ii = 0; ii < N-1; ii++)
				{
					du(ii,p) = (u(ii+1,p) - u(ii,p))/(eta[ii+1] - eta[ii]);
				}
				du(N-1,p) = (u(N-1,p) - u(N-2,p))/(eta[N-1] - eta[N-2]);
			}
		}
		else
		{
			for(int ii = 0; ii < N; ii++)
			{
				for(int p = 0; p < P-1; p++)
				{
					du(ii,p) = (u(ii,p+1) - u(ii,p))/(eta[p+1] - eta[p]);
				}
				du(ii,P-1) = (u(ii,P-1) - u(ii,P-2))/(eta[P-1] - eta[P-2]);
			}
		}

		return du;
	}

	template<class T, class U> std::vector<T> diffxy(std::vector<U> eta, std::vector<T> u)
	{
		if(eta.empty() || u.empty())
		{
			throw KeyCppException("Error in diffxy()! Empty vector supplied!");
		}
		if(eta.size() != u.size())
		{
			throw KeyCppException("Error in diffxy()! Vector sizes are not compatible!");
		}
		int N = u.size();

		std::vector<T> du(N);
		for(int ii = 0; ii < N-1; ii++)
		{
			du[ii] = (u[ii+1] - u[ii])/(eta[ii+1] - eta[ii]);
		}
		du[N-1] = (u[N-1] - u[N-2])/(eta[N-1] - eta[N-2]);

		return du;
	}

	template<class T> std::vector<std::complex<double> > fft(std::vector<T> u, int N = -1)
	{
		if(u.empty())
		{
			throw KeyCppException("Error in fft()! Empty vector supplied!");
		}
		
		if(N < 0)
		{
			N = u.size();
		}
		
		kiss_fft_cpx *cx_in = new kiss_fft_cpx[N];
		kiss_fft_cpx *cx_out = new kiss_fft_cpx[N];

		std::vector<std::complex<double> > u_hat(N);

		for(int ii = 0; ii < N; ii++)
		{
			cx_in[ii].r = real((std::complex<double>)u[ii]);
			cx_in[ii].i = real((std::complex<double>)u[ii]);
		}

		kiss_fft_cfg cfg = kiss_fft_alloc(N,false,0,0);
		kiss_fft(cfg,cx_in,cx_out);

		for(int ii = 0; ii < N; ii++)
		{
			u_hat[ii] = std::complex<T>((T)cx_out[ii].r,(T)cx_out[ii].i)*2.0/(T)N;
		}

		free(cfg);
		delete cx_in;
		delete cx_out;

		return u_hat;
	}

	template<class T, class U, class F>
	matrix<T> ode45(F odeClass, std::vector<U> x_ode, std::vector<T> ICs, double abs_tol = 1.0e-10, double rel_tol = 1.0e-6)
	{
		if(x_ode.empty())
		{
			throw KeyCppException("Error in ode45()! Vector x_ode cannot be empty!");
		}
		if(ICs.empty())
		{
			throw KeyCppException("Error in ode45()! Must provide initial conditions!");
		}
		if(x_ode.size() < 2)
		{
			throw KeyCppException("Error in ode45()! Invalid vector x_ode!");
		}
		U x0 = x_ode[0];
		U xf = x_ode[x_ode.size()-1];
		U delta_x0 = x_ode[1] - x_ode[0];
		std::vector<std::vector<T> > y_temp;
		std::vector<U> x_temp;

		{
			using namespace boost::numeric::odeint;
			size_t steps = integrate_adaptive(make_controlled<runge_kutta_cash_karp54<std::vector<T> > >(abs_tol, rel_tol), odeClass, ICs, x0, xf, delta_x0, observe<std::vector<T>,U>(y_temp,x_temp));
		}
		matrix<T> y = interp1(x_temp, y_temp, x_ode);
		return y;
	}
	
	inline void set(Figure &h, std::string property, double val)
	{
		h.set(property,val);
	}
}

#endif
