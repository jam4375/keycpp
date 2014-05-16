// Spline.h -- Class to perform cubic spline interpolation

#ifndef SPLINE_H_
#define SPLINE_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include "Matrix.h"

namespace keycpp
{
	class SplineException : public std::runtime_error
	{
		public:
			SplineException(const std::string &msg) : std::runtime_error(msg){}
	};
	
	class Extrap
	{
		public:
			bool isString;
			bool isDouble;
			std::string extrap_string;
			double extrap_val;
			Extrap() : isString(false), isDouble(false), extrap_string(""), extrap_val() {}
			Extrap(const char * extrap) : isString(true), isDouble(false), extrap_string(extrap), extrap_val() {}
			Extrap(double extrap) : isString(false), isDouble(true), extrap_string(""), extrap_val(extrap) {}
	};

	template<class U, class T>
	class Spline
	{
	private:
		int n;
		
		matrix<U,2> x;
		matrix<T,2> y;

		matrix<T,2> s{};
		matrix<T,2> h{};
		matrix<T,2> f{};

		matrix<T,2> tri{}; // tridiagonal array
		matrix<T,2> a{};
		matrix<T,2> b{};
		matrix<T,2> c{};
		matrix<T,2> d{};
	
		int find_spline();  // create the tridiagonal matrix
		int tridiagonal();  // solve the tridiagonal matrix and find coefficients
		int sort_xy(); // Sort x and y values by x
		
		Extrap extrap;
	
	public:
		Spline(int N, matrix<U,2> X, matrix<T,2> Y, Extrap extrap_in);
		virtual ~Spline();
	
		int compute_spline() {find_spline(); return tridiagonal();};
		T J(U);  // compute interpolation based on cubic Spline
	};


	template<class U, class T> Spline<U,T>::Spline(int N, matrix<U,2> X, matrix<T,2> Y, Extrap extrap_in) : n(N), x(X), y(Y), extrap(extrap_in)
	{
		if(N != X.length() || N != Y.length())
		{
			std::cout << "ERROR!! Problem with size of inputed matrices to Spline!\n";
		}
		if(!X.isVec() || !Y.isVec())
		{
			std::cout << "ERROR!! Problem with size of inputed matrices to Spline!\n";
		}

		s = matrix<T,2>(n);
		h = matrix<T,2>(n-1);
		f = matrix<T,2>(n-1);
		a = matrix<T,2>(n-1);
		b = matrix<T,2>(n-1);
		c = matrix<T,2>(n-1);
		d = matrix<T,2>(n-1);
	}

	/** \brief Spline destructor, nothing needs to be done. */
	template<class U, class T> Spline<U,T>::~Spline()
	{
	}

	template<class U, class T> int Spline<U,T>::find_spline()
	{
		// Sort x and y values:
		sort_xy();

		// Create tridiagonal matrix

		// compute each h
		for(int i = 0; i < n-1; i++)
		{
			h(i) = x(i+1) - x(i);
		}

		// Compute each f
		for(int i = 0; i < n-1; i++)
		{
			f(i) = (y(i+1) - y(i))/h(i);
		}

		//allocate memory for tridiagonal matrix
		tri = matrix<T,2>(n-2,4);

		// add the right side
		for(int i = 0; i < n-2; i++)
		{
			tri(i,3) = 6.0*(f(i+1)-f(i));
		}

		//fill two side bands
		tri(0,0) = 0;
		tri(n-3,2) = 0;
		for(int i = 1; i < n-2; i ++)
		{
			tri(i,0) = h(i);
			tri(i-1,2) = h(i);
		}

		//fill middle band
		for(int i = 0; i < n-2; i++)
		{
			tri(i,1) = 2.0*(h(i)+h(i+1));
		}

		//set natural conditions
		s(0) = 0.0;
		s(n-1) = 0.0;

		return 1;
	}

	template<class U, class T> int Spline<U,T>::tridiagonal()
	{
		// solve tridiagonal matrix
		for(int i = 1; i < n-2; i++)
		{
			tri(i,0) = tri(i,0)/tri(i-1,1);
			tri(i,1) = tri(i,1) - tri(i,0)*tri(i-1,2);
			tri(i,3) = tri(i,3) - tri(i,0)*tri(i-1,3);
		}

		tri(n-3,3) = tri(n-3,3)/tri(n-3,1);
		for(int i = n-4; i >=0; i--)
		{
			tri(i,3) = (tri(i,3) - tri(i,2)*tri(i+1,3))/tri(i,1);
		}

		// find each curvature
		// Use natural spline boundary conditions:
		s(0) = 0;
		s(n-1) = 0;
		for(int i = 1; i < n-1; i ++)
		{
			s(i) = tri(i-1,3);
		}

		// compute each coefficient
		for(int i = 0; i < n-1; i++)
		{
			a(i) = (s(i+1) - s(i))/(6.0*h(i));
			b(i) = s(i)/2.0;
			c(i) = (y(i+1) - y(i))/h(i) - (2.0*h(i)*s(i) + h(i)*s(i+1))/6.0;
			d(i) = y(i);
		}

		return 1;
	}

	template<class U, class T> int Spline<U,T>::sort_xy()
	{
		bool swapped = true;
		U temp;
		int temp_i;
		matrix<int,2> index(x.length());
		for(int ii = 0; ii < x.length(); ii++)
		{
			index(ii) = ii;
		}
		while(swapped)
		{     
			swapped = false;
			for(int ii = 1; ii < x.length(); ii++)
			{
				if((x(ii-1)) > (x(ii)))
				{
					temp = x(ii-1);
					x(ii-1) = x(ii);
					x(ii) = temp;
					temp_i = index(ii-1);
					index(ii-1) = index(ii);
					index(ii) = temp_i;
					swapped = true;
				}
			}
		}
		matrix<T,2> y_temp(y.length());
		for(int ii = 0; ii < y.length(); ii++)
		{
			y_temp(ii) = y(index(ii));
		}
		for(int ii = 0; ii < y.length(); ii++)
		{
			y(ii) = y_temp(ii);
		}

		return 1;
	}

	template<class U, class T> T Spline<U,T>::J(U x_interp)
	{
		// find which segment to use
		int index = 0;
		for(int i = 0; i < n-1; i++)
		{
			if(x_interp > x(i) && x_interp <= x(i+1))
			{
				index = i;
			}
		}
		if(x_interp <= x(0))
		{
			if(x_interp < x(0))
			{
				if(extrap.isString)
				{
					if(extrap.extrap_string.compare("extrap") == 0)
					{
						index = 0;
					}
					else
					{
						throw SplineException("ERROR! Could not interpolate!! Unknown string in interp1!");
						return nan("");
					}
				}
				else if(extrap.isDouble)
				{
					return extrap.extrap_val;
				}
				else
				{
					throw SplineException("ERROR! Could not interpolate!!");
					return nan("");
				}
			}
		}
		if(x_interp >= x(n-1))
		{
			if(x_interp > x(n-1))
			{
				if(extrap.isString)
				{
					if(extrap.extrap_string.compare("extrap") == 0)
					{
						index = n-2;
					}
					else
					{
						throw SplineException("ERROR! Could not interpolate!! Unknown string in interp1!");
						return nan("");
					}
				}
				else if(extrap.isDouble)
				{
					return extrap.extrap_val;
				}
				else
				{
					throw SplineException("ERROR! Could not interpolate!!");
					return nan("");
				}
			}
		}

		// return interpolated value
		return (a(index)*pow(x_interp-x(index),3.0) + b(index)*pow(x_interp-x(index),2.0) + c(index)*(x_interp-x(index)) + d(index));
	}
}

#endif
