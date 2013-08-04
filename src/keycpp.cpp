#include <stdlib.h>
#include <ctime>
#include "keycpp.h"
using namespace std;

namespace keycpp
{
	double rand_limits(double a, double b)
	{
		return (a + (double)rand()/((double)RAND_MAX)*(b-a));
	}

	complex<double> interp1(double *x, complex<double> **y, double x_interp, int index, int N)
	{
		for(int ii = 0; ii < (N-1); ii++)
		{
			if(x_interp == x[ii])
			{
				return y[ii][index];
			}
			else if((x_interp > x[ii] && x_interp < x[ii+1]) || (x_interp < x[ii] && x_interp > x[ii+1]))
			{
				return (y[ii][index] + (x_interp - x[ii])*(y[ii+1][index] - y[ii][index])/(x[ii+1] - x[ii]));
			}
		}

		if(x_interp == x[N-1])
		{
			return y[N-1][index];
		}

		cout << "ERROR! Could not interpolate!!\n";
		return complex<double> (0,0);
	}


	/** \brief Generalized complex-valued eigenvalue solver using LAPACK function call. 
	 *  
	 *  This function returns the eigenvalues(lambda) of the complex-valued generalized
	 *  eigenvalue problem: Ax_r = lambda*Bx_r or x_l^T*A = lambda*x_l^T*B. The eigenvalues
	 *  are returned by default. To return the right or left eigenvectors, supply the
	 *  function with a complex<double> matrix object in the 3rd or 4th parameters, respectively.
	 */
	vector<complex<double> > eig(matrix<complex<double> > A, matrix<complex<double> > B, matrix<complex<double> > *vr_return, matrix<complex<double> > *vl_return)
	{
		int n, lda, ldb, ldvl, ldvr, lwork, info;
		n = lda = ldb = A.size(1);
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

		complex<double> *a = new complex<double>[n*n];
		complex<double> *b = new complex<double>[n*n];
		complex<double> *vr = new complex<double>[n*n];
		complex<double> *vl = new complex<double>[n*n];
		complex<double> *alpha = new complex<double>[n];
		complex<double> *beta = new complex<double>[n];
		complex<double> *work = new complex<double>[lwork];
		double *rwork = new double[8*n];
		for(int ii = 0; ii < n; ii++)
		{
			for(int jj = 0; jj < n; jj++)
			{
				a[ii*n + jj] = A(jj,ii);
				b[ii*n + jj] = B(jj,ii);
			}
		}

		zggev_(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alpha, beta, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);

		vector<complex<double> > lambda(n);
		for(int ii = 0; ii < n; ii++)
		{
			lambda[ii] = alpha[ii]/beta[ii];
		}
		if(jobvr == 'V')
		{
			(*vr_return) = matrix<complex<double> >(n,n);
			for(int ii = 0; ii < n; ii++)
			{
				for(int jj = 0; jj < n; jj++)
				{
					(*vr_return)(jj,ii) = vr[ii*n + jj];
				}
			}
		}
		if(jobvl == 'V')
		{
			(*vl_return) = matrix<complex<double> >(n,n);
			for(int ii = 0; ii < n; ii++)
			{
				for(int jj = 0; jj < n; jj++)
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
}

