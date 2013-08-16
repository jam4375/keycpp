#include <stdlib.h>
#include <ctime>
#include "keycpp.h"
using namespace std;

namespace keycpp
{	
    /** \brief Returns a random double between 0 and 1.0
	 */
	double rand()
	{
		return ((double)std::rand()/((double)RAND_MAX));
	}
	
	/** \brief Returns an N x N matrix of random doubles between 0 and 1.0
	 */
	matrix<double> rand(int N)
	{
	    matrix<double> A(N,N);
	    for(int ii = 0; ii < N; ii++)
	    {
	        for(int jj = 0; jj < N; jj++)
	        {
	            A(ii,jj) = ((double)std::rand()/((double)RAND_MAX));
	        }
        }
	    
		return A;
	}
	
	/** \brief Returns an M x N matrix of random doubles between 0 and 1.0
	 */
	matrix<double> rand(int M, int N)
	{
	    matrix<double> A(M,N);
	    for(int ii = 0; ii < M; ii++)
	    {
	        for(int jj = 0; jj < N; jj++)
	        {
	            A(ii,jj) = ((double)std::rand()/((double)RAND_MAX));
	        }
        }
	    
		return A;
	}

	/** \brief Generalized complex-valued eigenvalue solver using LAPACK function call. 
	 *  
	 *  This function returns the eigenvalues(lambda) of the complex-valued generalized
	 *  eigenvalue problem: Ax_r = lambda*Bx_r or x_l^T*A = lambda*x_l^T*B. The eigenvalues
	 *  are returned by default. To return the right or left eigenvectors, supply the
	 *  function with a complex<double> matrix object in the 3rd or 4th parameters, respectively.
	 */
    vector<complex<double> > eig(const matrix<complex<double> > A, const matrix<complex<double> > B, matrix<complex<double> > *vr_return, matrix<complex<double> > *vl_return)
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
	
	/** \brief Complex-valued eigenvalue solver using LAPACK function call. 
	 *  
	 *  This function returns the eigenvalues(lambda) of the complex-valued
	 *  eigenvalue problem: Ax_r = lambda*x_r or x_l^T*A = lambda*x_l^T. The eigenvalues
	 *  are returned by default. To return the right or left eigenvectors, supply the
	 *  function with a complex<double> matrix object in the 2nd or 3rd parameters, respectively.
	 */
    vector<complex<double> > eig(const matrix<complex<double> > A, matrix<complex<double> > *vr_return, matrix<complex<double> > *vl_return)
	{
		int n, lda, ldb, ldvl, ldvr, lwork, info;
		n = lda = ldb = A.size(1);
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

		complex<double> *a = new complex<double>[n*n];
		complex<double> *vr = new complex<double>[n*n];
		complex<double> *vl = new complex<double>[n*n];
		complex<double> *w = new complex<double>[n];
		complex<double> *work = new complex<double>[lwork];
		double *rwork = new double[lwork];
		for(int ii = 0; ii < n; ii++)
		{
			for(int jj = 0; jj < n; jj++)
			{
				a[ii*n + jj] = A(jj,ii);
			}
		}

	    zgeev_(&jobvl, &jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);

		vector<complex<double> > lambda(n);
		for(int ii = 0; ii < n; ii++)
		{
			lambda[ii] = w[ii];
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
	 *  function with a complex<double> matrix object in the 2nd or 3rd parameters, respectively.
	 */
    vector<complex<double> > eig(const matrix<double> A, matrix<complex<double> > *vr_return, matrix<complex<double> > *vl_return)
	{
	    matrix<complex<double>> B(A.size(1),A.size(2));
	    for(int ii = 0; ii < B.size(1); ii++)
	    {
	        for(int jj = 0; jj < B.size(2); jj++)
	        {
	            B(ii,jj) = A(ii,jj);
	        }
	    }
	    return eig(B, vr_return, vl_return);
	}
	
	double rcond(const matrix<double> A)
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
        for(int ii = 0; ii < A.size(2); ii++)
        {
            for(int jj = 0; jj < A.size(1); jj++)
            {
                x[ii*A.size(1) + jj] = A(jj,ii);
            }
        }
        n = A.size(1);
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
	
	double rcond(const matrix<std::complex<double>> A)
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
        for(int ii = 0; ii < A.size(2); ii++)
        {
            for(int jj = 0; jj < A.size(1); jj++)
            {
                x[ii*A.size(1) + jj] = A(jj,ii);
            }
        }
        n = A.size(1);
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
	
	std::vector<complex<double>> linsolve(const matrix<complex<double>>& A_in,
	                                      const std::vector<complex<double>>& b_in)
	{
		if(b_in.empty() || A_in.size(1) <= 0 || A_in.size(2) <= 0)
		{
			throw KeyCppException("Error in linsolve()! Empty matrix or vector supplied!");
		}
		if(A_in.size(2) != b_in.size())
		{
			throw KeyCppException("Error in linsolve()! Matrix and vector sizes are not compatible!");
		}
		
		int n = b_in.size();
		int m = A_in.size(2);
		int nrhs = 1;

        int info, lda;
        double anorm, rcond;
        
        int *iw = new int[A_in.size(1)];
        std::complex<double> *w1 = new std::complex<double>[A_in.size(1)*A_in.size(2) + 64];
        double *w2 = new double[A_in.size(1)*A_in.size(2) + 64];
        std::complex<double> *A = new std::complex<double>[A_in.size(1)*A_in.size(2)];
        for(int ii = 0; ii < A_in.size(2); ii++)
        {
            for(int jj = 0; jj < A_in.size(1); jj++)
            {
                A[ii*A_in.size(1) + jj] = A_in(jj,ii);
            }
        }
        lda = n;

        /* Computes the norm of A */
        anorm = zlange_("1", &n, &n, A, &lda, w2);

        /* Modifies A in place with a LU decomposition */
        zgetrf_(&n, &n, A, &lda, iw, &info);
        if(info != 0)
        {
            if(info > 0)
            {
                cerr << "Warning: Matrix is singular. Results may be inaccurate.\n";
            }
            else
            {
                throw KeyCppException("Unknown error in linsolve()!");
            }
        }

        /* Computes the reciprocal norm */
        zgecon_("1", &n, A, &lda, &anorm, &rcond, w1, w2, &info);
        if(info != 0)
        {
            throw KeyCppException("Unknown error in linsolve()!");
        }
        
        if(rcond < 1e-15)
        {
            cerr << "Warning: Matrix is close to singular or badly scaled. Results may be inaccurate.\nrcond = " << rcond << endl;
        }
        
        complex<double> *b = new complex<double>[n];
        for(int ii = 0; ii < n; ii++)
        {
            b[ii] = b_in[ii];
        }
        zgetrs_("N", &m, &nrhs, A, &lda, iw, b, &n, &info);
        if(info != 0)
        {
            throw KeyCppException("Unknown error in linsolve()!");
        }
        
        vector<complex<double>> x_out(n);
        for(int ii = 0; ii < n; ii++)
        {
            x_out[ii] = b[ii];
        }
        
        delete [] iw;
        delete [] w1;
        delete [] w2;
        delete [] A;
        delete [] b;
        
		return x_out;
	}
	
	
	std::vector<double> linsolve(const matrix<double>& A_in,
	                                      const std::vector<double>& b_in)
	{
		if(b_in.empty() || A_in.size(1) <= 0 || A_in.size(2) <= 0)
		{
			throw KeyCppException("Error in linsolve()! Empty matrix or vector supplied!");
		}
		if(A_in.size(2) != b_in.size())
		{
			throw KeyCppException("Error in linsolve()! Matrix and vector sizes are not compatible!");
		}
		
		int n = b_in.size();
		int m = A_in.size(2);
		int nrhs = 1;

        int info, lda;
        double anorm, rcond;
        
        int *iw = new int[A_in.size(1)];
        double *w1 = new double[A_in.size(1)*A_in.size(2) + 64];
        int *w2 = new int[A_in.size(1)*A_in.size(2) + 64];
        double *A = new double[A_in.size(1)*A_in.size(2)];
        for(int ii = 0; ii < A_in.size(2); ii++)
        {
            for(int jj = 0; jj < A_in.size(1); jj++)
            {
                A[ii*A_in.size(1) + jj] = A_in(jj,ii);
            }
        }
        lda = n;

        /* Computes the norm of A */
        anorm = dlange_("1", &n, &n, A, &lda, w1);

        /* Modifies A in place with a LU decomposition */
        dgetrf_(&n, &n, A, &lda, iw, &info);
        if(info != 0)
        {
            if(info > 0)
            {
                cerr << "Warning: Matrix is singular. Results may be inaccurate.\n";
            }
            else
            {
                throw KeyCppException("Unknown error in linsolve()!");
            }
        }

        /* Computes the reciprocal norm */
        dgecon_("1", &n, A, &lda, &anorm, &rcond, w1, w2, &info);
        if(info != 0)
        {
            throw KeyCppException("Unknown error in linsolve()!");
        }
        
        if(rcond < 1e-15)
        {
            cerr << "Warning: Matrix is close to singular or badly scaled. Results may be inaccurate.\nrcond = " << rcond << endl;
        }
        
        double *b = new double[n];
        for(int ii = 0; ii < n; ii++)
        {
            b[ii] = b_in[ii];
        }
        
        dgetrs_("N", &m, &nrhs, A, &lda, iw, b, &n, &info);
        if(info != 0)
        {
            throw KeyCppException("Unknown error in linsolve()!");
        }
        
        vector<double> x_out(n);
        for(int ii = 0; ii < n; ii++)
        {
            x_out[ii] = b[ii];
        }
        
        delete [] iw;
        delete [] w1;
        delete [] w2;
        delete [] A;
        delete [] b;
        
		return x_out;
	}
	
	matrix<double> inv(const matrix<double>& A_in)
	{
	    if(A_in.size(1) <= 0 || A_in.size(2) <= 0)
		{
			throw KeyCppException("Error in inv()! Empty matrix supplied!");
		}
		if(A_in.size(1) != A_in.size(2))
		{
		    throw KeyCppException("Error in inv()! Matrix must be square!");
		}
		
		int n = A_in.size(1);
		int m = A_in.size(2);

        int info, lda;
        double anorm, rcond;
        
        int *iw = new int[A_in.size(1)];
        int lwork = A_in.size(1)*A_in.size(2) + 64;
        double *w1 = new double[lwork];
        int *w2 = new int[lwork];
        double *A = new double[A_in.size(1)*A_in.size(2)];
        for(int ii = 0; ii < A_in.size(2); ii++)
        {
            for(int jj = 0; jj < A_in.size(1); jj++)
            {
                A[ii*A_in.size(1) + jj] = A_in(jj,ii);
            }
        }
        lda = n;

        /* Computes the norm of A */
        anorm = dlange_("1", &n, &n, A, &lda, w1);

        /* Modifies A in place with a LU decomposition */
        dgetrf_(&n, &n, A, &lda, iw, &info);
        if(info != 0)
        {
            throw KeyCppException("Unknown error in inv()!");
        }

        /* Computes the reciprocal norm */
        dgecon_("1", &n, A, &lda, &anorm, &rcond, w1, w2, &info);
        if(info != 0)
        {
            throw KeyCppException("Unknown error in inv()!");
        }
        
        if(rcond < 1e-15)
        {
            cerr << "Warning: Matrix is close to singular or badly scaled. Results may be inaccurate.\nrcond = " << rcond << endl;
        }

        dgetri_(&n,A,&n,iw,w1,&lwork,&info);

        matrix<double> A_out(A_in.size(1),A_in.size(2));
        for(int ii = 0; ii < n; ii++)
        {
            for(int jj = 0; jj < n; jj++)
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
    
    matrix<complex<double>> inv(const matrix<complex<double>>& A_in)
	{
		if(A_in.size(1) <= 0 || A_in.size(2) <= 0)
		{
			throw KeyCppException("Error in inv()! Empty matrix supplied!");
		}
		if(A_in.size(1) != A_in.size(2))
		{
			throw KeyCppException("Error in inv()! Matrix must be square!");
		}
		
		int n = A_in.size(1);
		int m = A_in.size(2);
		int nrhs = 1;

        int info, lda;
        double anorm, rcond;
        
        int *iw = new int[A_in.size(1)];
        int lwork = A_in.size(1)*A_in.size(2) + 64;
        std::complex<double> *w1 = new std::complex<double>[lwork];
        double *w2 = new double[lwork];
        std::complex<double> *A = new std::complex<double>[A_in.size(1)*A_in.size(2)];
        for(int ii = 0; ii < A_in.size(2); ii++)
        {
            for(int jj = 0; jj < A_in.size(1); jj++)
            {
                A[ii*A_in.size(1) + jj] = A_in(jj,ii);
            }
        }
        lda = n;

        /* Computes the norm of A */
        anorm = zlange_("1", &n, &n, A, &lda, w2);

        /* Modifies A in place with a LU decomposition */
        zgetrf_(&n, &n, A, &lda, iw, &info);
        if(info != 0)
        {
            throw KeyCppException("Unknown error in inv()!");
        }

        /* Computes the reciprocal norm */
        zgecon_("1", &n, A, &lda, &anorm, &rcond, w1, w2, &info);
        if(info != 0)
        {
            throw KeyCppException("Unknown error in inv()!");
        }
        
        if(rcond < 1e-15)
        {
            cerr << "Warning: Matrix is close to singular or badly scaled. Results may be inaccurate.\nrcond = " << rcond << endl;
        }
        
        zgetri_(&n,A,&n,iw,w1,&lwork,&info);

        matrix<complex<double>> A_out(A_in.size(1),A_in.size(2));
        for(int ii = 0; ii < n; ii++)
        {
            for(int jj = 0; jj < n; jj++)
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
	
	double norm(const matrix<double> A_in, std::string method)
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
	        int m = A_in.size(1);
	        int n = A_in.size(2);

            int info, lda;
            int lwork = A_in.size(1)*A_in.size(2) + 64;
            double *w1 = new double[lwork];
            double *A = new double[A_in.size(1)*A_in.size(2)];
            for(int ii = 0; ii < A_in.size(2); ii++)
            {
                for(int jj = 0; jj < A_in.size(1); jj++)
                {
                    A[ii*A_in.size(1) + jj] = A_in(jj,ii);
                }
            }
            lda = m;
            anorm = dlange_(method.c_str(), &m, &n, A, &lda, w1);
            
            delete [] w1;
            delete [] A;
        }
        
        return anorm;
	}
	
	/** \brief Computes the singular value decomposition of matrix A_in.
	 */
	SVD_type<double,double> svd(matrix<double> A_in, std::string method)
	{
	    if(A_in.size(1) <= 0 || A_in.size(2) <= 0)
		{
			throw KeyCppException("Error in svd()! Empty matrix supplied!");
		}
		std::transform(method.begin(), method.end(), method.begin(), ::tolower);
		
	    int m = A_in.size(1);
	    int n = A_in.size(2);

        int info, lda, ldu, ldvt;
        double anorm;
        
        int col_u;
        
        SVD_type<double,double> out;
        int lwork = A_in.size(1)*A_in.size(2) + 64;
        double *work = new double[lwork];
        double *A = new double[A_in.size(1)*A_in.size(2)];
        double *U;
        double *S;
        double *VT;
        int min_dim;
        if(m > n)
        {
            min_dim = n;
        }
        else
        {
            min_dim = m;
        }
        S = new double[min_dim];
        for(int ii = 0; ii < A_in.size(2); ii++)
        {
            for(int jj = 0; jj < A_in.size(1); jj++)
            {
                A[ii*A_in.size(1) + jj] = A_in(jj,ii);
            }
        }
        lda = m;
        ldu = m;
        string jobu, jobvt;
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
		
		dgesvd_(jobu.c_str(), jobvt.c_str(), &m, &n, A, &lda, S, U, &ldu, VT, &ldvt, work, &lwork, &info);
                 
        if(info != 0)
        {
            throw KeyCppException("Unknown error in SVD()!");
        }
        
        out.U = matrix<double>(m,col_u);
        out.V = matrix<double>(n,ldvt);
        
        for(int ii = 0; ii < min_dim; ii++)
        {
            out.S(ii,ii) = S[ii];
        }
        for(int ii = 0; ii < m; ii++)
        {
            for(int jj = 0; jj < col_u; jj++)
            {
                out.U(ii,jj) = U[jj*m + ii];
            }
        }
        
        for(int ii = 0; ii < n; ii++)
        {
            for(int jj = 0; jj < ldvt; jj++)
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
	
	double norm(const matrix<complex<double>> A_in, std::string method)
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
	        int m = A_in.size(1);
	        int n = A_in.size(2);

            int info, lda;
            int lwork = A_in.size(1)*A_in.size(2) + 64;
            double *w1 = new double[lwork];
            complex<double> *A = new complex<double>[A_in.size(1)*A_in.size(2)];
            for(int ii = 0; ii < A_in.size(2); ii++)
            {
                for(int jj = 0; jj < A_in.size(1); jj++)
                {
                    A[ii*A_in.size(1) + jj] = A_in(jj,ii);
                }
            }
            lda = m;
            anorm = zlange_(method.c_str(), &m, &n, A, &lda, w1);
            
            delete [] w1;
            delete [] A;
        }
        
        return anorm;
	}
	
	/** \brief Computes the singular value decomposition of matrix A_in.
	 */
	SVD_type<complex<double>, double> svd(matrix<complex<double>> A_in, std::string method)
	{
	    if(A_in.size(1) <= 0 || A_in.size(2) <= 0)
		{
			throw KeyCppException("Error in svd()! Empty matrix supplied!");
		}
		std::transform(method.begin(), method.end(), method.begin(), ::tolower);
		
	    int m = A_in.size(1);
	    int n = A_in.size(2);

        int info, lda, ldu, ldvt;
        double anorm;
        
        int col_u;
        
        SVD_type<std::complex<double>,double> out;
        int lwork = A_in.size(1)*A_in.size(2) + 64;
        double *rwork = new double[lwork];
        complex<double> *work = new complex<double>[lwork];
        complex<double> *A = new complex<double>[A_in.size(1)*A_in.size(2)];
        complex<double> *U;
        double *S;
        complex<double> *VT;
        int min_dim;
        if(m > n)
        {
            min_dim = n;
        }
        else
        {
            min_dim = m;
        }
        S = new double[min_dim];
        for(int ii = 0; ii < A_in.size(2); ii++)
        {
            for(int jj = 0; jj < A_in.size(1); jj++)
            {
                A[ii*A_in.size(1) + jj] = A_in(jj,ii);
            }
        }
        lda = m;
        ldu = m;
        string jobu, jobvt;
        if(method.compare("0") == 0)
        {
            jobu = "S";
            jobvt = "A";
            U = new complex<double>[m*min_dim];
            col_u = min_dim;
            VT = new complex<double>[n*n];
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
            U = new complex<double>[m*min_dim];
            col_u = min_dim;
            VT = new complex<double>[n*min_dim];
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
            U = new complex<double>[m*m];
            col_u = m;
            VT = new complex<double>[n*n];
            ldvt = n;
            out.S = matrix<double>(m,n);
        }
        else
        {
            throw KeyCppException("Unknown argument to svd()!");
        }
		
		zgesvd_(jobu.c_str(), jobvt.c_str(), &m, &n, A, &lda, S, U, &ldu, VT, &ldvt, work, &lwork, rwork, &info);
                 
        if(info != 0)
        {
            throw KeyCppException("Unknown error in SVD()!");
        }
        
        out.U = matrix<complex<double>>(m,col_u);
        out.V = matrix<complex<double>>(n,ldvt);
        
        for(int ii = 0; ii < min_dim; ii++)
        {
            out.S(ii,ii) = S[ii];
        }
        for(int ii = 0; ii < m; ii++)
        {
            for(int jj = 0; jj < col_u; jj++)
            {
                out.U(ii,jj) = U[jj*m + ii];
            }
        }
        
        for(int ii = 0; ii < n; ii++)
        {
            for(int jj = 0; jj < ldvt; jj++)
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

