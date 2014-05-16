#include <cmath>
#include <complex>
#include <vector>
#include <sstream>
#include "Matrix.h"

namespace keycpp
{
	class eigs_opt
	{
	    public:
	        int maxit{300};
	        matrix<std::complex<double>> v0{};
	};
	
	class znaupdException : public std::runtime_error
	{
		public:
		    int code;
			znaupdException(int ii) : std::runtime_error(""), code(ii), m_what("")
			{
                switch(code)
                {
                    case 1:
                        m_what = "\nMaximum number of iterations taken. All possible eigenvalues of OP has been found. IPARAM(5) returns the number of wanted converged Ritz values.";
                        break;
                    case 2:
                        m_what = "\nNo longer an informational error. Deprecated starting with release 2 of ARPACK.";
                        break;
                    case 3:
                        m_what = "\nNo shifts could be applied during a cycle of the Implicitly restarted Arnoldi iteration. One possibility is to increase the size of NCV relative to NEV.";
                        break;
                    case -1:
                        m_what = "\nN must be positive.";
                        break;
                    case -2:
                        m_what = "\nNEV must be positive.";
                        break;
                    case -3:
                        m_what = "\nNCV-NEV >= 2 and less than or equal to N.";
                        break;
                    case -4:
                        m_what = "\nThe maximum number of Arnoldi update iteration must be greater than zero.";
                        break;
                    case -5:
                        m_what = "\nWHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'";
                        break;
                    case -6:
                        m_what = "\nBMAT must be one of 'I' or 'G'.";
                        break;
                    case -7:
                        m_what = "\nLength of private work array is not sufficient.";
                        break;
                    case -8:
                        m_what = "\nError return from LAPACK eigenvalue calculation;";
                        break;
                    case -9:
                        m_what = "\nStarting vector is zero.";
                        break;
                    case -10:
                        m_what = "\nIPARAM(7) must be 1,2,3.";
                        break;
                    case -11:
                        m_what = "\nIPARAM(7) = 1 and BMAT = 'G' are incompatable.";
                        break;
                    case -12:
                        m_what = "\nIPARAM(1) must be equal to 0 or 1.";
                        break;
                    case -9999:
                        m_what = "\nCould not build an Arnoldi factorization. User input error highly likely.  Please check actual array dimensions and layout. IPARAM(5) returns the size of the current Arnoldi factorization.";
                        break;
                    default:
                        std::stringstream sstm;
                        sstm << "\nUnknown error in znaupd with info = " << code;
                        m_what = sstm.str();
                        break;
                }
			}
			
            virtual ~znaupdException() throw () {}
			
	        const char* what() const throw()
	        {
		        return m_what.c_str();
	        }
	        
	    private:
		    std::string m_what;
	};
	
	class zneupdException : public std::runtime_error
	{
		public:
		    int code;
			zneupdException(int ii) : std::runtime_error(""), code(ii), m_what("")
			{
                switch(code)
                {
                    case 1:
                        m_what = "\nThe Schur form computed by LAPACK routine csheqr could not be reordered by LAPACK routine ztrsen. Re-enter subroutine ZNEUPD with IPARAM(5)=NCV and increase the size of the array D to have dimension at least dimension NCV and allocate at least NCV columns for Z. NOTE: Not necessary if Z and V share the same space. Please notify the authors if this error occurs.";
                        break;
                    case -1:
                        m_what = "\nN must be positive.";
                        break;
                    case -2:
                        m_what = "\nNEV must be positive.";
                        break;
                    case -3:
                        m_what = "\nNCV-NEV >= 2 and less than or equal to N.";
                        break;
                    case -5:
                        m_what = "\nWHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'";
                        break;
                    case -6:
                        m_what = "\nBMAT must be one of 'I' or 'G'.";
                        break;
                    case -7:
                        m_what = "\nLength of private work WORKL array is not sufficient.";
                        break;
                    case -8:
                        m_what = "\nError return from LAPACK eigenvalue calculation. This should never happened.";
                        break;
                    case -9:
                        m_what = "\nError return from calculation of eigenvectors. Informational error from LAPACK routine ztrevc.";
                        break;
                    case -10:
                        m_what = "\nIPARAM(7) must be 1,2,3";
                        break;
                    case -11:
                        m_what = "\nIPARAM(7) = 1 and BMAT = 'G' are incompatible.";
                        break;
                    case -12:
                        m_what = "\nHOWMNY = 'S' not yet implemented";
                        break;
                    case -13:
                        m_what = "\nHOWMNY must be one of 'A' or 'P' if RVEC = .true.";
                        break;
                    case -14:
                        m_what = "\nZNAUPD did not find any eigenvalues to sufficient accuracy.";
                        break;
                    default:
                        std::stringstream sstm;
                        sstm << "\nUnknown error in zneupd with info = " << code;
                        m_what = sstm.str();
                        break;
                }
			}
			
            virtual ~zneupdException() throw () {}
			
	        const char* what() const throw()
	        {
		        return m_what.c_str();
	        }
	        
	    private:
		    std::string m_what;
	};

    extern "C" void znaupd_(int *ido, char *bmat, int *n, char *which,
			    int *nev, double *tol, std::complex<double> *resid,
			    int *ncv, std::complex<double> *v, int *ldv,
			    int *iparam, int *ipntr, std::complex<double> *workd,
			    std::complex<double> *workl, int *lworkl,
			    double *rwork, int *info);

    extern "C" void zneupd_(int *rvec, char *All, int *select,
			    std::complex<double> *d, std::complex<double> *v, int *ldv,
			    std::complex<double> *sigma, std::complex<double> *workev, char *bmat,
			    int *n, char *which, int *nev, double *tol,
			    std::complex<double> *resid, int *ncv,
			    std::complex<double> *, int *, int *iparam,
			    int *ipntr, std::complex<double> *workd,
			    std::complex<double> *workl, int *lworkl,
			    double *rwork, int *ierr);
			    

    inline void mv(int n, std::complex<double> *in, std::complex<double> *out, const matrix<std::complex<double>> &A)
    {
	    int m = n;
        char TRANS = 'N';
        std::complex<double> ALPHA = 1.0;
        int LDA = m;
        int INCX = 1;
        std::complex<double> BETA = 0.0;
        int INCb = 1;
        
        zgemv_(&TRANS, &m, &n, &ALPHA, &A.mData[0], &LDA, &in[0],&INCX, &BETA, &out[0], &INCb);
    }

    inline void mv_special(int n, std::complex<double> *in, std::complex<double> *out, const matrix<std::complex<double>> &A, matrix<std::complex<double>> &Y, int *iw)
    {
	    int m = n;
        char TRANS = 'N';
        std::complex<double> ALPHA = 1.0;
        int LDA = m;
        int INCX = 1;
        std::complex<double> BETA = 0.0;
        int INCb = 1;
        
        zgemv_(&TRANS, &m, &n, &ALPHA, &A.mData[0], &LDA, &in[0],&INCX, &BETA, &out[0], &INCb);
        
        int nrhs = 1;
        int mm = Y.size(2);
        int nn = Y.size(1);
        int lda = nn;
        
        int info2;

        zgetrs_("N", &mm, &nrhs, &Y.mData[0], &lda, iw, &out[0], &nn, &info2);
        
        if(info2 != 0)
        {
            throw KeyCppException("Unknown error in linsolve()!");
        }
    }

    template<int type>
    void znaupd(int n, int nev, matrix<std::complex<double>> &Evals, std::string which, const matrix<std::complex<double>,2,type> &A, eigs_opt *opt_base)
    {
        eigs_opt opt;
        
        if(opt_base == NULL)
        {
            opt = eigs_opt();
        }
        else
        {
            opt = *opt_base;
        }
        
        int ido = 0;

        char bmat[2] = "I";
        std::vector<char> which_char(which.c_str(), which.c_str() + which.size() + 1u);

        double tol = 0.0;

        std::complex<double> *resid;
        resid = new std::complex<double>[n];

        int ncv = 2*nev + 1;
        if(ncv>=n)
        {
            ncv = n;
        }

        std::complex<double> *v;
        int ldv = n;
        v = new std::complex<double>[ldv*ncv];

        int *iparam;
        iparam = new int[11];
        iparam[0] = 1;
        iparam[2] = opt.maxit;
        iparam[6] = 1;

        int *ipntr;
        ipntr = new int[14];

        std::complex<double> *workd;
        workd = new std::complex<double>[3*n];

        int lworkl = 3*ncv*ncv + 5*ncv;
        std::complex<double> *workl;
        workl = new std::complex<double>[lworkl];

        double *rwork;
        rwork = new double[ncv];
        

        int info = 0;
        if(opt.v0.numel() == (size_t)n)
        {
            info = 1;
            for(int ii = 0; ii < n; ii++)
            {
                resid[ii] = opt.v0(ii);
            }
        }

        int rvec = 0;

        int *select;
        select = new int[ncv];
        std::complex<double> *d;
        d = new std::complex<double>[ncv];
        std::complex<double> sigma;

        std::complex<double> *workev;
        workev = new std::complex<double>[3*ncv];

        int ierr;

        do
        {
            znaupd_(&ido, bmat, &n, &which_char[0], &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);

            if((ido==1)||(ido==-1))
            {
                mv(n, workd+ipntr[0]-1, workd+ipntr[1]-1,A);
            }
        } while((ido==1)||(ido==-1));


        if(info<0)
        {
            if(info == -5)
            {
                throw KeyCppException("Invalid sigma provided to eigs!");
            }
            std::cout << "Error with znaupd, info = " << info << "\n";
            std::cout << "Check documentation in znaupd\n\n";
        }
        else
        {
            char howmny[] = "All";
            zneupd_(&rvec, howmny, select, d, v, &ldv, &sigma, workev,
            bmat, &n, &which_char[0], &nev, &tol, resid, &ncv, v, &ldv,
            iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);

            if(ierr!=0)
            {
                std::cout << "Error with zneupd, info = " << ierr << "\n";
                std::cout << "Check the documentation of zneupd.\n\n";
            }
            else if(info==1)
            {
                std::cout << "Maximum number of iterations reached.\n\n";
            }else if(info==3)
            {
                std::cout << "No shifts could be applied during implicit\n";
                std::cout << "Arnoldi update, try increasing NCV.\n\n";
            }


            Evals = matrix<std::complex<double>>(nev,nev);
            for(int i=0; i<nev; i++)
            {
                Evals(i,i) = d[nev-1-i];
            }

            // Sort the eigenvalues by real part

            std::complex<double> temp;
            for(int i=0; i<nev; i++)
            {
                for(int j=i; j<nev; j++)
                {
                    if(Evals(j,j).real() > Evals(i,i).real())
                    {
                        temp = Evals(j,j);
                        Evals(j,j) = Evals(i,i);
                        Evals(i,i) = temp;
                    }
                }
            }
            
            opt.v0 = matrix<std::complex<double>>(n);
            for(int ii = 0; ii < n; ii++)
            {
                opt.v0(ii) = resid[ii];
            }
            
            delete [] rwork;
            delete [] workev;
            delete [] resid;
            delete [] v;
            delete [] iparam;
            delete [] ipntr;
            delete [] workd;
            delete [] workl;
            delete [] select;
            delete [] d;
        }
    }

    template<int type>
    void znaupd(int n, int nev, matrix<std::complex<double>> &Evals, matrix<std::complex<double>> &Evecs, std::string which, const matrix<std::complex<double>,2,type> &A, eigs_opt *opt_base)
    {
        eigs_opt opt;
        
        if(opt_base == NULL)
        {
            opt = eigs_opt();
        }
        else
        {
            opt = *opt_base;
        }
        
        int ido = 0;
        char bmat[2] = "I";
        //char which[3] = "SM";
        std::vector<char> which_char(which.c_str(), which.c_str() + which.size() + 1u);
        double tol = 0.0;
        std::complex<double> *resid;
        resid = new std::complex<double>[n];
        int ncv = 2*nev + 1;
        if (ncv>n)
        {
            ncv = n;
        }
        std::complex<double> *v;
        int ldv = n;
        v = new std::complex<double>[ldv*ncv];

        int *iparam;
        iparam = new int[11];
        iparam[0] = 1;
        iparam[2] = opt.maxit;
        iparam[6] = 1;
        int *ipntr;
        ipntr = new int[14];
        std::complex<double> *workd;
        workd = new std::complex<double>[3*n];
        int lworkl = 3*ncv*ncv + 5*ncv;
        std::complex<double> *workl;
        workl = new std::complex<double>[lworkl];
        double *rwork;
        rwork = new double[ncv];
        int info = 0;
        if(opt.v0.numel() == (size_t)n)
        {
            info = 1;
            for(int ii = 0; ii < n; ii++)
            {
                resid[ii] = opt.v0(ii);
            }
        }
        int rvec = 1;  // Changed from above
        int *select;
        select = new int[ncv];
        std::complex<double> *d;
        d = new std::complex<double>[ncv];
        std::complex<double> sigma;
        std::complex<double> *workev;
        workev = new std::complex<double>[3*ncv];
        int ierr;

        do
        {
            znaupd_(&ido, bmat, &n, &which_char[0], &nev, &tol, resid, 
            &ncv, v, &ldv, iparam, ipntr, workd, workl,
            &lworkl, rwork, &info);

            if((ido==1)||(ido==-1))
            {
                mv(n, workd+ipntr[0]-1, workd+ipntr[1]-1,A);
            }
        } while((ido==1)||(ido==-1));

        if(info<0)
        {
            if(info == -5)
            {
                throw KeyCppException("Invalid sigma provided to eigs!");
            }
            std::cout << "Error with znaupd, info = " << info << "\n";
            std::cout << "Check documentation in znaupd\n\n";
        }
        else
        {
            char howmny[] = "All";
            zneupd_(&rvec, howmny, select, d, v, &ldv, &sigma, workev,
            bmat, &n, &which_char[0], &nev, &tol, resid, &ncv, v, &ldv,
            iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);

            if(ierr!=0)
            {
                std::cout << "Error with zneupd, info = " << ierr << "\n";
                std::cout << "Check the documentation of zneupd.\n\n";
            }
            else if(info==1)
            {
                std::cout << "Maximum number of iterations reached.\n\n";
            }
            else if(info==3)
            {
                std::cout << "No shifts could be applied during implicit\n";
                std::cout << "Arnoldi update, try increasing NCV.\n\n";
            }

            Evals = matrix<std::complex<double>>(nev,nev);
            for(int i=0; i<nev; i++)
            {
                Evals(i,i) = d[i];
            }
            
            Evecs = matrix<std::complex<double>>(n,nev);
            for(int i=0; i<nev; i++)
            {
                for(int j=0; j<n; j++)
                {
                    Evecs(j,i) = v[i*n+j];
                }
            }

            std::complex<double> temp;
            for(int i=0; i<nev; i++)
            {
                for(int j=i; j<nev; j++)
                {
                    if(Evals(j,j).real() > Evals(i,i).real())
                    {
                        temp = Evals(j,j);
                        Evals(j,j) = Evals(i,i);
                        Evals(i,i) = temp;
                        for(int k=0; k<n; k++)
                        {
                            temp = Evecs(k,i);
                            Evecs(k,i) = Evecs(k,j);
                            Evecs(k,j) = temp;
                        }
                    }
                }
            }
            
            opt.v0 = matrix<std::complex<double>>(n);
            for(int ii = 0; ii < n; ii++)
            {
                opt.v0(ii) = resid[ii];
            }

            delete [] rwork;
            delete [] workev;
            delete [] resid;
            delete [] v;
            delete [] iparam;
            delete [] ipntr;
            delete [] workd;
            delete [] workl;
            delete [] select;
            delete [] d;
        }
    }

    inline matrix<std::complex<double>,2> eigs(const matrix<std::complex<double>,2,DENSE_MATRIX> &A, size_t k = 6, std::string sigma = "LM", matrix<std::complex<double>,2> *vr_return = NULL, eigs_opt *opt = NULL)
    {
        if(A.empty())
        {
            throw KeyCppException("Empty matrix in eigs!");
        }
        if(A.size(1) != A.size(2))
        {
            throw KeyCppException("Non-square matrix in eigs!");
        }
        
        if((A.size(1) - k) < 2)
        {
            if(k > A.size(1))
            {
                throw KeyCppException("Too many eigenvalues requested!");
            }
            else
            {
                return eig(A,vr_return);
            }
        }
        
        matrix<std::complex<double>,2> eig_val;
        if(vr_return == NULL)
        {
            int n = A.size(1);
            znaupd(n, k, eig_val, sigma, A, opt);
        }
        else
        {
            int n = A.size(1);
            znaupd(n, k, eig_val, *vr_return, sigma, A, opt);
        }
        return eig_val;
    }

    inline void znaupd_shift_invert(int n, int nev, matrix<std::complex<double>> &Evals, std::complex<double> sigma, const matrix<std::complex<double>> &A, eigs_opt *opt_base)
    {
        eigs_opt opt;
        
        if(opt_base == NULL)
        {
            opt = eigs_opt();
        }
        else
        {
            opt = *opt_base;
        }
        
        int ido = 0;

        char bmat[2] = "I";
        std::string which = "LM";
        std::vector<char> which_char(which.c_str(), which.c_str() + which.size() + 1u);
        
        int *iw = new int[A.size(1)];
        matrix<std::complex<double>> Y = lu(A-sigma*eye<std::complex<double>>(n),iw);

        double tol = 0.0; /* Sets the tolerance; tol<=0 specifies 
        machine precision */

        std::complex<double> *resid;
        resid = new std::complex<double>[n];

        int ncv = 2*nev + 1; /* The largest number of basis vectors that will
        be used in the Implicitly Restarted Arnoldi
        Process.  Work per major iteration is
        proportional to N*NCV*NCV. */
        if(ncv>=n)
        {
            ncv = n;
        }

        std::complex<double> *v;
        int ldv = n;
        v = new std::complex<double>[ldv*ncv];

        int *iparam;
        iparam = new int[11]; /* An array used to pass information to the routines
        about their functional modes. */
        iparam[0] = 1;   // Specifies the shift strategy (1->exact)
        iparam[2] = opt.maxit; // Maximum number of iterations
        iparam[6] = 3;   /* Sets the mode of znaupd.
        1 is exact shifting,
        2 is user-supplied shifts,
        3 is shift-invert mode,
        4 is buckling mode,
        5 is Cayley mode. */

        int *ipntr;
        ipntr = new int[14]; /* Indicates the locations in the work array workd
        where the input and output vectors in the
        callback routine are located. */

        std::complex<double> *workd;
        workd = new std::complex<double>[3*n];

        int lworkl = 3*ncv*ncv + 5*ncv; /* Length of the workl array */
        std::complex<double> *workl;
        workl = new std::complex<double>[lworkl];

        double *rwork;
        rwork = new double[ncv];

        int info = 0; /* Passes convergence information out of the iteration
        routine. */
        if(opt.v0.numel() == (size_t)n)
        {
            info = 1;
            for(int ii = 0; ii < n; ii++)
            {
                resid[ii] = opt.v0(ii);
            }
        }

        int rvec = 0; /* Specifies that eigenvectors should not be calculated */

        int *select;
        select = new int[ncv];
        std::complex<double> *d;
        d = new std::complex<double>[ncv]; /* This vector will return the
         eigenvalues from the second routine,
         dseupd. */

        std::complex<double> *workev;
        workev = new std::complex<double>[3*ncv]; // I don't know what this is used for

        int ierr;

        /* Here we enter the main loop where the calculations are
        performed.  The communication parameter ido tells us when
        the desired tolerance is reached, and at that point we exit
        and extract the solutions. */

        do
        {
            znaupd_(&ido, bmat, &n, &which_char[0], &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);

            if((ido == -1) || (ido == 1))
            {
                int nrhs = 1;
                int mm = Y.size(2);
                int nn = Y.size(1);
                int lda = nn;
                
                int info2;
                
                for(int ii = 0; ii < n; ii++)
                {
                    workd[ipntr[1] - 1 + ii] = workd[ipntr[0] - 1 + ii];
                }
            
                zgetrs_("N", &mm, &nrhs, &Y.mData[0], &lda, iw, &workd[ipntr[1] - 1], &nn, &info2);
                
                if(info2 != 0)
                {
                    throw KeyCppException("Unknown error in linsolve()!");
                }
            }
        } while((ido == 1) || (ido == -1));

        /* From those results, the eigenvalues and vectors are
        extracted. */

        if(info < 0)
        {
            if(info == -5)
            {
                throw KeyCppException("Invalid sigma provided to eigs!");
            }
            std::cout << "Error with znaupd, info = " << info << "\n";
            std::cout << "Check documentation in znaupd\n\n";
        }
        else
        {
            char howmny[] = "All";
            zneupd_(&rvec, howmny, select, d, v, &ldv, &sigma, workev,
            bmat, &n, &which_char[0], &nev, &tol, resid, &ncv, v, &ldv,
            iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);

            if(ierr!=0)
            {
                std::cout << "Error with zneupd, info = " << ierr << "\n";
                std::cout << "Check the documentation of zneupd.\n\n";
            }
            else if(info==1)
            {
                std::cout << "Maximum number of iterations reached.\n\n";
            }else if(info==3)
            {
                std::cout << "No shifts could be applied during implicit\n";
                std::cout << "Arnoldi update, try increasing NCV.\n\n";
            }

            /* Before exiting, we copy the solution information over to
            the arrays of the calling program, then clean up the
            memory used by this routine.  For some reason, when I
            don't find the eigenvectors I need to reverse the order of
            the values. */

            Evals = matrix<std::complex<double>>(nev,nev);
            for(int i=0; i<nev; i++)
            {
                Evals(i,i) = d[nev-1-i];
            }

            // Sort the energies by real part

            std::complex<double> temp;
            for(int i=0; i<nev; i++)
            {
                for(int j=i; j<nev; j++)
                {
                    if(Evals(j,j).real() > Evals(i,i).real())
                    {
                        temp = Evals(j,j);
                        Evals(j,j) = Evals(i,i);
                        Evals(i,i) = temp;
                    }
                }
            }
            
            opt.v0 = matrix<std::complex<double>>(n);
            for(int ii = 0; ii < n; ii++)
            {
                opt.v0(ii) = resid[ii];
            }

            delete [] iw;
            delete [] rwork;
            delete [] workev;
            delete [] resid;
            delete [] v;
            delete [] iparam;
            delete [] ipntr;
            delete [] workd;
            delete [] workl;
            delete [] select;
            delete [] d;
        }
    }

    inline void znaupd_shift_invert(int n, int nev, matrix<std::complex<double>> &Evals, matrix<std::complex<double>> &Evecs, std::complex<double> sigma, const matrix<std::complex<double>> &A, eigs_opt *opt_base)
    {
        eigs_opt opt;
        
        if(opt_base == NULL)
        {
            opt = eigs_opt();
        }
        else
        {
            opt = *opt_base;
        }
        
        int ido = 0;
        char bmat[2] = "I";
        std::string which = "LM";
        std::vector<char> which_char(which.c_str(), which.c_str() + which.size() + 1u);
        
        int *iw = new int[A.size(1)];
        matrix<std::complex<double>> Y = lu(A-sigma*eye<std::complex<double>>(n),iw);
        
        double tol = 0.0;
        std::complex<double> *resid;
        resid = new std::complex<double>[n];
        int ncv = 2*nev + 1;
        if(ncv>n)
        {
            ncv = n;
        }
        std::complex<double> *v;
        int ldv = n;
        v = new std::complex<double>[ldv*ncv];

        int *iparam;
        iparam = new int[11];
        iparam[0] = 1;
        iparam[2] = opt.maxit;
        iparam[6] = 3;
        int *ipntr;
        ipntr = new int[14];
        std::complex<double> *workd;
        workd = new std::complex<double>[3*n];
        int lworkl = 3*ncv*ncv + 5*ncv;
        std::complex<double> *workl;
        workl = new std::complex<double>[lworkl];
        double *rwork;
        rwork = new double[ncv];
        int info = 0;
        if(opt.v0.numel() == (size_t)n)
        {
            info = 1;
            for(int ii = 0; ii < n; ii++)
            {
                resid[ii] = opt.v0(ii);
            }
        }
        int rvec = 1;  // Changed from above
        int *select;
        select = new int[ncv];
        std::complex<double> *d;
        d = new std::complex<double>[ncv];
        std::complex<double> *workev;
        workev = new std::complex<double>[3*ncv];
        int ierr;

        do
        {
            znaupd_(&ido, bmat, &n, &which_char[0], &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);

            if((ido == -1) || (ido == 1))
            {
                int nrhs = 1;
                int mm = Y.size(2);
                int nn = Y.size(1);
                int lda = nn;
                
                int info2;
                
                for(int ii = 0; ii < n; ii++)
                {
                    workd[ipntr[1] - 1 + ii] = workd[ipntr[0] - 1 + ii];
                }
            
                zgetrs_("N", &mm, &nrhs, &Y.mData[0], &lda, iw, &workd[ipntr[1] - 1], &nn, &info2);
                
                if(info2 != 0)
                {
                    throw KeyCppException("Unknown error in linsolve()!");
                }
            }
        } while((ido == 1) || (ido == -1));

        if(info<0)
        {
            if(info == -5)
            {
                throw KeyCppException("Invalid sigma provided to eigs!");
            }
            std::cout << "Error with znaupd, info = " << info << "\n";
            std::cout << "Check documentation in znaupd\n\n";
        }
        else
        {
            char howmny[] = "All";
            zneupd_(&rvec, howmny, select, d, v, &ldv, &sigma, workev,
            bmat, &n, &which_char[0], &nev, &tol, resid, &ncv, v, &ldv,
            iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);

            if(ierr!=0)
            {
                std::cout << "Error with zneupd, info = " << ierr << "\n";
                std::cout << "Check the documentation of zneupd.\n\n";
            }
            else if(info==1)
            {
                std::cout << "Maximum number of iterations reached.\n\n";
            }
            else if(info==3)
            {
                std::cout << "No shifts could be applied during implicit\n";
                std::cout << "Arnoldi update, try increasing NCV.\n\n";
            }

            Evals = matrix<std::complex<double>>(nev,nev);
            for(int i=0; i<nev; i++)
            {
                Evals(i,i) = d[i];
            }
            
            Evecs = matrix<std::complex<double>>(n,nev);
            for(int i=0; i<nev; i++)
            {
                for(int j=0; j<n; j++)
                {
                    Evecs(j,i) = v[i*n+j];
                }
            }

            std::complex<double> temp;
            for(int i=0; i<nev; i++)
            {
                for(int j=i; j<nev; j++)
                {
                    if(Evals(j,j).real() > Evals(i,i).real())
                    {
                        temp = Evals(j,j);
                        Evals(j,j) = Evals(i,i);
                        Evals(i,i) = temp;
                        for(int k=0; k<n; k++)
                        {
                            temp = Evecs(k,i);
                            Evecs(k,i) = Evecs(k,j);
                            Evecs(k,j) = temp;
                        }
                    }
                }
            }
            
            opt.v0 = matrix<std::complex<double>>(n);
            for(int ii = 0; ii < n; ii++)
            {
                opt.v0(ii) = resid[ii];
            }

            delete [] iw;
            delete [] rwork;
            delete [] workev;
            delete [] resid;
            delete [] v;
            delete [] iparam;
            delete [] ipntr;
            delete [] workd;
            delete [] workl;
            delete [] select;
            delete [] d;
        }
    }


    inline matrix<std::complex<double>,2> eigs(const matrix<std::complex<double>,2> &A, size_t k, std::complex<double> sigma, matrix<std::complex<double>,2> *vr_return = NULL, eigs_opt *opt = NULL)
    {
        if(A.empty())
        {
            throw KeyCppException("Empty matrix in eigs!");
        }
        if(A.size(1) != A.size(2))
        {
            throw KeyCppException("Non-square matrix in eigs!");
        }
        
        if((A.size(1) - k) < 2)
        {
            if(k > A.size(1))
            {
                throw KeyCppException("Too many eigenvalues requested!");
            }
            else
            {
                return eig(A,vr_return);
            }
        }
        
        matrix<std::complex<double>,2> eig_val;
        if(vr_return == NULL)
        {
            int n = A.size(1);
            znaupd_shift_invert(n, k, eig_val, sigma, A, opt);
        }
        else
        {
            int n = A.size(1);
            znaupd_shift_invert(n, k, eig_val, *vr_return, sigma, A, opt);
        }
        return eig_val;
    }

    inline void znaupd(int n, int nev, matrix<std::complex<double>> &Evals, std::string which, const matrix<std::complex<double>> &A, const matrix<std::complex<double>> &B, eigs_opt *opt_base)
    {
        eigs_opt opt;
        
        if(opt_base == NULL)
        {
            opt = eigs_opt();
        }
        else
        {
            opt = *opt_base;
        }
        
        int ido = 0;

        char bmat[2] = "I";
        std::vector<char> which_char(which.c_str(), which.c_str() + which.size() + 1u);
        
        int *iw = new int[A.size(1)];
        matrix<std::complex<double>> Y = lu(B,iw);

        double tol = 0.0;

        std::complex<double> *resid;
        resid = new std::complex<double>[n];

        int ncv = 2*nev + 1;
        if(ncv>=n)
        {
            ncv = n;
        }

        std::complex<double> *v;
        int ldv = n;
        v = new std::complex<double>[ldv*ncv];

        int *iparam;
        iparam = new int[11];
        iparam[0] = 1;
        iparam[2] = opt.maxit;
        iparam[6] = 1;

        int *ipntr;
        ipntr = new int[14];

        std::complex<double> *workd;
        workd = new std::complex<double>[3*n];

        int lworkl = 3*ncv*ncv + 5*ncv;
        std::complex<double> *workl;
        workl = new std::complex<double>[lworkl];

        double *rwork;
        rwork = new double[ncv];

        int info = 0;
        if(opt.v0.numel() == (size_t)n)
        {
            info = 1;
            for(int ii = 0; ii < n; ii++)
            {
                resid[ii] = opt.v0(ii);
            }
        }

        int rvec = 0;

        int *select;
        select = new int[ncv];
        std::complex<double> *d;
        d = new std::complex<double>[ncv];
        std::complex<double> sigma;

        std::complex<double> *workev;
        workev = new std::complex<double>[3*ncv];

        int ierr;

        do
        {
            znaupd_(&ido, bmat, &n, &which_char[0], &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);

            
            if((ido == 1) || (ido == -1))
            {
                mv_special(n, workd+ipntr[0]-1, workd+ipntr[1]-1,A,Y,iw);
            }
        } while((ido == 1) || (ido == -1));


        if(info<0)
        {
            if(info == -5)
            {
                throw KeyCppException("Invalid sigma provided to eigs!");
            }
            std::cout << "Error with znaupd, info = " << info << "\n";
            std::cout << "Check documentation in znaupd\n\n";
        }
        else
        {
            char howmny[] = "All";
            zneupd_(&rvec, howmny, select, d, v, &ldv, &sigma, workev,
            bmat, &n, &which_char[0], &nev, &tol, resid, &ncv, v, &ldv,
            iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);

            if(ierr!=0)
            {
                std::cout << "Error with zneupd, info = " << ierr << "\n";
                std::cout << "Check the documentation of zneupd.\n\n";
            }
            else if(info==1)
            {
                std::cout << "Maximum number of iterations reached.\n\n";
            }else if(info==3)
            {
                std::cout << "No shifts could be applied during implicit\n";
                std::cout << "Arnoldi update, try increasing NCV.\n\n";
            }


            Evals = matrix<std::complex<double>>(nev,nev);
            for(int i=0; i<nev; i++)
            {
                Evals(i,i) = d[nev-1-i];
            }

            // Sort the eigenvalues by real part

            std::complex<double> temp;
            for(int i=0; i<nev; i++)
            {
                for(int j=i; j<nev; j++)
                {
                    if(Evals(j,j).real() > Evals(i,i).real())
                    {
                        temp = Evals(j,j);
                        Evals(j,j) = Evals(i,i);
                        Evals(i,i) = temp;
                    }
                }
            }
            
            opt.v0 = matrix<std::complex<double>>(n);
            for(int ii = 0; ii < n; ii++)
            {
                opt.v0(ii) = resid[ii];
            }

            delete [] iw;
            delete [] rwork;
            delete [] workev;
            delete [] resid;
            delete [] v;
            delete [] iparam;
            delete [] ipntr;
            delete [] workd;
            delete [] workl;
            delete [] select;
            delete [] d;
        }
    }

    inline void znaupd(int n, int nev, matrix<std::complex<double>> &Evals, matrix<std::complex<double>> &Evecs, std::string which, const matrix<std::complex<double>> &A, const matrix<std::complex<double>> &B, eigs_opt *opt_base)
    {
        eigs_opt opt;
        
        if(opt_base == NULL)
        {
            opt = eigs_opt();
        }
        else
        {
            opt = *opt_base;
        }
        
        int ido = 0;
        char bmat[2] = "I";
        //char which[3] = "SM";
        std::vector<char> which_char(which.c_str(), which.c_str() + which.size() + 1u);
        
        int *iw = new int[A.size(1)];
        matrix<std::complex<double>> Y = lu(B,iw);
        
        double tol = 0.0;
        std::complex<double> *resid;
        resid = new std::complex<double>[n];
        int ncv = 2*nev + 1;
        if (ncv>n)
        {
            ncv = n;
        }
        std::complex<double> *v;
        int ldv = n;
        v = new std::complex<double>[ldv*ncv];

        int *iparam;
        iparam = new int[11];
        iparam[0] = 1;
        iparam[2] = opt.maxit;
        iparam[6] = 1;
        int *ipntr;
        ipntr = new int[14];
        std::complex<double> *workd;
        workd = new std::complex<double>[3*n];
        int lworkl = 3*ncv*ncv + 5*ncv;
        std::complex<double> *workl;
        workl = new std::complex<double>[lworkl];
        double *rwork;
        rwork = new double[ncv];
        int info = 0;
        if(opt.v0.numel() == (size_t)n)
        {
            info = 1;
            for(int ii = 0; ii < n; ii++)
            {
                resid[ii] = opt.v0(ii);
            }
        }
        int rvec = 1;  // Changed from above
        int *select;
        select = new int[ncv];
        std::complex<double> *d;
        d = new std::complex<double>[ncv];
        std::complex<double> sigma;
        std::complex<double> *workev;
        workev = new std::complex<double>[3*ncv];
        int ierr;

        do
        {
            znaupd_(&ido, bmat, &n, &which_char[0], &nev, &tol, resid, 
            &ncv, v, &ldv, iparam, ipntr, workd, workl,
            &lworkl, rwork, &info);

            if((ido == 1) || (ido == -1))
            {
                mv_special(n, workd+ipntr[0]-1, workd+ipntr[1]-1,A,Y,iw);
            }
        } while((ido==1)||(ido==-1));

        if(info<0)
        {
            if(info == -5)
            {
                throw KeyCppException("Invalid sigma provided to eigs!");
            }
            std::cout << "Error with znaupd, info = " << info << "\n";
            std::cout << "Check documentation in znaupd\n\n";
        }
        else
        {
            char howmny[] = "All";
            zneupd_(&rvec, howmny, select, d, v, &ldv, &sigma, workev,
            bmat, &n, &which_char[0], &nev, &tol, resid, &ncv, v, &ldv,
            iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);

            if(ierr!=0)
            {
                std::cout << "Error with zneupd, info = " << ierr << "\n";
                std::cout << "Check the documentation of zneupd.\n\n";
            }
            else if(info==1)
            {
                std::cout << "Maximum number of iterations reached.\n\n";
            }
            else if(info==3)
            {
                std::cout << "No shifts could be applied during implicit\n";
                std::cout << "Arnoldi update, try increasing NCV.\n\n";
            }

            Evals = matrix<std::complex<double>>(nev,nev);
            for(int i=0; i<nev; i++)
            {
                Evals(i,i) = d[i];
            }
            
            Evecs = matrix<std::complex<double>>(n,nev);
            for(int i=0; i<nev; i++)
            {
                for(int j=0; j<n; j++)
                {
                    Evecs(j,i) = v[i*n+j];
                }
            }

            std::complex<double> temp;
            for(int i=0; i<nev; i++)
            {
                for(int j=i; j<nev; j++)
                {
                    if(Evals(j,j).real() > Evals(i,i).real())
                    {
                        temp = Evals(j,j);
                        Evals(j,j) = Evals(i,i);
                        Evals(i,i) = temp;
                        for(int k=0; k<n; k++)
                        {
                            temp = Evecs(k,i);
                            Evecs(k,i) = Evecs(k,j);
                            Evecs(k,j) = temp;
                        }
                    }
                }
            }
            
            opt.v0 = matrix<std::complex<double>>(n);
            for(int ii = 0; ii < n; ii++)
            {
                opt.v0(ii) = resid[ii];
            }

            delete [] iw;
            delete [] rwork;
            delete [] workev;
            delete [] resid;
            delete [] v;
            delete [] iparam;
            delete [] ipntr;
            delete [] workd;
            delete [] workl;
            delete [] select;
            delete [] d;
        }
    }

    inline matrix<std::complex<double>,2> eigs(const matrix<std::complex<double>,2> &A, const matrix<std::complex<double>,2> &B, size_t k = 6, std::string sigma = "LM", matrix<std::complex<double>,2> *vr_return = NULL, eigs_opt *opt = NULL)
    {
        if(A.empty() || B.empty())
        {
            throw KeyCppException("Empty matrix in eigs!");
        }
        if(A.size(1) != A.size(2) || B.size(1) != B.size(2))
        {
            throw KeyCppException("Non-square matrix in eigs!");
        }
        if(A.size(1) != B.size(1))
        {
            throw KeyCppException("Invalid matrices in eigs!");
        }
        
        if((A.size(1) - k) < 2)
        {
            if(k > A.size(1))
            {
                throw KeyCppException("Too many eigenvalues requested!");
            }
            else
            {
                return eig(A,B,vr_return);
            }
        }
        
        matrix<std::complex<double>,2> eig_val;
        if(vr_return == NULL)
        {
            int n = A.size(1);
            znaupd(n, k, eig_val, sigma, A, B, opt);
        }
        else
        {
            int n = A.size(1);
            znaupd(n, k, eig_val, *vr_return, sigma, A, B, opt);
        }
        return eig_val;
    }

    inline void znaupd_shift_invert(int n, int nev, matrix<std::complex<double>> &Evals, std::complex<double> sigma, const matrix<std::complex<double>> &A, const matrix<std::complex<double>> &B, eigs_opt *opt_base)
    {
        eigs_opt opt;
        
        if(opt_base == NULL)
        {
            opt = eigs_opt();
        }
        else
        {
            opt = *opt_base;
        }
        
        int ido = 0;

        char bmat[2] = "I";
        std::string which = "LM";
        std::vector<char> which_char(which.c_str(), which.c_str() + which.size() + 1u);
        
        int *iw = new int[A.size(1)];
        matrix<std::complex<double>> Y = lu(A - sigma*B,iw);

        double tol = 0.0;

        std::complex<double> *resid;
        resid = new std::complex<double>[n];

        int ncv = 2*nev + 1;
        if(ncv > n)
        {
            ncv = n;
        }

        std::complex<double> *v;
        int ldv = n;
        v = new std::complex<double>[ldv*ncv];

        int *iparam;
        iparam = new int[11];
        iparam[0] = 1;
        iparam[2] = opt.maxit;
        iparam[6] = 1;

        int *ipntr;
        ipntr = new int[14];

        std::complex<double> *workd;
        workd = new std::complex<double>[3*n];

        int lworkl = 3*ncv*ncv + 5*ncv;
        std::complex<double> *workl;
        workl = new std::complex<double>[lworkl];

        double *rwork;
        rwork = new double[ncv];

        int info = 0;
        if(opt.v0.numel() == (size_t)n)
        {
            info = 1;
            for(int ii = 0; ii < n; ii++)
            {
                resid[ii] = opt.v0(ii);
            }
        }

        int rvec = 0;

        int *select;
        select = new int[ncv];
        std::complex<double> *d;
        d = new std::complex<double>[ncv];

        std::complex<double> *workev;
        workev = new std::complex<double>[3*ncv];
        std::complex<double> sigma_temp;

        int ierr;
        do
        {
            znaupd_(&ido, bmat, &n, &which_char[0], &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);

            
            if((ido == 1) || (ido == -1))
            {
                mv_special(n, workd+ipntr[0]-1, workd+ipntr[1]-1,B,Y,iw);
            }
        } while((ido == 1) || (ido == -1));

        if(info != 0)
        {
            throw znaupdException(info);
        }
        else
        {
            char howmny[] = "All";
            zneupd_(&rvec, howmny, select, d, v, &ldv, &sigma_temp, workev,
            bmat, &n, &which_char[0], &nev, &tol, resid, &ncv, v, &ldv,
            iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);

            if(ierr!=0)
            {
                throw zneupdException(ierr);
            }


            Evals = matrix<std::complex<double>>(nev,nev);
            for(int i=0; i<nev; i++)
            {
                Evals(i,i) = sigma + 1.0/d[nev-1-i];
            }

            // Sort the eigenvalues by real part

            std::complex<double> temp;
            for(int i=0; i<nev; i++)
            {
                for(int j=i; j<nev; j++)
                {
                    if(Evals(j,j).real() > Evals(i,i).real())
                    {
                        temp = Evals(j,j);
                        Evals(j,j) = Evals(i,i);
                        Evals(i,i) = temp;
                    }
                }
            }
            
            opt.v0 = matrix<std::complex<double>>(n);
            for(int ii = 0; ii < n; ii++)
            {
                opt.v0(ii) = resid[ii];
            }

            delete [] iw;
            delete [] rwork;
            delete [] workev;
            delete [] resid;
            delete [] v;
            delete [] iparam;
            delete [] ipntr;
            delete [] workd;
            delete [] workl;
            delete [] select;
            delete [] d;
        }
    }

    inline void znaupd_shift_invert(int n, int nev, matrix<std::complex<double>> &Evals, matrix<std::complex<double>> &Evecs, std::complex<double> sigma, const matrix<std::complex<double>> &A, const matrix<std::complex<double>> &B, eigs_opt *opt_base)
    {
        eigs_opt opt;
        
        if(opt_base == NULL)
        {
            opt = eigs_opt();
        }
        else
        {
            opt = *opt_base;
        }
        
        int ido = 0;
        char bmat[2] = "I";
        std::string which = "LM";
        std::vector<char> which_char(which.c_str(), which.c_str() + which.size() + 1u);
        
        int *iw = new int[A.size(1)];
        matrix<std::complex<double>> Y = lu(A - sigma*B,iw);
        
        double tol = 0.0;
        std::complex<double> *resid;
        resid = new std::complex<double>[n];
        int ncv = 2*nev + 1;
        if(ncv > n)
        {
            ncv = n;
        }
        std::complex<double> *v;
        int ldv = n;
        v = new std::complex<double>[ldv*ncv];
        int *iparam;
        iparam = new int[11];
        
        iparam[0] = 1;
        iparam[2] = opt.maxit;
        iparam[6] = 1;
        int *ipntr;
        ipntr = new int[14];
        std::complex<double> *workd;
        workd = new std::complex<double>[3*n];
        int lworkl = 3*ncv*ncv + 5*ncv;
        std::complex<double> *workl;
        workl = new std::complex<double>[lworkl];
        double *rwork;
        rwork = new double[ncv];
        int info = 0;
        if(opt.v0.numel() == (size_t)n)
        {
            info = 1;
            for(int ii = 0; ii < n; ii++)
            {
                resid[ii] = opt.v0(ii);
            }
        }
        int rvec = 1;  // Changed from above
        int *select;
        select = new int[ncv];
        std::complex<double> *d;
        std::complex<double> sigma_temp;
        d = new std::complex<double>[ncv];
        std::complex<double> *workev;
        workev = new std::complex<double>[3*ncv];
        int ierr;

        do
        {
            znaupd_(&ido, bmat, &n, &which_char[0], &nev, &tol, resid, 
            &ncv, v, &ldv, iparam, ipntr, workd, workl,
            &lworkl, rwork, &info);

            if((ido == 1) || (ido == -1))
            {
                mv_special(n, workd+ipntr[0]-1, workd+ipntr[1]-1,B,Y,iw);
            }
        } while((ido==1)||(ido==-1));

        if(info != 0)
        {
            throw znaupdException(info);
        }
        else
        {
            char howmny[] = "All";
            zneupd_(&rvec, howmny, select, d, v, &ldv, &sigma_temp, workev,
            bmat, &n, &which_char[0], &nev, &tol, resid, &ncv, v, &ldv,
            iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);

            if(ierr!=0)
            {
                throw zneupdException(ierr);
            }

            Evals = matrix<std::complex<double>>(nev,nev);
            for(int i=0; i<nev; i++)
            {
                Evals(i,i) = sigma + 1.0/d[i];
            }
            
            Evecs = matrix<std::complex<double>>(n,nev);
            for(int i=0; i<nev; i++)
            {
                for(int j=0; j<n; j++)
                {
                    Evecs(j,i) = v[i*n+j];
                }
            }

            std::complex<double> temp;
            for(int i=0; i<nev; i++)
            {
                for(int j=i; j<nev; j++)
                {
                    if(Evals(j,j).real() > Evals(i,i).real())
                    {
                        temp = Evals(j,j);
                        Evals(j,j) = Evals(i,i);
                        Evals(i,i) = temp;
                        for(int k=0; k<n; k++)
                        {
                            temp = Evecs(k,i);
                            Evecs(k,i) = Evecs(k,j);
                            Evecs(k,j) = temp;
                        }
                    }
                }
            }
            
            opt.v0 = matrix<std::complex<double>>(n);
            for(int ii = 0; ii < n; ii++)
            {
                opt.v0(ii) = resid[ii];
            }
            
            delete [] iw;
            delete [] rwork;
            delete [] workev;
            delete [] resid;
            delete [] v;
            delete [] iparam;
            delete [] ipntr;
            delete [] workd;
            delete [] workl;
            delete [] select;
            delete [] d;
        }
    }


    inline matrix<std::complex<double>,2> eigs(const matrix<std::complex<double>,2> &A, const matrix<std::complex<double>,2> &B, size_t k, std::complex<double> sigma, matrix<std::complex<double>,2> *vr_return = NULL, eigs_opt *opt = NULL)
    {
        if(A.empty() || B.empty())
        {
            throw KeyCppException("Empty matrix in eigs!");
        }
        if(A.size(1) != A.size(2) || B.size(1) != B.size(2))
        {
            throw KeyCppException("Non-square matrix in eigs!");
        }
        if(A.size(1) != B.size(1) || A.size(2) != B.size(2))
        {
            throw KeyCppException("Matrix sizes not compatible in eigs!");
        }
        
        if((A.size(1) - k) < 2)
        {
            if(k > A.size(1))
            {
                throw KeyCppException("Too many eigenvalues requested!");
            }
            else
            {
                return eig(A,B,vr_return);
            }
        }
        
        matrix<std::complex<double>,2> eig_val;
        if(vr_return == NULL)
        {
            int n = A.size(1);
            znaupd_shift_invert(n, k, eig_val, sigma, A, B, opt);
        }
        else
        {
            int n = A.size(1);
            znaupd_shift_invert(n, k, eig_val, *vr_return, sigma, A, B, opt);
        }
        return eig_val;
    }
}
