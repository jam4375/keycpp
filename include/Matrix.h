// matrix.h -- matrix class

#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
#include <stdexcept>
#include <vector>
#include <stdarg.h>
#include "vector_k.h"

namespace keycpp
{
    extern "C"{
	/** \brief This provides a C interface to BLAS's double matrix-matrix multiplication function. */
	void dgemm_(const char *TRANSA, const char *TRANSB, const int *M, const int *N,
	            const int *K, double *ALPHA, const double *A, const int *LDA, const double *B,
	            const int *LDB, double *BETA, double *C, const int *LDC);
	            
	/** \brief This provides a C interface to BLAS's complex double matrix-matrix multiplication function. */
	void zgemm_(const char *TRANSA, const char *TRANSB, const int *M, const int *N,
	            const int *K, std::complex<double> *ALPHA, const std::complex<double> *A, const int *LDA, const std::complex<double> *B,
	            const int *LDB, std::complex<double> *BETA, std::complex<double> *C, const int *LDC);
	            
	/** \brief This provides a C interface to BLAS's double matrix-vector multiplication function. */
	void dgemv_(const char * TRANS, const int *M, const int *N, const double *ALPHA, const double *A, const int *LDA, const double *X,
	            const int *INCX, const double *BETA, double *Y, const int *INCY);
	            
	/** \brief This provides a C interface to BLAS's complex double matrix-vector multiplication function. */
	void zgemv_(const char * TRANS, const int *M, const int *N, const std::complex<double> *ALPHA,
	            const std::complex<double> *A, const int *LDA, const std::complex<double> *X,
	            const int *INCX, const std::complex<double> *BETA, std::complex<double> *Y, const int *INCY);
	}

	class MatrixException : public std::runtime_error
	{
	    public:
		MatrixException(const std::string &msg) : std::runtime_error(msg){}
	};
	
	template<class T, size_t dim>
	class matrix
	{
	public:
        typedef typename keycpp::vector_k<T>::iterator iterator;
        typedef typename keycpp::vector_k<T>::const_iterator const_iterator;
	
		matrix();
		matrix(const matrix<T,dim> &A);
		matrix(bool, vector_k<T>& mat);
		matrix(const size_t &d1, const size_t &d2 = 0, ...);
		matrix(const vector_k<vector_k<T>>& mat);
		matrix(const vector_k<T>& mat);
		matrix(const std::initializer_list<std::initializer_list<T>>& lst);
		matrix(const std::initializer_list<T>& lst);
		T& operator()(const size_t &i, const size_t &j, const size_t &k, ...);
		T operator()(const size_t &i, const size_t &j, const size_t &k, ...) const;
		T& operator()(const size_t &i, const size_t &j);
		T operator()(const size_t &i, const size_t &j) const;
		T& operator()(const size_t &i);
		T operator()(const size_t &i) const;
		/*vector_k<T> operator*(const vector_k<T> &x) const;
		matrix<T,dim> operator*(const matrix<T,dim-1> &x) const;
		matrix<T,dim> operator*(const matrix<T,dim+1> &x) const;*/
		template<class U>
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> operator*(const matrix<U,dim> &B) const;
		//template<class U>
		//matrix<decltype(std::declval<T>()*std::declval<U>()),dim> operator+(const matrix<U,dim> &B) const;
		template<class U>
		matrix<T,dim>& operator+=(const matrix<U,dim> &B);
		template<class U>
		matrix<T,dim>& operator-=(const matrix<U,dim> &B);
		template<class U>
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> operator-(const matrix<U,dim> &B) const;
		bool operator!=(const matrix<T,dim> &B) const;
		bool operator==(const matrix<T,dim> &B) const;
        template<class U>
		matrix<T,dim>& operator=(const matrix<U,dim> &v);
		size_t size(const size_t &n) const;
		void resize(const std::array<size_t,dim> &pSize);
		void resize(const matrix<size_t,2> &pSize);
		void resize(const size_t &pSize);
		bool empty() const;
        iterator begin() {return mData.begin();};
        iterator end() {return mData.end();};
        const_iterator begin() const {return mData.begin();};
        const_iterator end() const {return mData.end();};
		matrix<T,2> row(const size_t &i);
		matrix<T,2> row(const size_t &i) const;
		matrix<T,2> col(const size_t &j) const;
		matrix<T,2> col(const size_t &j);
		void reserve(const size_t &N);
		vector_k<T> mData;
		std::array<size_t,dim> mSize;
        size_t length() const
        {
            return *std::max_element(mSize.begin(), mSize.end());
        };
        size_t numel() const
        {
            return mData.size();
        };
        
        bool isVec() const
        {
            if(dim != 2)
            {
                return false;
            }
            if(mSize[0] == 1 || mSize[1] == 1)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
		
        operator vector_k<T>()
        {
            static_assert(dim == 2,"Conversion from matrix to vector_k requires matrix with dimension 2.");
            if(mSize[0] != 1)
            {
                throw MatrixException("Cannot convert to vector_k(). Number of rows must equal 1.");
            }
            vector_k<T> v1(mSize[1]);
            for(size_t ii = 0; ii < mSize[1]; ii++)
            {
                v1[ii] = mData[ii];
            }
            return v1;
        };
        
        operator T()
        {
            if(mData.size() != 1)
            {
                throw MatrixException("Cannot convert to T. Number of elements must be 1.");
            }

            return mData[0];
        };
	};

	template<class T, size_t dim>
	matrix<T,dim>::matrix() : mData()
	{
	    static_assert(dim >= 2,"Minimum matrix dimension is 2.");
	    for(size_t ii = 0; ii < dim; ii++)
	    {
	        mSize[ii] = 0;
	    }
	}

	template<class T, size_t dim>
	matrix<T,dim>::matrix(const matrix<T,dim> &A) : mData()
	{
	    static_assert(dim >= 2,"Minimum matrix dimension is 2.");
	    for(size_t ii = 0; ii < dim; ii++)
	    {
	        mSize[ii] = A.mSize[ii];
	    }
        mData.resize(A.numel());
	    for(size_t ii = 0; ii < A.numel(); ii++)
	    {
	        mData[ii] = A.mData[ii];
	    }
	}

	template<class T, size_t dim>
	matrix<T,dim>::matrix(bool, vector_k<T>& mat) : mData(&mat[0],mat.size(),mat.get_inc())
	{
	    static_assert(dim >= 2,"Minimum matrix dimension is 2.");
	    static_assert(dim == 2,"This constructor is a specialization for matrices of dimension 2.");
		if(mat.empty())
		{
			throw MatrixException("Cannot assign empty vector to a matrix object!");
		}
		
		mSize[0] = 1;
		mSize[1] = mat.size();
	}

    // First parameter is a trick to not override the default constructor...
	template<class T, size_t dim>
	matrix<T,dim>::matrix(const size_t & d1, const size_t & d2, ...) : mData()
	{
	    static_assert(dim >= 2,"Minimum matrix dimension is 2.");
        
        size_t total_size;
        if(d2 == 0)
        {
            mSize[0] = 1;
            mSize[1] = d1;
            total_size = d1;
        }
		else
		{
		    mSize[0] = d1;
		    mSize[1] = d2;
            total_size = d1*d2;
		}
		
		if(dim > 2)
		{
            va_list ap;
            va_start(ap, d2);
            for(size_t ii = 2; ii < dim; ii++)
            {
                mSize[ii] = va_arg(ap, size_t);
                if(mSize[ii] < 1)
		        {
			        throw MatrixException("Tried to create matrix with invalid size.");
		        }
		        total_size *= mSize[ii];
            }
            va_end(ap);
        }
        mData.resize(total_size);
	}

	template<class T, size_t dim>
	matrix<T,dim>::matrix(const vector_k<vector_k<T> >& mat) : mData(mat.size()*mat[0].size())
	{
	    static_assert(dim >= 2,"Minimum matrix dimension is 2.");
	    static_assert(dim == 2,"This constructor is a specialization for matrices of dimension 2.");
		if(mat.empty())
		{
			throw MatrixException("Cannot assign empty vector of vectors to a matrix object!");
		}
		else if(mat[0].empty())
		{
			throw MatrixException("Cannot assign empty vector of vectors to a matrix object!");
		}
		
		mSize[0] = mat.size();
		mSize[1] = mat[0].size();
		
		for(size_t ii = 0; ii < mSize[0]; ii++)
		{
			for(size_t jj = 0; jj < mSize[1]; jj++)
			{
				mData[jj*mSize[0] + ii] = mat[ii][jj];
			}
		}
	}

	template<class T, size_t dim>
	matrix<T,dim>::matrix(const vector_k<T>& mat) : mData(mat.size())
	{
	    static_assert(dim >= 2,"Minimum matrix dimension is 2.");
	    static_assert(dim == 2,"This constructor is a specialization for matrices of dimension 2.");
		if(mat.empty())
		{
			throw MatrixException("Cannot assign empty vector to a matrix object!");
		}
		
		mSize[0] = 1;
		mSize[1] = mat.size();
		
		for(size_t ii = 0; ii < mSize[1]; ii++)
		{
			mData[ii] = mat[ii];
		}
	}

	template<class T, size_t dim>
	matrix<T,dim>::matrix(const std::initializer_list<std::initializer_list<T> >& lst) : matrix(lst.size(), lst.size() ? lst.begin()->size() : 0)
	{
	    static_assert(dim >= 2,"Minimum matrix dimension is 2.");
	    static_assert(dim == 2,"This constructor is a specialization for matrices of dimension 2.");
		if(lst.size() <= 0 || lst.begin()->size() <= 0)
		{
			throw MatrixException("Cannot assign empty initializer list to a matrix object!");
		}
		size_t ii = 0, jj = 0;
		for(const auto& l : lst)
		{
			for(const auto& v : l)
			{
				mData[jj*mSize[0] + ii] = v;
				jj++;
			}
			jj = 0;
			ii++;
		}
	}

	template<class T, size_t dim>
	matrix<T,dim>::matrix(const std::initializer_list<T>& lst) : matrix(lst.size())
	{
	    static_assert(dim >= 2,"Minimum matrix dimension is 2.");
	    static_assert(dim == 2,"This constructor is a specialization for matrices of dimension 2.");
		if(lst.size() <= 0)
		{
			throw MatrixException("Cannot assign empty initializer list to a matrix object!");
		}
		size_t ii = 0;
		for(const auto& l : lst)
		{
			mData[ii] = l;
			ii++;
		}
	}

	template<class T, size_t dim>
	T& matrix<T,dim>::operator()(const size_t &i, const size_t &j, const size_t &k, ...)
	{
	    static_assert(dim >= 3,"This function is a specialization for matrices of dimension greater than or equal to 3.");
		if(mData.empty())
		{
			throw MatrixException("Cannot access member of empty matrix!");
		}
		if(i > (mSize[0]-1) || j > (mSize[1]-1) || k > (mSize[2]-1))
		{
			throw MatrixException("Tried to access invalid matrix member!");
		}
		
		std::array<size_t,dim> index_all;
		index_all[0] = i;
		index_all[1] = j;
		index_all[2] = k;
		
		if(dim > 3)
		{
		    va_list ap;
            va_start(ap, i);
            for(size_t ii = 3; ii < dim; ii++)
            {
                index_all[ii] = va_arg(ap, size_t);
                if(index_all[ii] > (mSize[ii]-1))
		        {
			        throw MatrixException("Tried to access invalid matrix member!");
		        }
            }
            va_end(ap);
        }
        size_t index = index_all[0];
        
        for(size_t ii = 1; ii < dim; ii++)
        {
            size_t prod = 1;
            for(size_t jj = 0; jj < ii; jj++)
            {
                prod *= mSize[jj];
            }
            index += index_all[ii]*prod;
        }
		
		return mData[index];
	}

	template<class T, size_t dim>
	T matrix<T,dim>::operator()(const size_t &i, const size_t &j, const size_t &k, ...) const
	{
	    static_assert(dim >= 3,"This function is a specialization for matrices of dimension greater than or equal to 3.");
		if(mData.empty())
		{
			throw MatrixException("Cannot access member of empty matrix!");
		}
		if(i > (mSize[0]-1) || j > (mSize[1]-1) || k > (mSize[2]-1))
		{
			throw MatrixException("Tried to access invalid matrix member!");
		}
		
		std::array<size_t,dim> index_all;
		index_all[0] = i;
		index_all[1] = j;
		index_all[2] = k;
		
		if(dim > 3)
		{
		    va_list ap;
            va_start(ap, j);
            for(size_t ii = 3; ii < dim; ii++)
            {
                index_all[ii] = va_arg(ap, size_t);
                if(index_all[ii] > (mSize[ii]-1))
		        {
			        throw MatrixException("Tried to access invalid matrix member!");
		        }
            }
            va_end(ap);
        }
        size_t index = index_all[0];
        
        for(size_t ii = 1; ii < dim; ii++)
        {
            size_t prod = 1;
            for(size_t jj = 0; jj < ii; jj++)
            {
                prod *= mSize[jj];
            }
            index += index_all[ii]*prod;
        }
		
		return mData[index];
	}

	template<class T, size_t dim>
	T& matrix<T,dim>::operator()(const size_t &i, const size_t &j)
	{
	    static_assert(dim == 2,"This function is a specialization for matrices of dimension 2.");
		if(i > (mSize[0]-1) || j > (mSize[1]-1))
		{
			throw MatrixException("Tried to access invalid matrix member!");
		}
		
		return mData[j*mSize[0] + i];
	}

	template<class T, size_t dim>
	T matrix<T,dim>::operator()(const size_t &i, const size_t &j) const
	{
	    static_assert(dim == 2,"This function is a specialization for matrices of dimension 2.");
		if(i > (mSize[0]-1) || j > (mSize[1]-1))
		{
			throw MatrixException("Tried to access invalid matrix member!");
		}
		
		return mData[j*mSize[0] + i];
	}

	template<class T,size_t dim>
	T& matrix<T,dim>::operator()(const size_t &i)
	{
		if(i > mData.size()-1)
		{
			throw MatrixException("Tried to access invalid matrix member!");
		}
		
	    return mData[i];
	}

	template<class T, size_t dim>
	T matrix<T,dim>::operator()(const size_t &i) const
	{
		if(i > mData.size()-1)
		{
			throw MatrixException("Tried to access invalid matrix member!");
		}
		
	    return mData[i];
	}

	/*template<class T, size_t dim>
	vector_k<T> matrix<T,dim>::operator*(const vector_k<T> &x) const
	{
	    static_assert(dim == 2,"Matrix-Vector multiplication is for matrices of dimension 2 only.");
		if(x.empty())
		{
			throw MatrixException("Vector `x` cannot be empty in matrix-vector multiplication!");
		}
		if(x.size() != mSize[1])
		{
			throw MatrixException("Matrix and vector dimensions are not compatible in matrix-vector multiplication!");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		vector_k<T> b(mSize[0]);
		for(size_t ii = 0; ii < mSize[0]; ii++)
		{
			b[ii] = 0.0;
			for(size_t jj = 0; jj < mSize[1]; jj++)
			{
				b[ii] += mData[jj*mSize[0] + ii]*x[jj];
			}
		}
		return b;
	}*/

	/*template<class T, size_t dim>
	matrix<T,dim> matrix<T,dim>::operator*(const matrix<T,dim-1> &x) const
	{
	    static_assert(dim == 2,"Matrix-Vector multiplication is for matrices of dimension 2 only.");
		if(x.empty())
		{
			throw MatrixException("Vector `x` cannot be empty in matrix-vector multiplication!");
		}
		if(mSize[1] != 1)
		{
			throw MatrixException("Matrix and vector dimensions are not compatible in matrix-vector multiplication!");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		matrix<T,2> b(mSize[0],x.size(1));
		for(size_t ii = 0; ii < mSize[0]; ii++)
		{
		    for(size_t jj = 0; jj < x.size(1); jj++)
		    {
			    b(ii,jj) = mData[ii]*x(jj);
			}
		}
		return b;
	}*/

	/*template<class T, size_t dim>
	matrix<T,dim> matrix<T,dim>::operator*(const matrix<T,dim+1> &x) const
	{
	    static_assert(dim == 1,"This function requires dimension 1.");
		if(x.empty())
		{
			throw MatrixException("Vector `x` cannot be empty in matrix-vector multiplication!");
		}
		if(mSize[0] != x.size(1))
		{
			throw MatrixException("Matrix and vector dimensions are not compatible in matrix-vector multiplication!");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		matrix<T,1> b(x.size(2));
		for(size_t ii = 0; ii < x.size(2); ii++)
		{
		    for(size_t jj = 0; jj < x.size(1); jj++)
		    {
			    b(ii) += mData[jj]*x(jj,ii);
			}
		}
		return b;
	}*/
	
	/*template<>
	inline vector_k<double> matrix<double,2>::operator*(const vector_k<double> &x) const
	{
		if(x.empty())
		{
			throw MatrixException("Vector `x` cannot be empty in matrix-vector multiplication!");
		}
		if(x.size() != mSize[1])
		{
			throw MatrixException("Matrix and vector dimensions are not compatible in matrix-vector multiplication!");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		vector_k<double> b(mSize[0]);
		int m = (int)mSize[0];
		int n = (int)mSize[1];
        char TRANS = 'N';
        double ALPHA = 1.0;
        int LDA = m;
        int INCX = (int)x.get_inc();
        double BETA = 0.0;
        int INCb = (int)b.get_inc();

        dgemv_(&TRANS, &m, &n, &ALPHA, &mData[0], &LDA, &x[0],&INCX, &BETA, &b[0], &INCb);
        
        return b;
    }*/
	
	/*template<>
	inline vector_k<std::complex<double>> matrix<std::complex<double>,2>::operator*(const vector_k<std::complex<double>> &x) const
	{
		if(x.empty())
		{
			throw MatrixException("Vector `x` cannot be empty in matrix-vector multiplication!");
		}
		if(x.size() != mSize[1])
		{
			throw MatrixException("Matrix and vector dimensions are not compatible in matrix-vector multiplication!");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		vector_k<std::complex<double>> b(mSize[0]);
		int m = (int)mSize[0];
		int n = (int)mSize[1];
        char TRANS = 'N';
        std::complex<double> ALPHA = 1.0;
        int LDA = m;
        int INCX = (int)x.get_inc();
        std::complex<double> BETA = 0.0;
        int INCb = (int)b.get_inc();
        
        zgemv_(&TRANS, &m, &n, &ALPHA, &mData[0], &LDA, &x[0],&INCX, &BETA, &b[0], &INCb);
        
        return b;
    }*/

	template<class T, size_t dim>
	template<class U>
	matrix<decltype(std::declval<T>()*std::declval<U>()),dim> matrix<T,dim>::operator*(const matrix<U,dim> &B) const
	{
	    static_assert(dim == 2,"Matrix multiplication not supported for dimensions greater than 2.");
		if(B.size(1) <= 0 || B.size(2) <= 0)
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(mSize[1] != B.size(1))
		{
			throw MatrixException("Matrix dimensions are not compatible in matrix-matrix multiplication!");
		}
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> C(mSize[0],B.size(2));
		for(size_t kk = 0; kk < B.size(2); kk++)
		{
			for(size_t ii = 0; ii < mSize[0]; ii++)
			{
				C(ii,kk) = 0.0;
				for(size_t jj = 0; jj < mSize[1]; jj++)
				{
					C(ii,kk) += mData[jj*mSize[0] + ii]*B(jj,kk);
				}
			}
		}
		return C;
	}
	
	template<>
	template<>
	inline matrix<double,2> matrix<double,2>::operator*(const matrix<double,2> &B) const
	{
		if(B.size(1) <= 0 || B.size(2) <= 0)
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(mSize[1] != B.size(1))
		{
			throw MatrixException("Matrix dimensions are not compatible in matrix-matrix multiplication!");
		}
		matrix<double> C(mSize[0],B.size(2));
		int m = (int)mSize[0];
		int k = (int)mSize[1];
		int n = (int)B.size(2);
        char TRANS = 'N';
        double ALPHA = 1.0;
        int LDA = m;
        int LDB = k;
        double BETA = 0.0;
        int LDC = m;

        dgemm_(&TRANS, &TRANS, &m, &n, &k, &ALPHA, &mData[0], &LDA, &B.mData[0], &LDB, &BETA, &C.mData[0], &LDC);
        
        return C;
    }
	
	template<>
	template<>
	inline matrix<std::complex<double>,2> matrix<std::complex<double>,2>::operator*(const matrix<std::complex<double>,2> &B) const
	{
		if(B.size(1) <= 0 || B.size(2) <= 0)
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(mSize[1] != B.size(1))
		{
			throw MatrixException("Matrix dimensions are not compatible in matrix-matrix multiplication!");
		}
		matrix<std::complex<double>> C(mSize[0],B.size(2));
		int m = (int)mSize[0];
		int k = (int)mSize[1];
		int n = (int)B.size(2);
        char TRANS = 'N';
        std::complex<double> ALPHA = 1.0;
        int LDA = m;
        int LDB = k;
        std::complex<double> BETA = 0.0;
        int LDC = m;

        zgemm_(&TRANS, &TRANS, &m, &n, &k, &ALPHA, &mData[0], &LDA, &B.mData[0], &LDB, &BETA, &C.mData[0], &LDC);
        
        return C;
    }

	/*template<class T, size_t dim>
	template<class U>
	matrix<decltype(std::declval<T>()*std::declval<U>()),dim> matrix<T,dim>::operator+(const matrix<U,dim> &B) const
	{
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		for(int ii = 0; ii < dim; ii++)
		{
		    if(B.size(ii+1) <= 0)
		    {
		        throw MatrixException("Cannot perform operation on empty matrix!");
		    }
		    if(mSize[ii] != B.size(ii+1))
		    {
		        throw MatrixException("Matrix dimensions are not compatible in matrix-matrix addition!");
		    }
		}
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> C;
		C.resize(mSize);
		for(size_t ii = 0; ii < C.mData.size(); ii++)
		{
			C.mData[ii] = this->mData[ii] + B.mData[ii];
		}
		return C;
	}*/
	
	template<class T, size_t dim>
	template<class U>
	matrix<T,dim>& matrix<T,dim>::operator+=(const matrix<U,dim> &B)
	{
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		for(int ii = 0; ii < dim; ii++)
		{
		    if(B.size(ii+1) <= 0)
		    {
		        throw MatrixException("Cannot perform operation on empty matrix!");
		    }
		    if(mSize[ii] != B.size(ii+1))
		    {
		        throw MatrixException("Matrix dimensions are not compatible in matrix-matrix addition!");
		    }
		}
		
		for(size_t ii = 0; ii < this->mData.size(); ii++)
		{
			this->mData[ii] += B.mData[ii];
		}
		return *this;
	}
	
	template<class T, size_t dim>
	template<class U>
	matrix<T,dim>& matrix<T,dim>::operator-=(const matrix<U,dim> &B)
	{
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		for(int ii = 0; ii < dim; ii++)
		{
		    if(B.size(ii+1) <= 0)
		    {
		        throw MatrixException("Cannot perform operation on empty matrix!");
		    }
		    if(mSize[ii] != B.size(ii+1))
		    {
		        throw MatrixException("Matrix dimensions are not compatible in matrix-matrix addition!");
		    }
		}
		
		for(size_t ii = 0; ii < this->mData.size(); ii++)
		{
			this->mData[ii] -= B.mData[ii];
		}
		return *this;
	}

	template<class T, size_t dim>
	template<class U>
	matrix<decltype(std::declval<T>()*std::declval<U>()),dim> matrix<T,dim>::operator-(const matrix<U,dim> &B) const
	{
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		for(int ii = 0; ii < dim; ii++)
		{
		    if(B.size(ii+1) <= 0)
		    {
		        throw MatrixException("Cannot perform operation on empty matrix!");
		    }
		    if(mSize[ii] != B.size(ii+1))
		    {
		        throw MatrixException("Matrix dimensions are not compatible in matrix-matrix addition!");
		    }
		}
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> C;
		C.resize(mSize);
		for(size_t ii = 0; ii < C.mData.size(); ii++)
		{
			C.mData[ii] = this->mData[ii] - B.mData[ii];
		}
		return C;
	}

	template<class T, size_t dim>
	bool matrix<T,dim>::operator!=(const matrix<T,dim> &B) const
	{
		if(mData.empty() || B.empty())
		{
			return true;
		}
		for(int ii = 0; ii < dim; ii++)
		{
		    if(mSize[ii] != B.size(ii+1))
		    {
		        return true;
		    }
		}
		for(size_t ii = 0; ii < mData.size(); ii++)
		{
			if(this->mData[ii] != B.mData[ii])
			{
			    return true;
			}
		}
		return false;
	}

	template<class T, size_t dim>
	bool matrix<T,dim>::operator==(const matrix<T,dim> &B) const
	{
		if(mData.empty() || B.empty())
		{
			return false;
		}
		for(int ii = 0; ii < dim; ii++)
		{
		    if(mSize[ii] != B.size(ii+1))
		    {
		        return false;
		    }
		}
		for(size_t ii = 0; ii < mData.size(); ii++)
		{
			if(this->mData[ii] != B.mData[ii])
			{
			    return false;
			}
		}
		return true;
	}

    template<class T, size_t dim>
    template<class U>
    matrix<T,dim>& matrix<T,dim>::operator=(const matrix<U,dim> &v)
    {
        if(mSize != v.mSize)
        {
            resize(v.mSize);
        }
        for(size_t ii = 0; ii < v.mData.size(); ii++)
        {
            mData[ii] = v.mData[ii];
        }
        
        return *this;
    }

	template<class T, size_t dim>
	size_t matrix<T,dim>::size(const size_t &n) const
	{
		if(n > dim || n == 0)
		{
			throw MatrixException("Invalid dimension in size().");
		}
		return mSize[n-1];
	}
	
	template<class T, size_t dim>
	void matrix<T,dim>::resize(const std::array<size_t,dim> &pSize)
	{
		matrix<T,dim> B = *this;
		std::array<size_t,dim> bakSize = mSize;
		
	    mSize = pSize;
	    size_t total_size = 1;
	    for(size_t ii = 0; ii < dim; ii++)
	    {
	        total_size *= mSize[ii];
	    }
	    
	    for(size_t ii = 0; ii < mData.size(); ii++)
	    {
	        mData[ii] = T();
	    }
        mData.resize(total_size);
        
		for(size_t index_old = 0; index_old < B.mData.size(); index_old++)
		{
		    size_t temp = index_old;
		    std::array<size_t, dim> index_all;
		    bool cont = false;
            for(size_t ii = dim; ii-- > 0; )
            {
                size_t prod = 1;
                for(size_t jj = 0; jj < ii; jj++)
                {
                    prod *= bakSize[jj];
                }
                index_all[ii] = (temp/prod);
                temp -= index_all[ii]*prod;
                if(index_all[ii] >= mSize[ii])
                {
                    cont = true;
                }
            }
            if(cont)
            {
                continue;
            }
            
            size_t index_new = index_all[0];
            
            for(size_t ii = 1; ii < dim; ii++)
            {
                size_t prod = 1;
                for(size_t jj = 0; jj < ii; jj++)
                {
                    prod *= mSize[jj];
                }
                index_new += index_all[ii]*prod;
            }
            
            mData[index_new] = B.mData[index_old];
        }
	}
	
	template<class T, size_t dim>
	void matrix<T,dim>::resize(const matrix<size_t> &pSize)
	{
	    if(pSize.size(2) != dim)
	    {
	        throw MatrixException("Invalid size provided to resize in Matrix.h!");
	    }
	
		matrix<T,dim> B = *this;
		std::array<size_t,dim> bakSize = mSize;
		
	    size_t total_size = 1;
	    for(size_t ii = 0; ii < dim; ii++)
	    {
	        mSize[ii] = pSize(ii);
	        total_size *= mSize[ii];
	    }
	    
	    for(size_t ii = 0; ii < mData.size(); ii++)
	    {
	        mData[ii] = T();
	    }
        mData.resize(total_size);
        
		for(size_t index_old = 0; index_old < B.mData.size(); index_old++)
		{
		    size_t temp = index_old;
		    std::array<size_t, dim> index_all;
		    bool cont = false;
            for(size_t ii = dim; ii-- > 0; )
            {
                size_t prod = 1;
                for(size_t jj = 0; jj < ii; jj++)
                {
                    prod *= bakSize[jj];
                }
                index_all[ii] = (temp/prod);
                temp -= index_all[ii]*prod;
                if(index_all[ii] >= mSize[ii])
                {
                    cont = true;
                }
            }
            if(cont)
            {
                continue;
            }
            
            size_t index_new = index_all[0];
            
            for(size_t ii = 1; ii < dim; ii++)
            {
                size_t prod = 1;
                for(size_t jj = 0; jj < ii; jj++)
                {
                    prod *= mSize[jj];
                }
                index_new += index_all[ii]*prod;
            }
            
            mData[index_new] = B.mData[index_old];
        }
	}
       
    template<class T, size_t dim>
    void matrix<T,dim>::resize(const size_t &pSize)
    {
        static_assert(dim == 2,"This version of resize is not supported for dimensions greater than 2.");
        if(!this->isVec())
        {
            mSize[0] = pSize;
        }
        else
        {
            if(mSize[1] == 1)
            {
                mSize[0] = pSize;
            }
            else if(mSize[0] == 1)
            {
                mSize[1] = pSize;
            }
        }
        
        mData.resize(pSize);
    }

	
	template<class T, size_t dim>
	bool matrix<T,dim>::empty() const
	{
		for(size_t ii = 0; ii < dim; ii++)
		{
		    if(mSize[ii] == 0)
		    {
			    return true;
			}
		}
		
		return false;
	}
	
	template<class T, size_t dim>
	matrix<T,2> matrix<T,dim>::row(const size_t &n)
	{
	    static_assert(dim == 2,"This function is only available for matrices of dimension 2.");
		if(n > mSize[0])
		{
			throw MatrixException("Invalid row index in row().");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot use method row() on empty matrix!");
		}
		vector_k<T> temp(&mData[n],mSize[1],mSize[0]);
		
		matrix<T,2> Row(true, temp);
		return Row;
	}

	template<class T, size_t dim>
	matrix<T,2> matrix<T,dim>::col(const size_t &n)
	{
	    static_assert(dim == 2,"This function is only available for matrices of dimension 2.");
		if(n > mSize[1])
		{
			throw MatrixException("Invalid column index in col().");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot use method col() on empty matrix!");
		}
		
		vector_k<T> temp(&mData[n*mSize[0]],mSize[0],1);
		
		matrix<T,2> Col(true, temp);
		return Col;
	}

	template<class T, size_t dim>
	matrix<T,2> matrix<T,dim>::row(const size_t &n) const
	{
	    static_assert(dim == 2,"This function is only available for matrices of dimension 2.");
		if(n > mSize[0] || n < 0)
		{
			throw MatrixException("Invalid row index in row().");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot use method row() on empty matrix!");
		}
		matrix<T,2> Row(mSize[1]);

		for(size_t jj = 0; jj < mSize[1]; jj++)
		{
			Row(jj) = mData[jj*mSize[0] + n];
		}

		return Row;
	}

	template<class T, size_t dim>
	matrix<T,2> matrix<T,dim>::col(const size_t &n) const
	{
	    static_assert(dim == 2,"This function is only available for matrices of dimension 2.");
		if(n > mSize[1])
		{
			throw MatrixException("Invalid column index in col().");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot use method col() on empty matrix!");
		}
		matrix<T,2> Col(mSize[0]);

		for(size_t ii = 0; ii < mSize[0]; ii++)
		{
			Col(ii) = mData[n*mSize[0] + ii];
		}

		return Col;
	}
	
	template<class T, size_t dim>
	void matrix<T,dim>::reserve(const size_t &N)
	{
		mData.reserve(N);
	}


    // Matrix Operations:
    
	template<class T, class U, size_t dim> matrix<decltype(std::declval<T>()*std::declval<U>()),dim> operator+(const matrix<T,dim>& A, const U& a)
	{
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> C;
		C.resize(size(A));
		for(size_t ii = 0; ii < C.numel(); ii++)
		{
			C(ii) = A(ii) + a;
		}
		return C;
	}

	template<class T, class U, size_t dim> matrix<decltype(std::declval<T>()*std::declval<U>()),dim> operator+(const U& a, const matrix<T,dim>& A)
	{
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> C;
		C.resize(size(A));
		for(size_t ii = 0; ii < C.numel(); ii++)
		{
			C(ii) = A(ii) + a;
		}
		return C;
	}
	
	template<class T, class U, size_t dim>
	matrix<decltype(std::declval<T>()*std::declval<U>()),dim> operator+(const matrix<U,dim> &A, const matrix<T,dim> &B)
	{
		if(A.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		for(int ii = 0; ii < dim; ii++)
		{
		    if(B.size(ii+1) <= 0)
		    {
		        throw MatrixException("Cannot perform operation on empty matrix!");
		    }
		    if(A.size(ii+1) != B.size(ii+1))
		    {
		        throw MatrixException("Matrix dimensions are not compatible in matrix-matrix addition!");
		    }
		}
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> C;
		C.resize(size(A));
		for(size_t ii = 0; ii < C.numel(); ii++)
		{
			C(ii) = A(ii) + B(ii);
		}
		return C;
	}

	template<class T, class U,size_t dim> matrix<decltype(std::declval<T>() - std::declval<U>()),dim> operator-(const matrix<T,dim>& A, const U& a)
	{
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> result;
		result.resize(size(A));
		for(size_t ii = 0; ii < result.numel(); ii++)
		{
			result(ii) = A(ii) - a;
		}
		return result;
	}

	template<class T, class U,size_t dim> matrix<decltype(std::declval<U>() - std::declval<T>()),dim> operator-(const U& a, const matrix<T,dim>& A)
	{
		matrix<decltype(std::declval<U>() - std::declval<T>()),dim> result;
		result.resize(size(A));
		for(size_t ii = 0; ii < result.numel(); ii++)
		{
			result(ii) = a - A(ii);
		}
		return result;
	}
	
	template<class T, class U, size_t dim>
	matrix<decltype(std::declval<U>() - std::declval<T>()),dim> operator-(const matrix<U,dim> &A, const matrix<T,dim> &B)
	{
		if(A.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		for(int ii = 0; ii < dim; ii++)
		{
		    if(B.size(ii+1) <= 0)
		    {
		        throw MatrixException("Cannot perform operation on empty matrix!");
		    }
		    if(A.size(ii+1) != B.size(ii+1))
		    {
		        throw MatrixException("Matrix dimensions are not compatible in matrix-matrix addition!");
		    }
		}
		matrix<decltype(std::declval<U>() - std::declval<T>()),dim> C;
		C.resize(size(A));
		for(size_t ii = 0; ii < C.numel(); ii++)
		{
			C(ii) = A(ii) - B(ii);
		}
		return C;
	}

	template<class T, class U,size_t dim> matrix<decltype(std::declval<T>()*std::declval<U>()),dim> operator*(const T& a, const matrix<U,dim>& A)
	{
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> B;
		B.resize(size(A));
		for(size_t ii = 0; ii < B.numel(); ii++)
		{
			B(ii) = a*A(ii);
		}
		return B;
	}
	
	template<class T, class U,size_t dim> matrix<decltype(std::declval<T>()*std::declval<U>()),dim> operator*(const matrix<U,dim>& A, const T& a)
	{
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> B;
		B.resize(size(A));
		for(size_t ii = 0; ii < B.numel(); ii++)
		{
			B(ii) = a*A(ii);
		}
		return B;
	}
	
	template<class T,size_t dim> matrix<T,dim> operator-(const matrix<T,dim>& A)
	{
		matrix<T,dim> B;
		B.resize(size(A));
		for(size_t ii = 0; ii < B.numel(); ii++)
		{
			B(ii) = -A(ii);
		}
		return B;
	}
	
	template<class T,size_t dim> matrix<T,dim> operator+(const matrix<T,dim>& A)
	{
		matrix<T,dim> B(A.size(1),A.size(2));
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
		B.resize(size(A));
		
		for(size_t ii = 0; ii < A.numel(); ii++)
		{
			B(ii) = A(ii)/a;
		}
		return B;
	}
	
	template<class T, class U,size_t dim> matrix<decltype(std::declval<T>()*std::declval<U>()),dim> operator/(const U& a, const matrix<T,dim>& A)
	{
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim> B;
		B.resize(size(A));
		
		for(size_t ii = 0; ii < A.numel(); ii++)
		{
			B(ii) = a/A(ii);
		}
		return B;
	}
	
	template<class T, class U,size_t dim> matrix<decltype(std::declval<T>()*std::declval<U>()),dim> operator/(const matrix<U,dim>& A, const matrix<T,dim>& B)
	{
	    if(A.numel() == 1)
	    {
		    matrix<decltype(std::declval<T>()*std::declval<U>()),dim> C;
		    C.resize(size(B));
		
		    for(size_t ii = 0; ii < C.numel(); ii++)
		    {
			    C(ii) = A(0)/B(ii);
		    }
		    return C;
	    }
	    else if(B.numel() == 1)
	    {
		    matrix<decltype(std::declval<T>()*std::declval<U>()),dim> C;
		    C.resize(size(A));
		
		    for(size_t ii = 0; ii < C.numel(); ii++)
		    {
			    C(ii) = A(ii)/B(0);
		    }
		    return C;
	    }
	    else
	    {
	        throw MatrixException("Cannot divide two matrices. Use rdivide or ldivide instead.");
	        
	        matrix<decltype(std::declval<T>()*std::declval<U>()),dim> C;
	        return C;
	    }
	}
	
	template<class T,size_t dim> T mat2num(const matrix<T,dim>& A)
	{
		if(A.numel() == 1)
		{
		    return A(0);
		}
		else
		{
		    throw MatrixException("Error in mat2num! Matrix must have only 1 element!");
		    return T();
		}
	}
}

#endif
