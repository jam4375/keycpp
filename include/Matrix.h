// Matrix.h -- matrix class

#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
#include <stdexcept>
#include <vector>
#include <stdarg.h>
#include "vector_k.h"

#define DENSE_MATRIX 0
#define SPARSE_MATRIX 1

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

    class eigs_opt;
	
	std::string sprintf(const std::string &fmt, ...);
    template<class T> void disp(const T &x, std::ostream& outStream = std::cout);
	template<class T> matrix<T,2> max(const matrix<T,2> &A);
	template<class T, size_t dim> matrix<size_t,2> size(const matrix<T,dim> &A);

	class MatrixException : public std::runtime_error
	{
	public:
		MatrixException(const std::string &msg) : std::runtime_error(msg){}
	};
	
	class span
	{
    private:
        size_t m_first, m_last;
        int m_inc;
        bool all = false;
        
	public:
        span() : m_first(), m_last(), m_inc(1), all(true) {};
        span(size_t first, size_t last) : m_first(first), m_last(last), m_inc(1) {};
        span(size_t first, int increment, size_t last) : m_first(first), m_last(last), m_inc(increment) {};
        
        size_t get_first() const {return m_first;};
        size_t get_last() const {return m_last;};
        int get_inc() const {return m_inc;};
        bool get_all() const {return all;};
	};

	
	template<class T, size_t dim>
	class matrix<T,dim,DENSE_MATRIX>
	{
	public:
        typedef typename keycpp::vector_k<T>::iterator iterator;
        typedef typename keycpp::vector_k<T>::const_iterator const_iterator;
	
		matrix();
		matrix(const matrix<T,dim> &A);
		matrix(bool, vector_k<T>& mat, size_t nrows, size_t ncols);
		matrix(bool, vector_k<T>& mat, size_t nrows, size_t ncols, const vector_k<span> &span_all);
		matrix(const size_t &d1, const size_t &d2 = 0, ...);
		matrix(const vector_k<vector_k<T>>& mat);
		matrix(const vector_k<T>& mat);
		matrix(const std::initializer_list<std::initializer_list<T>>& lst);
		matrix(const std::initializer_list<T>& lst);
		T& operator()(const size_t &i, const size_t &j, const size_t &k, ...);
		const T& operator()(const size_t &i, const size_t &j, const size_t &k, ...) const;
		T& operator()(const size_t &i, const size_t &j);
		const T& operator()(const size_t &i, const size_t &j) const;
		matrix<T,2> operator()(const size_t &i, const span &cols);
		matrix<T,2> operator()(const span &rows, const size_t &j);
		matrix<T,2> operator()(const span &rows, const span &cols);
		matrix<T,2> operator()(const size_t &i, const span &cols) const;
		matrix<T,2> operator()(const span &rows, const size_t &j) const;
		matrix<T,2> operator()(const span &rows, const span &cols) const;
		T& operator()(const size_t &i);
		const T& operator()(const size_t &i) const;
		template<class U>
        matrix<T,dim>& operator+=(const matrix<U,dim> &B);
		template<class U>
        matrix<T,dim>& operator-=(const matrix<U,dim> &B);
		bool operator!=(const matrix<T,dim> &B) const;
		bool operator==(const matrix<T,dim> &B) const;
		template<class U>
		matrix<decltype(std::declval<T>()*std::declval<U>()),dim>& trans_mult(const matrix<U,dim> &v) const;
        template<class U>
		matrix<T,dim>& operator=(const matrix<U,dim> &v);
        template<class U>
		matrix<T,dim>& operator=(const vector_k<U> &v);
		matrix<T,dim>& operator=(const T &a);
		template<class U>
		matrix<T,dim>& operator=(const U &a);
		matrix<T,dim>& operator=(const matrix<T,dim> &v);
		matrix<T,dim>& operator=(const vector_k<T> &v);
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
        size_t length() const
        {
            if(!submat)
            {
                return *std::max_element(mSize.begin(), mSize.end());
            }
            else
            {
                size_t temp = 0;
                for(size_t ii = 0; ii < dim; ii++)
                {
                    if(this->size(ii+1) > temp)
                    {
                        temp = this->size(ii+1);
                    }
                }
                return temp;
            }
        };
        size_t numel() const
        {
            if(!submat)
            {
                return mData.size();
            }
            else
            {
                size_t temp = 1;
                for(size_t ii = 0; ii < dim; ii++)
                {
                    temp *= this->size(ii+1);
                }
                return temp;
            }
        };
        
        bool isVec() const
        {
            if(dim != 2)
            {
                return false;
            }
            if(this->size(1) == 1 || this->size(2) == 1)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        
        size_t get_inc() const {return mData.get_inc();};
        
        // Friend functions
        template<class V, class U>
        friend matrix<decltype(std::declval<V>()*std::declval<U>()),2> operator*(const matrix<V,2> &A, const matrix<U,2> &B);
        friend double norm(const matrix<std::complex<double>> &A_in, std::string method);
        friend double norm(const matrix<double> &A_in, std::string method);
        friend void mv_special(int n, std::complex<double> *in, std::complex<double> *out, const matrix<std::complex<double>> &A, matrix<std::complex<double>> &Y, int *iw);
        friend void mv(int n, std::complex<double> *in, std::complex<double> *out, const matrix<std::complex<double>> &A);
        friend void znaupd_shift_invert(int n, int nev, matrix<std::complex<double>> &Evals, std::complex<double> sigma, const matrix<std::complex<double>> &A, eigs_opt *opt);
        friend void znaupd_shift_invert(int n, int nev, matrix<std::complex<double>> &Evals, matrix<std::complex<double>> &Evecs, std::complex<double> sigma, const matrix<std::complex<double>> &A, eigs_opt *opt);
		
        operator vector_k<T>()
        {
            static_assert(dim == 2,"Conversion from matrix to vector_k requires matrix with dimension 2.");
            if(this->size(1) != 1)
            {
                throw MatrixException("Cannot convert to vector_k(). Number of rows must equal 1.");
            }
            vector_k<T> v1(this->size(2));
            for(size_t ii = 0; ii < this->size(2); ii++)
            {
                v1[ii] = this->operator()(ii);
            }
            return v1;
        };
        
        operator T()
        {
            if(this->numel() != 1)
            {
                throw MatrixException("Cannot convert to T. Number of elements must be 1.");
            }

            return this->operator()(0);
        };
        
        template<class U>
        operator matrix<U>()
        {
            matrix<U> B(this->size(1),this->size(2));
            for(size_t ii = 0; ii < this->numel(); ii++)
            {
                B(ii) = this->operator()(ii);
            }

            return B;
        }
        
        bool get_submat() const
        {
            return submat;
        }
        
        std::array<size_t,dim> get_mSize() const
        {
            return mSize;
        }
        
        vector_k<T> get_mData() const
        {
            return mData;
        }
        
        vector_k<span> get_mSpan() const
        {
            return mSpan;
        }
        
    private:
		vector_k<T> mData;
		std::array<size_t,dim> mSize;
		vector_k<span> mSpan;
		bool submat = false;
	};

	template<class T, size_t dim>
	matrix<T,dim,DENSE_MATRIX>::matrix() : mData()
	{
	    static_assert(dim >= 2,"Minimum matrix dimension is 2.");
	    for(size_t ii = 0; ii < dim; ii++)
	    {
	        mSize[ii] = 0;
	    }
	}

	template<class T, size_t dim>
	matrix<T,dim,DENSE_MATRIX>::matrix(const matrix<T,dim> &A) : mData()
	{
	    static_assert(dim >= 2,"Minimum matrix dimension is 2.");
	    for(size_t ii = 0; ii < dim; ii++)
	    {
	        mSize[ii] = A.size(ii+1);
	    }
        mData.resize(A.numel());
	    for(size_t ii = 0; ii < A.numel(); ii++)
	    {
	        mData[ii] = A(ii);
	    }
	}

	template<class T, size_t dim>
	matrix<T,dim,DENSE_MATRIX>::matrix(bool, vector_k<T>& mat, size_t nrows, size_t ncols) : mData(&mat[0],mat.size(),mat.get_inc())
	{
	    static_assert(dim >= 2,"Minimum matrix dimension is 2.");
	    static_assert(dim == 2,"This constructor is a specialization for matrices of dimension 2.");
		if(mat.empty())
		{
			throw MatrixException("Cannot assign empty vector to a matrix object!");
		}
		
		mSize[0] = nrows;
		mSize[1] = ncols;
	}

	template<class T, size_t dim>
	matrix<T,dim,DENSE_MATRIX>::matrix(bool, vector_k<T>& mat, size_t nrows, size_t ncols, const vector_k<span> &span_all) : mData(&mat[0],mat.size(),mat.get_inc()), mSpan(span_all), submat(true)
	{
	    static_assert(dim >= 2,"Minimum matrix dimension is 2.");
	    static_assert(dim == 2,"This constructor is a specialization for matrices of dimension 2.");
		if(mat.empty())
		{
			throw MatrixException("Cannot assign empty vector to a matrix object!");
		}
		
		mSize[0] = nrows;
		mSize[1] = ncols;
	}

    // First parameter is a trick to not override the default constructor...
	template<class T, size_t dim>
	matrix<T,dim,DENSE_MATRIX>::matrix(const size_t & d1, const size_t & d2, ...) : mData()
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
	matrix<T,dim,DENSE_MATRIX>::matrix(const vector_k<vector_k<T> >& mat) : mData(mat.size()*mat[0].size())
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
	matrix<T,dim,DENSE_MATRIX>::matrix(const vector_k<T>& mat) : mData(mat.size())
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
	matrix<T,dim,DENSE_MATRIX>::matrix(const std::initializer_list<std::initializer_list<T> >& lst) : matrix(lst.size(), lst.size() ? lst.begin()->size() : 0)
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
	matrix<T,dim,DENSE_MATRIX>::matrix(const std::initializer_list<T>& lst) : matrix(lst.size())
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
	T& matrix<T,dim,DENSE_MATRIX>::operator()(const size_t &i, const size_t &j, const size_t &k, ...)
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
	const T& matrix<T,dim,DENSE_MATRIX>::operator()(const size_t &i, const size_t &j, const size_t &k, ...) const
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
	T& matrix<T,dim,DENSE_MATRIX>::operator()(const size_t &i, const size_t &j)
	{
	    static_assert(dim == 2,"This function is a specialization for matrices of dimension 2.");
		if(i > (this->size(1)-1) || j > (this->size(2)-1))
		{
			throw MatrixException("Tried to access invalid matrix member!");
		}
		if(!submat)
		{
		    return mData[j*mSize[0] + i];
		}
		else
		{
	        size_t r1, c1;
	        int r_inc, c_inc;
	        
	        if(mSpan[0].get_all())
	        {
	            r1 = 0;
	            r_inc = 1;
	        }
	        else
	        {
	            r1 = mSpan[0].get_first();
	            r_inc = mSpan[0].get_inc();
	        }
	        if(mSpan[1].get_all())
	        {
	            c1 = 0;
	            c_inc = 1;
	        }
	        else
	        {
	            c1 = mSpan[1].get_first();
	            c_inc = mSpan[1].get_inc();
	        }
	        
	        size_t i2 = r1 + i*r_inc;
	        size_t j2 = c1 + j*c_inc;
		    if(i2 > (mSize[0]-1) || j2 > (mSize[1]-1))
		    {
			    throw MatrixException("Tried to access invalid matrix member!");
		    }
		
		    return mData[j2*mSize[0] + i2];
		}
	}

	template<class T, size_t dim>
	const T& matrix<T,dim,DENSE_MATRIX>::operator()(const size_t &i, const size_t &j) const
	{
	    static_assert(dim == 2,"This function is a specialization for matrices of dimension 2.");
		if(i > (this->size(1)-1) || j > (this->size(2)-1))
		{
			throw MatrixException("Tried to access invalid matrix member!");
		}
		if(!submat)
		{
		    return mData[j*mSize[0] + i];
		}
		else
		{
	        size_t r1, c1;
	        int r_inc, c_inc;
	        
	        if(mSpan[0].get_all())
	        {
	            r1 = 0;
	            r_inc = 1;
	        }
	        else
	        {
	            r1 = mSpan[0].get_first();
	            r_inc = mSpan[0].get_inc();
	        }
	        if(mSpan[1].get_all())
	        {
	            c1 = 0;
	            c_inc = 1;
	        }
	        else
	        {
	            c1 = mSpan[1].get_first();
	            c_inc = mSpan[1].get_inc();
	        }
	        
	        size_t i2 = r1 + i*r_inc;
	        size_t j2 = c1 + j*c_inc;
		    if(i2 > (mSize[0]-1) || j2 > (mSize[1]-1))
		    {
			    throw MatrixException("Tried to access invalid matrix member!");
		    }
		
		    return mData[j2*mSize[0] + i2];
		}
	}
	
	template<class T, size_t dim>
	matrix<T,2> matrix<T,dim,DENSE_MATRIX>::operator()(const size_t &i, const span &cols)
	{
	    static_assert(dim == 2,"This function is a specialization for matrices of dimension 2.");
	    
		vector_k<T> temp(&mData[0],mData.size(),1);
		
		vector_k<span> span_all(dim);
		span_all[0] = span(i,i);
		span_all[1] = cols;
		
		matrix<T,2> submat_return(true, temp, mSize[0], mSize[1], span_all);
		
		return submat_return;
	}
	
	template<class T, size_t dim>
	matrix<T,2> matrix<T,dim,DENSE_MATRIX>::operator()(const span &rows, const size_t &j)
	{
	    static_assert(dim == 2,"This function is a specialization for matrices of dimension 2.");
	    
		vector_k<T> temp(&mData[0],mData.size(),1);
		
		vector_k<span> span_all(dim);
		span_all[0] = rows;
		span_all[1] = span(j,j);
		
		matrix<T,2> submat_return(true, temp, mSize[0], mSize[1], span_all);
		
		return submat_return;
	}
	
	template<class T, size_t dim>
	matrix<T,2> matrix<T,dim,DENSE_MATRIX>::operator()(const span &rows, const span &cols)
	{
	    static_assert(dim == 2,"This function is a specialization for matrices of dimension 2.");
	    
		vector_k<T> temp(&mData[0],mData.size(),1);
		
		vector_k<span> span_all(dim);
		span_all[0] = rows;
		span_all[1] = cols;
		
		matrix<T,2> submat_return(true, temp, mSize[0], mSize[1], span_all);
		
		return submat_return;
	}
	
	template<class T, size_t dim>
	matrix<T,2> matrix<T,dim,DENSE_MATRIX>::operator()(const size_t &i, const span &cols) const
	{
	    static_assert(dim == 2,"This function is a specialization for matrices of dimension 2.");
	    
	    size_t c1, c2;
	    int c_inc;
	    
	    int N_j;
	    if(cols.get_all())
	    {
	        N_j = mSize[1];
	        c1 = 0;
	        c2 = N_j-1;
	        c_inc = 1;
	    }
	    else
	    {
	        c1 = cols.get_first();
	        c2 = cols.get_last();
	        c_inc = cols.get_inc();
	        N_j = (c2 - c1)/c_inc + 1;
	    }
	    if(N_j < 0)
	    {
	        N_j = 0;
	        matrix<T,2> B;
	        return B;
	    }
	    matrix<size_t> j(N_j);
	    for(int jj = 0; jj < N_j; jj++)
	    {
	        j(jj) = c1 + jj*c_inc;
	    }
		if(i > (mSize[0]-1) || max(j) > (mSize[1]-1))
		{
			throw MatrixException("Tried to access invalid matrix member!");
		}
	    
	    matrix<T,2> B(1,N_j);
        for(size_t jj = 0; jj < (size_t)N_j; jj++)
        {
            B(jj) = mData[j(jj)*mSize[0] + i];
        }
		
		return B;
	}
	
	template<class T, size_t dim>
	matrix<T,2> matrix<T,dim,DENSE_MATRIX>::operator()(const span &rows, const size_t &j) const
	{
	    static_assert(dim == 2,"This function is a specialization for matrices of dimension 2.");
	    
	    size_t r1, r2;
	    int r_inc;
	    
	    int N_i;
	    if(rows.get_all())
	    {
	        N_i = mSize[0];
	        r1 = 0;
	        r2 = N_i-1;
	        r_inc = 1;
	    }
	    else
	    {
	        r1 = rows.get_first();
	        r2 = rows.get_last();
	        r_inc = rows.get_inc();
	        N_i = (r2 - r1)/r_inc + 1;
	    }
	    if(N_i < 0)
	    {
	        N_i = 0;
	        matrix<T,2> B;
	        return B;
	    }
	    matrix<size_t> i(N_i);
	    for(int ii = 0; ii < N_i; ii++)
	    {
	        i(ii) = r1 + ii*r_inc;
	    }
		if(max(i) > (mSize[0]-1) || j > (mSize[1]-1))
		{
			throw MatrixException("Tried to access invalid matrix member!");
		}
	    
	    matrix<T,2> B(N_i,1);
	    for(size_t ii = 0; ii < (size_t)N_i; ii++)
	    {
	        B(ii) = mData[j*mSize[0] + i(ii)];
	    }
		
		return B;
	}
	
	template<class T, size_t dim>
	matrix<T,2> matrix<T,dim,DENSE_MATRIX>::operator()(const span &rows, const span &cols) const
	{
	    static_assert(dim == 2,"This function is a specialization for matrices of dimension 2.");
	    
	    size_t r1, r2, c1, c2;
	    int r_inc, c_inc;
	    
	    int N_i;
	    if(rows.get_all())
	    {
	        N_i = mSize[0];
	        r1 = 0;
	        r2 = N_i-1;
	        r_inc = 1;
	    }
	    else
	    {
	        r1 = rows.get_first();
	        r2 = rows.get_last();
	        r_inc = rows.get_inc();
	        N_i = (r2 - r1)/r_inc + 1;
	    }
	    if(N_i < 0)
	    {
	        N_i = 0;
	        matrix<T,2> B;
	        return B;
	    }
	    int N_j;
	    if(cols.get_all())
	    {
	        N_j = mSize[1];
	        c1 = 0;
	        c2 = N_j-1;
	        c_inc = 1;
	    }
	    else
	    {
	        c1 = cols.get_first();
	        c2 = cols.get_last();
	        c_inc = cols.get_inc();
	        N_j = (c2 - c1)/c_inc + 1;
	    }
	    if(N_j < 0)
	    {
	        N_j = 0;
	        matrix<T,2> B;
	        return B;
	    }
	    matrix<size_t> i(N_i), j(N_j);
	    for(int ii = 0; ii < N_i; ii++)
	    {
	        i(ii) = r1 + ii*r_inc;
	    }
	    for(int jj = 0; jj < N_j; jj++)
	    {
	        j(jj) = c1 + jj*c_inc;
	    }
		if(max(i) > (mSize[0]-1) || max(j) > (mSize[1]-1))
		{
			throw MatrixException("Tried to access invalid matrix member!");
		}
	    
	    matrix<T,2> B(N_i,N_j);
	    for(size_t ii = 0; ii < (size_t)N_i; ii++)
	    {
	        for(size_t jj = 0; jj < (size_t)N_j; jj++)
	        {
	            B(ii,jj) = mData[j(jj)*mSize[0] + i(ii)];
	        }
	    }
		
		return B;
	}

	template<class T,size_t dim>
	T& matrix<T,dim,DENSE_MATRIX>::operator()(const size_t &i)
	{
	    if(i > this->numel()-1)
	    {
		    throw MatrixException("Tried to access invalid matrix member!");
	    }
	    if(!submat)
	    {
	        return mData[i];
	    }
	    else
	    {
	        size_t ii,jj;
	        ii = i % this->size(1);
	        jj = i/this->size(1);
	        return this->operator()(ii,jj);
	    }
	}

	template<class T, size_t dim>
	const T& matrix<T,dim,DENSE_MATRIX>::operator()(const size_t &i) const
	{
	    if(i > this->numel()-1)
	    {
		    throw MatrixException("Tried to access invalid matrix member!");
	    }
	    if(!submat)
	    {
	        return mData[i];
	    }
	    else
	    {
	        size_t ii,jj;
	        ii = i % this->size(1);
	        jj = i/this->size(1);
	        return this->operator()(ii,jj);
	    }
	}
	
	template<class T, size_t dim>
	template<class U>
	matrix<T,dim>& matrix<T,dim,DENSE_MATRIX>::operator+=(const matrix<U,dim> &B)
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
			this->mData[ii] += B(ii);
		}
		return *this;
	}
	
	template<class T, size_t dim>
	template<class U>
	matrix<T,dim>& matrix<T,dim,DENSE_MATRIX>::operator-=(const matrix<U,dim> &B)
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
			this->mData[ii] -= B(ii);
		}
		return *this;
	}

	template<class T, size_t dim>
	bool matrix<T,dim,DENSE_MATRIX>::operator==(const matrix<T,dim> &B) const
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
			if(this->mData[ii] != B(ii))
			{
			    return false;
			}
		}
		return true;
	}
	
	
	template<class T, size_t dim>
	template<class U>
	matrix<decltype(std::declval<T>()*std::declval<U>()),dim>& matrix<T,dim,DENSE_MATRIX>::trans_mult(const matrix<U,dim> &B) const
	{
	    static_assert(dim == 2, "This method is only supported for matrices of dimension 2.");
		if(this->empty() || B.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(this->size(1) != B.size(1))
		{
			throw MatrixException("Matrix dimensions are not compatible in matrix-matrix multiplication!");
		}
		matrix<decltype(std::declval<T>()*std::declval<U>()),2> C(this->size(2),B.size(2));
		for(size_t kk = 0; kk < B.size(2); kk++)
		{
			for(size_t ii = 0; ii < this->size(2); ii++)
			{
				C(ii,kk) = 0.0;
				for(size_t jj = 0; jj < this->size(1); jj++)
				{
					C(ii,kk) += this->operator()(jj,ii)*B(jj,kk);
				}
			}
		}
		return C;
	}

	template<class T, size_t dim>
	bool matrix<T,dim,DENSE_MATRIX>::operator!=(const matrix<T,dim> &B) const
	{
		return !(*this == B);
	}

    template<class T, size_t dim>
    template<class U>
    matrix<T,dim>& matrix<T,dim,DENSE_MATRIX>::operator=(const matrix<U,dim,DENSE_MATRIX> &v)
    {
        if(this->empty())
        {
            matrix<size_t> temp(dim);
            for(size_t ii = 0; ii < dim; ii++)
            {
                temp(ii) = v.size(ii+1);
            }
            resize(temp);
        }
        if(!submat && !v.get_submat())
        {
            if(mSize != v.get_mSize() && !(numel() == v.numel() && isVec() == v.isVec()))
            {
                resize(v.get_mSize());
            }
            for(size_t ii = 0; ii < v.numel(); ii++)
            {
                mData[ii] = v(ii);
            }
        }
        else
        {
            for(size_t ii = 0; ii < dim; ii++)
            {
                if(this->size(ii+1) != v.size(ii+1))
                {
                    throw MatrixException("Size of submatrix is not compatible in assignment!");
                }
            }
            for(size_t ii = 0; ii < v.numel(); ii++)
            {
                this->operator()(ii) = v(ii);
            }
        }
        
        return *this;
    }

    template<class T, size_t dim>
    template<class U>
    matrix<T,dim>& matrix<T,dim,DENSE_MATRIX>::operator=(const vector_k<U> &v)
    {
        if(this->empty())
        {
            resize(v.size());
        }
        if(!submat)
        {
            if(numel() != v.size())
            {
                resize(v.size());
            }
            for(size_t ii = 0; ii < v.size(); ii++)
            {
                mData[ii] = v[ii];
            }
        }
        else
        {
            if(this->numel() != v.size())
            {
                throw MatrixException("Size of submatrix is not compatible in assignment!");
            }
            
            for(size_t ii = 0; ii < v.size(); ii++)
            {
                this->operator()(ii) = v[ii];
            }
        }
        
        return *this;
    }
    
    template<class T, size_t dim>
    matrix<T,dim>& matrix<T,dim,DENSE_MATRIX>::operator=(const matrix<T,dim,DENSE_MATRIX> &v)
    {
        if(this == &v)
        {
            return *this;
        }
        if(this->empty())
        {
            matrix<size_t> temp(dim);
            for(size_t ii = 0; ii < dim; ii++)
            {
                temp(ii) = v.size(ii+1);
            }
            resize(temp);
        }
        if(!submat && !v.get_submat())
        {
            if(mSize != v.get_mSize() && !(numel() == v.numel() && isVec() == v.isVec()))
            {
                resize(v.get_mSize());
            }
            for(size_t ii = 0; ii < v.numel(); ii++)
            {
                mData[ii] = v(ii);
            }
        }
        else
        {
            for(size_t ii = 0; ii < dim; ii++)
            {
                if(this->size(ii+1) != v.size(ii+1))
                {
                    throw MatrixException("Size of submatrix is not compatible in assignment!");
                }
            }
            for(size_t ii = 0; ii < v.numel(); ii++)
            {
                this->operator()(ii) = v(ii);
            }
        }
        
        return *this;
    }
    
    template<class T, size_t dim>
    matrix<T,dim>& matrix<T,dim,DENSE_MATRIX>::operator=(const vector_k<T> &v)
    {
        if(this->empty())
        {
            resize(v.size());
        }
        if(!submat)
        {
            if(numel() != v.size())
            {
                resize(v.size());
            }
            for(size_t ii = 0; ii < v.size(); ii++)
            {
                mData[ii] = v[ii];
            }
        }
        else
        {
            if(this->numel() != v.size())
            {
                throw MatrixException("Size of submatrix is not compatible in assignment!");
            }
            
            for(size_t ii = 0; ii < v.size(); ii++)
            {
                this->operator()(ii) = v[ii];
            }
        }
        
        return *this;
    }
    
    template<class T, size_t dim>
    template<class U>
    matrix<T,dim>& matrix<T,dim,DENSE_MATRIX>::operator=(const U &a)
    {
        if(this->empty())
        {
            matrix<size_t> temp(dim);
            for(size_t ii = 0; ii < dim; ii++)
            {
                temp(ii) = 1;
            }
            resize(temp);
        }
        if(!submat)
        {
            for(size_t ii = 0; ii < numel(); ii++)
            {
                mData[ii] = a;
            }
        }
        else
        {
            for(size_t ii = 0; ii < numel(); ii++)
            {
                this->operator()(ii) = a;
            }
        }
        
        return *this;
    }
    
    template<class T, size_t dim>
    matrix<T,dim>& matrix<T,dim,DENSE_MATRIX>::operator=(const T &a)
    {
        if(this->empty())
        {
            matrix<size_t> temp(dim);
            for(size_t ii = 0; ii < dim; ii++)
            {
                temp(ii) = 1;
            }
            resize(temp);
        }
        if(!submat)
        {
            for(size_t ii = 0; ii < numel(); ii++)
            {
                mData[ii] = a;
            }
        }
        else
        {
            for(size_t ii = 0; ii < numel(); ii++)
            {
                this->operator()(ii) = a;
            }
        }
        
        return *this;
    }

	template<class T, size_t dim>
	size_t matrix<T,dim,DENSE_MATRIX>::size(const size_t &n) const
	{
	    if(n > dim || n == 0)
	    {
		    throw MatrixException("Invalid dimension in size().");
	    }
	    if(!submat)
	    {
		    return mSize[n-1];
		}
		else
		{
		    int N_i;
		    if(mSpan[n-1].get_all())
	        {
	            N_i = mSize[n-1];
	        }
	        else
	        {
	            int r1 = mSpan[n-1].get_first();
	            int r2 = mSpan[n-1].get_last();
	            int r_inc = mSpan[n-1].get_inc();
	            N_i = (r2 - r1)/r_inc + 1;
	        }
	        if(N_i < 0)
	        {
	            N_i = 0;
	        }
	        return ((size_t)N_i);
		}
	}
	
	template<class T, size_t dim>
	void matrix<T,dim,DENSE_MATRIX>::resize(const std::array<size_t,dim> &pSize)
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
        
		for(size_t index_old = 0; index_old < B.numel(); index_old++)
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
            
            mData[index_new] = B(index_old);
        }
	}
	
	template<class T, size_t dim>
	void matrix<T,dim,DENSE_MATRIX>::resize(const matrix<size_t> &pSize)
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
        
		for(size_t index_old = 0; index_old < B.numel(); index_old++)
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
            
            mData[index_new] = B(index_old);
        }
	}
       
    template<class T, size_t dim>
    void matrix<T,dim,DENSE_MATRIX>::resize(const size_t &pSize)
    {
        static_assert(dim == 2,"This version of resize is not supported for dimensions greater than 2.");
        if(!this->isVec())
        {
            mSize[0] = pSize;
            if(mSize[1] == 0)
            {
                mSize[1] = 1;
            }
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
	bool matrix<T,dim,DENSE_MATRIX>::empty() const
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
	matrix<T,2> matrix<T,dim,DENSE_MATRIX>::row(const size_t &n)
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
		
		matrix<T,2> Row(true, temp, 1, mSize[1]);
		return Row;
	}

	template<class T, size_t dim>
	matrix<T,2> matrix<T,dim,DENSE_MATRIX>::col(const size_t &n)
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
		
		matrix<T,2> Col(true, temp, mSize[0], 1);
		
		return Col;
	}

	template<class T, size_t dim>
	matrix<T,2> matrix<T,dim,DENSE_MATRIX>::row(const size_t &n) const
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
	matrix<T,2> matrix<T,dim,DENSE_MATRIX>::col(const size_t &n) const
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
		matrix<T,2> Col(mSize[0],1);

		for(size_t ii = 0; ii < mSize[0]; ii++)
		{
			Col(ii) = mData[n*mSize[0] + ii];
		}

		return Col;
	}
	
	template<class T, size_t dim>
	void matrix<T,dim,DENSE_MATRIX>::reserve(const size_t &N)
	{
		mData.reserve(N);
	}
	
    // Matrix Operations:

	template<class V,class U>
	matrix<decltype(std::declval<V>()*std::declval<U>()),2> operator*(const matrix<V,2> &A, const matrix<U,2> &B)
	{
		if(A.empty() || B.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(A.size(2) != B.size(1))
		{
			throw MatrixException("Matrix dimensions are not compatible in matrix-matrix multiplication!");
		}
		matrix<decltype(std::declval<V>()*std::declval<U>()),2> C(A.size(1),B.size(2));
		for(size_t kk = 0; kk < B.size(2); kk++)
		{
			for(size_t ii = 0; ii < A.size(1); ii++)
			{
				C(ii,kk) = 0.0;
				for(size_t jj = 0; jj < A.size(2); jj++)
				{
					C(ii,kk) += A(ii,jj)*B(jj,kk);
				}
			}
		}
		return C;
	}
	
	template<>
	inline matrix<double,2> operator*(const matrix<double,2> &A, const matrix<double,2> &B)
	{
		if(A.empty() || B.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(A.size(2) != B.size(1))
		{
			throw MatrixException("Matrix dimensions are not compatible in matrix-matrix multiplication!");
		}
		
		matrix<double> C(A.size(1),B.size(2));
		
		if(A.get_inc() == 1 && B.get_inc() == 1 && !A.get_submat() && !B.get_submat())
		{
		    int m = (int)A.size(1);
		    int k = (int)A.size(2);
		    int n = (int)B.size(2);
            char TRANS = 'N';
            double ALPHA = 1.0;
            int LDA = m;
            int LDB = k;
            double BETA = 0.0;
            int LDC = m;

            dgemm_(&TRANS, &TRANS, &m, &n, &k, &ALPHA, &A.mData[0], &LDA, &B.mData[0], &LDB, &BETA, &C.mData[0], &LDC);
        }
        else
        {
            for(size_t kk = 0; kk < B.size(2); kk++)
		    {
			    for(size_t ii = 0; ii < A.size(1); ii++)
			    {
				    C(ii,kk) = 0.0;
				    for(size_t jj = 0; jj < A.size(2); jj++)
				    {
					    C(ii,kk) += A(ii,jj)*B(jj,kk);
				    }
			    }
		    }
        }
		return C;
    }
	
	template<>
	inline matrix<std::complex<double>,2> operator*(const matrix<std::complex<double>,2> &A, const matrix<std::complex<double>,2> &B)
	{
		if(A.empty() || B.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(A.size(2) != B.size(1))
		{
			throw MatrixException("Matrix dimensions are not compatible in matrix-matrix multiplication!");
		}
		
		matrix<std::complex<double>> C(A.size(1),B.size(2));
		if(A.get_inc() == 1 && B.get_inc() == 1 && !A.get_submat() && !B.get_submat())
		{
		    int m = (int)A.size(1);
		    int k = (int)A.size(2);
		    int n = (int)B.size(2);
            char TRANS = 'N';
            std::complex<double> ALPHA = 1.0;
            int LDA = m;
            int LDB = k;
            std::complex<double> BETA = 0.0;
            int LDC = m;

            zgemm_(&TRANS, &TRANS, &m, &n, &k, &ALPHA, &A.mData[0], &LDA, &B.mData[0], &LDB, &BETA, &C.mData[0], &LDC);
        }
        else
        {
            for(size_t kk = 0; kk < B.size(2); kk++)
		    {
			    for(size_t ii = 0; ii < A.size(1); ii++)
			    {
				    C(ii,kk) = 0.0;
				    for(size_t jj = 0; jj < A.size(2); jj++)
				    {
					    C(ii,kk) += A(ii,jj)*B(jj,kk);
				    }
			    }
		    }
        }
		return C;
    }
    
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
		if(A.empty() || B.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		for(int ii = 0; ii < dim; ii++)
		{
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
		if(A.empty() || B.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		for(int ii = 0; ii < dim; ii++)
		{
		    if(A.size(ii+1) != B.size(ii+1))
		    {
		        throw MatrixException("Matrix dimensions are not compatible in matrix-matrix subtraction!");
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
		matrix<T,dim> B;
		B.resize(size(A));
		for(size_t ii = 0; ii < B.numel(); ii++)
		{
		    B(ii) = A(ii);
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
                out << A(ii) << std::endl;
            }
        }
        
        return out;
    }
}

#endif
