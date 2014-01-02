// SparseMatrix.h -- matrix class

#ifndef SPARSE_MATRIX_H_
#define SPARSE_MATRIX_H_

#include <iostream>
#include <stdexcept>
#include <vector>
#include <stdarg.h>
#include "vector_k.h"
#include "Eigen/Sparse"

#define DENSE_MATRIX 0
#define SPARSE_MATRIX 1

namespace keycpp
{
    template<class T> void disp(const T &x, std::ostream& outStream = std::cout);
    template<class T> struct Sort_Matrix;
	template<class T> Sort_Matrix<T> sort(const matrix<T> &A, const size_t &dim = 2, std::string method = "ascend");
	
    template<class T>
    matrix<T,2,SPARSE_MATRIX> sparse(const matrix<T,2,0> &A)
    {
        matrix<T,2,SPARSE_MATRIX> B(A.size(1),A.size(2));
        for(size_t ii = 0; ii < A.size(1); ii++)
        {
            for(size_t jj = 0; jj < A.size(2); jj++)
            {
                if(abs(A(ii,jj)) > eps)
                {
                    B(ii,jj) = A(ii,jj);
                }
            }
        }
        B.update();
        return B;
    }

	class SparseMatrixException : public std::runtime_error
	{
	    public:
		SparseMatrixException(const std::string &msg) : std::runtime_error(msg){}
	};
	
	template<class T>
	class matrix<T,2,SPARSE_MATRIX>
	{
	public:	
		matrix();
		matrix(const matrix<T,2,SPARSE_MATRIX> &A);
		matrix(const size_t &d1, const size_t &d2);
		matrix(const std::initializer_list<std::initializer_list<T> >& lst);
		matrix(const std::initializer_list<T>& lst);
		void update() const;
		T& operator()(const size_t &i, const size_t &j);
		T operator()(const size_t &i, const size_t &j) const;
		T& operator()(const size_t &i);
		T operator()(const size_t &i) const;
		matrix<T,2,SPARSE_MATRIX> operator()(const size_t &i, const span &cols) const;
		matrix<T,2,SPARSE_MATRIX> operator()(const span &rows, const size_t &j) const;
		matrix<T,2,SPARSE_MATRIX> operator()(const span &rows, const span &cols) const;
		matrix<T,2,SPARSE_MATRIX> row(const size_t &i) const;
		matrix<T,2,SPARSE_MATRIX> col(const size_t &j) const;
		bool operator!=(const matrix<T,2,SPARSE_MATRIX> &B) const;
		bool operator==(const matrix<T,2,SPARSE_MATRIX> &B) const;
        template<class U>
		matrix<T,2,SPARSE_MATRIX>& operator=(const matrix<U,2,SPARSE_MATRIX> &v);
		size_t size(const size_t &n) const;
		bool empty() const;
		void reserve(const size_t &N);
        size_t length() const
        {
            return *std::max_element(mSize.begin(), mSize.end());
        };
        size_t numel() const
        {
            return mSize[0]*mSize[1];
        };
        
        bool isVec() const
        {
            if(mSize[0] == 1 || mSize[1] == 1)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        
        template<class V, class U>
	    friend matrix<decltype(std::declval<V>()*std::declval<U>()),2,0> operator*(const matrix<V,2,SPARSE_MATRIX> &A, const matrix<U,2,0> &x);
	    template<class V, class U>
	    friend matrix<decltype(std::declval<V>()*std::declval<U>()),2,SPARSE_MATRIX> operator*(const matrix<V,2,SPARSE_MATRIX> &A, const matrix<U,2,SPARSE_MATRIX> &B);
	    template<class V>
        friend matrix<V,2,0> full(const matrix<V,2,SPARSE_MATRIX> &A);
	    template<class V, class U>
	    friend matrix<decltype(std::declval<V>() + std::declval<U>()),2,SPARSE_MATRIX> operator+(const matrix<V,2,SPARSE_MATRIX> &A, const matrix<U,2,SPARSE_MATRIX> &B);
	    template<class V, class U>
	    friend matrix<decltype(std::declval<V>() + std::declval<U>()),2,SPARSE_MATRIX> operator-(const matrix<V,2,SPARSE_MATRIX> &A, const matrix<U,2,SPARSE_MATRIX> &B);
	    template<class V>
	    friend matrix<V,2,SPARSE_MATRIX> transpose(const matrix<V,2,SPARSE_MATRIX> &A);
	    friend void mv(int n, std::complex<double> *in, std::complex<double> *out, const matrix<std::complex<double>,2,SPARSE_MATRIX> &A);
	    friend void znaupd(int n, int nev, matrix<std::complex<double>> &Evals, std::string which, const matrix<std::complex<double>,2,SPARSE_MATRIX> &A, const matrix<std::complex<double>,2,SPARSE_MATRIX> &B);
	    friend void znaupd(int n, int nev, matrix<std::complex<double>> &Evals, matrix<std::complex<double>> &Evecs, std::string which, const matrix<std::complex<double>,2,SPARSE_MATRIX> &A, const matrix<std::complex<double>,2,SPARSE_MATRIX> &B);
        
        operator T()
        {
            if(mData.size() != 1)
            {
                throw MatrixException("Cannot convert to T. Number of elements must be 1.");
            }

            return mData[0];
        };
        
        operator matrix<T>()
        {
            return full(*this);
        };
        
    private:
		mutable vector_k<T> mData;
		mutable vector_k<size_t> rowInd;
		mutable vector_k<size_t> colPtr;
		std::array<size_t,2> mSize;
		
		mutable vector_k<Eigen::Triplet<T>> data_unsorted;
		mutable bool m_update = false;
	};

	template<class T>
	matrix<T,2,SPARSE_MATRIX>::matrix() : mData(), rowInd(), colPtr()
	{
	    for(size_t ii = 0; ii < 2; ii++)
	    {
	        mSize[ii] = 0;
	    }
	}

	template<class T>
	matrix<T,2,SPARSE_MATRIX>::matrix(const matrix<T,2,SPARSE_MATRIX> &A) : mData(), rowInd(), colPtr()
	{
	    for(size_t ii = 0; ii < 2; ii++)
	    {
	        mSize[ii] = A.mSize[ii];
	    }
        mData.resize(A.mData.size());
        rowInd.resize(A.mData.size());
	    for(size_t ii = 0; ii < A.mData.size(); ii++)
	    {
	        mData[ii] = A.mData[ii];
	        rowInd[ii] = A.rowInd[ii];
	    }
	    colPtr.resize(mSize[1]+1);
		for(size_t ii = 0; ii < mSize[1]+1; ii++)
		{
			colPtr[ii] = A.colPtr[ii];
		}
	}

	template<class T>
	matrix<T,2,SPARSE_MATRIX>::matrix(const size_t & d1, const size_t & d2) : mData(), rowInd(), colPtr()
	{
		mSize[0] = d1;
		mSize[1] = d2;
		
		colPtr.resize(mSize[1]+1);
	}

	template<class T>
	matrix<T,2,SPARSE_MATRIX>::matrix(const std::initializer_list<std::initializer_list<T> >& lst)
	{
		*this = sparse(matrix<T,2,0>(lst));
	}

	template<class T>
	matrix<T,2,SPARSE_MATRIX>::matrix(const std::initializer_list<T>& lst)
	{
		*this = sparse(matrix<T,2,0>(lst));
	}
	
	template<class T>
	int compare(const void * a, const void * b)
	{
	    Eigen::Triplet<T> A = *(Eigen::Triplet<T> *)a;
	    Eigen::Triplet<T> B = *(Eigen::Triplet<T> *)b;
	    
	    if(A.col() < B.col())
	    {
	        return -1;
	    }
	    else if(A.col() == B.col() && A.row() < B.row())
	    {
	        return -1;
	    }
	    else
	    {
	        return 1;
	    }
	}
	
	template<class T>
	void matrix<T,2,SPARSE_MATRIX>::update() const
	{
	    if(!m_update)
	    {
	        return;
	    }
	    
		mData.resize(data_unsorted.size());
		rowInd.resize(data_unsorted.size());
		
		std::qsort(&data_unsorted[0], data_unsorted.size(), sizeof(Eigen::Triplet<T>), compare<T>);
		
	    for(size_t jj = 0; jj < data_unsorted.size(); jj++)
	    {
	        mData[jj] = data_unsorted[jj].value();
	        rowInd[jj] = data_unsorted[jj].row();
	    }
	    colPtr[0] = 0;
	    for(size_t jj = 1; jj < colPtr.size()-1; jj++)
	    {
	        for(size_t ii = colPtr[jj-1]+1; ii < data_unsorted.size(); ii++)
	        {
                if(data_unsorted[ii].col() > data_unsorted[ii-1].col())
                {
                    colPtr[jj] = ii;
                    break;
                }
            }
	    }
	    colPtr[colPtr.size()-1] = data_unsorted.size();
	    m_update = false;
	    
		return;
	}

	template<class T>
	T& matrix<T,2,SPARSE_MATRIX>::operator()(const size_t &i, const size_t &j)
	{
		if(i > (mSize[0]-1) || j > (mSize[1]-1))
		{
			throw MatrixException("Tried to access invalid matrix member!");
		}
		
		for(size_t ii = colPtr[j]; ii < colPtr[j+1]; ii++)
		{
		    if(rowInd[ii] == i)
		    {
		        return mData[ii];
		    }
		}
		
		for(size_t ii = 0; ii < data_unsorted.size(); ii++)
		{
		    if(data_unsorted[ii].row() == i && data_unsorted[ii].col() == j)
		    {
		        return data_unsorted[ii].value();
		    }
		}
		
		m_update = true;
		
		Eigen::Triplet<T> temp;
		temp.row() = i;
		temp.col() = j;
		temp.value() = T();
		
		data_unsorted.push_back(temp);
		
		return data_unsorted[data_unsorted.size()-1].value();
	}

	template<class T>
	T matrix<T,2,SPARSE_MATRIX>::operator()(const size_t &i, const size_t &j) const
	{
		if(i > (mSize[0]-1) || j > (mSize[1]-1))
		{
			throw MatrixException("Tried to access invalid matrix member!");
		}
		
		for(size_t ii = colPtr[j]; ii < colPtr[j+1]; ii++)
		{
		    if(rowInd[ii] == i)
		    {
		        return mData[ii];
		    }
		}
		
		for(size_t ii = 0; ii < data_unsorted.size(); ii++)
		{
		    if(data_unsorted[ii].row() == i && data_unsorted[ii].col() == j)
		    {
		        return data_unsorted[ii].value();
		    }
		}
		
		return T();
	}
	
	template<class T>
	matrix<T,2,SPARSE_MATRIX> matrix<T,2,SPARSE_MATRIX>::operator()(const size_t &i, const span &cols) const
	{
	    update();
	    size_t c1, c2;
	    int c_inc;
	    
	    int N_j;
	    if(cols.get_all())
	    {
	        N_j = mSize[0];
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
	        matrix<T,2,SPARSE_MATRIX> B;
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
	    
	    matrix<T,2,SPARSE_MATRIX> B(1,N_j);
        for(size_t n = 0; n < N_j; n++)
        {
		    for(size_t ii = colPtr[j(n)]; ii < colPtr[j(n)+1]; ii++)
		    {
		        if(rowInd[ii] == i)
		        {
			        B(j(n)) = mData[ii];
			    }
		    }
		}
		B.update();
		
		return B;
	}
	
	template<class T>
	matrix<T,2,SPARSE_MATRIX> matrix<T,2,SPARSE_MATRIX>::operator()(const span &rows, const size_t &j) const
	{
	    update();
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
	        matrix<T,2,SPARSE_MATRIX> B;
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
	    
	    matrix<T,2,SPARSE_MATRIX> B(N_i,1);
	    for(size_t n = 0; n < N_i; n++)
	    {
	        for(size_t ii = colPtr[j]; ii < colPtr[j+1]; ii++)
	        {
	            if(i(n) == rowInd[ii])
	            {
		            B(rowInd[ii]) = mData[ii];
		        }
	        }
	    }
		B.update();
		
		return B;
	}
	
	template<class T>
	matrix<T,2,SPARSE_MATRIX> matrix<T,2,SPARSE_MATRIX>::operator()(const span &rows, const span &cols) const
	{
	    update();
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
	        matrix<T,2,SPARSE_MATRIX> B;
	        return B;
	    }
	    int N_j;
	    if(cols.get_all())
	    {
	        N_j = mSize[0];
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
	        matrix<T,2,SPARSE_MATRIX> B;
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
	    
	    matrix<T,2,SPARSE_MATRIX> B(N_i,N_j);
	    for(size_t n = 0; n < N_i; n++)
	    {
	        for(size_t jj = 0; jj < N_j; jj++)
	        {
	            for(size_t ii = colPtr[j(jj)]; ii < colPtr[j(jj)+1]; ii++)
	            {
	                if(i(n) == rowInd[ii])
	                {
		                B(rowInd[ii],j(jj)) = mData[ii];
		            }
	            }
	        }
	    }
		B.update();
		
		return B;
	}

	template<class T>
	T& matrix<T,2,SPARSE_MATRIX>::operator()(const size_t &i)
	{
		return this->operator()(i % mSize[0], i/mSize[0]);
	}

	template<class T>
	T matrix<T,2,SPARSE_MATRIX>::operator()(const size_t &i) const
	{
		return this->operator()(i % mSize[0], i/mSize[0]);
	}

	template<class T>
	matrix<T,2,SPARSE_MATRIX> matrix<T,2,SPARSE_MATRIX>::row(const size_t &n) const
	{
	    update();
		if(n > mSize[0])
		{
			throw MatrixException("Invalid row index in row().");
		}
		matrix<T,2,SPARSE_MATRIX> Row(1,mSize[1]);
        
        for(size_t ii = 0; ii < mSize[1]; ii++)
        {
		    for(size_t jj = colPtr[ii]; jj < colPtr[ii+1]; jj++)
		    {
		        if(rowInd[jj] == n)
		        {
			        Row(ii) = mData[jj];
			    }
		    }
		}
		Row.update();

		return Row;
	}

	template<class T>
	matrix<T,2,SPARSE_MATRIX> matrix<T,2,SPARSE_MATRIX>::col(const size_t &n) const
	{
	    update();
		if(n > mSize[1])
		{
			throw MatrixException("Invalid column index in col().");
		}
		matrix<T,2,SPARSE_MATRIX> Col(mSize[0],1);

		for(size_t ii = colPtr[n]; ii < colPtr[n+1]; ii++)
		{
			Col(rowInd[ii]) = mData[ii];
		}
		Col.update();

		return Col;
	}

	template<class T>
	bool matrix<T,2,SPARSE_MATRIX>::operator!=(const matrix<T,2,SPARSE_MATRIX> &B) const
	{
	    update();
		if(mData.empty() || B.empty())
		{
			return true;
		}
		for(int ii = 0; ii < 2; ii++)
		{
		    if(mSize[ii] != B.size(ii+1))
		    {
		        return true;
		    }
		}
		for(size_t ii = 0; ii < mData.size(); ii++)
		{
			if(this->mData[ii] != B.mData[ii] || this->rowInd[ii] != B.rowInd[ii])
			{
			    return true;
			}
		}
		for(size_t ii = 0; ii < mSize[1]+1; ii++)
		{
			if(this->colPtr[ii] != B.colPtr[ii])
			{
			    return true;
			}
		}
		return false;
	}

	template<class T>
	bool matrix<T,2,SPARSE_MATRIX>::operator==(const matrix<T,2,SPARSE_MATRIX> &B) const
	{
	    update();
		if(mData.empty() || B.empty())
		{
			return false;
		}
		for(int ii = 0; ii < 2; ii++)
		{
		    if(mSize[ii] != B.size(ii+1))
		    {
		        return false;
		    }
		}
		for(size_t ii = 0; ii < mData.size(); ii++)
		{
			if(this->mData[ii] != B.mData[ii] || this->rowInd[ii] != B.rowInd[ii])
			{
			    return false;
			}
		}
		for(size_t ii = 0; ii < mSize[1]+1; ii++)
		{
			if(this->colPtr[ii] != B.colPtr[ii])
			{
			    return false;
			}
		}
		return true;
	}

    template<class T>
    template<class U>
    matrix<T,2,SPARSE_MATRIX>& matrix<T,2,SPARSE_MATRIX>::operator=(const matrix<U,2,SPARSE_MATRIX> &v)
    {
        v.update();
        for(size_t ii = 0; ii < 2; ii++)
	    {
	        mSize[ii] = v.mSize[ii];
	    }
        mData.resize(v.mData.size());
        rowInd.resize(v.mData.size());
        for(size_t ii = 0; ii < v.mData.size(); ii++)
        {
            mData[ii] = v.mData[ii];
            rowInd[ii] = v.rowInd[ii];
        }
	    colPtr.resize(mSize[1]+1);
		for(size_t ii = 0; ii < mSize[1]+1; ii++)
		{
			colPtr[ii] = v.colPtr[ii];
		}
        
        return *this;
    }

	template<class T>
	size_t matrix<T,2,SPARSE_MATRIX>::size(const size_t &n) const
	{
		if(n > 2 || n == 0)
		{
			throw MatrixException("Invalid dimension in size().");
		}
		return mSize[n-1];
	}
	
	template<class T>
	bool matrix<T,2,SPARSE_MATRIX>::empty() const
	{
		for(size_t ii = 0; ii < 2; ii++)
		{
		    if(mSize[ii] == 0)
		    {
			    return true;
			}
		}
		
		return false;
	}
	
	template<class T>
	void matrix<T,2,SPARSE_MATRIX>::reserve(const size_t &N)
	{
		mData.reserve(N);
	}
	
	template<class T, class U>
	matrix<decltype(std::declval<T>()*std::declval<U>()),2,0> operator*(const matrix<T,2,SPARSE_MATRIX> &A, const matrix<U,2,0> &B)
	{
	    A.update();
	    if(A.size(2) != B.size(1))
	    {
	        throw SparseMatrixException("Matrix and vector have incompatible sizes!");
	    }
	    
	    matrix<decltype(std::declval<T>()*std::declval<U>()),2,0> C(A.size(1),B.size(2));
	    for(int kk = 0; kk < B.size(2); kk++)
	    {
	        for(int jj = 0; jj < B.size(1); jj++)
	        {
	            for(int ii = A.colPtr[jj]; ii < A.colPtr[jj+1]; ii++)
	            {
	                C(A.rowInd[ii],kk) += A.mData[ii]*B(jj,kk);
	            }
	        }
	    }
	    return C;
	}
	
	template<class V, class U>
	matrix<decltype(std::declval<V>() + std::declval<U>()),2,SPARSE_MATRIX> operator+(const matrix<V,2,SPARSE_MATRIX> &A, const matrix<U,2,SPARSE_MATRIX> &B)
	{
	    A.update();
	    B.update();
	    if(A.size(1) != B.size(1) && A.size(2) != B.size(2))
	    {
	        throw SparseMatrixException("Matrices have incompatible sizes!");
	    }
	    matrix<decltype(std::declval<V>() + std::declval<U>()),2,SPARSE_MATRIX> C = A;
	
        for(size_t n = 0; n < B.size(2); n++)
        {
		    for(size_t ii = B.colPtr[n]; ii < B.colPtr[n+1]; ii++)
		    {
			    C(B.rowInd[ii],n) += B.mData[ii];
		    }
		}
		C.update();
		
	    return C;
	}
	
	template<class V, class U>
	matrix<decltype(std::declval<V>() + std::declval<U>()),2,SPARSE_MATRIX> operator-(const matrix<V,2,SPARSE_MATRIX> &A, const matrix<U,2,SPARSE_MATRIX> &B)
	{
	    A.update();
	    B.update();
	    if(A.size(1) != B.size(1) && A.size(2) != B.size(2))
	    {
	        throw SparseMatrixException("Matrices have incompatible sizes!");
	    }
	    matrix<decltype(std::declval<V>() + std::declval<U>()),2,SPARSE_MATRIX> C = A;
	
        for(size_t n = 0; n < B.size(2); n++)
        {
		    for(size_t ii = B.colPtr[n]; ii < B.colPtr[n+1]; ii++)
		    {
			    C(B.rowInd[ii],n) -= B.mData[ii];
		    }
		}
		C.update();
	    return C;
	}
	
	template<class V>
	matrix<V,2,SPARSE_MATRIX> transpose(const matrix<V,2,SPARSE_MATRIX> &A)
	{
	    A.update();
	    matrix<V,2,SPARSE_MATRIX> C(A.size(2),A.size(1));
	
        for(size_t n = 0; n < A.size(2); n++)
        {
		    for(size_t ii = A.colPtr[n]; ii < A.colPtr[n+1]; ii++)
		    {
			    C(n,A.rowInd[ii]) = A.mData[ii];
		    }
		}
		C.update();
	    return C;
	}
	
	template<class T, class U>
	matrix<decltype(std::declval<T>()*std::declval<U>()),2,SPARSE_MATRIX> operator*(const matrix<T,2,SPARSE_MATRIX> &A, const matrix<U,2,SPARSE_MATRIX> &B)
	{
	    A.update();
	    B.update();
	    if(A.size(2) != B.size(1))
	    {
	        throw SparseMatrixException("Matrix and vector have incompatible sizes!");
	    }
	    
	    matrix<decltype(std::declval<T>()*std::declval<U>()),2,SPARSE_MATRIX> C(A.size(1),B.size(2));
	    auto B_temp = transpose(B);
	    for(int n = 0; n < A.size(2); n++)
	    {
            for(int jj = B_temp.colPtr[n]; jj < B_temp.colPtr[n+1]; jj++)
            {
	            for(int ii = A.colPtr[n]; ii < A.colPtr[n+1]; ii++)
	            {
	                C(A.rowInd[ii],B_temp.rowInd[jj]) += A.mData[ii]*B_temp.mData[jj];
	            }
	        }
	    }
		C.update();
	    
	    return C;
	}
	
	template<class T> T mat2num(const matrix<T,2,SPARSE_MATRIX>& A)
	{
	    A.update();
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
	
	template<class T>
	std::ostream& operator<<(std::ostream &out, const matrix<T,2,SPARSE_MATRIX> &A)
	{
	    A.update();
        for(size_t ii = 0; ii < A.size(1); ii++)
        {
            for(size_t jj = 0; jj < A.size(2); jj++)
            {
                out << A(ii,jj) << " ";
            }
            out << std::endl;
        }
        
        return out;
    }
    
    template<class T>
    matrix<T,2,0> full(const matrix<T,2,SPARSE_MATRIX> &A)
    {
        A.update();
        matrix<T,2,0> B(A.size(1),A.size(2));
        
        for(size_t n = 0; n < A.size(2); n++)
        {
		    for(size_t ii = A.colPtr[n]; ii < A.colPtr[n+1]; ii++)
		    {
			    B(A.rowInd[ii],n) = A.mData[ii];
		    }
		}
		
        return B;
    }
}

#endif
