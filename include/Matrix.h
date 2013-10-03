// matrix.h -- matrix class

#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
#include <stdexcept>
#include <vector>
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

	struct matrix_size_type
	{
		long unsigned int rows;
		long unsigned int cols;
	};

	class MatrixException : public std::runtime_error
	{
	    public:
		MatrixException(const std::string &msg) : std::runtime_error(msg){}
	};
	
	template<class T>
	class matrix
	{
	public:
		matrix();
		matrix(const long unsigned int &rows, const long unsigned int &cols);
		matrix(const vector_k<vector_k<T>>& mat);
		matrix(const std::initializer_list<std::initializer_list<T>>& lst);
		T& operator()(const long unsigned int &i, const long unsigned int &j);
		T operator()(const long unsigned int &i, const long unsigned int &j) const;
		vector_k<T> operator*(const vector_k<T> &x) const;
		matrix<T> operator*(const matrix<T> &B) const;
		matrix<T> operator+(const matrix<T> &B) const;
		matrix<T>& operator+=(const matrix<T> &B);
		matrix<T> operator-(const matrix<T> &B) const;
		long unsigned int size(const long unsigned int &n) const;
		bool empty() const;
		vector_k<T> row(const long unsigned int &i);
		void setLastRow(const vector_k<T> &row);
		void addLastRow(const vector_k<T> &row);
		vector_k<T> getRow(const long unsigned int &i) const;
		vector_k<T> getLastRow() const;
		vector_k<T> getCol(const long unsigned int &j) const;
		vector_k<T> col(const long unsigned int &j);
		void reserve(const long unsigned int &N);
		vector_k<T> mData;

	private:
		long unsigned int mRows;
		long unsigned int mCols;
	};

	template<class T>
	matrix<T>::matrix() : mData(), mRows(0), mCols(0)
	{
	}

	template<class T>
	matrix<T>::matrix(const long unsigned int &rows, const long unsigned int &cols)
	: mData(rows * cols),
	  mRows(rows),
	  mCols(cols)
	{
		if(rows <= 0 || cols <= 0)
		{
			throw MatrixException("Tried to create matrix with invalid size. Number of rows and columns must be greater than zero!");
		}
	}

	template<class T>
	matrix<T>::matrix(const vector_k<vector_k<T> >& mat) : mData(mat[0].size()*mat.size()), mRows(mat.size()), mCols(mat[0].size())
	{
		if(mat.empty())
		{
			throw MatrixException("Cannot assign empty vector of vectors to a matrix object!");
		}
		else if(mat[0].empty())
		{
			throw MatrixException("Cannot assign empty vector of vectors to a matrix object!");
		}
		
		for(long unsigned int ii = 0; ii < mRows; ii++)
		{
			for(long unsigned int jj = 0; jj < mCols; jj++)
			{
				mData[jj*mRows + ii] = mat[ii][jj];
			}
		}
	}

	template<class T>
	matrix<T>::matrix(const std::initializer_list<std::initializer_list<T> >& lst) : matrix(lst.size(), lst.size() ? lst.begin()->size() : 0)
	{
		if(lst.size() <= 0 || lst.begin()->size() <= 0)
		{
			throw MatrixException("Cannot assign empty initializer list to a matrix object!");
		}
		long unsigned int ii = 0, jj = 0;
		for(const auto& l : lst)
		{
			for(const auto& v : l)
			{
				mData[jj*mRows + ii] = v;
				jj++;
			}
			jj = 0;
			ii++;
		}
	}

	template<class T>
	T& matrix<T>::operator()(const long unsigned int &i, const long unsigned int &j)
	{
		if(mData.empty())
		{
			throw MatrixException("Cannot access member of empty matrix!");
		}
		if(i > (mRows-1) || j > (mCols-1))
		{
			throw MatrixException("Tried to access invalid matrix member!");
		}
		return mData[j*mRows + i];
	}

	template<class T>
	T matrix<T>::operator()(const long unsigned int &i, const long unsigned int &j) const
	{
		if(mData.empty())
		{
			throw MatrixException("Cannot access member of empty matrix!");
		}
		if(i > (mRows-1) || j > (mCols-1))
		{
			throw MatrixException("Tried to access invalid matrix member!");
		}
		return mData[j*mRows + i];
	}

	template<class T>
	vector_k<T> matrix<T>::operator*(const vector_k<T> &x) const
	{
		if(x.empty())
		{
			throw MatrixException("Vector `x` cannot be empty in matrix-vector multiplication!");
		}
		if(x.size() != mCols)
		{
			throw MatrixException("Matrix and vector dimensions are not compatible in matrix-vector multiplication!");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		vector_k<T> b(mRows);
		for(long unsigned int ii = 0; ii < mRows; ii++)
		{
			b[ii] = 0.0;
			for(long unsigned int jj = 0; jj < mCols; jj++)
			{
				b[ii] += mData[jj*mRows + ii]*x[jj];
			}
		}
		return b;
	}
	
	template<>
	inline vector_k<double> matrix<double>::operator*(const vector_k<double> &x) const
	{
		if(x.empty())
		{
			throw MatrixException("Vector `x` cannot be empty in matrix-vector multiplication!");
		}
		if(x.size() != mCols)
		{
			throw MatrixException("Matrix and vector dimensions are not compatible in matrix-vector multiplication!");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		vector_k<double> b(mRows);
		int m = (int)mRows;
		int n = (int)mCols;
        char TRANS = 'N';
        double ALPHA = 1.0;
        int LDA = m;
        int INCX = 1;
        double BETA = 0.0;
        int INCY = 1;

        dgemv_(&TRANS, &m, &n, &ALPHA, &mData[0], &LDA, &x[0],&INCX, &BETA, &b[0], &INCY);
        
        return b;
    }
	
	template<>
	inline vector_k<std::complex<double>> matrix<std::complex<double>>::operator*(const vector_k<std::complex<double>> &x) const
	{
		if(x.empty())
		{
			throw MatrixException("Vector `x` cannot be empty in matrix-vector multiplication!");
		}
		if(x.size() != mCols)
		{
			throw MatrixException("Matrix and vector dimensions are not compatible in matrix-vector multiplication!");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		vector_k<std::complex<double>> b(mRows);
		int m = (int)mRows;
		int n = (int)mCols;
        char TRANS = 'N';
        std::complex<double> ALPHA = 1.0;
        int LDA = m;
        int INCX = 1;
        std::complex<double> BETA = 0.0;
        int INCY = 1;

        zgemv_(&TRANS, &m, &n, &ALPHA, &mData[0], &LDA, &x[0],&INCX, &BETA, &b[0], &INCY);
        
        return b;
    }

	template<class T>
	matrix<T> matrix<T>::operator*(const matrix<T> &B) const
	{
		if(B.size(1) <= 0 || B.size(2) <= 0)
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(mCols != B.size(1))
		{
			throw MatrixException("Matrix dimensions are not compatible in matrix-matrix multiplication!");
		}
		matrix<T> C(mRows,B.size(2));
		for(long unsigned int kk = 0; kk < B.size(2); kk++)
		{
			for(long unsigned int ii = 0; ii < mRows; ii++)
			{
				C(ii,kk) = 0.0;
				for(long unsigned int jj = 0; jj < mCols; jj++)
				{
					C(ii,kk) += mData[jj*mRows + ii]*B(jj,kk);
				}
			}
		}
		return C;
	}
	
	template<>
	inline matrix<double> matrix<double>::operator*(const matrix<double> &B) const
	{
		if(B.size(1) <= 0 || B.size(2) <= 0)
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(mCols != B.size(1))
		{
			throw MatrixException("Matrix dimensions are not compatible in matrix-matrix multiplication!");
		}
		matrix<double> C(mRows,B.size(2));
		int m = (int)mRows;
		int k = (int)mCols;
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
	inline matrix<std::complex<double>> matrix<std::complex<double>>::operator*(const matrix<std::complex<double>> &B) const
	{
		if(B.size(1) <= 0 || B.size(2) <= 0)
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(mCols != B.size(1))
		{
			throw MatrixException("Matrix dimensions are not compatible in matrix-matrix multiplication!");
		}
		matrix<std::complex<double>> C(mRows,B.size(2));
		int n = (int)mRows;
		int k = (int)mCols;
		int m = (int)B.size(2);
        char TRANS = 'N';
        std::complex<double> ALPHA = 1.0;
        int LDA = n;
        int LDB = k;
        std::complex<double> BETA = 0.0;
        int LDC = n;

        zgemm_(&TRANS, &TRANS, &m, &n, &k, &ALPHA, &mData[0], &LDA, &B.mData[0], &LDB, &BETA, &C.mData[0], &LDC);
        
        return C;
    }

	template<class T>
	matrix<T> matrix<T>::operator+(const matrix<T> &B) const
	{
		if(B.size(1) <= 0 || B.size(2) <= 0)
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(mRows != B.size(1) || mCols != B.size(2))
		{
			throw MatrixException("Matrix dimensions are not compatible in matrix-matrix addition!");
		}
		matrix<T> C(mRows,mCols);
		for(long unsigned int ii = 0; ii < mRows; ii++)
		{
			for(long unsigned int jj = 0; jj < mCols; jj++)
			{
				C(ii,jj) = mData[jj*mRows + ii] + B(ii,jj);
			}
		}
		return C;
	}
	
	template<class T>
	matrix<T>& matrix<T>::operator+=(const matrix<T> &B)
	{
		if(B.size(1) <= 0 || B.size(2) <= 0)
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(mRows != B.size(1) || mCols != B.size(2))
		{
			throw MatrixException("Matrix dimensions are not compatible in matrix-matrix addition!");
		}
		for(long unsigned int ii = 0; ii < mRows; ii++)
		{
			for(long unsigned int jj = 0; jj < mCols; jj++)
			{
				mData[jj*mRows + ii] += B(ii,jj);
			}
		}
		return *this;
	}

	template<class T>
	matrix<T> matrix<T>::operator-(const matrix<T> &B) const
	{
		if(B.size(1) <= 0 || B.size(2) <= 0)
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot perform operation on empty matrix!");
		}
		if(mRows != B.size(1) || mCols != B.size(2))
		{
			throw MatrixException("Matrix dimensions are not compatible in matrix-matrix subtraction!");
		}
		matrix<T> C(mRows,mCols);
		for(long unsigned int ii = 0; ii < mRows; ii++)
		{
			for(long unsigned int jj = 0; jj < mCols; jj++)
			{
				C(ii,jj) = mData[jj*mRows + ii] - B(ii,jj);
			}
		}
		return C;
	}

	template<class T>
	long unsigned int matrix<T>::size(const long unsigned int &n) const
	{
		if(n == 1)
		{
			return mRows;
		}
		else if(n == 2)
		{
			return mCols;
		}
		else
		{
			throw MatrixException("Invalid dimension in size().");
		}
	}

	template<class T>
	bool matrix<T>::empty() const
	{
		if(mRows == 0 || mCols == 0)
		{
			return true;
		}
		else
		{
		    return false;
		}
	}
	
	template<class T>
	vector_k<T> matrix<T>::row(const long unsigned int &n)
	{
		if(n > mRows)
		{
			throw MatrixException("Invalid row index in row().");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot use method row() on empty matrix!");
		}
		
		vector_k<T> Row(&mData[n],mCols,mRows);
		return Row;
	}

	template<class T>
	void matrix<T>::setLastRow(const vector_k<T> &p_row)
	{
		if(p_row.size() != mCols)
		{
			throw MatrixException("Vector dimension is incompatible with matrix dimension in setLastRow().");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot use method setLastRow() on empty matrix!");
		}

		for(long unsigned int jj = 0; jj < mCols; jj++)
		{
			mData[jj*mRows + (mRows-1)] = p_row[jj];
		}
	}

	template<class T>
	void matrix<T>::addLastRow(const vector_k<T> &p_row)
	{
		if(p_row.size() != mCols)
		{
			throw MatrixException("Vector dimension is incompatible with matrix dimension in addLastRow().");
		}

		mRows++;
		vector_k<T> temp = mData;
		mData.resize(mRows*mCols);
		for(long unsigned int jj = 0; jj < mCols; jj++)
		{
		    for(long unsigned int ii = 0; ii < (mRows-1); ii++)
		    {
		        mData[jj*mRows + ii] = temp[jj*(mRows-1) + ii];
		    }
		}
		for(long unsigned int jj = 0; jj < mCols; jj++)
		{
			mData[jj*mRows + (mRows-1)] = p_row[jj];
		}
	}

	template<class T>
	vector_k<T> matrix<T>::col(const long unsigned int &n)
	{
		if(n > mCols || n < 0)
		{
			throw MatrixException("Invalid column index in col().");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot use method col() on empty matrix!");
		}
		
		vector_k<T> Col(&mData[n*mRows],mRows,1);
		return Col;
	}

	template<class T>
	vector_k<T> matrix<T>::getRow(const long unsigned int &n) const
	{
		if(n > mRows || n < 0)
		{
			throw MatrixException("Invalid row index in getRow().");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot use method getRow() on empty matrix!");
		}
		vector_k<T> Row(mCols);

		for(long unsigned int jj = 0; jj < mCols; jj++)
		{
			Row[jj] = mData[jj*mRows + n];
		}

		return Row;
	}

	template<class T>
	vector_k<T> matrix<T>::getLastRow() const
	{
		if(mData.empty())
		{
			throw MatrixException("Cannot use method getLastRow() on empty matrix!");
		}
		vector_k<T> Row(mCols);

		for(long unsigned int jj = 0; jj < mCols; jj++)
		{
			Row[jj] = mData[jj*mRows + (mRows-1)];
		}

		return Row;
	}

	template<class T>
	vector_k<T> matrix<T>::getCol(const long unsigned int &n) const
	{
		if(n > mCols)
		{
			throw MatrixException("Invalid column index in getCol().");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot use method getCol() on empty matrix!");
		}
		vector_k<T> Col(mRows);

		for(long unsigned int ii = 0; ii < mRows; ii++)
		{
			Col[ii] = mData[n*mRows + ii];
		}

		return Col;
	}
	
	template<class T>
	void matrix<T>::reserve(const long unsigned int &N)
	{
		mData.reserve(N);
	}
}
#endif
