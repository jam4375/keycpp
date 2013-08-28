// matrix.h -- matrix class

#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
#include <stdexcept>
#include <vector>

namespace keycpp
{
    extern "C"{
	/** \brief This provides a C interface to BLAS's double matrix-matrix multiplication function. */
	void dgemm_(const char *TRANSA, const char *TRANSB, const int *M, const int *N,
	            const int *K, double *ALPHA, double *A, const int *LDA, double *B,
	            const int *LDB, double *BETA, double *C, const int *LDC);
	            
	/** \brief This provides a C interface to BLAS's complex double matrix-matrix multiplication function. */
	void zgemm_(const char *TRANSA, const char *TRANSB, const int *M, const int *N,
	            const int *K, std::complex<double> *ALPHA, std::complex<double> *A, const int *LDA, std::complex<double> *B,
	            const int *LDB, std::complex<double> *BETA, std::complex<double> *C, const int *LDC);
	            
	/** \brief This provides a C interface to BLAS's double matrix-vector multiplication function. */
	void dgemv_(const char * TRANS, const int *M, const int *N, const double *ALPHA, const double *A, const int *LDA, double *X,
	            const int *INCX, const double *BETA, double *Y, const int *INCY);
	            
	/** \brief This provides a C interface to BLAS's complex double matrix-vector multiplication function. */
	void zgemv_(const char * TRANS, const int *M, const int *N, const std::complex<double> *ALPHA,
	            const std::complex<double> *A, const int *LDA, std::complex<double> *X,
	            const int *INCX, const std::complex<double> *BETA, std::complex<double> *Y, const int *INCY);
	}

	struct matrix_size_type
	{
		int rows;
		int cols;
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
		matrix(const int &rows, const int &cols);
		matrix(const std::vector<std::vector<T>>& mat);
		matrix(const std::initializer_list<std::initializer_list<T>>& lst);
		T& operator()(const int &i, const int &j);
		T operator()(const int &i, const int &j) const;
		std::vector<T> operator*(const std::vector<T> &x) const;
		matrix<T> operator*(const matrix<T> &B) const;
		matrix<T> operator+(const matrix<T> &B) const;
		matrix<T>& operator+=(const matrix<T> &B);
		matrix<T> operator-(const matrix<T> &B) const;
		int size(const int &n) const;
		bool empty() const;
		int setRow(const std::vector<T> &row, const int &i);
		int setLastRow(const std::vector<T> &row);
		int addLastRow(const std::vector<T> &row);
		int setCol(const std::vector<T> &col, const int &j);
		std::vector<T> getRow(const int &i) const;
		std::vector<T> getLastRow() const;
		std::vector<T> getCol(const int &j) const;
		int reserve(const int &N);

	private:
		int mRows;
		int mCols;
		std::vector<T> mData;
	};

	template<class T>
	matrix<T>::matrix()
	{
	}

	template<class T>
	matrix<T>::matrix(const int &rows, const int &cols)
	: mRows(rows),
	  mCols(cols),
	  mData(rows * cols)
	{
		if(rows <= 0 || cols <= 0)
		{
			throw MatrixException("Tried to create matrix with invalid size. Number of rows and columns must be greater than zero!");
		}
	}

	template<class T>
	matrix<T>::matrix(const std::vector<std::vector<T> >& mat) : mRows(mat.size()), mCols(mat[0].size()), mData(mRows * mCols)
	{
		if(mat.empty())
		{
			throw MatrixException("Cannot assign empty vector of vectors to a matrix object!");
		}
		else if(mat[0].empty())
		{
			throw MatrixException("Cannot assign empty vector of vectors to a matrix object!");
		}
		for(int ii = 0; ii < mRows; ii++)
		{
			for(int jj = 0; jj < mCols; jj++)
			{
				mData[ii*mCols + jj] = mat[ii][jj];
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
		int ii = 0, jj = 0;
		for(const auto& l : lst)
		{
			for(const auto& v : l)
			{
				mData[ii*mCols + jj] = v;
				jj++;
			}
			jj = 0;
			ii++;
		}
	}

	template<class T>
	T& matrix<T>::operator()(const int &i, const int &j)
	{
		if(mData.empty())
		{
			throw MatrixException("Cannot access member of empty matrix!");
		}
		if(i < 0 || i > (mRows-1) || j < 0 || j > (mCols-1))
		{
			throw MatrixException("Tried to access invalid matrix member!");
		}
		return mData[i*mCols + j];
	}

	template<class T>
	T matrix<T>::operator()(const int &i, const int &j) const
	{
		if(mData.empty())
		{
			throw MatrixException("Cannot access member of empty matrix!");
		}
		if(i < 0 || i > (mRows-1) || j < 0 || j > (mCols-1))
		{
			throw MatrixException("Tried to access invalid matrix member!");
		}
		return mData[i*mCols + j];
	}

	template<class T>
	std::vector<T> matrix<T>::operator*(const std::vector<T> &x) const
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
		std::vector<T> b(mRows);
		for(int ii = 0; ii < mRows; ii++)
		{
			b[ii] = 0.0;
			for(int jj = 0; jj < mCols; jj++)
			{
				b[ii] += mData[ii * mCols + jj]*x[jj];
			}
		}
		return b;
	}
	
	template<>
	inline std::vector<double> matrix<double>::operator*(const std::vector<double> &x) const
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
		std::vector<double> b(mRows);
		int m = mRows;
		int n = mCols;
		double *A = new double[m*n];
		double *X = new double[n];
		for(int ii = 0; ii < n; ii++)
		{
			for(int jj = 0; jj < m; jj++)
			{
				A[ii*m + jj] = mData[jj*mCols + ii];
			}
		}
		for(int ii = 0; ii < n; ii++)
		{
			X[ii] = x[ii];
		}
        char TRANS = 'N';
        double ALPHA = 1.0;
        int LDA = m;
        int INCX = 1;
        double BETA = 0.0;
        double *y = new double[m];
        int INCY = 1;

        dgemv_(&TRANS, &m, &n, &ALPHA, A, &LDA, X,&INCX, &BETA, y, &INCY);
        
        for(int ii = 0; ii < m; ii++)
        {
            b[ii] = y[ii];
        }
        
        delete [] A;
        delete [] X;
        delete [] y;
        
        return b;
    }
	
	template<>
	inline std::vector<std::complex<double>> matrix<std::complex<double>>::operator*(const std::vector<std::complex<double>> &x) const
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
		std::vector<std::complex<double>> b(mRows);
		int m = mRows;
		int n = mCols;
		std::complex<double> *A = new std::complex<double>[m*n];
		std::complex<double> *X = new std::complex<double>[n];
		for(int ii = 0; ii < n; ii++)
		{
			for(int jj = 0; jj < m; jj++)
			{
				A[ii*m + jj] = mData[jj*mCols + ii];
			}
		}
		for(int ii = 0; ii < n; ii++)
		{
			X[ii] = x[ii];
		}
        char TRANS = 'N';
        std::complex<double> ALPHA = 1.0;
        int LDA = m;
        int INCX = 1;
        std::complex<double> BETA = 0.0;
        std::complex<double> *y = new std::complex<double>[m];
        int INCY = 1;

        zgemv_(&TRANS, &m, &n, &ALPHA, A, &LDA, X,&INCX, &BETA, y, &INCY);
        
        for(int ii = 0; ii < m; ii++)
        {
            b[ii] = y[ii];
        }
        
        delete [] A;
        delete [] X;
        delete [] y;
        
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
		for(int kk = 0; kk < B.size(2); kk++)
		{
			for(int ii = 0; ii < mRows; ii++)
			{
				C(ii,kk) = 0.0;
				for(int jj = 0; jj < mCols; jj++)
				{
					C(ii,kk) += mData[ii * mCols + jj]*B(jj,kk);
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
		int m = mRows;
		int k = mCols;
		int n = B.size(2);
		double *a = new double[m*k];
		double *b = new double[k*n];
		for(int ii = 0; ii < k; ii++)
		{
			for(int jj = 0; jj < m; jj++)
			{
				a[ii*m + jj] = mData[jj*mCols + ii];
			}
		}
		for(int ii = 0; ii < n; ii++)
		{
			for(int jj = 0; jj < k; jj++)
			{
				b[ii*k + jj] = B(jj,ii);
			}
		}
        char TRANS = 'N';
        double ALPHA = 1.0;
        int LDA = m;
        int LDB = k;
        double BETA = 0.0;
        double *c = new double[n*m];
        int LDC = m;

        dgemm_(&TRANS, &TRANS, &m, &n, &k, &ALPHA, a, &LDA, b, &LDB, &BETA, c, &LDC);
        
        for(int ii = 0; ii < n; ii++)
        {
            for(int jj = 0; jj < m; jj++)
            {
                C(jj,ii) = c[ii*m + jj];
            }
        }
        
        delete [] a;
        delete [] b;
        delete [] c;
        
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
		int n = mRows;
		int k = mCols;
		int m = B.size(2);
		std::complex<double> *a = new std::complex<double>[n*k];
		std::complex<double> *b = new std::complex<double>[k*m];
		for(int ii = 0; ii < k; ii++)
		{
			for(int jj = 0; jj < n; jj++)
			{
				a[ii*n + jj] = mData[jj*mCols + ii];
			}
		}
		for(int ii = 0; ii < m; ii++)
		{
			for(int jj = 0; jj < k; jj++)
			{
				b[ii*k + jj] = B(jj,ii);
			}
		}
        char TRANS = 'N';
        std::complex<double> ALPHA = 1.0;
        int LDA = n;
        int LDB = k;
        std::complex<double> BETA = 0.0;
        std::complex<double> *c = new std::complex<double>[n*m];
        int LDC = n;

        zgemm_(&TRANS, &TRANS, &m, &n, &k, &ALPHA, a, &LDA, b, &LDB, &BETA, c, &LDC);
        
        for(int ii = 0; ii < m; ii++)
        {
            for(int jj = 0; jj < n; jj++)
            {
                C(jj,ii) = c[ii*n + jj];
            }
        }
        
        delete [] a;
        delete [] b;
        delete [] c;
        
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
		for(int ii = 0; ii < mRows; ii++)
		{
			for(int jj = 0; jj < mCols; jj++)
			{
				C(ii,jj) = mData[ii * mCols + jj] + B(ii,jj);
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
		for(int ii = 0; ii < mRows; ii++)
		{
			for(int jj = 0; jj < mCols; jj++)
			{
				mData[ii * mCols + jj] += B(ii,jj);
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
		for(int ii = 0; ii < mRows; ii++)
		{
			for(int jj = 0; jj < mCols; jj++)
			{
				C(ii,jj) = mData[ii * mCols + jj] - B(ii,jj);
			}
		}
		return C;
	}

	template<class T>
	int matrix<T>::size(const int &n) const
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
	int matrix<T>::setRow(const std::vector<T> &row, const int &n)
	{
		if(row.size() != mCols)
		{
			throw MatrixException("Vector dimension is incompatible with matrix dimension in setRow().");
		}
		if(n > mRows || n < 0)
		{
			throw MatrixException("Invalid row index in setRow().");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot use method setRow() on empty matrix!");
		}

		for(int jj = 0; jj < mCols; jj++)
		{
			mData[n * mCols + jj] = row[jj];
		}
		return 0;
	}

	template<class T>
	int matrix<T>::setLastRow(const std::vector<T> &row)
	{
		if(row.size() != mCols)
		{
			throw MatrixException("Vector dimension is incompatible with matrix dimension in setLastRow().");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot use method setLastRow() on empty matrix!");
		}

		for(int jj = 0; jj < mCols; jj++)
		{
			mData[(mRows-1) * mCols + jj] = row[jj];
		}
		return 0;
	}

	template<class T>
	int matrix<T>::addLastRow(const std::vector<T> &row)
	{
		if(row.size() != mCols)
		{
			throw MatrixException("Vector dimension is incompatible with matrix dimension in addLastRow().");
		}

		mRows++;
		mData.resize(mRows*mCols);
		for(int jj = 0; jj < mCols; jj++)
		{
			mData[(mRows-1) * mCols + jj] = row[jj];
		}
		return 0;
	}

	template<class T>
	int matrix<T>::setCol(const std::vector<T> &col, const int &n)
	{
		if(col.size() != mRows)
		{
			throw MatrixException("Vector dimension is incompatible with matrix dimension in setCol().");
		}
		if(n > mCols || n < 0)
		{
			throw MatrixException("Invalid column index in setCol().");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot use method setCol() on empty matrix!");
		}
		for(int ii = 0; ii < mRows; ii++)
		{
			mData[ii * mCols + n] = col[ii];
		}
		return 0;
	}

	template<class T>
	std::vector<T> matrix<T>::getRow(const int &n) const
	{
		if(n > mRows || n < 0)
		{
			throw MatrixException("Invalid row index in getRow().");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot use method getRow() on empty matrix!");
		}
		std::vector<T> Row(mCols);

		for(int jj = 0; jj < mCols; jj++)
		{
			Row[jj] = mData[n * mCols + jj];
		}

		return Row;
	}

	template<class T>
	std::vector<T> matrix<T>::getLastRow() const
	{
		if(mData.empty())
		{
			throw MatrixException("Cannot use method getLastRow() on empty matrix!");
		}
		std::vector<T> Row(mCols);

		for(int jj = 0; jj < mCols; jj++)
		{
			Row[jj] = mData[(mRows-1) * mCols + jj];
		}

		return Row;
	}

	template<class T>
	std::vector<T> matrix<T>::getCol(const int &n) const
	{
		if(n > mCols || n < 0)
		{
			throw MatrixException("Invalid column index in getCol().");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot use method getCol() on empty matrix!");
		}
		std::vector<T> Col(mRows);

		for(int ii = 0; ii < mRows; ii++)
		{
			Col[ii] = mData[ii * mCols + n];
		}

		return Col;
	}

	template<class T>
	int matrix<T>::reserve(const int &N)
	{
		if(N < 0)
		{
			throw MatrixException("Invalid reserve size!");
		}

		mData.reserve(N);

		return 0;
	}
}
#endif
