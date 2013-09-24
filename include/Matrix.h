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
		int rows;
		int cols;
	};

	class MatrixException : public std::runtime_error
	{
	    public:
		MatrixException(const std::string &msg) : std::runtime_error(msg){}
	};
	
	template<class T>
	class vector_ref
	{
	public:
	    vector_ref(T *pData, int pN, int pinc) {Data = pData; N = pN; inc = pinc;};
	    T *Data;
	    int N;
	    int inc;
		vector_ref<T>& operator=(const std::vector<T> &v1)
		{
		    if(v1.size() != N)
		    {
		        throw MatrixException("Vectors must be of same length for assignment with vector_ref!");
		    }
		    for(int ii = 0; ii < N; ii++)
		    {
		        Data[ii*inc] = v1[ii];
		    }
		    return *this;
		};
		operator std::vector<T>()
		{
		    std::vector<T> v1(N);
		    for(int ii = 0; ii < N; ii++)
		    {
		        v1[ii] = Data[ii*inc];
		    }
		    return v1;
		};
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
		vector_ref<T> row(const int &i);
		int setLastRow(const std::vector<T> &row);
		int addLastRow(const std::vector<T> &row);
		std::vector<T> getRow(const int &i) const;
		std::vector<T> getLastRow() const;
		std::vector<T> getCol(const int &j) const;
		vector_ref<T> col(const int &j);
		int reserve(const int &N);
		std::vector<T> mData;

	private:
		int mRows;
		int mCols;
	};

	template<class T>
	matrix<T>::matrix() : mRows(0), mCols(0)
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
		int ii = 0, jj = 0;
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
		return mData[j*mRows + i];
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
		return mData[j*mRows + i];
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
				b[ii] += mData[jj*mRows + ii]*x[jj];
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
        char TRANS = 'N';
        double ALPHA = 1.0;
        int LDA = m;
        int INCX = 1;
        double BETA = 0.0;
        double *y = new double[m];
        int INCY = 1;

        dgemv_(&TRANS, &m, &n, &ALPHA, &mData[0], &LDA, &x[0],&INCX, &BETA, &b[0], &INCY);
        
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
        char TRANS = 'N';
        std::complex<double> ALPHA = 1.0;
        int LDA = m;
        int INCX = 1;
        std::complex<double> BETA = 0.0;
        std::complex<double> *y = new std::complex<double>[m];
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
		for(int kk = 0; kk < B.size(2); kk++)
		{
			for(int ii = 0; ii < mRows; ii++)
			{
				C(ii,kk) = 0.0;
				for(int jj = 0; jj < mCols; jj++)
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
		int m = mRows;
		int k = mCols;
		int n = B.size(2);
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
		int n = mRows;
		int k = mCols;
		int m = B.size(2);
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
		for(int ii = 0; ii < mRows; ii++)
		{
			for(int jj = 0; jj < mCols; jj++)
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
		for(int ii = 0; ii < mRows; ii++)
		{
			for(int jj = 0; jj < mCols; jj++)
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
		for(int ii = 0; ii < mRows; ii++)
		{
			for(int jj = 0; jj < mCols; jj++)
			{
				C(ii,jj) = mData[jj*mRows + ii] - B(ii,jj);
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
	vector_ref<T> matrix<T>::row(const int &n)
	{
		if(n > mRows || n < 0)
		{
			throw MatrixException("Invalid row index in setRow().");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot use method setRow() on empty matrix!");
		}
		
		vector_ref<T> Row(&mData[n],mCols,mRows);
		return Row;
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
			mData[jj*mRows + (mRows-1)] = row[jj];
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
		std::vector<T> temp = mData;
		mData.resize(mRows*mCols);
		for(int jj = 0; jj < mCols; jj++)
		{
		    for(int ii = 0; ii < (mRows-1); ii++)
		    {
		        mData[jj*mRows + ii] = temp[jj*(mRows-1) + ii];
		    }
		}
		for(int jj = 0; jj < mCols; jj++)
		{
			mData[jj*mRows + (mRows-1)] = row[jj];
		}
		return 0;
	}

	template<class T>
	vector_ref<T> matrix<T>::col(const int &n)
	{
		if(n > mCols || n < 0)
		{
			throw MatrixException("Invalid column index in setCol().");
		}
		if(mData.empty())
		{
			throw MatrixException("Cannot use method setCol() on empty matrix!");
		}
		
		vector_ref<T> Col(&mData[n*mRows],mRows,1);
		return Col;
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
			Row[jj] = mData[jj*mRows + n];
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
			Row[jj] = mData[jj*mRows + (mRows-1)];
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
			Col[ii] = mData[n*mRows + ii];
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
