// matrix.h -- matrix class

#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
#include <stdexcept>
#include <vector>

namespace keycpp
{

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
		matrix(int rows, int cols);
		matrix(const std::vector<std::vector<T>>& mat);
		matrix(const std::initializer_list<std::initializer_list<T>>& lst);
		T& operator()(int i, int j);
		T operator()(int i, int j) const;
		std::vector<T> operator*(const std::vector<T> x) const;
		matrix<T> operator*(const matrix<T> B) const;
		matrix<T> operator+(const matrix<T> B) const;
		matrix<T> operator-(const matrix<T> B) const;
		int size(int n) const;
		bool empty() const;
		int setRow(std::vector<T> row, int i);
		int setLastRow(std::vector<T> row);
		int addLastRow(std::vector<T> row);
		int setCol(std::vector<T> col, int j);
		std::vector<T> getRow(int i);
		std::vector<T> getLastRow();
		std::vector<T> getCol(int j);
		int reserve(int N);

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
	matrix<T>::matrix(int rows, int cols)
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
	T& matrix<T>::operator()(int i, int j)
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
	T matrix<T>::operator()(int i, int j) const
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
	std::vector<T> matrix<T>::operator*(const std::vector<T> x) const
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

	template<class T>
	matrix<T> matrix<T>::operator*(const matrix<T> B) const
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


	template<class T>
	matrix<T> matrix<T>::operator+(const matrix<T> B) const
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
	matrix<T> matrix<T>::operator-(const matrix<T> B) const
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
	int matrix<T>::size(int n) const
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
	int matrix<T>::setRow(std::vector<T> row, int n)
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
	int matrix<T>::setLastRow(std::vector<T> row)
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
	int matrix<T>::addLastRow(std::vector<T> row)
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
	int matrix<T>::setCol(std::vector<T> col, int n)
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
	std::vector<T> matrix<T>::getRow(int n)
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
	std::vector<T> matrix<T>::getLastRow()
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
	std::vector<T> matrix<T>::getCol(int n)
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
	int matrix<T>::reserve(int N)
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
