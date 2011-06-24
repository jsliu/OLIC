/*
*
* Template Numerical Toolkit (TNT)
*
* Mathematical and Computational Sciences Division
* National Institute of Technology,
* Gaithersburg, MD USA
*
*
* This software was developed at the National Institute of Standards and
* Technology (NIST) by employees of the Federal Government in the course
* of their official duties. Pursuant to title 17 Section 105 of the
* United States Code, this software is not subject to copyright protection
* and is in the public domain. NIST assumes no responsibility whatsoever for
* its use by other parties, and makes no guarantees, expressed or implied,
* about its quality, reliability, or any other characteristic.
*
*/


#ifndef TNT_SPARSE_MATRIX_H
#define TNT_SPARSE_MATRIX_H

#include <vector>
#include "tnt_vector.h"
#include "tnt_sparse_vector.h"

namespace TNT
{

//namespace Linear_Algebra
//{


/**
	Read-only view of a sparse matrix in compressed-row storage
	format.  Neither array elements (nonzeros) nor sparsity
	structure can be modified.  If modifications are required,
	create a new view.

	<p>
	Index values begin at 1.

	<p>
	<b>Storage requirements:</b> An (m x n) matrix with
	nz nonzeros requires no more than  ((T+2I)*nz)
	bytes, where T is the size of data elements and
	I is the size of integer subscripts.
	

*/
template <class T>
class Sparse_Matrix
{
  private:

		// compressed row storage M rows of <value, col_index> pairs
		//
		std::vector< Sparse_Vector<T> > S_;    

		// diagonals are stored seprately
		//
		std::vector<T> diag_;

  	int num_rows_;        // number of rows
  	int num_cols_;        // number of cols
		int num_nonzeros_; 		// number of nonzeros

													// Used only in multi-step constructions.  This
													// allows one to build a sparse matrix with 
													// multiple calls to insert().
													//
		int internal_state_;	// 0 if closed (no more inserts) , 1 if open;



public:


	Sparse_Matrix(Subscript M, Subscript N):
								S_(M), diag_(M < N ? M : N, T(0)), 
								num_rows_(M), num_cols_(N), num_nonzeros_(0),
								internal_state_(1){};

			

	Sparse_Matrix(Subscript M, Subscript N, Subscript nz, const T* val,
				const Subscript *r, const Subscript *c):
					S_(M),
					diag_( M < N ? M : N, T(0) ),
					num_rows_(M), 
					num_cols_(N), 
					num_nonzeros_(0),
					internal_state_(1)
					{
						insert(nz, val, r, c);	
						close();
					};
   

	int is_closed() { return internal_state_; }

	void   insert(const T& val, Subscript i, Subscript j)
	{
		  if (internal_state_ == 0) return;

			i--;
			j--;

			if (i==j)
				diag_[i] = val;
			else
				S_[i].insert(val, j);

			num_nonzeros_++;
	}

	void   insert(Subscript nz, const T*  val, const Subscript *i, 
					const Subscript *j)
	{
		  if (internal_state_ == 0) return;

			for (int count=0; count<nz; count++)
			{
				insert(val[count], i[count], j[count]);	
			}
	}


	
	void 		close()
	{
			internal_state_ = 0;
	}

  inline   int    num_rows() const {return num_rows_;}
  inline   int    num_cols() const {return num_cols_;}
  int          num_nonzeros() const {return num_nonzeros_;}


	Vector<T> diag() const
	{
		return Vector<T>(diag_.size(),  &diag_[0]);
	}

	Vector<T> mult(const Vector<T> &x) const
	{
		int M = num_rows();
		Vector<T> y(M);
		for (int i=0; i<M; i++)
		{
			y[i] = dot_product(S_[i], x);
			y[i] += diag_[i] * x[i];
		}
		return y;
	}


	inline double norm() const
	{
		T sum(0.0);
		for (int i=0; i<num_rows_; i++)
		{
			// sum += S_[i].norm()* S_[i].norm();
			for (typename Sparse_Vector<T>::const_iterator p = S_[i].begin();
						p < S_[i].end(); p++ )
			{
				// sum += S_[i].norm();
				sum += S_[i].value(p) * S_[i].value(p);
			}
		}
		return sqrt(sum);
	}

			
};


template <class T>
inline Vector<T> operator*(const Sparse_Matrix<T> &S, const Vector<T> &x) 
{
	return S.mult(x);
}

template <class T>
inline double norm(const Sparse_Matrix<T> &S)
{
	return S.norm();
}


//}  /* namspace TNT::Linear_Algebra */
}	 /* namspace TNT */

#endif
