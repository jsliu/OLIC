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


#ifndef TNT_SPARSE_VECTOR_H
#define TNT_SPARSE_VECTOR_H

#include <vector>
#include "tnt_vector.h"

namespace TNT
{
//namespace Linear_Algebra
//{


/**

		Sparse Vector.

	<p>
	Index values begin at 1.  Thus S(3) = 6.5 means that the third, not
	the fourth element, is set to 6.5.

	<p>
	S.value(i) is the ith nonzero in the list, and S.index(i) is the 1-based
	index value, i.e.  S.value(1) is the first nonzero of S.

	

*/
template <class T>
class Sparse_Vector
{
	/* internally, everything is zero based */

  private:

		// (value, index) pairs
		//
		std::vector<std::pair<T, Subscript>  > s_;    

		// diagonals are stored seprately
		//

  	int dim_;        			// dimension
		int num_nonzeros_; 		// number of nonzeros




public:
		typedef typename std::vector<std::pair<T, Subscript> >::const_iterator 
																															const_iterator;


		const_iterator begin() const { return s_.begin(); }
		const_iterator end() const { return s_.end(); }

		inline const T& value(Subscript i) const { return s_[i-1].first; }
		inline const T& value(const_iterator p) const {return  p->first; }

		inline Subscript index(Subscript i) const {return  s_[i-1].second; }
		inline Subscript index(const_iterator p) const {return  p->second; }
	  
    inline T dot_product(const Vector<T> x) const
    {
        T sum(0);

        for ( const_iterator p = s_.begin(); p < s_.end(); p++)
        {
            sum += value(p) *  x(index(p));
        }
       return sum;
   	}

	Sparse_Vector() : s_(), dim_(0), num_nonzeros_(0) {}
	Sparse_Vector(Subscript N): dim_(N), num_nonzeros_(0) {}

			

	Sparse_Vector(Subscript N, Subscript nz, const T* Val, const Subscript *I):
					dim_(N),
					num_nonzeros_(0)
					{
						insert(nz, Val, I);	
					}
   


	void   insert(const T& val, Subscript i)
	{

			s_.push_back( std::make_pair(val,i) );
			num_nonzeros_++;
	}

	void   insert(Subscript nz, const T* Val, const Subscript *I)
	{

			for (int count=0; count<nz; count++)
			{
				insert(Val[count], I[count]);	
			}
	}


	

  inline   int    dim() const 	{return dim_;}
  int          num_nonzeros() const {return num_nonzeros_;}




	inline double norm() const
	{
		T sum(0.0);

		for (const_iterator p = s_.begin(); p < s_.end(); p++)
		{
			sum += value(p) * value(p); 
		}

		return sqrt(sum);
	}

			
};

template <class T>
inline const T& value(typename Sparse_Vector<T>::const_iterator p) 
{ 
	return p->first; 
}

template <class T>
inline const T& index(typename Sparse_Vector<T>::const_iterator p) 
{ 
	return p->second; 
}


template <class T>
inline T dot_product(const Sparse_Vector<T> &s, const Vector<T> &x)
{
	return s.dot_product(x);
}

template <class T>
inline T dot_product(const Vector<T> &x, const Sparse_Vector<T> &s)
{
	return s.dot_product(x);
}

template <class T>
inline double norm( const Sparse_Vector<T> & s)
{
	return s.norm();
}

//}  /* namspace TNT::Linear_Algebra */
}	 /* namspace TNT */

#endif
