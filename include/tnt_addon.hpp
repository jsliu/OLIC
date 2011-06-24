/****************************************************
 *
 *  LU decomposition for band diagonal matrix
 *  The algorithm comes from Numerical Recipes (C version)
 *  I modified it to make it fit for C++
 *
 ****************************************************/
#ifndef _TNT_ADD_ON_H
#define _TNT_ADD_ON_H


#include "tnt/tnt_matrix.h"
#include "tnt/tnt_linalg.h"

using namespace TNT;
using namespace Linear_Algebra;

template <class Real>
class BANDLU
{
public:
		BANDLU(const Matrix<Real>&, const int, const int, const int, const int);
		~BANDLU();
		Matrix<Real> getL();
		Matrix<Real> getU();
		Vector<Real> solve(const Vector<Real>&);
		Matrix<Real> solve(const Matrix<Real>&);
		Real det();

private:
		Matrix<Real> L_;                  // lower triangular matrix
		Matrix<Real> U_;                  // upper triangular matrix
		Vector<int> piv;                  // record the permutation order
		int m, n, pivsign;
		int m_left, m_right;

		LU<Real> *lu;                     // lu decomposition for non banddiagonal matrix
		inline void SWAP(Real, Real);
};


template <class Real>
BANDLU<Real>::BANDLU(const Matrix<Real> &A, const int m1, const int m2, const int row, const int col)
		:piv(row),m(row),n(col),m_left(m1),m_right(m2)
/****************************************************************************************************
 *    
 *    calculate the LU decomposition of sub(A), where A is a banddiagonal matrix
 *    m1 is the subdiagonal rows and m2 in the superdiagonal rows;
 *    row and col are those corresponding to the submatrix of A
 *
 ***************************************************************************************************/
{
		const double TINY = 1.0e-20;
		assert(m==n);
		assert(m<=A.num_rows());

		int mm = m_left+m_right+1;
		if(mm>=n)
		{
				Matrix<double> Atmp(m,n);
				for(int i=1;i<=m;i++)
				{
						Atmp(i,i) = A(i,i);	
                        // e.g. if m==1, i-m_left will be 0, which is not allowed
						if(m>m_left)
						{				
								if(i>m_left && i<=m-m_right)
								{
										Atmp(i,i-m_left) = A(i,i-m_left);
										Atmp(i,i+m_right) = A(i,i+m_right);
								}
			
								if(i==1)
										Atmp(i,i+m_right) = A(i,i+m_right);
								if(i==m)
										Atmp(i,i-m_left) = A(i,i-m_left);
						}											   
				}
				lu = new LU<Real>(Atmp);
		}
		else 
		{
				Matrix<Real> Ltmp(m,m_left);
 
				// store the matrix compactly as a m x mm matrix 
				Matrix<Real> Utmp(m,mm);
				for(int i=1;i<=m;i++)
						Utmp(i,m_left+1) = A(i,i);

				// dist is the distance before the column [m_left+1]
				for(int dist=1;dist<=m_left;dist++)
						for(int i=dist+1;i<=m;i++)
								Utmp(i,m_left+1-dist) = A(i,i-dist);
				// dist is the distance after the column [m_left+1]

				for(int dist=1;dist<=m_right;dist++)
						for(int i=1;i<=m-dist;i++)
								Utmp(i,m_left+1+dist) = A(i+dist,i);

				int l = m_left;
				pivsign = 1;

				for(int i=1;i<m;i++)
						piv(i) = i;
  
				// Rearrange the storage a bit
				for(int i=1;i<=m_left;i++)
				{
						for(int j=m_left+2-i;j<=mm;j++)
								Utmp(i,j-l) = Utmp(i,j);
						l--;
						for(int j=mm-l;j<=mm;j++)
								Utmp(i,j) = 0.0;
				}

				l = m_left;
				// For each row
				for(int k=1;k<=m;k++)
				{
						Real dum = Utmp(k,1);
						int i = k;
						if(l<n)
								l++;

						// Find the pivot
						for(int j=k+1;j<=l;j++)
								if(abs(Utmp(j,1))>abs(dum))
								{
										dum = Utmp(j,1);
										i = j;
								}
		
						// Record the permutation order
						piv(k) = i;

						// Proceed with tiny pivot
						if(dum==0.0)
								Utmp(k,1) = TINY;
						if(i!=k)
						{
								pivsign = -pivsign;
								for(int j=1;j<=mm;j++)
										SWAP(Utmp(k,j),Utmp(i,j));
						}

						// Do the elimination
						for(i=k+1;i<=l;i++)
						{
								dum = Utmp(i,1)/Utmp(k,1);
								Ltmp(k,i-k) = dum;
								for(int j=2;j<=mm;j++)
										Utmp(i,j-1) = Utmp(i,j)-dum*Utmp(k,j);
								Utmp(i,mm) = 0.0;
						}
				}
				L_ = Ltmp;
				U_ = Utmp;
		}  // if else
}


template <class Real>
BANDLU<Real>:: ~BANDLU()
{
		int mm = m_left+m_right+1;
		if(mm>=n)
				delete lu;
}


template <class Real>
Matrix<Real> BANDLU<Real>::getL()
{
		int lrow = L_.num_rows();
		Matrix<Real> Lmat(lrow,lrow);

		for(int i=1;i<=lrow;i++)
				Lmat(i,i) = 1.0;

		for(int j=1;j<=m_left;j++)
				for(int i=1;i<=lrow-j;i++)
						Lmat(i+j,i) = L_(i,j);
		
		return Lmat;     
}


template <class Real>
Matrix<Real> BANDLU<Real>::getU()
{
		int mm = m_left+m_right+1;
		int urow = U_.num_rows();
		Matrix<Real> Umat(urow,urow);

		for(int j=1;j<=mm;j++)
				for(int i=1;i<=urow-j+1;i++)
						Umat(i,i+j-1) = U_(i,j);
   
		return Umat;    
}


template <class Real>
Vector<Real> BANDLU<Real>::solve(const Vector<Real> &B)
{
		assert(B.size()>=m);

		int mm = m_left+m_right+1;
		int l = m_left;
		Vector<Real> Xvec(m);
		for(int i=1;i<=m;i++)
				Xvec(i) = B(i);

		// Forward substitution, unscrambling the permuted rows
		for(int k=1;k<=n;k++)
		{
				int i = piv(k);
				if(i!=k)
						SWAP(Xvec(k),Xvec(i));
				if(l<n)
						l++;
				for(i=k+1;i<=l;i++)
						Xvec(i) -= L_(k,i-k)*Xvec(k);
		}

		l = 1;
		// Back substitution
		for(int i=n;i>=1;i--)
		{
				Real dum = Xvec(i);
				for(int k=2;k<=l;k++)
						dum -= U_(i,k)*Xvec(k+i-1);
				Xvec(i) = dum/U_(i,1);
				if(l<mm)
						l++;
		}
	
		return Xvec;					
}


template <class Real>
Matrix<Real> BANDLU<Real>::solve(const Matrix<Real> &B)
{
		assert(B.num_rows()>=m);

		Vector<Real> LUcolj(m);
		int nrow = n;
		int ncol = m;
		Matrix<Real> Xmat(nrow,ncol);

		int mm = m_left+m_right+1;
		if(mm>=n)
		{
				Matrix<Real> id(n,n);
				for(int i=1;i<=n;i++)
						id(i,i) = B(i,i);
				Xmat = (*lu).solve(id);
		}
		else 
		{
				for(int j=1;j<=n;j++)
				{
						for(int i=1;i<=m;i++)
								LUcolj(i) = B(i,j);
			 
						Vector<Real> vec_tmp = solve(LUcolj);
						for(int i=1;i<=m;i++)
								Xmat(i,j) = vec_tmp(i);
				}
		}
		return Xmat;
}


template <class Real>
Real BANDLU<Real>::det()
{
		Real d;
		int mm = m_left+m_right+1;
		if(mm>=n)
				d = (*lu).det();		
		else 
		{
				if (m != n) 
						return Real(0);
      
				d = Real(pivsign);
				for (int j = 1; j <= n; j++)
						d *= U_(j,1);
		}
		return d;
}


template <class Real>
inline void BANDLU<Real>::SWAP(Real a, Real b)
{
		Real tmp = a;
		a = b;
		b = tmp;
}


/***********************************************************************
 *
 *                      Inverse a tridiagonal matrix
 *
 *  To solve the equation AX = B
 * 
 **********************************************************************/
template <class Real>
class TRIDAG
{
public:
		TRIDAG(const Matrix<Real>&, int, int);
		TRIDAG(const Vector<Real>&, const Vector<Real>&, const Vector<Real>&, int, int); 
		~TRIDAG();
		Vector<Real> solve(const Vector<Real>&);
		Matrix<Real> solve(const Matrix<Real>&);
		Matrix<Real> inverse();                                   // matrix inversion

private:
		Vector<Real> *diag;
		Vector<Real> *subdiag;
		Vector<Real> *supdiag;
		int m,n;                                     // m==n the rows and columns of A
};

template <class Real>
TRIDAG<Real>::TRIDAG(const Matrix<Real> &A, int row, int col)
		:m(row),n(col)
{
		diag = new Vector<Real>(n);
		subdiag = new Vector<Real>(n-1);
		supdiag = new Vector<Real>(n-1);

		for(int i=1;i<n;i++)
		{
				(*diag)(i) = A(i,i);
				(*subdiag)(i) = A(i+1,i);
				(*supdiag)(i) = A(i,i+1);
		}
		(*diag)(n) = A(n,n);
}


template <class Real>
TRIDAG<Real>::TRIDAG(const Vector<Real> &Diag, const Vector<Real> &SubDiag, const Vector<Real> &SupDiag, int row, int col)
		:m(row),n(col)
{
		diag = new Vector<Real>(n);
		subdiag = new Vector<Real>(n-1);
		supdiag = new Vector<Real>(n-1);

		for(int i=1;i<n;i++)
		{
				(*diag)(i) = Diag(i);
				(*subdiag)(i) = SubDiag(i);
				(*supdiag)(i) = SupDiag(i);
		}
		(*diag)(n) = Diag(n);
}

 
template <class Real>
TRIDAG<Real>::~TRIDAG()
{
		delete diag;
		delete subdiag;
		delete supdiag;
}

template <class Real>
Vector<Real> TRIDAG<Real>::solve(const Vector<Real> &B)
{
		Vector<Real> Xvec(n);
		Vector<Real> gam(n);
		if((*diag)(1)== Real (0))
				cerr<<"Error 1 in TRIDAG.\n";

		double bet = (*diag)(1);
		Xvec(1) = B(1)/bet;

		// Decomposition and forward substitution
		for(int j=2;j<=n;j++)
		{
				gam(j) = (*supdiag)(j-1)/bet;
				bet = (*diag)(j)-(*subdiag)(j-1)*gam(j);
				if(bet==Real (0))
						cerr<<"Error 2 in TRIDAG.\n";
				Xvec(j) = (B(j)-(*subdiag)(j-1)*Xvec(j-1))/bet;
		}

		// Backward substitution
		for(int j=n-1;j>=1;j--)
				Xvec(j) -= gam(j+1)*Xvec(j+1);

		return Xvec;
}

template <class Real>
Matrix<Real> TRIDAG<Real>::solve(const Matrix<Real> &B)
{
		Matrix<Real> Xmat(n,m);
		for(int j=1;j<=n;j++)
		{
				Vector<Real> Bcolj(m);
				for(int i=1;i<=m;i++)
						Bcolj(i) = B(i,j);
			 
				Vector<Real> vec_tmp = solve(Bcolj);
				for(int i=1;i<=m;i++)
						Xmat(i,j) = vec_tmp(i);
		}

		return Xmat;
}

template<class Real>
Matrix<Real> TRIDAG<Real>::inverse()
{
		Matrix<Real> Xmat(n,m);
		for(int j=1;j<=n;j++)
		{
				Vector<Real> Bcolj(m);
				Bcolj(j) = Real (1);
			 
				Vector<Real> vec_tmp = solve(Bcolj);
				for(int i=1;i<=m;i++)
						Xmat(i,j) = vec_tmp(i);
		}

		return Xmat;
}


#endif  // _TNT_ADD_ON_H
