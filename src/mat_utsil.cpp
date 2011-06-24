#include "mat_utsil.hpp"

Matrix<double> sub_mult(const Matrix<double> &A, const Matrix<double> &B, int left, int right)
/*****************************************************************************
 *
 *   Calculate a dot product of A x sub(B)
 *   A:M x N matrix
 *   B:N x K matrix, will subtract some columns to match A'columns
 *
 ****************************************************************************/
{
		int M = A.num_rows();
		int N = A.num_cols();
		int K = B.num_cols();

		assert(K>left+right);
		Matrix<double> tmp(M,K-left-right);

		for(int i=0;i<M;i++)
		{
				for(int k=left;k<K-right;k++)
				{
						double sum = 0.0;
						for (int j=0;j<N;j++)
								sum = sum + A[i][j] * B[j][k];
	 	
						tmp[i][k-left] = sum;
				}
		}

	return tmp;
}


Matrix<double> sub_mult2(const Matrix<double> &A, const Matrix<double> &B, int left, int right)
/*****************************************************************************
 *
 *   Calculate a dot product of sub(A) x B
 *   A:M x N matrix, will subtract some columns to match B'columns
 *   B:P x K matrix
 *
 ****************************************************************************/
{
		int M = A.num_rows();
		int N = A.num_cols();
		int P = B.num_rows();
		int K = B.num_cols();

		assert(N-P==left+right);
		Matrix<double> tmp(M,K);

		for(int i=0;i<M;i++)
		{
				for(int k=0;k<K;k++)
				{
						double sum = 0.0;
						for (int j=left;j<N-right;j++)
								sum = sum + A[i][j-left] * B[j-left][k];
	 	
						tmp[i][k] = sum;
				}
		}

		return tmp;
}


Matrix<double> transpose_submult(const Matrix<double> &A, const Matrix<double> &B, int left, int right)
/************************************************************************
 *
 *  Calculate a dot product of t(sub(A))*B
 *  A: M x N matrix, will subtract some rows to match B'columns
 *  B: M x K matrix
 *
 ***********************************************************************/
{
		int M = A.num_cols();
		int N = A.num_rows();
		int K = B.num_cols();

		assert(M>left+right);
		Matrix<double> tmp(M-left-right,K);
		double sum;

        for(int i=left; i<M-right; i++) 
				for(int k=0; k<K; k++)
				{
						sum = 0;
						for (int j=0; j<N; j++)  
								sum = sum +  A[j][i] * B[j][k];

						tmp[i-left][k] = sum;
				}

		return tmp;
}


Matrix<double> transpose_submult2(const Matrix<double> &A, const Matrix<double> &B)
/*******************************************************************************************
 *
 *  Calculate a dot product of t(A)*sub(B)
 *  A: M x N matrix
 *  B: P x Q matrix
 *  subtract rows and columns from the top and right in B to match the dimension of A 
 *
 *******************************************************************************************/
{
		int M = A.num_cols();
		int N = A.num_rows();
		int P = B.num_rows();
		int Q = B.num_cols();
		assert(P>N && Q>M);

		Matrix<double> tmp(M,N);
		double sum;
 
        for(int i=0; i<M; i++) 
				for(int k=0; k<N; k++)
				{
						sum = 0;
						for (int j=0; j<N; j++)  
								sum = sum +  A[j][i] * B[j+Q-N][k+Q-N];
						tmp[i][k] = sum;
				}

		return tmp;
}


Matrix<double> transpose_submult3(const Matrix<double> &A, const Matrix<double> &B, int left, int right)
/************************************************************************
 *
 *  Calculate a dot product of t(sub(A))*sub(B)
 *  A: M x N matrix
 *  B: M x N matrix
 *
 ***********************************************************************/
{
		int M = A.num_cols();
		int N = A.num_rows();

		assert(M>left+right);
		Matrix<double> tmp(M-left-right,M-left-right);
		double sum;

        for(int i=left; i<M-right; i++) 
				for(int k=left; k<M-right; k++)
				{
						sum = 0;
						for (int j=0; j<N; j++)  
								sum = sum +  A[j][i] * B[j][k];

						tmp[i-left][k-left] = sum;
				}

		return tmp;
}


Matrix<double> transpose_submult4(const Matrix<double> &A, const Matrix<double> &B, int left, int right)
/************************************************************************
 *
 *  Calculate a dot product of A*t(sub(B))
 *  A: M x N matrix
 *  B: P x K matrix, will subtract some rows to match A's columns
 *
 ************************************************************************/
{
		int M = A.num_rows();
		int N = A.num_cols();
		int P = B.num_cols();
		int K = B.num_rows();

		assert(P-N==left+right);
		Matrix<double> tmp(M,K);
		double sum;

        for(int i=0; i<M; i++) 
				for(int k=0; k<K; k++)
				{
						sum = 0;
						for (int j=left; j<P-right; j++)  
								sum = sum +  A[i][j-left] * B[k][j-left];

						tmp[i][k] = sum;
				}

		return tmp;
}


Matrix<double> transpose_mult2(const Matrix<double> &A, const Matrix<double> &B)
/************************************************************************
 *
 *  Calculate a dot product of A*t(B)
 *  A: M x N matrix
 *  B: K x N matrix
 *
 ************************************************************************/
{
    assert(A.num_cols() == B.num_cols());
    int M = A.num_rows();
    int N = A.num_cols();
    int K = B.num_rows();

    Matrix<double> tmp(M,K);
    double sum;

    for (int i=0; i<M; i++) 
			for (int k=0; k<K; k++)
			{
					sum = 0;
					for (int j=0; j<N; j++) 
							sum = sum +  A[i][j] * B[k][j];

					tmp[i][k] = sum;
			}
   
    return tmp;
}

Matrix<double> sub_minus(const Matrix<double> &A, const Matrix<double> &B)
/************************************************************************
 *
 *  Calculate sub(A)-B
 *  A: M x N matrix, 
 *  will subtract the rows and columns from the botton row and right column
 *  B: P x Q matrix
 *  sub(A): P x Q matrix
 *
 ***********************************************************************/
{
		int P = B.num_rows();
		int Q = B.num_cols();

		Matrix<double> tmp(P,Q);
        for (int i=0; i<P; i++)
				for (int j=0; j<Q; j++)
						tmp[i][j] = A[i][j] - B[i][j];
 
		return tmp;

}


Matrix<double> sub_minus2(const Matrix<double> &A, const Matrix<double> &B)
/************************************************************************
 *
 *  Calculate sub(A)-B
 *  A: M x N matrix, 
 *  will subtract the rows and columns from the top row and left column
 *  B: P x Q matrix
 *  sub(A): P x Q matrix
 *
 ***********************************************************************/
{
		int M = A.num_rows();
		int N = A.num_cols();
		int P = B.num_rows();
		int Q = B.num_cols();

	    assert(M>P && N>Q);
		int row_dist = M-P;
		int col_dist = N-Q;

		Matrix<double> tmp(P,Q);
        for (int i=0; i<P; i++)
				for (int j=0; j<Q; j++)
						tmp[i][j] = A[i+row_dist][j+col_dist] - B[i][j];

		return tmp;

}

