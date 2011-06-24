#include <iostream>
using namespace std;

// calculate Ax=b, where x,b are vectors, A is a tridiagonal matrix
double* solve_vec(double *diag, double *subdiag, double *supdiag, double *b, int size)
{
		double gam[size];
		double inv[size];
	 
		if(diag[0]==0.0)
				cerr<<"Error 1 in tridag\n";
		double bet;
		inv[0] = b[0]/(bet=(double) diag[0]);

		// Decomposition and forward substitution
		for(int i=1;i<size;i++)
		{
				gam[i] = supdiag[i-1]/bet;
				bet = diag[i]-subdiag[i-1]*gam[i];
				if(bet==0.0)
						cerr<<"Error 2 in tridag\n";
				inv[i] = (b[i]-subdiag[i-1]*inv[i-1])/bet;
		}

		// Bacward substitution
		for(int i=size-2;i>=0;i--)
				inv[i] -= gam[i+1]*inv[i+1];

		double *tmp = inv;
		return tmp;
}


extern "C" {
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

		// inverse the tridiagonal matrix
		SEXP inv_tridag(SEXP Diag, SEXP SubDiag, SEXP SupDiag)
		{
				int len = length(Diag);
				double diag[len],subdiag[len-1],supdiag[len-1];

				Diag = coerceVector(Diag,REALSXP);
				SubDiag = coerceVector(SubDiag,REALSXP);
				SupDiag = coerceVector(SupDiag,REALSXP);

				for(int i=0;i<len-1;i++)
				{
						diag[i] = REAL(Diag)[i];
						subdiag[i] = REAL(SubDiag)[i];
						supdiag[i] = REAL(SubDiag)[i];
				}
				diag[len-1] = REAL(Diag)[len-1];

				SEXP Rinv;
				PROTECT(Rinv = allocMatrix(REALSXP,len,len));

				for(int j=0;j<len;j++)
				{
						double Bvec[len];
						for(int i=0;i<len;i++)
								Bvec[i] = 0.0;
						Bvec[j] = 1.0;

						double *Mcolj = solve_vec(diag,subdiag,supdiag,Bvec,len);
						for(int i=0;i<len;i++)
								REAL(Rinv)[j+i*len] = Mcolj[i];
				}
				UNPROTECT(1);
				return Rinv;
		}

		// calculate t(A) x B
		SEXP transpose_mult(SEXP MatA, SEXP MatB)
		{
				MatA = coerceVector(MatA,REALSXP);
				MatB = coerceVector(MatB,REALSXP);

				SEXP RdimA = getAttrib(MatA, R_DimSymbol);
				SEXP RdimB = getAttrib(MatB, R_DimSymbol);

				int AI = INTEGER(RdimA)[0];
				int AJ = INTEGER(RdimA)[1];
				int BI = INTEGER(RdimB)[0];
				int BJ = INTEGER(RdimB)[1];

				if(AI!=BI)
						error("The two matrices don't match.");

				SEXP MatAB;
				PROTECT(MatAB = allocMatrix(REALSXP,AJ,BJ));
				for(int i=0;i<AJ;i++)
						for(int j=0;j<BJ;j++)
						{
								double sum = 0.0;
								for(int k=0;k<AI;k++)
										sum += REAL(MatA)[k+i*AI] * REAL(MatB)[k+j*BI];
								
								REAL(MatAB)[i+j*AJ] = sum;
						}
				UNPROTECT(1);
				return MatAB;
		}

}
