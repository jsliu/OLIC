#ifndef _MAT_UTSIL_H
#define _MAT_UTSIL_H

#include "tnt/tnt_matrix.h"
using namespace TNT;

Matrix<double> sub_mult(const Matrix<double>&, const Matrix<double>&, int , int);
Matrix<double> sub_mult2(const Matrix<double>&, const Matrix<double>&, int, int);
Matrix<double> transpose_submult(const Matrix<double>&, const Matrix<double>&, int ,int);
Matrix<double> transpose_submult2(const Matrix<double>&, const Matrix<double>&);
Matrix<double> transpose_submult3(const Matrix<double>&, const Matrix<double>&, int, int);
Matrix<double> transpose_submult4(const Matrix<double>&, const Matrix<double>&, int, int);
Matrix<double> transpose_mult2(const Matrix<double> &A, const Matrix<double> &B);
Matrix<double> sub_minus(const Matrix<double>&, const Matrix<double>&);
Matrix<double> sub_minus2(const Matrix<double>&, const Matrix<double> &);


#endif
