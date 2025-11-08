#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

double** createMatrix(int m, int n);
void freeMat(double** mat, int m, int n);
double** transpose(double** mat, int m, int n);
void setCol(double** mat, double* vec, int n, int col);
void MatVec(double** mat, double* vec, double* res, int m, int n);
void addVec(double* vec1, double* vec2, double* res, int n);
void scaleVec(double* vec, double scale, int n);
void copyVec(double* dest, double* src, int n);
double innerProduct(double* vec1, double* vec2, int n);
double findNorm(double* vec, int n);
void normalize(double* vec, int n);
void reorthogonalize(double* x, double** Q, int len, int numCols);
double** matMul(double** a, double** b, int m, int k, int n);
double frobErr(double** mat1, double** mat2, int m, int n);
double frob(double** mat1, int m, int n);

#endif
