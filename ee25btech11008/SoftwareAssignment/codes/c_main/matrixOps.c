#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../c_libs/matrixOps.h"

double** createMatrix(int m, int n) {
    double** mat = (double**) malloc(m*sizeof(double*));
    for (int i=0; i<m; i++) {
        double* row = (double*) calloc(n, sizeof(double));
        mat[i] = row;
    }
    return mat;
}

void setCol(double** mat, double* vec, int n, int col) {
    for (int i=0; i<n; i++) {
        mat[i][col] = vec[i];
    }
}

double** transpose(double** mat, int m, int n) {
    double** matT = createMatrix(n, m);
    for (int i=0; i<m; i++) {
        setCol(matT, mat[i], n, i);
    }
    return matT;
}

void addVec(double* vec1, double* vec2, double* res, int n) {
    for (int i=0; i<n; i++) {
        res[i] = vec1[i] + vec2[i];
    }
}

void scaleVec(double* vec, double scale, int n) {
    for (int i=0; i<n; i++) {
        vec[i] *= scale;
    }
}

double innerProduct(double* vec1, double* vec2, int n) {
    double ans = 0;
    for (int i=0; i<n; i++) {
        ans += vec1[i]*vec2[i];
    }
    return ans;
}

void MatVec(double** mat, double* vec, double* res, int m, int n) {
    for (int i=0; i<m; i++) {
        res[i] = innerProduct(mat[i], vec, n);
    }
}

void copyVec(double* dest, double* src, int n) {
    for (int i=0; i<n; i++) {
        dest[i] = src[i];
    }
}

double findNorm(double* vec, int n) {
    double squareSum = 0;
    for (int i=0; i<n; i++) {
        squareSum += vec[i]*vec[i];
    }
    return sqrt(squareSum);
}

// Add norm = 0 case later
void normalize(double* vec, int n) {
    double norm = findNorm(vec, n);
    if (norm<1e-12) return;
    for (int i=0; i<n; i++) {
        vec[i] /= norm;
    }
}

void freeMat(double** mat, int m, int n) {
    for (int i=0; i<m; i++) {
        free(mat[i]);
    }
    free(mat);
}

void reorthogonalize(double* x, double** Q, int len, int numCols) {
    for (int j = 0; j < numCols; j++) {
        double dot = 0.0;
        for (int i = 0; i < len; i++) {
            dot += x[i] * Q[i][j];
        }
        for (int i = 0; i < len; i++) {
            x[i] -= dot * Q[i][j];
        }
    }
}

double** matMul(double** a, double** b, int m, int k, int n) {
    double** res = (double**)malloc(m * sizeof(double*));
    for (int i = 0; i < m; i++) {
        res[i] = (double*)calloc(n, sizeof(double));
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int x = 0; x < k; x++) {
                res[i][j] += a[i][x] * b[x][j];
            }
        }
    }

    return res;
}

double frobErr(double** mat1, double** mat2, int m, int n) {
    double squareSum = 0;
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            squareSum += (mat1[i][j]-mat2[i][j])*(mat1[i][j]-mat2[i][j]);
        }
    }
    return sqrt(squareSum);
}

double frob(double** mat1, int m, int n) {
    double squareSum = 0;
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
            squareSum += (mat1[i][j])*(mat1[i][j]);
        }
    }
    return sqrt(squareSum);
}