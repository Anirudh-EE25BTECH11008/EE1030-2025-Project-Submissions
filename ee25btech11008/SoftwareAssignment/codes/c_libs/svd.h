#ifndef SVD_H
#define SVD_H

#include "matrixOps.h"

double wilkinsonShift(double** bk, int k);
void rightGivens(double** matrix, int rows, int col1, int col2, double c, double s);
void leftGivens(double** matrix, int cols, int row1, int row2, double c, double s);
void svd(double** ub, double** bk, double** vb, int m, int n, int k);

#endif
