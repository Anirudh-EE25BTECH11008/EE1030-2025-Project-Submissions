#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../c_libs/matrixOps.h"
#include "../c_libs/svd.h"

double wilkinsonShift(double** bk, int k) {
    if (k > 2) {
        double a = bk[k-2][k-2]*bk[k-2][k-2] + bk[k-3][k-2]*bk[k-3][k-2];
        double b = bk[k-2][k-1]*bk[k-2][k-2];
        double c = bk[k-1][k-1]*bk[k-1][k-1] + bk[k-2][k-1]*bk[k-2][k-1];
        double d = (a-c)/2.0;
        double p = sqrt(d*d + b*b);
        if (d < 0) {
            p *= -1;
        }
        double shift = c - b*b/(d + p);
        return shift;
    } 
    return 0;
}

void rightGivens(double** matrix, int n, int col1, int col2, double c, double s) {
    for (int i=0; i<n; i++) {
        double temp1 = matrix[i][col1];
        double temp2 = matrix[i][col2];
        matrix[i][col1] = c*temp1 + s*temp2;
        matrix[i][col2] = -s*temp1 + c*temp2;
    }
}

void leftGivens(double** matrix, int n, int row1, int row2, double c, double s) {
    for (int j=0; j<n; j++) {
        double temp1 = matrix[row1][j];
        double temp2 = matrix[row2][j];
        matrix[row1][j] = c * temp1 + s*temp2;
        matrix[row2][j] = -s*temp1 + c*temp2;
    }
}

void svd(double** uk, double** bk, double** vk, int m, int n, int k) {
    double eps = 1e-12;
    int max_iter = 1000;
    
    while(max_iter--) {
        int p = 0;
        int q = 0;
        for (int i=0; i<k-1; i++) {
            if (fabs(bk[i][i+1]) < eps*(fabs(bk[i][i]) + fabs(bk[i+1][i+1]))) {
                bk[i][i+1] = 0;
                q = i+2;
            }
        }
        for (int i=q-2; i>0; i--) {
            if (fabs(bk[i-1][i]) < eps*(fabs(bk[i-1][i-1]) + fabs(bk[i][i]))) {
                p = i;
                break;
            }
        }
        if (q==0) {
            break;
        }
        
        double shift = wilkinsonShift(bk, q);
        double e1 = bk[p][p]*bk[p][p] - shift;
        double e2 = bk[p][p]*bk[p][p+1];
        
        for (int i=p; i<q-1; i++) {
            double r = sqrt(e1*e1 + e2*e2);
            if (r<eps) {
                break;
            }
            double c = e1/r;
            double s = e2/r;
            rightGivens(bk, k, i, i+1, c, s);
            rightGivens(vk, n, i, i+1, c, s);
            
            e1 = bk[i][i];
            e2 = bk[i+1][i];
            r = sqrt(e1*e1 + e2*e2);
            if (r<eps) {
                break;
            }
            c = e1/r;
            s = e2/r;
            leftGivens(bk, k, i, i+1, c, s);
            rightGivens(uk, m, i, i+1, c, s);
            
            if (i<q-2) {
                e1 = bk[i][i+1];
                e2 = bk[i][i+2];
            }
        }
    }
}