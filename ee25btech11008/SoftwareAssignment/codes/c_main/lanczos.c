#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../c_libs/matrixOps.h"
#include "../c_libs/lanczos.h"

// Add norm = 0 case later
void bidiagonalize(double** A, double** uk, double** vk, double** bk, int m, int n, int k) {
    // Initialising
    srand(42);
    double** At = transpose(A, m, n);
    double* v = calloc(n, sizeof(double));
    double* u = calloc(m, sizeof(double));
    for (int i=0; i<n; i++) {
        v[i] = (double)rand()/(RAND_MAX);
    }
    normalize(v, n);
    setCol(vk, v, n, 0);
    double* temp1 = (double*) calloc(m, sizeof(double));
    double* temp2 = (double*) calloc(n, sizeof(double));
    double* Av = (double*) calloc(m, sizeof(double));
    double* AtU = (double*) calloc(n, sizeof(double));

    // A * v`i = alpha`i * u`i + beta`i-1 * u`i-1
    // A^T * u`i = alpha`i * v`i + beta`i * v`i+1
    int currU = 0;
    int currV = 1;
    double alpha = 0;
    double beta = 0; 
    MatVec(A, v, Av, m, n);
    alpha = findNorm(Av, m);
    bk[currV-1][currU] = alpha;
    scaleVec(Av, 1/alpha, m);
    copyVec(u, Av, m);
    reorthogonalize(u, uk, m, currU);
    normalize(u, m);
    setCol(uk, u, m, currU);
    currU++;
    MatVec(At, u, AtU, n, m);
    scaleVec(v, -1*alpha, n);
    addVec(AtU, v, temp2, n);
    beta = findNorm(temp2, n);
    bk[currV-1][currU] = beta;
    scaleVec(temp2, 1/beta, n);
    copyVec(v, temp2, n);
    reorthogonalize(v, vk, n, currV);
    normalize(v, n);
    setCol(vk, v, n, currV);
    currV++;


    // The loop
    while(currU!=k) {
        MatVec(A, v, Av, m, n);
        scaleVec(u, -1*beta, m);
        addVec(Av, u, temp1, m);
        alpha = findNorm(temp1, m);
        if (alpha<1e-12) break;
        bk[currV-1][currU] = alpha;
        scaleVec(temp1, 1/alpha, m);
        copyVec(u, temp1, m);
        reorthogonalize(u, uk, m, currU);
        normalize(u, m);
        setCol(uk, u, m, currU);
        currU++;
        if (currV<k) {
            MatVec(At, u, AtU, n, m);
            scaleVec(v, -1*alpha, n);
            addVec(AtU, v, temp2, n);
            beta = findNorm(temp2, n);
            if (beta<1e-12) break;
            bk[currV-1][currU] = beta;
            scaleVec(temp2, 1/beta, n);
            copyVec(v, temp2, n);
            reorthogonalize(v, vk, n, currV); 
            normalize(v, n);
            setCol(vk, v, n, currV);
            currV++;
        }
    }

    // Free Memory
    free(temp1);
    free(temp2);
    free(v);
    free(Av);
    free(u);
    free(AtU);
    freeMat(At, n, m);
}