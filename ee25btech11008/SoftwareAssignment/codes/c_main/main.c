#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../c_libs/matrixOps.h"
#include "../c_libs/lanczos.h"
#include "../c_libs/svd.h"
#include "../c_libs/imgHandling.h"

double** do_it(double** A, int m, int n, int k) {
    double** uk = createMatrix(m, k);
    double** vk = createMatrix(n, k);
    double** bk = createMatrix(k, k);
    bidiagonalize(A, uk, vk, bk, m, n, k);
    svd(uk, bk, vk, m, n, k);
    double** vt = transpose(vk, n, k);
    double** res_t = matMul(uk, bk, m, k, k);
    double** res = matMul(res_t, vt, m, k, n);
    freeMat(res_t, m, k);
    freeMat(uk, m, k);
    freeMat(vk, n, k);
    freeMat(bk, k, k);
    freeMat(vt, k, n);
    return res;
}

int main() {
    // Initialisation and Defaults
    int m, n, k, mode, l;
    double **A, **AR, **AB, **AG;
    char filename1[100];
    char filename2[100];
    char ext[10];
    char input_path[256];
    char output_path1[256];
    char output_path2[256];
    char output_path3[256];
    char final_answer[256];
    int kArr[6] = {5, 20, 50, 100, 200, 400};
    l = 6;
    ImageData *img;
    strcpy(filename1, "globe.jpg");
    strcpy(filename2, "globe");
    mode = 0;

    // Taking Input - Interactive Mode
    // printf("Enter k: ");
    // scanf("%d", &k); 
    // l = 1;
    // kArr[0] = k
    // printf("Enter Input Filename: ");
    // scanf("%s", filename1);
    // printf("Enter Output Filename: ");
    // scanf("%s", filename2);

    // Processing
    snprintf(input_path, sizeof(input_path), "../../figs/inp_jpg_png/%s", filename1);
    snprintf(output_path1, sizeof(output_path1), "../../figs/inp_pgm_ppm/%s", filename2);
    snprintf(output_path2, sizeof(output_path2), "../../figs/out_pgm_ppm/%s", filename2);
    snprintf(output_path3, sizeof(output_path3), "../../figs/out_jpg_png/%s", filename2);
    printf("Compressing for k=%d ...\n", kArr[0]);
    img = loadImageSTB(input_path, output_path1, &mode);
    m = img->height;
    n = img->width;
    for (int i=0; i<l; i++) {
        k = kArr[i];
        if (k<fmin(m, n)) {
            if (i!=0) {
                printf("\nCompressing for k=%d ...\n", kArr[i]);
                img = loadImageSTB(input_path, output_path1, &mode);
            }
            if (mode==0) {
                A = img->gray;
                img->gray = do_it(A, m, n, k);
                double frob_err= frobErr(A, img->gray, m, n);
                double frob_A = frob(A, m, n);
                printf("Absolute Frobenius norm: %lf\n", frob_err);
                printf("Relative Frobenius norm: %lf\n", frob_err/frob_A);
                printf("Total Energy captued (in %%): %lf\n", (1 - (frob_err*frob_err)/(frob_A*frob_A))*100);
            } else if (mode==1) {
                AR = img->R;
                AG = img->G;
                AB = img->B;
                img->R = do_it(AR, m, n, k);
                img->G = do_it(AG, m, n, k);
                img->B = do_it(AB, m, n, k);
                double frob_err1 = frobErr(AR, img->R, m, n);
                double frob_A1 = frob(AR, m, n);
                double frob_err2 = frobErr(AG, img->G, m, n);
                double frob_A2 = frob(AG, m, n);
                double frob_err3 = frobErr(AB, img->B, m, n);
                double frob_A3 = frob(AB, m, n);
                double frob_err = sqrt(frob_err1*frob_err1 + frob_err2*frob_err2 + frob_err3*frob_err3);
                double frob_A = sqrt(frob_A1*frob_A1 + frob_A2*frob_A2 + frob_A3*frob_A3);
                printf("Absolute Frobenius norm: %lf\n", frob_err);
                printf("Relative Frobenius norm: %lf\n", frob_err/frob_A);
                printf("Total Energy captued (in %%): %lf\n", (1 - (frob_err*frob_err)/(frob_A*frob_A))*100);
            }
            saveCompressedImage(output_path2, output_path3, img, mode, k);
            freeImageData(img);
        } else {
            break;
        }
    }

    // Free Memory
    if (mode==0) {
        freeMat(A, m, n);
    } else if (mode==1) {
        freeMat(AR, m, n);
        freeMat(AG, m, n);
        freeMat(AB, m, n);
    }
}

