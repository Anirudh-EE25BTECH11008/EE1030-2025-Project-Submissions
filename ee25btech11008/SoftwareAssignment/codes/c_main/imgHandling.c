#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../c_libs/stb_image.h"
#include "../c_libs/stb_image_write.h"
#include "../c_libs/imgHandling.h"
#include "../c_libs/matrixOps.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static void ucharToDouble(double **mat, unsigned char *data, int w, int h, int c, int offset) {
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            mat[y][x] = (double)data[(y * w + x) * c + offset];
        }
    }
}

ImageData *loadImageSTB(const char *input_filename1, const char *input_filename2, int* channels) {
    int w, h, c;
    unsigned char *data = stbi_load(input_filename1, &w, &h, &c, 0);
    if (!data) {
        fprintf(stderr, "Error: cannot open %s\n", input_filename1);
        return NULL;
    }
    
    printf("Loaded %s - Width=%d, Height=%d, Channels=%d\n", input_filename1, w, h, c);

    ImageData *img = malloc(sizeof(ImageData));
    img->width = w;
    img->height = h;
    img->gray = img->R = img->G = img->B = NULL;

    if (c>2) {
        img->R = createMatrix(h, w);
        img->G = createMatrix(h, w);
        img->B = createMatrix(h, w);
        ucharToDouble(img->R, data, w, h, c, 0);
        ucharToDouble(img->G, data, w, h, c, 1);
        ucharToDouble(img->B, data, w, h, c, 2);
        int totalPixels = h * w;
        int colorPixels = 0;
        
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                double r = img->R[y][x], g = img->G[y][x], b = img->B[y][x];
                double diff1 = fabs(r - g);
                double diff2 = fabs(r - b);
                double maxPixelDiff = (diff1 > diff2 ? diff1 : diff2);
                if (maxPixelDiff > 2.0) {
                    colorPixels++;
                }
            }
        }

        double colorPercentage = (100.0 * colorPixels) / totalPixels;
        int isGrayscale = (colorPercentage < 2);
        
        if (isGrayscale) {
            img->gray = createMatrix(h, w);
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    img->gray[y][x] = img->G[y][x]; 
                }
            }
            freeMat(img->R, h, w);
            freeMat(img->G, h, w);
            freeMat(img->B, h, w);
            img->R = img->G = img->B = NULL;
            *channels = 0; 
        } else {
            *channels = 1; 
        }
    } else {
        img->gray = createMatrix(h, w);
        ucharToDouble(img->gray, data, w, h, c >= 1 ? c : 1, 0);
        *channels = 0;
    }

    stbi_image_free(data);

    printf("Writing output with mode=%d (%s)\n", *channels, (*channels == 0 ? "PGM" : "PPM"));
    char outpath[256];
    snprintf(outpath, sizeof(outpath), "%s.%s", input_filename2, (*channels == 0 ? "pgm" : "ppm"));
    if (*channels == 0) {
        writePGM(outpath, img->gray, w, h);
    } else {
        writePPM(outpath, img->R, img->G, img->B, w, h);
    }

    return img;
}

void writePGM(const char *filename, double **array, int w, int h) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) { 
        fprintf(stderr, "Cannot write %s\n", filename); 
        return; 
    }
    fprintf(fp, "P5\n%d %d\n255\n", w, h);
    for (int y = 0; y < h; y++) {
        unsigned char *row = malloc(w);
        for (int x = 0; x < w; x++) {
            double val = array[y][x];
            if (val < 0) val = 0;
            if (val > 255) val = 255;
            row[x] = (unsigned char)(val + 0.5);
        }
        fwrite(row, 1, w, fp);
        free(row);
    }
    fclose(fp);
}

void writePPM(const char *filename, double **R, double **G, double **B, int w, int h) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) { 
        fprintf(stderr, "Cannot write %s\n", filename); 
        return; 
    }
    fprintf(fp, "P6\n%d %d\n255\n", w, h);
    unsigned char *row = malloc(3 * w);
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            double r = R[y][x], g = G[y][x], b = B[y][x];
            if (r < 0) r = 0; if (r > 255) r = 255;
            if (g < 0) g = 0; if (g > 255) g = 255;
            if (b < 0) b = 0; if (b > 255) b = 255;
            row[3*x] = (unsigned char)(r + 0.5);
            row[3*x+1] = (unsigned char)(g + 0.5);
            row[3*x+2] = (unsigned char)(b + 0.5);
        }
        fwrite(row, 3, w, fp);
    }
    free(row);
    fclose(fp);
}

void writeJPEG(const char *filename, ImageData *img, int w, int h, int isColor) {
    unsigned char *data;
    int channels = isColor ? 3 : 1;
    data = malloc(w * h * channels);
    
    if (isColor) {
        // RGB image
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                double r = img->R[y][x], g = img->G[y][x], b = img->B[y][x];
                if (r < 0) r = 0; if (r > 255) r = 255;
                if (g < 0) g = 0; if (g > 255) g = 255;
                if (b < 0) b = 0; if (b > 255) b = 255;
                int idx = (y * w + x) * 3;
                data[idx] = (unsigned char)(r + 0.5);
                data[idx + 1] = (unsigned char)(g + 0.5);
                data[idx + 2] = (unsigned char)(b + 0.5);
            }
        }
    } else {
        // Grayscale image
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                double val = img->gray[y][x];
                if (val < 0) val = 0;
                if (val > 255) val = 255;
                data[y * w + x] = (unsigned char)(val + 0.5);
            }
        }
    }
    
    stbi_write_jpg(filename, w, h, channels, data, 90);
    free(data);
}

void saveCompressedImage(const char *base_filename1, const char *base_filename2, ImageData *img, int isColor, int k) {
    char outpath[256];
    snprintf(outpath, sizeof(outpath), "%s_%d.%s", base_filename1, k, (isColor ? "out.ppm" : "out.pgm"));
    if (isColor) {
        writePPM(outpath, img->R, img->G, img->B, img->width, img->height);
    } else {
        writePGM(outpath, img->gray, img->width, img->height);
    }
    char jpegpath[256];
    snprintf(jpegpath, sizeof(jpegpath), "%s_%d.out.jpg", base_filename2, k);
    writeJPEG(jpegpath, img, img->width, img->height, isColor);
}

void freeImageData(ImageData *img) {
    if (!img) return;
    if (img->gray) freeMat(img->gray, img->height, img->width);
    if (img->R) freeMat(img->R, img->height, img->width);
    if (img->G) freeMat(img->G, img->height, img->width);
    if (img->B) freeMat(img->B, img->height, img->width);
    free(img);
}