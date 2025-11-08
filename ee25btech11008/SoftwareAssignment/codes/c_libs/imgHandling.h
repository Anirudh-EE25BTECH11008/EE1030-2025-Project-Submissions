#ifndef IMG_HANDLING_H
#define IMG_HANDLING_H

typedef struct {
    int width;
    int height;
    double **gray;
    double **R;
    double **G;
    double **B;
} ImageData;

ImageData *loadImageSTB(const char *input_filename, const char *input_filename_final, int* channels);
void writePGM(const char *output_filename, double **array, int width, int height);
void writePPM(const char *output_filename, double **R, double **G, double **B, int width, int height);
void saveCompressedImage(const char *base_filename1, const char *base_filename2, ImageData *img, int isColor, int k);
void freeImageData(ImageData *img);
void writeJPEG(const char *filename, ImageData *img, int w, int h, int isColor);;

#endif
