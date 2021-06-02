#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dense_matrix.h"

void create_double_dense_matrix(double_dense_matrix *mat, int n_rows, int n_cols) {
    mat->n_rows = n_rows;
    mat->n_cols = n_cols;
    mat->data = malloc(n_rows * n_cols * sizeof(double));
}

void free_double_dense_matrix(double_dense_matrix* matrix) {
    free(matrix->data);
}

double_dense_matrix load_double_dense_matrix(const char *filename) {
    int n_rows, n_cols, i;
    int row, col;
    double value;
    double_dense_matrix loaded_matrix;

    // load file
    FILE *in_file = fopen(filename, "r");

    // check file was opened correctly
    if (in_file == NULL)
    {
        printf("Error! Could not open file\n");
        exit(-1);
    }

    fscanf(in_file, "%d %d", &n_rows, &n_cols);

    create_double_dense_matrix(&loaded_matrix, n_rows, n_cols);

    // read each element of the matrix
    for (i = 0; i < n_rows * n_cols; i++)
    {
        
        fscanf(in_file, "%d %d %lf", &row, &col, &value);

        loaded_matrix.data[row * n_cols + col] = value;
    }
    fclose(in_file);
    return loaded_matrix;
}

void util_print_dense(const double_dense_matrix mat){
    printf("dense [%dx%d]\n", mat.n_rows, mat.n_cols);
    for (int i=0; i<mat.n_cols*mat.n_rows; i++)
        printf("- (%d,%d): %lf\n", i/mat.n_cols, i%mat.n_cols, mat.data[i]);
    printf("\n");
}


void transpose_dense_matrix(double_dense_matrix *mat){
    int n_cols = mat->n_cols, n_rows=mat->n_rows;
    double *transposed = malloc(n_cols*n_rows*sizeof(double));

    for (int i=0; i<mat->n_cols*mat->n_rows; i++){
        int new_idx = i/n_cols   + (i%n_cols)*n_rows;
        transposed[new_idx] = mat->data[i];
    }
    free(mat->data);
    mat->data = transposed;
    mat->n_cols = n_rows;
    mat->n_rows = n_cols;
}