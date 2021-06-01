#include <stdio.h>
#include <stdlib.h>
#include "dense_matrix.h"

void create_double_dense_matrix(double_dense_matrix *mat, int n_rows, int n_cols) {
    mat->n_rows = n_rows;
    mat->n_cols = n_cols;
    mat->data = malloc(n_rows * n_cols * sizeof(double));
}

void clean_double_dense_matrix(double_dense_matrix* matrix) {
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
