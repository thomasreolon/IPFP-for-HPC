#ifndef __DENSE_MATRIX_H
#define __DENSE_MATRIX_H 1

typedef struct dense_matrix_struct
{
    int n_rows;
    int n_cols;
    double* data;
} double_dense_matrix;


void create_double_dense_matrix(double_dense_matrix *mat, int n_rows, int n_cols);

void clean_double_dense_matrix(double_dense_matrix* matrix);

double_dense_matrix load_double_dense_matrix(const char *filename);

#endif