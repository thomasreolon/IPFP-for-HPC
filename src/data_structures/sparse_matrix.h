#ifndef __SPARSE_MATRIX_H
#define __SPARSE_MATRIX_H 1

#include "dense_matrix.h"

typedef struct matrix_element_struct
{
    int row;
    int col;
    double val;
} matrix_element;

typedef struct sparse_matrix_struct
{
    int n_rows;
    int n_cols;
    int n_elements;
    matrix_element* data;
} double_sparse_matrix;

/**
 * Initilizes a sparse matrix allocating space for the data
 **/
void create_double_sparse_matrix(double_sparse_matrix *mat, int n_rows, int n_cols, int n_elements);

/**
 * Frees the allocated space for the matrix's data
 **/
void free_double_sparse_matrix(double_sparse_matrix *mat);


/**
 * Loads a sparse matrix from a file
 * 
 * example file
 * 3 5 2            --> n_rows n_cols n_elements
 * 0 0 100.0        --> element 1
 * 2 4 21.15        --> element 2
 **/
double_sparse_matrix load_double_sparse_matrix(const char *filename);


double_sparse_matrix clone_submatrix(const double_sparse_matrix original_submatrix);

void util_print_sparse(const double_sparse_matrix submatrix);

#endif