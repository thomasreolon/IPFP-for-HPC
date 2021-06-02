#ifndef __IPFP_FUNCTIONS_H
#define __IPFP_FUNCTIONS_H 1

#include "sparse_matrix.h"
#include "dense_matrix.h"


double_dense_matrix get_alphas_row(double_sparse_matrix *working_matrix, double_dense_matrix *local_cbg, int i,  int n_threads);


double_dense_matrix get_alphas_col(double_sparse_matrix *working_matrix, double_dense_matrix *local_poi, int i, int n_threads);

void multiply_cols_by_alphas(double_sparse_matrix *working_matrix, double_dense_matrix *alphas);

void multiply_rows_by_alphas(double_sparse_matrix *working_matrix, double_dense_matrix *alphas);

void save_result_matrix(const char *filename, matrix_element *elements, int n_rows, int n_cols, int n_elements);

double_dense_matrix sparse_to_dense(const double_sparse_matrix mat);

#endif