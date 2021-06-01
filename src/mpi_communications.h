#include "data_structures/sparse_matrix.h"
#include "data_structures/dense_matrix.h"


#ifndef __MPI_COMMUNICATIONS_H
#define __MPI_COMMUNICATIONS_H 1


double_sparse_matrix broadcast_sparse_matrix(double_sparse_matrix *mat, int rank);


double_sparse_matrix scatter_poi_matrix(double_sparse_matrix *mat, int n_rows, int tot_columns, int *counts, int* displacements);


double_dense_matrix scatter_cbg_marginals(double_dense_matrix* mat, int tot_rows, int n_cols, int *counts, int* displacements);


#endif