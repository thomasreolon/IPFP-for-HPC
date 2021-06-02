#include "data_structures/sparse_matrix.h"
#include "data_structures/dense_matrix.h"


#ifndef __MPI_COMMUNICATIONS_H
#define __MPI_COMMUNICATIONS_H 1

int get_num_hours(double_dense_matrix *poi, int rank, int argc, char *argv[]);


double_sparse_matrix broadcast_sparse_matrix(double_sparse_matrix *mat, int rank);


double_dense_matrix scatter_marginals(double_dense_matrix* mat, int tot_rows, int n_cols, int *counts, int* displacements);

void send_result_matrix_to_master(const double_sparse_matrix mat, int hour);

void receive_and_save_result_matrix(const char* path, const double_sparse_matrix aggregate_mat);

#endif