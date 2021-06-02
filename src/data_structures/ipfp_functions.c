#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include "ipfp_functions.h"


double_dense_matrix get_alphas_row(double_sparse_matrix *working_matrix, double_dense_matrix *local_cbg, int i){
    double_dense_matrix alphas, dense_poi;
    create_double_dense_matrix(&alphas, 1, local_cbg->n_cols);

    // sum along columns
    memset( alphas.data, 0, local_cbg->n_cols*sizeof(double) );
    for (int j=0; j<working_matrix->n_elements; j++){
        int idx = working_matrix->data[j].col;
        alphas.data[idx] += working_matrix->data[j].val;
    }
    
    // alphas
    for (int j=0; j<local_cbg->n_cols; j++){
        if (alphas.data[j]>0)
            alphas.data[j] = local_cbg->data[(i*alphas.n_cols)+j] / alphas.data[j];
    }

    return alphas;
}


double_dense_matrix get_alphas_col(double_sparse_matrix *working_matrix, double_sparse_matrix *local_poi, int i){
    double_dense_matrix alphas, dense_poi;
    create_double_dense_matrix(&alphas, local_poi->n_rows, 1);

    // sum along columns
    memset( alphas.data, 0, alphas.n_rows*sizeof(double) );
    for (int j=0; j<working_matrix->n_elements; j++){
        int idx = working_matrix->data[j].row;
        alphas.data[idx] += working_matrix->data[j].val;
    }

    // make poi in dense format to speed up calculations
    create_double_dense_matrix(&dense_poi, alphas.n_cols, 1);
    memset( dense_poi.data, 0, alphas.n_rows*sizeof(double) );
    for (int j=0; j<local_poi->n_elements; j++){
        if (local_poi->data[j].col==i)
            dense_poi.data[local_poi->data[j].row] = local_poi->data[j].val;
    }
    
    // alphas
    for (int j=0; j<alphas.n_rows; j++){
        if (alphas.data[j]>0)
            alphas.data[j] = dense_poi.data[j] / alphas.data[j];
    }

    free_double_dense_matrix(&dense_poi);

    return alphas;
}


void multiply_cols_by_alphas(double_sparse_matrix *working_matrix, double_dense_matrix *alphas){
    for (int i=0; i<working_matrix->n_elements; i++){
        double coeff = alphas->data[working_matrix->data[i].col];
        working_matrix->data[i].val *= coeff;
        if (working_matrix->data[i].val<1.)
            working_matrix->data[i].val = 1.;
    }
}
void multiply_rows_by_alphas(double_sparse_matrix *working_matrix, double_dense_matrix *alphas){
    for (int i=0; i<working_matrix->n_elements; i++){
        double coeff = alphas->data[working_matrix->data[i].row];
        working_matrix->data[i].val *= coeff;
        if (working_matrix->data[i].val<1.)
            working_matrix->data[i].val = 1.;
    }
}
