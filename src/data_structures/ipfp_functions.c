#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <float.h>
#include "ipfp_functions.h"


double_dense_matrix get_alphas_row(double_sparse_matrix *working_matrix, double_dense_matrix *local_cbg, int i){
    double_dense_matrix alphas;
    create_double_dense_matrix(&alphas, 1, local_cbg->n_cols);

    // sum along columns
    memset( alphas.data, 0, local_cbg->n_cols*sizeof(double) );
    for (int j=0; j<working_matrix->n_elements; j++){
        int idx = working_matrix->data[j].col;
        alphas.data[idx] += working_matrix->data[j].val;
    }
    
    // alphas
    for (int j=0; j<local_cbg->n_cols; j++){
        if (alphas.data[j]<FLT_EPSILON)
            alphas.data[j] = 1.;
        alphas.data[j] = local_cbg->data[(i*alphas.n_cols)+j] / alphas.data[j];
    }

    return alphas;
}

double_dense_matrix get_alphas_col(double_sparse_matrix *working_matrix, double_dense_matrix *local_poi, int i){
    double_dense_matrix alphas;
    create_double_dense_matrix(&alphas, local_poi->n_rows, 1);

    // sum along columns
    memset( alphas.data, 0, alphas.n_rows*sizeof(double) );
    for (int j=0; j<working_matrix->n_elements; j++){
        int idx = working_matrix->data[j].row;
        alphas.data[idx] += working_matrix->data[j].val;
    }
    
    // alphas
    for (int j=0; j<alphas.n_rows; j++){
        if (alphas.data[j]<FLT_EPSILON)
            alphas.data[j] = 1.;
        alphas.data[j] = local_poi->data[(j*alphas.n_cols)+i] / alphas.data[j];
    }

    return alphas;
}



void multiply_cols_by_alphas(double_sparse_matrix *working_matrix, double_dense_matrix *alphas){
    for (int i=0; i<working_matrix->n_elements; i++){
        double coeff = alphas->data[working_matrix->data[i].col];
        working_matrix->data[i].val *= coeff;
    }
}
void multiply_rows_by_alphas(double_sparse_matrix *working_matrix, double_dense_matrix *alphas){
    for (int i=0; i<working_matrix->n_elements; i++){
        double coeff = alphas->data[working_matrix->data[i].row];
        working_matrix->data[i].val *= coeff;
    }
}


void save_result_matrix(const char *filename, matrix_element *elements, int n_rows, int n_cols, int n_elements)
{
    FILE *out_file = fopen(filename, "w+");

    printf("writing %s\n", filename);

    // check file was opened correctly
    if (out_file == NULL)
    {
        printf("Error! Could not open file\n");
        exit(-1);
    }

    // write the first line
    fprintf(out_file, "%d %d %d\n", n_rows, n_cols, n_elements);

    int i;
    // write elements
    for (i = 0; i < n_elements; i++)
    {
        fprintf(out_file, "%d %d %e\n", elements[i].row, elements[i].col, elements[i].val);
    }
}

double_dense_matrix sparse_to_dense(const double_sparse_matrix mat){
    double_dense_matrix dense;

    create_double_dense_matrix(&dense, mat.n_rows, mat.n_cols);
    memset(dense.data, 0, sizeof(double)*mat.n_rows*mat.n_cols);
    for(int i=0; i<mat.n_elements; i++){
        int idx = mat.data[i].row * mat.n_cols + mat.data[i].col;
        dense.data[idx] = mat.data[i].val;
    }
    return dense;
}
