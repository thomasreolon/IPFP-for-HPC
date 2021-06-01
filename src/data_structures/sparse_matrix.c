#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include "sparse_matrix.h"


void create_double_sparse_matrix(double_sparse_matrix *mat, int n_rows, int n_cols, int n_elements){
    mat->n_cols = n_cols;
    mat->n_rows = n_rows;
    mat->n_elements = n_elements;
    mat->data = malloc(sizeof(matrix_element)*n_elements);
}

void free_double_sparse_matrix(double_sparse_matrix *mat){
    free(mat->data);
}


double_sparse_matrix load_double_sparse_matrix(const char *filename)
{
    int n_rows, n_cols, n_elements, i;
    double_sparse_matrix loaded_matrix;

    // load file
    FILE *in_file = fopen(filename, "r");

    // check file was opened correctly
    if (in_file == NULL)
    {
        printf("Error! Could not open file %s\n", filename);
        MPI_Abort( MPI_COMM_WORLD , -1);
        MPI_Finalize();
        exit(-1);
    }

    fscanf(in_file, "%d %d %d", &n_rows, &n_cols, &n_elements);

    create_double_sparse_matrix(&loaded_matrix, n_rows, n_cols, n_elements);

    // read each element of the matrix
    for (i = 0; i < n_elements; i++)
    {
        int row;
        int col;
        double value;
        fscanf(in_file, "%d %d %lf", &row, &col, &value);

        loaded_matrix.data[i].row = row;
        loaded_matrix.data[i].col = col;
        loaded_matrix.data[i].val = value;
    }
    fclose(in_file);

    return loaded_matrix;
}

double_sparse_matrix clone_submatrix(const double_sparse_matrix original_submatrix){
    double_sparse_matrix clone;
    memcpy(&clone, &original_submatrix, sizeof(double_sparse_matrix)); // copy n_rows, n_cols, ...
    clone.data = malloc(sizeof(matrix_element)*clone.n_elements);
    memcpy(clone.data, original_submatrix.data, sizeof(matrix_element)*clone.n_elements); // copy the data
    return clone;
}