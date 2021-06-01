#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <stdio.h>
#include "data_structures/sparse_matrix.h"
#include "data_structures/dense_matrix.h"

void build_mpi_tuple(MPI_Datatype *mpi_tuple)
{
    matrix_element element;
    // array of structure member sizes
    int blocklengths[3] = {1, 1, 1};

    // structure member types
    MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_DOUBLE};

    // offset of structure members
    MPI_Aint offsets[3];
    MPI_Aint points[3];
    MPI_Get_address(&element.row, points);
    MPI_Get_address(&element.col, points + 1);
    MPI_Get_address(&element.val, points + 2);
    offsets[0] = 0;
    offsets[1] = points[1] - points[0];
    offsets[2] = points[2] - points[0];

    // create mpi struct
    MPI_Type_create_struct(3, blocklengths, offsets, types, mpi_tuple);
    MPI_Type_commit(mpi_tuple);
}

double_sparse_matrix broadcast_sparse_matrix(double_sparse_matrix *mat, int rank){
    int infos[3];
    MPI_Datatype mpi_tuple;

    infos[0] = mat->n_rows;
    infos[1] = mat->n_cols;
    infos[2] = mat->n_elements;
    
    MPI_Bcast( infos , 3 , MPI_INT , 0 , MPI_COMM_WORLD);

    if (rank!=0){
        create_double_sparse_matrix(mat, infos[0], infos[1], infos[2]);
    }

    build_mpi_tuple(&mpi_tuple);
    MPI_Bcast( mat->data , mat->n_elements , mpi_tuple , 0 , MPI_COMM_WORLD);
    MPI_Type_free(&mpi_tuple);
}

double_sparse_matrix scatter_poi_matrix(double_sparse_matrix *mat, int n_rows, int tot_columns, int *counts, int* displacements){
    int count_elements_per_hour[tot_columns], offsets[tot_columns];
    int world_size, rank;
    MPI_Datatype mpi_tuple;
    double_sparse_matrix local_mat;
    matrix_element *buffer;


    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    world_size--;

    if (rank==0){
        // to have contiguos columns
        memset( count_elements_per_hour, 0, tot_columns*sizeof(int) );
        buffer = malloc(mat->n_elements*sizeof(matrix_element));
        for (int i=0; i<mat->n_elements; i++)
            count_elements_per_hour[mat->data[i].col]++;
        offsets[0] =  count_elements_per_hour[0];
        for (int i=1; i<tot_columns; i++)
            offsets[i] = offsets[i-1] + count_elements_per_hour[i];

        for (int i=0; i<mat->n_elements; i++){
            int index = --offsets[mat->data[i].col];
            buffer[index] = mat->data[i];
        }
    }
    
    /// counts and displacements for scatterv
    int hour_per_proc = tot_columns / world_size;
    memset( counts, 0, (world_size+1)*sizeof(int) );
    MPI_Bcast( count_elements_per_hour , tot_columns , MPI_INT , 0 , MPI_COMM_WORLD);
    MPI_Bcast( offsets , tot_columns , MPI_INT , 0 , MPI_COMM_WORLD);
    counts[0]=0; displacements[0] = 0;
    for (int i=0; i<world_size; i++){
        for (int j=i*hour_per_proc; j<(i+1)*hour_per_proc; j++){
            counts[i+1] += count_elements_per_hour[j];
        }
        displacements[i+1] = offsets[i*hour_per_proc];
    }

    /// last process will do the last hours
    for (int j=world_size*hour_per_proc; j<tot_columns; j++){
        counts[world_size] += count_elements_per_hour[j];
    }

    // scatter hours between processes
    hour_per_proc = hour_per_proc + (rank==world_size)? tot_columns % world_size : 0;
    build_mpi_tuple(&mpi_tuple);
    create_double_sparse_matrix(&local_mat, n_rows, hour_per_proc, counts[rank]);
    MPI_Scatterv(buffer, counts, displacements, mpi_tuple, local_mat.data, counts[rank], mpi_tuple, 0, MPI_COMM_WORLD);

    MPI_Type_free(&mpi_tuple);
    if (rank==0) free(buffer);
    return local_mat;
}


double_dense_matrix scatter_cbg_marginals(double_dense_matrix* mat, int tot_rows, int n_cols, int *counts, int* displacements){
    double_dense_matrix local_mat;
    int world_size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    world_size--; // process 0 does only I/O on files
    
    // compute counts / offsets
    int hour_per_proc = tot_rows / world_size;
    counts[0]=0; displacements[0] = 0;
    for (int i=0; i<world_size; i++){
        counts[i+1] = n_cols * hour_per_proc;
        displacements[i+1] = i*hour_per_proc * n_cols;
    }
    counts[world_size] += (tot_rows % world_size)*n_cols;

    // allocate space to receive data
    hour_per_proc += (rank==world_size)? (tot_rows % world_size) : 0;
    create_double_dense_matrix(&local_mat, hour_per_proc, n_cols);

    // send data
    MPI_Scatterv(mat->data, counts, displacements, MPI_DOUBLE, local_mat.data, counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    return local_mat;
}