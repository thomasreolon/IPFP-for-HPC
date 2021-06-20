#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <stdio.h>
#include "data_structures/sparse_matrix.h"
#include "data_structures/dense_matrix.h"
#include "data_structures/ipfp_functions.h"


int get_num_hours(double_dense_matrix *cbg, int rank, int argc, char *argv[]){
    int num_hours = cbg->n_rows;
    if (argc==7 && atoi(argv[6])>0)
        num_hours = atoi(argv[6]);            // user specified the number of hours
    else
        MPI_Bcast( &num_hours , 1 , MPI_INT , 0 , MPI_COMM_WORLD);
    return num_hours;
}


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

double_dense_matrix scatter_poi_marginals(double_dense_matrix* mat, int tot_rows, int n_cols, int *counts, int* displacements){
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



double_dense_matrix scatter_marginals(double_dense_matrix* mat, int tot_rows, int n_cols, int *counts, int* displacements){
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


void send_result_matrix_to_master(const double_sparse_matrix mat, int hour){
    int infos[3];
    MPI_Datatype mpi_tuple;

    mat.data[mat.n_elements].row = hour;  //--> we use the attribute row to append the info "which hour is this matrix"

    build_mpi_tuple(&mpi_tuple);
    MPI_Send( mat.data, mat.n_elements+1, mpi_tuple, 0, 0, MPI_COMM_WORLD);
    MPI_Type_free(&mpi_tuple);
}

void receive_and_save_result_matrix(const char* path, const double_sparse_matrix aggregate_mat){
    int infos[3];
    char filename[1024];
    MPI_Datatype mpi_tuple;
    MPI_Status status;
    matrix_element *buffer = malloc((aggregate_mat.n_elements+1)*sizeof(matrix_element));

    build_mpi_tuple(&mpi_tuple);
    MPI_Recv(buffer, aggregate_mat.n_elements+1, mpi_tuple, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
    if (status.MPI_ERROR != MPI_SUCCESS){
        printf("error while receivind data from %d\n", status.MPI_SOURCE);
    }
    MPI_Type_free(&mpi_tuple);

    int hour = buffer[aggregate_mat.n_elements].row; // we used the trik of saving the hour in the last element
    sprintf(filename, "%s/result_%d.txt", path, hour);
    save_result_matrix(filename, buffer, aggregate_mat.n_rows, aggregate_mat.n_cols, aggregate_mat.n_elements);

    free(buffer);
}

