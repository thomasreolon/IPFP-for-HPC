#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "data_structures/sparse_matrix.h"
#include "data_structures/dense_matrix.h"
#include "data_structures/ipfp_functions.h"
#include "mpi_communications.h"

#define NUM_ITERATIONS 5

int main(int argc, char** argv) {
    if (argc < 6) {
        printf("Usage: distributedIPFP [aggregate_visit_matrix] [week_poi_marginals] [week_cbg_marginals] [output_dir] [n_threads] <num_hour>\n");
        exit(1);
    }

    // init MPI
	int world_size, world_rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // read input data from files.txt
    double_sparse_matrix aggregate_visit_matrix;
    double_sparse_matrix tmp_poi;
    double_dense_matrix poi_marginals_matrix;
    double_dense_matrix cbg_marginals_matrix;
    if (world_rank == 0) {
        // load matrices from files
        aggregate_visit_matrix = load_double_sparse_matrix(argv[1]);
        cbg_marginals_matrix   = load_double_dense_matrix(argv[3]);

        tmp_poi                = load_double_sparse_matrix(argv[2]);
        poi_marginals_matrix   = sparse_to_dense(tmp_poi);
        transpose_dense_matrix(&poi_marginals_matrix);
        free_double_sparse_matrix(&tmp_poi);
    }
    // how many output matrices
    const int n_threads = atoi(argv[5]);
    const int num_hours = get_num_hours(&cbg_marginals_matrix, world_rank, argc, argv);
    int i, j; // used in for loops

    int poi_counts[world_size], poi_displacements[world_size];
    int cbg_counts[world_size], cbg_displacements[world_size];
    double_dense_matrix local_poi, local_cbg;


    // broadcast matrix
    broadcast_sparse_matrix(&aggregate_visit_matrix, world_rank);

    // scatter cbg marginals
    local_cbg = scatter_marginals(&cbg_marginals_matrix, num_hours, aggregate_visit_matrix.n_cols, cbg_counts, cbg_displacements);

    // scatter poi marginals across processes
    local_poi = scatter_marginals(&poi_marginals_matrix, num_hours, aggregate_visit_matrix.n_rows, poi_counts, poi_displacements);
    transpose_dense_matrix(&local_poi);

    if (world_rank==0){
        free_double_dense_matrix(&poi_marginals_matrix);
        free_double_dense_matrix(&cbg_marginals_matrix);

        // write matrix when received
        for (i=0; i<num_hours; i++)
            receive_and_save_result_matrix(argv[4], aggregate_visit_matrix);
            
    }else{
        // elaborate matrix
# pragma omp parallel for num_threads(n_threads)
        for (i=0; i<local_poi.n_cols; i++){
            double_sparse_matrix working_matrix = clone_submatrix(aggregate_visit_matrix);

            // hour  =  offset wrt the other processes         + actual processed column
            int hour = (num_hours / (world_size-1))*(world_rank-1) + i;

            //// IPFP algorithm
            for (j=0; j<NUM_ITERATIONS; j++){
                if (j%2==0){
                    // alphas = cbg_marginals /  sum_along_columns
                    double_dense_matrix alphas = get_alphas_row(&working_matrix, &local_cbg, i, n_threads);
                    multiply_cols_by_alphas(&working_matrix, &alphas);
                    free_double_dense_matrix(&alphas);
                }else{
                    // alphas = poi_marginals /  sum_along_rows
                    double_dense_matrix alphas = get_alphas_col(&working_matrix, &local_poi, i, n_threads);  // hour->the values of the poi are not shifted 
                    multiply_rows_by_alphas(&working_matrix, &alphas);
                    free_double_dense_matrix(&alphas);
                }
            }
            send_result_matrix_to_master(working_matrix, hour);

            // free resources
            free_double_sparse_matrix(&working_matrix);
        }
        free_double_dense_matrix(&local_poi);
        free_double_dense_matrix(&local_cbg);
    }
    free_double_sparse_matrix(&aggregate_visit_matrix);


    MPI_Finalize();
    return 0;
}