#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "data_structures/sparse_matrix.h"
#include "data_structures/dense_matrix.h"
#include "data_structures/ipfp_functions.h"
#include "mpi_communications.h"

#define NUM_ITERATIONS 1
#define NUM_HOURS 3

int main(int argc, char** argv) {
    if (argc != 5) {
        printf("Usage: distributedIPFP [aggregate_visit_matrix] [week_poi_marginals] [week_cbg_marginals] [output_dir]\n");
        exit(1);
    }

    // init MPI
	int world_size, world_rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // read input data from files.txt
    double_sparse_matrix aggregate_visit_matrix;
    double_sparse_matrix poi_marginals_matrix;
    double_dense_matrix cbg_marginals_matrix;
    if (world_rank == 0) {
        // load matrices from files
        aggregate_visit_matrix = load_double_sparse_matrix(argv[1]);
        poi_marginals_matrix = load_double_sparse_matrix(argv[2]);
        cbg_marginals_matrix = load_double_dense_matrix(argv[3]);
    }

    int poi_counts[world_size], poi_displacements[world_size];
    int cbg_counts[world_size], cbg_displacements[world_size];
    double_sparse_matrix local_poi;
    double_dense_matrix local_cbg;

    // broadcast matrix
    broadcast_sparse_matrix(&aggregate_visit_matrix, world_rank);

    // scatter poi marginals across processes
    local_poi = scatter_poi_matrix(&poi_marginals_matrix, aggregate_visit_matrix.n_rows, NUM_HOURS, poi_counts, poi_displacements);

    // scatter cbg marginals
    local_cbg = scatter_cbg_marginals(&cbg_marginals_matrix, NUM_HOURS, aggregate_visit_matrix.n_cols, cbg_counts, cbg_displacements);

    if (world_rank==0){
        // write matrix when received
    }else{
        // elaborate matrix
        for (int i=0; i<local_poi.n_cols; i++){
            double_sparse_matrix working_matrix = clone_submatrix(aggregate_visit_matrix);

            for (int j=0; j<NUM_ITERATIONS; j++){
                if (j%2==1){
                    // alphas = cbg_marginals /  sum_along_columns
                    double_dense_matrix alphas = get_alphas_row(&working_matrix, &local_cbg, i);
                    multiply_cols_by_alphas(&working_matrix, &alphas);
                }else{
                    printf("--------------------\n");
                    util_print_sparse(working_matrix);
                    // alphas = poi_marginals /  sum_along_rows
                    double_dense_matrix alphas = get_alphas_col(&working_matrix, &local_poi, i);
                    util_print_dense(alphas);
                    multiply_rows_by_alphas(&working_matrix, &alphas);
                    util_print_sparse(working_matrix);
                }
            }

        }

    }


    MPI_Finalize();
    return 0;
}