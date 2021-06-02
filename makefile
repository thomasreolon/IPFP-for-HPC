
N_PROC  = 2
N_HOURS = -1


help:
		echo "use 'make run' to see a proof of concept\n you have downloaded the dataset you can use 'make run_real'"
		echo "you can set N_PROC to the number of processes to use and N_HOURS to compute only the first X hours (useful with run_real)\n"
		echo "example:   make run N_PROC=4 N_HOURS=3"

data/agg.txt:
		echo "2 2 2\n0 0 50\n0 1 20" > ./data/agg.txt
		echo "3 2\n0 0 100\n0 1 40\n1 0 200\n1 1 80\n2 0 1000\n2 1 400" > ./data/cbg.txt
		echo "2 3 3\n0 0 140\n0 1 280\n0 2 1400" > ./data/poi.txt

run: data/agg.txt
		mpicc ./src/data_structures/dense_matrix.c ./src/data_structures/sparse_matrix.c ./src/data_structures/ipfp_functions.c ./src/mpi_communications.c ./src/main.c
		mpirun -n $(N_PROC) ./a.out ./data/agg.txt ./data/poi.txt ./data/cbg.txt ./data $(N_HOURS)


# to run this script you need the dataset 
run_real: data/aggregate_visit_matrix.txt
		mpicc ./src/data_structures/dense_matrix.c ./src/data_structures/sparse_matrix.c ./src/data_structures/ipfp_functions.c ./src/mpi_communications.c ./src/main.c
		mpirun -n $(N_PROC) ./a.out ./data/aggregate_visit_matrix.txt ./data/poi_marginals_2020_03_02.txt ./data/cbg_marginals_2020_03_02.txt ./data $(N_HOURS)



.PHONY: run run_real help
.SILENT: run run_real data/agg.txt help