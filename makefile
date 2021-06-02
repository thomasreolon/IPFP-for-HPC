
run: 
		mpicc ./src/data_structures/dense_matrix.c ./src/data_structures/sparse_matrix.c ./src/mpi_communications.c ./src/main.c
		mpirun -n 2 ./a.out ./data/agg.txt ./data/poi.txt ./data/cbg.txt ./data

run3: 
		mpicc ./src/data_structures/dense_matrix.c ./src/data_structures/sparse_matrix.c ./src/data_structures/ipfp_functions.c ./src/mpi_communications.c ./src/main.c
		mpirun -n 3 ./a.out ./data/agg.txt ./data/poi.txt ./data/cbg.txt ./data


#target "clean" non è un file!
.PHONY: run
.SILENT: run