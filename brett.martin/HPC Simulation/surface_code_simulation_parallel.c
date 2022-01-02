/*	2d Lt Brett Martin
 *	Advisor: Dr. Laurence Merkle
 *	CSCE 656 - Parallel and Distributed Processing Algorithms
 *	Term project: parallel surface code simulation
 **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>


int main (int argc, char** argv) {

	//Inputs:
	int depth = 3;		// Or code distance
	int p_x_error;
	int p_z_error;
	int p_m_error; 		// eventually

	//Output:
	//list of records in which each possible combination of data qubit errors is associated with the resulting ancilla qubit values (may be probabilistic), along with the probability of the combination of data qubit errors (and measurement errors?)
	FILE *f = fopen("/tmp/brett/output.txt", "a");	// Change directory

	int num_data_qubits = depth ** 2;
	int data_qubit_x_error[ depth + 2][ depth + 2 ];
	int data_qubit_z_error[ depth + 2  ][ depth + 2 ];
	int ancilla_qubit_value[ depth + 1 ][ depth + 1 ];

	// MPI vars (NOTE: Per the Navy DSRC intro guide, the standard queue has a maximum cores per job of 8,168 cores)
	int iproc, nproc;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	// Parallel block decomposition variables
	int total_num_iter = 4 ** num_data_qubits;			// Total number of outer loop iterations
	int block_size = floor(total_num_iter/nproc);		// Number of iterations per process (block)

	if (total_num_iter%nproc > 0) { block_size += 1; }	// Add 1 if blocks don't divide evenly

	int iter_first = iproc * block_size;
	int iter_last = iter_first + block_size;

	for ( int i = iter_first; i < iter_last; i++ ) {

		// Exclude other processes outside of block
		if (i >= total_num_iter || (iproc >= iter_first && iproc < iter_last)) { continue; }

		int errors = i;
		double probability = 1.0;

		// initialize data_qubit_x_error and data_qubit_z_error to be full of FALSE (or 0) values
		memset(data_qubit_x_error, 0, sizeof(data_qubit_x_error[0][0]) * (depth + 2) * (depth + 2));
		memset(data_qubit_z_error, 0, sizeof(data_qubit_z_error[0][0]) * (depth + 2) * (depth + 2));
		memset(final_qubit_x_error, 0, sizeof(final_qubit_x_error[0][0]) * (depth + 2) * (depth + 2));
		memset(final_qubit_z_error, 0, sizeof(final_qubit_z_error[0][0]) * (depth + 2) * (depth + 2));

		for ( int j = 0; j < num_data_qubits; j++ ) {
			int this_error = errors % 2;
			data_qubit_x_error[ ( j / depth ) + 1][ ( j % depth ) + 1 ] = this_error;
			probability *= this_error ? p_x_error : (1.0 - p_x_error);
			errors /= 2;

			this_error = errors % 2;
			data_qubit_z_error[ ( j / depth ) + 1][ ( j % depth ) + 1 ] = this_error;
			errors /= 2;
			probability *= this_error ? p_z_error : (1.0 - p_z_error);
		}

		for ( int j = 0; j < ( depth + 1 ) ** 2 - 1; j++ ) {
			ancilla_qubit_value[ j / depth ][ j % depth ] = 
				  data_qubit_x_error[   j / depth       ][   j % depth       ]
				^ data_qubit_x_error[   j / depth       ][ ( j % depth ) + 1 ]
				^ data_qubit_x_error[ ( j / depth ) + 1 ][   j % depth       ]
				^ data_qubit_x_error[ ( j / depth ) + 1 ][ ( j % depth ) + 1 ];
			j++;
			ancilla_qubit_value[ j / depth ][ j % depth ] = 
				  data_qubit_z_error[   j / depth       ][   j % depth       ]
				^ data_qubit_z_error[   j / depth       ][ ( j % depth ) + 1 ]
				^ data_qubit_z_error[ ( j / depth ) + 1 ][   j % depth       ]
				^ data_qubit_z_error[ ( j / depth ) + 1 ][ ( j % depth ) + 1 ];
		}

		// Output the probability, actual data_qubit_x_error, data_qubit_z_error, and ancilla_qubit_value
		fwrite(probability, sizeof(double), sizeof(probability), f);
		fwrite(final_qubit_x_error, sizeof(int), sizeof(final_qubit_x_error), f);
		fwrite(final_qubit_z_error, sizeof(int), sizeof(final_qubit_z_error), f);
		fwrite(final_ancilla_qubit_value, sizeof(int), sizeof(final_ancilla_qubit_value), f);

	}
	
	MPI_Finalize();
	fclose(f);

	return 0;
}