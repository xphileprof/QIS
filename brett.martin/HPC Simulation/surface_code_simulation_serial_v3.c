/*	2d Lt Brett Martin
 *	Advisor: Dr. Laurence Merkle
 *	CSCE 656 - Parallel and Distributed Processing Algorithms
 *	Term project: parallel surface code simulation
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <inttypes.h>
#define TRI_SIZE 50


// Struct used in calculating and sorting all possible probabilities along with the number of samples per probability and the respective number of x and z errors
struct s_prob {
    int n_x_err, n_z_err, n_total_err;
    float probability, cumulative;
    unsigned long count;    // TODO: Change to double to have higher range
};


// Helper function passed into the quicksort algorithm to sort the struct s_prob by the probability values
int q_comp( const void *p1, const void *p2 ) {
    float a = ( (struct s_prob *) p1 )->probability;
    float b = ( (struct s_prob *) p2 )->probability;
    if (a < b) { return -1; }
    else if (a > b) { return 1; }
    return 0;
}


// Begin main function
int main ( int argc, char** argv ) {

	//Inputs:
	const int depth = 3;		// Or code distance
    const int num_samples = 1000;
	const float p_x_error = 0.05;
	const float p_z_error = 0.05;

	//Output:
	//list of records in which each possible combination of data qubit errors is associated with the resulting ancilla qubit values (may be probabilistic), along with the probability of the combination of data qubit errors (and measurement errors?)
	FILE *f = fopen("D:/Documents/AFIT/Research/test.csv", "w+");

    // Variables used in the generation of sample data
    int num_data_qubits = depth * depth;
	int data_qubit_x_error[ depth + 2][ depth + 2 ];
	int data_qubit_z_error[ depth + 2  ][ depth + 2 ];
	int ancilla_qubit_value[ depth + 1 ][ depth + 1 ];

    // Parallel block decomposition variables
    // TODO: Consider removing this
	uint64_t total_num_iter = pow(4, num_data_qubits);			// Total number of outer loop iterations

    // Generate pascal's triangle for use in calculating binomial coefficients ( time O(n^2) )
    // When referencing (n k) or "n choose k", just use: pascal_triangle[n][k]
    const int tri_size = 50;
    unsigned long pascal_triangle[ TRI_SIZE ][ TRI_SIZE ];

    for (int row = 0; row < tri_size; row++) {  // This generates a ton of empty space, but it should be okay with such negligible size
        for (int col = 0; col <= row; col++)
        {
            if (row == col || col == 0) {
                pascal_triangle[row][col] = (unsigned long) 1;
            } else {
                pascal_triangle[row][col] = (unsigned long) pascal_triangle[row-1][col] + pascal_triangle[row-1][col-1];
            }
        }
    }

    // Create CSV headers
    int ancilla = 0;
    for (int k=1; k < depth; k++) {
        if (k%2 == 0) {
            fprintf(f, "Z%d,", ancilla);
		    ancilla += 1;
        }
        for (int j=depth-1; j > 0; j--) {
            if (k == 1 && j%2 == 0) {
                fprintf(f, "X%d,", ancilla);
		        ancilla += 1;
            } else if (k == (depth - 1) && j%2 == 1) {
                fprintf(f, "X%d,", ancilla);
		        ancilla += 1;
            }
            if ((j+k)%2 == 0) {
                fprintf(f, "X%d,", ancilla);
		        ancilla += 1;
            } else {
                fprintf(f, "Z%d,", ancilla);
		        ancilla += 1;
            }
        }
        if (k%2 == 1) {
            fprintf(f, "Z%d,", ancilla);
		    ancilla += 1;
        }
    }
    fprintf(f, "Labels\n");

    // Determine the number of unique probabilities, the number of samples generated per probability
    // TODO: Comment out the following print statements after you're done debugging
    struct s_prob prob_values[num_data_qubits * num_data_qubits];
    memset(prob_values, 0, sizeof(prob_values));

    int prob_count = 0;
    for (int i = 0; i < num_data_qubits; i++) {
        for (int j = 0; j < num_data_qubits; j++) {
            prob_values[prob_count].n_total_err = prob_count;
            prob_values[prob_count].n_x_err = i;
            prob_values[prob_count].n_z_err = j;
            prob_values[prob_count].probability = (float) pow((1.0 - p_x_error), num_data_qubits - i) * pow((1.0 - p_z_error), num_data_qubits - j) * pow(p_x_error, i) * pow(p_z_error, j);
            prob_values[prob_count].count = (unsigned long) pascal_triangle[num_data_qubits][i] * pascal_triangle[num_data_qubits][j];
            printf("At index %d: prob: %e, count: %lu, errors: %d x and %d z\n", prob_values[prob_count].n_total_err, prob_values[prob_count].probability, prob_values[prob_count].count, prob_values[prob_count].n_x_err, prob_values[prob_count].n_z_err);
            prob_count++;
        }
    }

    printf("\n");

    // Sort by ascending probability and calculate cumulative probs
    qsort((void*)prob_values, num_data_qubits * num_data_qubits, sizeof(prob_values[0]), q_comp);

    prob_count = 0;
    float cum_prob = 0.0;
    for (int i = 0; i < num_data_qubits; i++) {
        for (int j = 0; j < num_data_qubits; j++) {
            cum_prob += (float)(prob_values[prob_count].probability * prob_values[prob_count].count);
            prob_values[prob_count].cumulative = cum_prob;
            printf("Code with %d errors (%d x and %d z): prob: %e, count: %lu, cumulative prob: %e\n", prob_values[prob_count].n_total_err, prob_values[prob_count].n_x_err, prob_values[prob_count].n_z_err, prob_values[prob_count].probability, prob_values[prob_count].count, prob_values[prob_count].cumulative);
            prob_count++;
        } 
    }

    // TODO: Determine which samples will be written to file here:
    //      - Sample from the above struct array to find out how many samples to take for each probability
    //      - For each probability, randomly select a configuration of the respective number of x and z errors          (Is this where I use the Fisher Yates shuffle, or is it on the prior line?)
    //      - Write that configuration to a bit string that can be used in lieu of the 'i' variable for the first nested for loop to function properly
    //      - Alternatively, instead of writing to the data qubit arrays using the existing method, write to them using the names of the selected qubits
    //      - Just ask Dr. Merkle. I'm not entirely sure what would work best here. I think just writing out the bit string to replace 'i' would be the best, but it would also require me to decompose the problem differently when I rebase the code to run in parallel
    int sample_count = 0;

    /* Alternative plan: 
     *  - Determine how many samples from each probability range I want to generate (again, ask Dr. Merkle about this)
     *  - For each probability, randomly select where the corresponding number of x and z qubits will end up on the surface code
     *  - Generate bit string with the surface code arrangement in the previous step where every even bit contains an x error and every odd bit contains a z error ( O(n) )
     *  - Add this to a 2d array of size [num_samples][2]
     *  - Decompose the list into processes and assign each entry a process (block decomposition)
     *  - Once in the while loop, iterate over the bit strings belonging to that process (exit the loop if the bit string is not assigned to this process)
     *  - Everything else runs as normal, but need to add communications to send samples from each process to the root process
     */

    int block_size;

    // Begin main sim outer loop
	while ( sample_count < block_size ) {

        // TODO: Figure out new errors value based on which sample to use
		int errors = i;
		double probability = 1.0;

        // Calculate number of errors and determine number of errors present
        unsigned int num_x_error_bits = 0;
        unsigned int num_z_error_bits = 0;
        unsigned int error_bit_string = errors;

        // TODO: use this to determine how many samples to take from iterations that have a certain number of errors
        // Alternatively, write the error bitstream prior to running the outer loop and change to a while loop to reduce runtime
        while(error_bit_string) {
            num_x_error_bits += error_bit_string & 1;   // Every even bit will contain possible x error
            error_bit_string >>= 1;
            num_z_error_bits += error_bit_string & 1;   // Every odd bit will contain possible z error
            error_bit_string >>= 1;
        }

		// initialize data_qubit_x_error and data_qubit_z_error to be full of FALSE (or 0) values
		memset(data_qubit_x_error, 0, sizeof(data_qubit_x_error));
		memset(data_qubit_z_error, 0, sizeof(data_qubit_z_error));
		memset(ancilla_qubit_value, 0, sizeof(ancilla_qubit_value));

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

		for ( int j = 0; j < pow( (depth + 1 ), 2) - 1; j++ ) {
            if ((j/(depth+1))%2 == 0) {
                ancilla_qubit_value[ j / (depth+1) ][ j % (depth+1) ] = 
                    data_qubit_x_error[   j / (depth+1)       ][   j % (depth+1)       ]
                    ^ data_qubit_x_error[   j / (depth+1)       ][ ( j % (depth+1) ) + 1 ]
                    ^ data_qubit_x_error[ ( j / (depth+1) ) + 1 ][   j % (depth+1)       ]
                    ^ data_qubit_x_error[ ( j / (depth+1) ) + 1 ][ ( j % (depth+1) ) + 1 ];
                j++;
                ancilla_qubit_value[ j / (depth+1) ][ j % (depth+1) ] = 
                    data_qubit_z_error[   j / (depth+1)       ][   j % (depth+1)       ]
                    ^ data_qubit_z_error[   j / (depth+1)       ][ ( j % (depth+1) ) + 1 ]
                    ^ data_qubit_z_error[ ( j / (depth+1) ) + 1 ][   j % (depth+1)       ]
                    ^ data_qubit_z_error[ ( j / (depth+1) ) + 1 ][ ( j % (depth+1) ) + 1 ];
            } else {
                ancilla_qubit_value[ j / (depth+1) ][ j % (depth+1) ] = 
                    data_qubit_z_error[   j / (depth+1)       ][   j % (depth+1)       ]
                    ^ data_qubit_z_error[   j / (depth+1)       ][ ( j % (depth+1) ) + 1 ]
                    ^ data_qubit_z_error[ ( j / (depth+1) ) + 1 ][   j % (depth+1)       ]
                    ^ data_qubit_z_error[ ( j / (depth+1) ) + 1 ][ ( j % (depth+1) ) + 1 ];
                j++;
                ancilla_qubit_value[ j / (depth+1) ][ j % (depth+1) ] = 
                    data_qubit_x_error[   j / (depth+1)       ][   j % (depth+1)       ]
                    ^ data_qubit_x_error[   j / (depth+1)       ][ ( j % (depth+1) ) + 1 ]
                    ^ data_qubit_x_error[ ( j / (depth+1) ) + 1 ][   j % (depth+1)       ]
                    ^ data_qubit_x_error[ ( j / (depth+1) ) + 1 ][ ( j % (depth+1) ) + 1 ];
            }
		}

        // Write out probability
        //fprintf(f, "%f,", probability);
        //probList[i][0] = i;
        //probList[i][1] = probability;

        // Output the ancilla qubit values in proper format
        int ancilla_value;
        for (int k=1; k < depth; k++) {
            if (k%2 == 0) {
                ancilla_value = (ancilla_qubit_value[depth][k] == 1) ? -1 : 1;
                fprintf(f, "%d,", ancilla_value);
            }
            for (int j=depth-1; j > 0; j--) {
                if (k == 1 && j%2 == 0) {
                    ancilla_value = (ancilla_qubit_value[j][k-1] == 1) ? -1 : 1;
                    fprintf(f, "%d,", ancilla_value);
                } else if (k == (depth - 1) && j%2 == 1) {
                    ancilla_value = (ancilla_qubit_value[j][k+1] == 1) ? -1 : 1;
                    fprintf(f, "%d,", ancilla_value);
                }
                ancilla_value = (ancilla_qubit_value[j][k] == 1) ? -1 : 1;
                fprintf(f, "%d,", ancilla_value);
            }
            if (k%2 == 1) {
                ancilla_value = (ancilla_qubit_value[0][k] == 1) ? -1 : 1;
                fprintf(f, "%d,", ancilla_value);
            }
        }

        // For printing label list:
        int max_char = 14 * (depth * depth) + 3;
        char label_list[max_char];
        strcpy(label_list, "\"[");
        char qubit_name[6];
        int first = 1;

		for (int k = 1; k < depth + 1; k++) {
            for (int j = depth; j > 0; j--) {
                if (data_qubit_x_error[j][k] == 1) {
                    if (first == 1) {
                        first = 0;
                    } else {
                        strcat(label_list, ", ");
                    }
                    sprintf(qubit_name, "'X%d%d'", (k-1), (3-j));
                    strcat(label_list, qubit_name);
                }
                if (data_qubit_z_error[j][k] == 1) {
                    if (first == 1) {
                        first = 0;
                    } else {
                        strcat(label_list, ", ");
                    }
                    sprintf(qubit_name, "'Z%d%d'", (k-1), (3-j));
                    strcat(label_list, qubit_name);
                }
            }
        }
        strcat(label_list, "]\"");
        fprintf(f, "%s\n", label_list);

	}

	return 0;
}