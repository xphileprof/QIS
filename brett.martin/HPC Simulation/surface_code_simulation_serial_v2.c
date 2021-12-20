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

struct s_prob {
    int n_x_err, n_z_err, n_total_err;
    float probability, cumulative;
    unsigned long count;    // TODO: Change to double to have higher range
};

int q_comp(const void *p1, const void *p2) {
    float a = ((struct s_prob *)p1)->probability;
    float b = ((struct s_prob *)p2)->probability;
    if (a < b) { return -1; }
    else if (a > b) { return 1; }
    return 0;
}

int main (int argc, char** argv) {

	//Inputs:
	const int depth = 3;		// Or code distance
    const int num_samples = 1000;
	const float p_x_error = 0.05;
	const float p_z_error = 0.05;

	//Output:
	//list of records in which each possible combination of data qubit errors is associated with the resulting ancilla qubit values (may be probabilistic), along with the probability of the combination of data qubit errors (and measurement errors?)
	FILE *f = fopen("D:/Documents/AFIT/Research/test.csv", "w+");

    // Generate pascal's triangle for use in calculating binomial coefficients: O(n^2)
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
    //fprintf(f, "Probability,");
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

	const int num_data_qubits = depth * depth;
	int data_qubit_x_error[ depth + 2][ depth + 2 ];
	int data_qubit_z_error[ depth + 2  ][ depth + 2 ];
	int ancilla_qubit_value[ depth + 1 ][ depth + 1 ];

	// Parallel block decomposition variables
	uint64_t total_num_iter = pow(4, num_data_qubits);			// Total number of outer loop iterations

    // Determine the number of unique probabilities, the number of samples generated per probability
    //float prob_values[(int) pow(num_data_qubits, 2)][2];
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

    // Sort by ascending probability
    // Also, include cumulative probs
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
    

    // Begin main sim outer loop
	for ( int i = 0; i < total_num_iter; i++ ) {

		int errors = i;
		double probability = 1.0;

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