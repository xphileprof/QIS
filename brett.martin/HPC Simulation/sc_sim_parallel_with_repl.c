/*	2d Lt Brett Martin
 *	Advisor: Dr. Laurence Merkle
 *	Surface code simulation - parallel random-sampled version without replacement
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <mpi.h>

#define TRI_SIZE 50


// Struct used in calculating and sorting all possible probabilities along with the number of samples per probability and the respective number of x and z errors
struct s_prob {
    int n_x_err, n_z_err, n_total_err;
    float probability, cumulative;
    unsigned long count, collected;
};


// Struct used to keep track of the x and z error syndrome for a given sample in the form of bit strings
struct error_syndrome {
    uint64_t x_syndrome, z_syndrome;
};


// Helper function passed into the quicksort algorithm to sort the struct s_prob by the probability values
int q_comp( const void *p1, const void *p2 ) {
    float a = ( (struct s_prob *) p1 )->probability;
    float b = ( (struct s_prob *) p2 )->probability;
    if (a < b) { return -1; }
    else if (a > b) { return 1; }
    return 0;
}


// Helper function to the Fisher Yates Shuffle function, swaps values at two indices of an array
void swap( int *first, int *second) {
    int buf = *first;
    *first = *second;
    *second = buf;
}


// Fisher Yates Shuffle for randomizing an array of indices
void fisher_yates( int array[], int size) {
    for (int i = size - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        swap(&array[i], &array[j]);
    }
}


// Begin main function
int main ( int argc, char** argv ) {

    // MPI vars (NOTE: Per the Navy DSRC intro guide, the standard queue has a maximum cores per job of 8,168 cores)
    int iproc, nproc;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // Execution time variables:
    double exec_time = 0.0;
    clock_t begin;
    if (iproc == 0) {
        begin = clock();
    }

    //Inputs:
    const int depth = 7;		// Or code distance
    const int num_samples = 1000;
    const float p_x_error = 0.05;
    const float p_z_error = 0.05;

    int max_char[] = {128, 640, 328};
    int i_max_char = max_char[(depth%3)];

    int block_size = floor(num_samples/nproc);
    if (num_samples%nproc > 0) { block_size++; }

    int iter_first = iproc * block_size;
    int iter_last = iter_first + block_size;

    MPI_Status status;
    MPI_File fh;

    const char* buf[200];

    // Set this to 1 if you want the sample set to include instances where no errors exist
    int include_no_err = 0;

    //Output:
    //list of records in which each possible combination of data qubit errors is associated with the resulting ancilla qubit values (may be probabilistic), along with the probability of the combination of data qubit errors (and measurement errors?)
    //FILE *f = fopen("D:/Documents/AFIT/Research/REPO/QIS/brett.martin/Neural Network/TEST/v3samples-d5-1000.csv", "w+");
    MPI_File_open(MPI_COMM_SELF, "testfile.csv", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, iproc * block_size * sizeof(buf), MPI_CHAR, MPI_CHAR, 'native', MPI_INFO_NULL);

    // Generate time-dependent seed for random number generation
    srand(time(NULL));

    // Variables used in the generation of sample data
    int num_data_qubits = depth * depth;
    int data_qubit_x_error[ depth + 2][ depth + 2 ];
    int data_qubit_z_error[ depth + 2  ][ depth + 2 ];
    int ancilla_qubit_value[ depth + 1 ][ depth + 1 ];

    // Generate pascal's triangle for use in calculating binomial coefficients ( time O(n^2), where n is the triangle size )
    // When referencing (n k) or "n choose k", just use: pascal_triangle[n][k]
    const int tri_size = 50;
    unsigned long pascal_triangle[ TRI_SIZE ][ TRI_SIZE ];

    for (int row = 0; row < tri_size; row++) {  // This generates a ton of empty space, but it should be okay with such negligible size
        for (int col = 0; col <= row; col++)
        {
            if (row == col || col == 0) { pascal_triangle[row][col] = (unsigned long) 1; }
            else { pascal_triangle[row][col] = (unsigned long) pascal_triangle[row-1][col] + pascal_triangle[row-1][col-1]; }
        }
    }

    // Create CSV headers
    if(iproc == 0) {
        char headers[i_max_char];
        // Create CSV headers
        char h_anc[4];
        int ancilla = 0;
        for (int k=1; k < depth; k++) {
            if (k%2 == 0) {
                sprintf(h_anc, "Z%d,", ancilla);
                strcat(headers, h_anc);
                ancilla++;
            }
            for (int j=depth-1; j > 0; j--) {
                if (k == 1 && j%2 == 0) {
                    sprintf(h_anc, "X%d,", ancilla);
                    strcat(headers, h_anc);
                    ancilla++;
                } else if (k == (depth - 1) && j%2 == 1) {
                    sprintf(h_anc, "X%d,", ancilla);
                    strcat(headers, h_anc);
                    ancilla++;
                }
                if ((j+k)%2 == 0) {
                    sprintf(h_anc, "X%d,", ancilla);
                    strcat(headers, h_anc);
                    ancilla++;
                } else {
                    sprintf(h_anc, "Z%d,", ancilla);
                    strcat(headers, h_anc);
                    ancilla++;
                }
            }
            if (k%2 == 1) {
                sprintf(h_anc, "Z%d,", ancilla);
                strcat(headers, h_anc);
                ancilla++;
            }
        }
        strcat(headers, "Labels");
        MPI_File_write(fh, headers, strlen(headers) * sizeof(char), MPI_CHAR, MPI_STATUS_IGNORE);
    }

    // Determine the number of unique probabilities, the number of samples generated per probability
    struct s_prob prob_values[(num_data_qubits + 1) * (num_data_qubits + 1)];
    memset(prob_values, 0, sizeof(prob_values));

    int prob_count = 0;
    for (int i = 0; i < num_data_qubits + 1; i++) {
        for (int j = 0; j < num_data_qubits + 1; j++) {
            prob_values[prob_count].n_total_err = i + j;
            prob_values[prob_count].n_x_err = i;
            prob_values[prob_count].n_z_err = j;
            prob_values[prob_count].probability = (float) pow((1.0 - p_x_error), num_data_qubits - i) * pow((1.0 - p_z_error), num_data_qubits - j) * pow(p_x_error, i) * pow(p_z_error, j);
            prob_values[prob_count].count = (unsigned long) pascal_triangle[num_data_qubits][i] * pascal_triangle[num_data_qubits][j];
            prob_count++;
        }
    }

    // Sort by ascending probability and calculate cumulative probs
    qsort((void*)prob_values, (num_data_qubits + 1) * (num_data_qubits + 1), sizeof(prob_values[0]), q_comp);

    prob_count = 0;
    float cum_prob = 0.0;
    for (int i = 0; i < num_data_qubits + 1; i++) {
        for (int j = 0; j < num_data_qubits + 1; j++) {
            cum_prob += (float)(prob_values[prob_count].probability * prob_values[prob_count].count);
            prob_values[prob_count].cumulative = cum_prob;
            prob_count++;
        }
    }

    // If all samples must contain at least one error, lower upper bound of random function to exclude that probability range
    float upper_prob_bound = include_no_err ? 1.0 : (prob_values[(num_data_qubits + 1) * (num_data_qubits + 1) - 2].cumulative - 0.01);

    // Begin main sim outer loop
    for (int i = 0; i < num_samples; i++) {

        // initialize data_qubit_x_error and data_qubit_z_error (and others) to be full of FALSE (or 0) values
        memset(data_qubit_x_error, 0, sizeof(data_qubit_x_error));
        memset(data_qubit_z_error, 0, sizeof(data_qubit_z_error));
        memset(ancilla_qubit_value, 0, sizeof(ancilla_qubit_value));

        // Generate random number between 0 and 1 (or upper bound)
        float rand_num = ((float) rand() / RAND_MAX) * upper_prob_bound;
        int sample_idx = 0;
        for (int j = 0; j < (num_data_qubits + 1) * (num_data_qubits + 1); j++) {
            // Index of that element is going to point to the entry in the struct array that will tell me how many of each error will occur, but only if all samples have not been generated
            if (prob_values[j].cumulative >= rand_num) {     // Check if the cumulative value is higher AND if there are enough samples to continue
                sample_idx = j;
                break;
            }
        }

        // Initialize array of length (num_data_qubits) where each element in the array is equal to its index
        int rand_array[num_data_qubits];
        for (int j = 0; j < num_data_qubits; j++) { rand_array[j] = j; }

        // Perform n_x_err and n_z_err iterations of the Fisher Yates Shuffle on the array (where the selected value is at index 0)
        // This will tell me the data qubits that have x errors
        int qubit_idx;
        for (int j = 0; j < prob_values[sample_idx].n_x_err; j++) {
            while (1) {    // Ensure the same error doesn't get written twice
                fisher_yates(rand_array, num_data_qubits);
                qubit_idx = rand_array[0];
                if (data_qubit_x_error[ (qubit_idx/depth) + 1 ][ (qubit_idx%depth) + 1 ] != 1) { break; }
            }
            data_qubit_x_error[ (qubit_idx/depth) + 1 ][ (qubit_idx%depth) + 1 ] = 1;
        }

        for (int j = 0; j < prob_values[sample_idx].n_z_err; j++) {
            while (1) {
                fisher_yates(rand_array, num_data_qubits);
                qubit_idx = rand_array[0];
                if (data_qubit_z_error[ (qubit_idx/depth) + 1 ][ (qubit_idx%depth) + 1 ] != 1) { break; }
            }
            data_qubit_z_error[ (qubit_idx/depth) + 1 ][ (qubit_idx%depth) + 1 ] = 1;
        }

        // Calculate ancilla qubit values
        for ( int j = 0; j < (depth + 1 ) * (depth + 1 ) - 1; j++ ) {
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

        char label_list[i_max_char];
        strcpy(label_list, "\n");
        char anc_name[3];

        // Output the ancilla qubit values in proper format
        int ancilla_value;
        for (int k=1; k < depth; k++) {
            if (k%2 == 0) {
                ancilla_value = (ancilla_qubit_value[depth][k] == 1) ? -1 : 1;
                sprintf(anc_name, "%d,", ancilla_value);
                strcat(label_list, anc_name);
            }
            for (int j=depth-1; j > 0; j--) {
                if (k == 1 && j%2 == 0) {
                    ancilla_value = (ancilla_qubit_value[j][k-1] == 1) ? -1 : 1;
                    sprintf(anc_name, "%d,", ancilla_value);
                    strcat(label_list, anc_name);
                } else if (k == (depth - 1) && j%2 == 1) {
                    ancilla_value = (ancilla_qubit_value[j][k+1] == 1) ? -1 : 1;
                    sprintf(anc_name, "%d,", ancilla_value);
                    strcat(label_list, anc_name);
                }
                ancilla_value = (ancilla_qubit_value[j][k] == 1) ? -1 : 1;
                sprintf(anc_name, "%d,", ancilla_value);
                strcat(label_list, anc_name);
            }
            if (k%2 == 1) {
                ancilla_value = (ancilla_qubit_value[0][k] == 1) ? -1 : 1;
                sprintf(anc_name, "%d,", ancilla_value);
                strcat(label_list, anc_name);
            }
        }

        // For printing label list:
        strcat(label_list, "\"[");
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
                    sprintf(qubit_name, "'X%d%d'", (k-1), (depth-j));
                    strcat(label_list, qubit_name);
                }
                if (data_qubit_z_error[j][k] == 1) {
                    if (first == 1) {
                        first = 0;
                    } else {
                        strcat(label_list, ", ");
                    }
                    sprintf(qubit_name, "'Z%d%d'", (k-1), (depth-j));
                    strcat(label_list, qubit_name);
                }
            }
        }
        strcat(label_list, "]\"");

        MPI_File_write(fh, label_list, strlen(label_list) * sizeof(char), MPI_CHAR, MPI_STATUS_IGNORE);

    }

    MPI_File_close(&fh);

    if (iproc == 0) {
        clock_t end = clock();
        exec_time += (double) (end - begin) / CLOCKS_PER_SEC;
        printf("Samples collected for %d, %d in %f seconds\n", depth, num_samples, exec_time);
    }

    MPI_Finalize();

    return 0;

}