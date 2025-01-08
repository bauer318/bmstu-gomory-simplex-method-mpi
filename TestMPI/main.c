#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_ITER 1000
#define EPSILON 1e-6
#define MASTER 0
#define FILE_NAME "result.txt"
#define COLS 300
#define ALL_ROWS 201
#define ALL_COLS 501
#define DATA_FILE_NAME "array_data.txt"


void print_matrix(double** matrix, int rows, int cols, FILE* file) {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            fprintf(file, "%8.3f ", matrix[i][j]);
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n---------------------------------------------------------------\n");
}

double* print_solution_and_get(double** tableau, double* basics, int rows, int cols, FILE* file) {
    fprintf(file, "Solution:\n");

    int basics_length = rows - 1;
    double* solution = (double*)calloc(basics_length, sizeof(double));
    if (solution == NULL) {
        fprintf(stderr, "Memory allocation failed for solution.\n");
        return NULL;
    }

    for (int i = 0; i < basics_length; i++) {
        solution[i] = tableau[i][cols - 1];
    }

    for (int i = 0; i < basics_length; i++) {
        int x_col = (int)basics[i];
        if (x_col != -1 && x_col < COLS) {
            fprintf(file, "x%d = %8.2f\n", x_col + 1, solution[i]);
        }
       
    }

    return solution;
}

// Allocate memory for a 2D array (matrix)
double** allocate_matrix(int rows, int cols) {
    double** matrix = malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        matrix[i] = malloc(cols * sizeof(double));
    }
    return matrix;
}

// Free memory of a 2D array
void free_matrix(double** matrix, int rows) {
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

int find_pivot_col(double** tableau, int rows, int cols) {
    int pivot_col = -1;
    double most_negative = 0;

    for (int col = 0; col < cols - 1; col++) {
        if (tableau[rows - 1][col] < most_negative) {
            most_negative = tableau[rows - 1][col];
            pivot_col = col;
        }
    }

    return pivot_col;
}


int is_integer(double value) {
    return fabs(value - round(value)) < 1e-6;
}

int exist_real_value(double** tableau, int rows, int cols, double* basics) {
    int exist_real_value = 0;

    for (int i = 0; i < rows - 1; i++) {
        int x_col = (int)basics[i];
        if (x_col < COLS) {
            if (!is_integer(tableau[i][cols - 1])) {
                exist_real_value = 1;
                break;
            }
        }
    }
    return exist_real_value;
}


int find_pivot_row(double** tableau, int rows, int cols, int pivot_col) {  
    int pivot_row = -1;
    double min_ratio = INFINITY;

    for (int row = 0; row < rows - 1; row ++) {
        if (tableau[row][pivot_col] > 0) {
            double ratio = tableau[row][cols - 1] / tableau[row][pivot_col];
            if (ratio < min_ratio && ratio > 0) {
                min_ratio = ratio;
                pivot_row = row;
            }
        }
    }

    return pivot_row;
}

int find_gomory_row_to_cut(double** tableau, int rows, int cols, double* basics) {
    int row_to_cut = -1;
    double max_fractional_part = 0.0;

    for (int i = 0; i < rows - 1; i++) {
        double value = tableau[i][cols - 1];
        double fractional_part = value - floor(value);

        if (fractional_part > max_fractional_part) {
            max_fractional_part = fractional_part;
            row_to_cut = i;
        }
        
    }

    return row_to_cut;
}

int find_gomory_column_to_add(double** tableau, int rows, int cols) {
    int gomory_column_to_add = -1;
    double min_value = INFINITY;

    for (int j = 0; j < cols - 2; j++) {
        double goromy_row_value = tableau[rows - 2][j];
        if (goromy_row_value != 0) {
            double last_row_gomory_row_rapport = tableau[rows - 1][j] / goromy_row_value;
            if (last_row_gomory_row_rapport <= min_value) {
                gomory_column_to_add = j;
                min_value = last_row_gomory_row_rapport;
            }
        }
    }
    return gomory_column_to_add;
}

double* extend_basics(double* basics, int old_cols, double init_value) {

    double* new_basics = (double*)malloc((old_cols + 1) * sizeof(double));
    if (!new_basics) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < old_cols; i++) {
        new_basics[i] = basics[i];
    }

    new_basics[old_cols] = init_value;

    return new_basics;
}

double try_to_convert_to_positive_zero(double current_value) {
    return round(current_value*1e9) == 0.0 ? 0.0 : current_value;
}

void pivot(double** tableau, int rows, int cols, int pivot_row, int pivot_col) {
    double pivot_value = tableau[pivot_row][pivot_col];

    for (int j = 0; j < cols; j++) {
        double value = tableau[pivot_row][j] / pivot_value;
        tableau[pivot_row][j] = try_to_convert_to_positive_zero(value);
    }
    for (int i = 0; i < rows; i++) {
        if (i != pivot_row) {
            double factor = tableau[i][pivot_col];
            for (int j = 0; j < cols; j++) {
                tableau[i][j] -= factor * tableau[pivot_row][j];
                tableau[i][j] = try_to_convert_to_positive_zero(tableau[i][j]);
            }
        }
    }
}
double** add_gomory_cut(double** tableau, int old_rows, int old_cols, int row_to_cut, int is_first_time) {
   
    int new_rows = old_rows + 1;
    int new_cols = old_cols + 1;
    int gomory_row = old_rows - 1;
    int gomory_col = old_cols - 1;

    double** result = allocate_matrix(new_rows, new_cols);

    for (int i = 0; i < new_rows; i++) {
        for (int j = 0; j < new_cols; j++) {
            int old_col = j == old_cols ? j - 1 : j;
            //Gomory row
            if (i == gomory_row) {
      
                if (j == gomory_col) {
                    //Intersection with gomory column
                    result[i][j] = 1;
                }
                else {
                    double value = tableau[row_to_cut][old_col];
                    double fractionalPart = value - floor(value);
                    result[i][j] = fractionalPart != 0 ? -1 * fractionalPart : fractionalPart;
                }
            }
            //Last row
            else if (i == old_rows) {
                if (j == gomory_col) {
                    //Intersection with gomory column 
                    result[i][j] = 0;
                }
                else {
                    int old_row = i - 1;

                    if (is_first_time) {
                        result[i][j] = tableau[old_row][old_col] != 0 ? -1 * tableau[old_row][old_col] : tableau[old_row][old_col];
                    }
                    else {
                        result[i][j] = tableau[old_row][old_col];
                    }
                }
            }
            else {
                if (j == gomory_col) {
                    result[i][j] = 0;
                }
                else {
                  
                    result[i][j] = tableau[i][old_col];
                }
            }
        }
    }
    
    return result;

}

void apply_gomory_cuts(double** tableau, int rows, int cols, double* basics, FILE* file) {
    int keep_apply_gomory_cut = exist_real_value(tableau, rows, cols, basics);
    int is_first_time = 1;
    int index = 0;
    int iteration = 0;
    
	if (keep_apply_gomory_cut) {
        fprintf(file, "\nApply Gomory\n");
	}
	int row_to_cut = -1;
	while (keep_apply_gomory_cut) {

        row_to_cut = find_gomory_row_to_cut(tableau, rows, cols, basics);
        if (row_to_cut == -1) {
            fprintf(file, "All solutions are integers.\n");
            break;
        }
        if (iteration == MAX_ITER) {
            fprintf(file, "\nNot found solution after %d iterations \n", MAX_ITER);
            fprintf(file, "-------------------------------------------------------------------\n");
            
            break;
        }
        fprintf(file, "Adding Gomory cut for row %d\n", row_to_cut);
        iteration++;

		tableau = add_gomory_cut(tableau, rows, cols, row_to_cut, is_first_time);
		rows++;
		cols++;
        int gomory_row = rows - 2;
        int gomory_col = find_gomory_column_to_add(tableau, rows, cols);

        basics = extend_basics(basics, rows - 2, gomory_col);

        pivot(tableau, rows, cols, gomory_row, gomory_col);
       
        fprintf(file, "-------------------------------------------------------------------\n");
        print_matrix(tableau, rows, cols, file);
		
		keep_apply_gomory_cut = exist_real_value(tableau, rows, cols,basics);
		is_first_time = 0;
		index++;
	}
    print_matrix(tableau, rows, cols, file);
    fprintf(file, "-------------------------------------------------------------------\n");
    double* solution = print_solution_and_get(tableau, basics, rows, cols, file);
	
        
}

void init_basic(double* basics, int length) {
    for (int i = 0; i < length; i++) {
        basics[i] = -1;
    }
}

int simplex_method(double** tableau, int rows, int cols,FILE* file) {
    double* basics = (double*)malloc((rows - 1) * sizeof(double));
    init_basic(basics, (rows - 1));
    
    while (1) {
        
        int pivot_col = find_pivot_col(tableau, rows, cols);
        if (pivot_col == -1) {
            
            double* solution = NULL;
            fprintf(file, "Optimal solution found\n");
            print_matrix(tableau, rows, cols, file);
            solution = print_solution_and_get(tableau, basics, rows, cols, file);

            apply_gomory_cuts(tableau, rows, cols, basics,file);
            
            return 1;
        }

        int pivot_row = find_pivot_row(tableau, rows, cols, pivot_col);
        double min_ratio = INFINITY;

        if (pivot_row == -1) {
            return 0;
        }

        pivot(tableau, rows, cols, pivot_row, pivot_col);
        basics[pivot_row] = pivot_col;
    }
}

void clear_file() {
    FILE* file = fopen(FILE_NAME, "w");
    if (file != NULL) {
        fclose(file);
    }
    else {
        fprintf(stderr, "Failed to clear file.\n");
    }
}

int read_array_from_file(const char* filename, double*** array) {
    
    *array = (double**)malloc(ALL_ROWS * sizeof(double*));
    if (*array == NULL) {
        fprintf(stderr, "Error allocating memory for tableau\n");
        return 1; 
    }

    for (int i = 0; i < ALL_ROWS; i++) {
        (*array)[i] = (double*)malloc(ALL_COLS * sizeof(double));
        if ((*array)[i] == NULL) {
            fprintf(stderr, "Error allocating memory for row %d of tableau\n", i);
            return 2;  
        }
    }

    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return 3; 
    }

    for (int i = 0; i < ALL_ROWS; i++) {
        for (int j = 0; j < ALL_COLS; j++) {
            int read_value = fscanf_s(file, "%lf", &(*array)[i][j]);
            if (read_value != 1) {
                fprintf(stderr, "Error reading value at [%d][%d] from file\n", i, j);
                fclose(file);
                return 4;
            }
        }
    }

    fclose(file);
    return 0;
}


int main(int argc, char* argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    FILE* file = fopen(FILE_NAME, "a");
    if (file == NULL) {
        fprintf(stderr, "Error opening file.\n");
        return -1;
    }
    double start_time = MPI_Wtime(); 
    clear_file();
    int rows = ALL_ROWS, cols = ALL_COLS;
    double** tableau = NULL;
    double* flat_tableau = NULL;

    if (rank == MASTER) {
		tableau = allocate_matrix(rows, cols);
        double** array = NULL;
      
        int read_file = read_array_from_file(DATA_FILE_NAME, &array);
        if (read_file != 0) {
            return 0;
        }
       
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                tableau[i][j] = array[i][j];
            }
        }

        flat_tableau = malloc(rows * cols * sizeof(double));
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                flat_tableau[i * cols + j] = tableau[i][j];
            }
        }
    }

    if (rank != MASTER) {
        flat_tableau = malloc(rows * cols * sizeof(double));
    }

    MPI_Bcast(flat_tableau, rows * cols, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    if (rank != MASTER) {
        tableau = allocate_matrix(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                tableau[i][j] = flat_tableau[i * cols + j];
            }
        }
    }

    if (rank == MASTER) {
        fprintf(file, "Initial\n");
        print_matrix(tableau, rows, cols,file);
    }

    int optimal = simplex_method(tableau, rows, cols,file);

    free_matrix(tableau, rows);
    free(flat_tableau);

    double end_time = MPI_Wtime();

    if (rank == MASTER) {
        fprintf(file,"Execution Time: %f seconds\n", end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}
