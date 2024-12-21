#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_ITER 100
#define EPSILON 1e-6
#define MASTER 0

void print_matrix(double** matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%8.3f ", matrix[i][j]);
        }
        printf("\n");
    }
}

double* print_solution_and_get(double** tableau, double* basics, int rows, int cols) {
    printf("Solution:\n");

    int basics_length = rows - 1;
    double* solution = (double*)calloc(basics_length, sizeof(double));
    if (solution == NULL) {
        fprintf(stderr, "Memory allocation failed for solution.\n");
        return NULL;
    }

    for (int i = 0; i < basics_length; i++) {
        int xCol = (int)basics[i];

        if (xCol >= 0 && xCol < cols) {
            solution[xCol] = tableau[i][cols - 1];
        }
        else {
            fprintf(stderr, "Warning: Invalid column index basics[%d] = %d\n", i, xCol);
        }
    }

    for (int i = 0; i < basics_length; i++) {
        printf("x%d = %.2f\n", i + 1, solution[i]);
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
    int rank, size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int pivotCol_local = -1;
    double mostNegative_local = 0;

    // Determine the range of columns this process will handle
    for (int col = rank; col < cols - 1; col += size) {
        if (tableau[rows - 1][col] < mostNegative_local) {
            mostNegative_local = tableau[rows - 1][col];
            pivotCol_local = col;
        }
    }

    // Prepare for the reduction operation
    struct {
        double value;
        int index;
    } local_result, global_result;

    local_result.value = mostNegative_local;
    local_result.index = pivotCol_local;

    // Reduce to find the most negative value and its index across all processes
    MPI_Allreduce(&local_result, &global_result, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);

    // Return the global pivot column index
    return global_result.index;
}


int is_integer(double value) {
    return fabs(value - round(value)) < 1e-6;
}

int exist_real_value(double** tableau, int rows, int cols) {
    int exist_real_value = 0;
    for (int i = 0; i < rows - 1; i++) {
        if (!is_integer(tableau[i][cols - 1])) {
            exist_real_value = 1;
            break;
        }
    }
    return exist_real_value;
}


int find_pivot_row(double** tableau, int rows, int cols, int pivot_col) {
    int rank, size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int pivotRow_local = -1;
    double minRatio_local = INFINITY;

    // Each process works on its subset of rows
    for (int row = rank; row < rows - 1; row += size) {
        if (tableau[row][pivot_col] > 0) {
            double ratio = tableau[row][cols - 1] / tableau[row][pivot_col];
            if (ratio < minRatio_local && ratio > 0) {
                minRatio_local = ratio;
                pivotRow_local = row;
            }
        }
    }

    // Prepare for the reduction operation
    struct {
        double value;
        int index;
    } local_result, global_result;

    local_result.value = minRatio_local;
    local_result.index = pivotRow_local;

    // Reduce to find the minimum ratio and its corresponding row across all processes
    MPI_Allreduce(&local_result, &global_result, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);

    // Return the global pivot row index
    return global_result.index;
}

int find_gomory_row_to_cut(double** tableau, int rows, int cols) {
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

void pivot(double** tableau, int rows, int cols, int pivot_row, int pivot_col) {
    double pivot_value = tableau[pivot_row][pivot_col];

    for (int j = 0; j < cols; j++) {
        tableau[pivot_row][j] /= pivot_value;
    }
    for (int i = 0; i < rows; i++) {
        if (i != pivot_row) {
            double factor = tableau[i][pivot_col];
            for (int j = 0; j < cols; j++) {
                tableau[i][j] -= factor * tableau[pivot_row][j];
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

	for (int i = 0; i < gomory_row; i++) {
		for (int j = 0; j < old_cols; j++) {
            if (j == gomory_col) {
                result[i][j] = 0.0;
            }
            else {
                result[i][j] = tableau[i][j];
            }
           
		}
	}

    for (int j = 0; j < old_cols; j++) {
        double value = tableau[row_to_cut][j];
        double fractionalPart = value - floor(value);
        if (j == gomory_col) {
            result[gomory_row][j] = 1;
        }
        else {
            result[gomory_row][j] = fractionalPart != 0 ? -1 * fractionalPart : fractionalPart;
        }
    }


    for (int j = 0; j <old_cols; j++) {
        double value = tableau[gomory_row][j];
        if (j == gomory_col) {
            result[old_rows][j] = 0.0;
        }
        else {
            if (is_first_time) {

                result[old_rows][j] = value != 0 ? -1 * value : value;
            }
            else {
                result[old_rows][j] = value;
            }
        }
    }

    //The last column of the extended tableau
    for (int i = 0; i < new_rows; i++) {
        if (i == gomory_row) {
            double value = tableau[row_to_cut][gomory_col];
            double fractionalPart = value - floor(value);
            result[i][old_cols] = fractionalPart != 0 ? -1 * fractionalPart : fractionalPart;
        }
        else {
            if (i < old_rows) {
                result[i][old_cols] = tableau[i][gomory_col];
            }
            else {
                result[i][old_cols] = is_first_time ? -1 * tableau[old_rows - 1][gomory_col] : tableau[old_rows - 1][gomory_col];
            }
        }
        
    }

    return result;

}
void apply_gomory_cuts(double** tableau, int rows, int cols, double* basics) {
    int keep_apply_gomory_cut = exist_real_value(tableau, rows, cols);
    int is_first_time = 1;

    if (keep_apply_gomory_cut) {
        printf("\nApply Gomory\n");
    }
   
    while (keep_apply_gomory_cut) {
        int row_to_cut = find_gomory_row_to_cut(tableau, rows, cols);
        if (row_to_cut == -1) {
            printf("All solutions are integers.\n");
            break;
        }
        printf("Adding Gomory cut for row %d\n", row_to_cut);
        tableau = add_gomory_cut(tableau, rows, cols,row_to_cut, is_first_time);
        rows++;
        cols++;
        print_matrix(tableau, rows, cols);
      
        int gomory_row = rows - 2;
        int gomory_col = find_gomory_column_to_add(tableau, rows, cols);

        basics = extend_basics(basics, rows - 2, gomory_col);

        pivot(tableau, rows, cols, gomory_row, gomory_col);
       
        double* solution = print_solution_and_get(tableau, basics, rows, cols);
        printf("-------------------------------------------------------------------\n");
        keep_apply_gomory_cut = exist_real_value(tableau, rows, cols);
        is_first_time = 0;
    }
    print_matrix(tableau, rows, cols);
    printf("-------------------------------------------------------------------\n");
}


// Perform Simplex method on the tableau
int simplex_method(double** tableau, int rows, int cols) {
    double* basics = (double*)malloc((rows - 1) * sizeof(double));
    int rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
  
    int basic_was_initialized = 0;
    while (1) {
        // Check for optimality
        int pivot_col = find_pivot_col(tableau, rows, cols);
        if (pivot_col == -1) {
            // Optimal solution found
            double* solution = NULL;
            if (rank == MASTER) {
                printf("Optimal solution found\n");
                print_matrix(tableau, rows, cols);
                solution = print_solution_and_get(tableau, basics, rows, cols);
                apply_gomory_cuts(tableau, rows, cols, basics);
                basic_was_initialized = 1;
            }
            
            return 1;
        }

        // Find pivot row
        int pivot_row = find_pivot_row(tableau, rows, cols, pivot_col);
        double min_ratio = INFINITY;

        if (pivot_row == -1) {
            // Unbounded solution
            return 0;
        }

        // Perform pivot operation
        pivot(tableau, rows, cols, pivot_row, pivot_col);

        basics[pivot_row] = pivot_col;
    }
}


int main(int argc, char* argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Define problem size and tableau (Master process initializes tableau)
    int rows = 3, cols = 5;
    double** tableau = NULL;
    double* flat_tableau = NULL;

    if (rank == MASTER) {

        // Maximize Z = 7x1 + 9x2
        // Subject to: 
        // -x1 + 3x2 <= 6
        // 7x1 +  x2 <= 35
        // Tableau format:
        // | Coefficients | RHS |

        tableau = allocate_matrix(rows, cols);
        double init_tableau[3][5] = {
            {-1, 3, 1, 0, 6},
            {7, 1, 0, 1, 35},
            {-7, -9, 0, 0, 0}
        };

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                tableau[i][j] = init_tableau[i][j];
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
        print_matrix(tableau, rows, cols);
    }

    // Perform Simplex
    int optimal = simplex_method(tableau, rows, cols);

    free_matrix(tableau, rows);
    free(flat_tableau);
    MPI_Finalize();
    return 0;
}
