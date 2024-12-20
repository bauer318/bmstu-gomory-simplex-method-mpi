#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_ITER 100
#define EPSILON 1e-6
#define MASTER 0
double* print_solution_and_get(double** tableau, double* basics, int rows, int cols, int need_extends_basic);

void print_matrix(double** matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%8.3f ", matrix[i][j]);
        }
        printf("\n");
    }
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

int find_pivot_col(double** tableau, int rows, int cols, MPI_Comm comm) {
    int rank, size;
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


int find_pivot_row(double** tableau, int rows, int cols, int pivot_col, MPI_Comm comm) {
    int rank, size;
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

void add_gomory_cut_m(double*** tableau, int* rows, int* cols, int gomory_row_to_cut, int is_first_time) {

    int new_rows = *rows + 1;
    int new_cols = *cols + 1;
    int gomory_column = *cols - 1;
    int gomory_row = *rows - 1;

    double** new_tableau = (double**)malloc(new_rows * sizeof(double*));
    for (int i = 0; i < new_rows; i++) {
        new_tableau[i] = (double*)malloc(new_cols * sizeof(double));
    }


    for (int i = 0; i < *rows - 1; i++) {
        for (int j = 0; j < *cols; j++) {
            if (j == gomory_column) {
                new_tableau[i][j] = 0.0;
            }
            else {
                new_tableau[i][j] = (*tableau)[i][j];
            }

        }
    }

    for (int j = 0; j < *cols; j++) {
        double value = (*tableau)[gomory_row_to_cut][j];
        double fractionalPart = value - floor(value);
        if (j == gomory_column) {
            new_tableau[gomory_row][j] = 1;
        }
        else {
            new_tableau[gomory_row][j] = fractionalPart != 0 ? -1 * fractionalPart : fractionalPart;
        }
    }

    for (int j = 0; j < *cols; j++) {
        double value = (*tableau)[*rows - 1][j];
        if (j == gomory_column) {
            new_tableau[*rows][j] = 0.0;
        }
        else {
            if (is_first_time) {

                new_tableau[*rows][j] = value != 0 ? -1 * value : value;
            }
            else {
                new_tableau[*rows][j] = value;
            }
        }
    }

    for (int i = 0; i < new_rows; i++) {
        if (i == gomory_row) {
            double value = (*tableau)[gomory_row_to_cut][gomory_column];
            double fractionalPart = value - floor(value);
            new_tableau[i][*cols] = fractionalPart != 0 ? -1 * fractionalPart : fractionalPart;
        }
        else {
            if (i == new_rows - 1 && is_first_time) {
                new_tableau[i][*cols] = -1 * (*tableau)[gomory_row][gomory_column];
            }
            else {
                new_tableau[i][*cols] = (*tableau)[i][gomory_column];
            }

        }
    }

    *tableau = new_tableau;
    *rows = new_rows;
    *cols = new_cols;

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

void apply_gomory_cuts(double** tableau, int rows, int cols, double* basics) {

    int keep_go = exist_real_value(tableau, rows, cols);
    int is_first_time = 1;
    int index = 0;

    while (keep_go) {
        int row_to_cut = find_gomory_row_to_cut(tableau, rows, cols);
        if (row_to_cut == -1) {
            printf("All solutions are integers.\n");
            break;
        }
        printf("Adding Gomory cut for row %d\n", row_to_cut);
        //add_gomory_cut_m(&tableau, &rows, &cols, row_to_cut, is_first_time);

        print_matrix(tableau, rows, cols);
        printf("-------------------------------------------------------------------\n");
        int gomory_row = rows - 2;
        int gomory_col = find_gomory_column_to_add(tableau, rows, cols);

        basics = extend_basics(basics, rows - 2, gomory_col);

        pivot(tableau, rows, cols, gomory_row, gomory_col);
        print_matrix(tableau, rows, cols);

        //double* solution = print_solution_and_get(tableau, basics, rows, cols, 1);
        keep_go = index <= 2;
        index++;
        is_first_time = 0;


    }
}


// Perform Simplex method on the tableau
int simplex_method(double** tableau, int rows, int cols, MPI_Comm comm) {
    double* basics = (double*)malloc((rows - 1) * sizeof(double));
    int rank;
    MPI_Comm_rank(comm, &rank);
    while (1) {
        // Check for optimality
        int pivot_col = find_pivot_col(tableau, rows, cols, comm);
        if (pivot_col == -1) {
            // Optimal solution found
            if (rank == MASTER) {
                double* solution = print_solution_and_get(tableau, basics, rows, cols, 0);
                apply_gomory_cuts(tableau, rows, cols, basics);
            }

            return 1;
        }

        // Find pivot row
        int pivot_row = find_pivot_row(tableau, rows, cols, pivot_col, comm);
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

double* print_solution_and_get(double** tableau, double* basics, int rows, int cols, int need_extends_basic) {
    printf("Solution:\n");

    int basicsLength = need_extends_basic ? rows + 1 : rows - 1;
    double* solution = (double*)calloc(basicsLength, sizeof(double));
    if (solution == NULL) {
        fprintf(stderr, "Memory allocation failed for solution.\n");
        return NULL;
    }

    for (int i = 0; i < basicsLength; i++) {
        int xCol = (int)basics[i];
        if (xCol >= 0 && xCol < cols) {
            solution[xCol] = tableau[i][cols - 1];
        }
        else {
            fprintf(stderr, "Warning: Invalid column index basics[%d] = %d\n", i, xCol);
        }
    }

    for (int i = 0; i < basicsLength; i++) {
        printf("x%d = %.2f\n", i + 1, solution[i]);
    }

    return solution;
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
    int optimal = simplex_method(tableau, rows, cols, MPI_COMM_WORLD);

    if (optimal == 1 && rank == MASTER) {
        print_matrix(tableau, rows, cols);
    }

    free_matrix(tableau, rows);
    free(flat_tableau);
    MPI_Finalize();
    return 0;
}
