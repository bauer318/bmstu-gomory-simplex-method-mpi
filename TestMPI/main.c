#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <Windows.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define MASTER 0
#define OUTPUT_FILE_NAME "matrix_output.txt"
#define COLS 300
#define ALL_ROWS 201
#define ALL_COLS 501
#define DATA_FILE_NAME "array_data.txt"

void clear_file() {
    FILE* file = fopen(OUTPUT_FILE_NAME, "w");
    if (file != NULL) {
        fclose(file);
    }
    else {
        fprintf(stderr, "Failed to clear file.\n");
    }
}

void write_matrix_ordered(double* local_matrix, int rows_per_process, int cols, int rank, int size, const char* title) {
    MPI_File fh;
    MPI_Status status;
    int extra_place = 200;

    int rc = MPI_File_open(MPI_COMM_WORLD, OUTPUT_FILE_NAME, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND, MPI_INFO_NULL, &fh);
    if (rc != MPI_SUCCESS) {
        fprintf(stderr, "Error opening file by process %d\n", rank);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    char* buffer = (char*)malloc(rows_per_process * cols * extra_place * sizeof(char));
    if (buffer == NULL) {
        fprintf(stderr, "Memory allocation failed on process %d\n", rank);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_MEM);
    }

    int pos = 0;

    if (rank == 0) {
        pos += snprintf(buffer + pos, rows_per_process * cols * extra_place - pos, "%s\n", title);
    }

    for (int i = 0; i < rows_per_process; i++) {
        for (int j = 0; j < cols; j++) {
            if (pos >= rows_per_process * cols * extra_place) {
                fprintf(stderr, "Buffer overflow detected on process %d\n", rank);
                MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
            }
            pos += snprintf(buffer + pos, rows_per_process * cols * extra_place - pos, "%8.3f ", local_matrix[i * cols + j]);
        }
        pos += snprintf(buffer + pos, rows_per_process * cols * extra_place - pos, "\n");
    }

    if (pos > 0) {
        MPI_File_write_ordered(fh, buffer, pos, MPI_CHAR, &status);
    }
    else {
        fprintf(stderr, "Invalid buffer size on process %d\n", rank);
        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_COUNT);
    }

    free(buffer);
    MPI_File_close(&fh);
}

int get_last_rank(int size) {
    return size - 1;
}

int is_last_rank(int rank, int size) {
    return rank == get_last_rank(size);
}

int* allocate_array(int cols) {
    int* allocated_array = malloc(cols * sizeof(int*));
    return allocated_array;
}

int get_global_row(int local_row, int cols, int rank, int* displs) {
    return (displs[rank] / cols) + local_row;
}

int get_local_row(int global_row, int cols, int rank, int* displs) {
    return global_row - (displs[rank] / cols);
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

void init_basic(double* basics, int length) {
    for (int i = 0; i < length; i++) {
        basics[i] = -1;
    }
}

void write_solution(double* local_matrix, double* basics, int local_rows, int rows, int cols, int rank, int size, int* displs) {
    MPI_File fh;
    MPI_Status status;
    int extra_place = 200;

    int rc = MPI_File_open(MPI_COMM_WORLD, OUTPUT_FILE_NAME, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND, MPI_INFO_NULL, &fh);
    if (rc != MPI_SUCCESS) {
        fprintf(stderr, "Error opening file by process %d\n", rank);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    // Так как базисные переменные не находятся в последней строке - строка целевой функции, то надо исключить ее
    // И последняя строка всегда будет в последным процессе.
    int basics_length = is_last_rank(rank, size) ? local_rows - 1 : local_rows;


    char* buffer = (char*)malloc(basics_length * extra_place * sizeof(char));
    int pos = 0;

    if (basics_length != 0) {
        pos += snprintf(buffer + pos, basics_length * extra_place - pos, "Solution from process %d\n", rank);
    }

    for (int i = 0; i < basics_length; i++) {
        int global_i = get_global_row(i, cols, rank, displs);
        int x_col = (int)basics[global_i];
        if (x_col != -1 && x_col < COLS) {
            pos += snprintf(buffer + pos, basics_length * extra_place - pos, "x%d = %8.2f\n", x_col + 1, local_matrix[i * cols + cols - 1]);
        }
    }

    MPI_File_write_ordered(fh, buffer, pos, MPI_CHAR, &status);

    free(buffer);
    MPI_File_close(&fh);
}

void write_simple_text(const char* text, int rank) {
    MPI_File fh;
    MPI_Status status;

    int rc = MPI_File_open(MPI_COMM_SELF, OUTPUT_FILE_NAME, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND, MPI_INFO_NULL, &fh);
    if (rc != MPI_SUCCESS) {
        fprintf(stderr, "Error opening file by process %d\n", rank);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    MPI_File_seek(fh, 0, MPI_SEEK_END);

    MPI_File_write(fh, text, strlen(text), MPI_CHAR, &status);

    MPI_File_close(&fh);
}


int is_integer(double value) {
    return fabs(value - round(value)) < 1e-6;
}

int exist_real_value(double* tableau, int rows, int cols, double* basics) {
    int exist_real_value = 0;

    for (int i = 0; i < rows - 1; i++) {
        int x_col = (int)basics[i];
        if (x_col < COLS) {
            if (!is_integer(tableau[i * cols + cols - 1])) {
                exist_real_value = 1;
                break;
            }
        }
    }
    return exist_real_value;
}


int find_global_pivot_col(double* local_matrix, int cols) {
    int pivot_col = -1;
    double most_negative = 0;

    for (int j = 0; j < cols - 1; j++) {
        if (local_matrix[j] < most_negative) {
            most_negative = local_matrix[j];
            pivot_col = j;
        }
    }

    return pivot_col;
}

int find_global_pivot_row(double* local_matrix, int local_rows, int cols, int pivot_col, int rank, int size, int* displs) {
    int pivot_row_local = -1;
    double min_ratio_local = INFINITY;
    int last_row = is_last_rank(rank, size) ? local_rows - 1 : local_rows;

    for (int row = 0; row < last_row; row++) {
        if (local_matrix[row * cols + pivot_col] > 0) {
            double ratio = local_matrix[row * cols + (cols - 1)] /
                local_matrix[row * cols + pivot_col];

            if (ratio < min_ratio_local && ratio > 0) {
                min_ratio_local = ratio;
                pivot_row_local = get_global_row(row, cols, rank, displs);
            }
        }
    }

    struct {
        double value;
        int index;
    } local_result, global_result;

    local_result.value = min_ratio_local;
    local_result.index = pivot_row_local;

    MPI_Allreduce(
        &local_result,
        &global_result,
        1,
        MPI_DOUBLE_INT,
        MPI_MINLOC,
        MPI_COMM_WORLD);

    return global_result.index;
}

double try_to_convert_to_positive_zero(double current_value) {
    return round(current_value * 1e9) == 0.0 ? 0.0 : current_value;
}

int is_local_matrix_contains_pivot_row(int start_global_row, int end_global_row, int pivot_row) {
    return start_global_row <= pivot_row && pivot_row <= end_global_row;
}

int find_rank_by_global_row_index(int global_row_index, int cols, int* displs) {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    for (int rank = 0; rank < size - 1; rank++) {
        int rank_max_row_index = (displs[rank + 1] / cols) - 1;
        if (global_row_index <= rank_max_row_index) {
            return rank;
        }
    }
    //The last rank 
    return size - 1;
}

void pivot(double* local_matrix, int local_rows, int cols, int pivot_row, int pivot_col, int rank, int* displs) {
    //Начальный индекс строки локальной матрицы по глобальной
    int start_global_row_index = get_global_row(0, cols, rank, displs);
    //Конечный индекс строки локальной матрицы по глобальной
    int end_global_row_index = get_global_row(local_rows - 1, cols, rank, displs);
    int pivot_local_row = -1;
    //Локальная свобобная строка
    double* pivoted_row = (double*)malloc(cols * sizeof(double));

    if (is_local_matrix_contains_pivot_row(start_global_row_index, end_global_row_index, pivot_row)) {
        pivot_local_row = get_local_row(pivot_row, cols, rank, displs);
        double pivot_value = local_matrix[pivot_local_row * cols + pivot_col];

        for (int j = 0; j < cols; j++) {
            double value = local_matrix[pivot_local_row * cols + j] / pivot_value;
            local_matrix[pivot_local_row * cols + j] = try_to_convert_to_positive_zero(value);
            pivoted_row[j] = local_matrix[pivot_local_row * cols + j];
        }

    }
    int pivoted_row_rank = find_rank_by_global_row_index(pivot_row, cols, displs);
    MPI_Bcast(&pivot_local_row, 1, MPI_INT, pivoted_row_rank, MPI_COMM_WORLD);
    MPI_Bcast(pivoted_row, cols, MPI_DOUBLE, pivoted_row_rank, MPI_COMM_WORLD);

    if (!is_local_matrix_contains_pivot_row(start_global_row_index, end_global_row_index, pivot_row)) {
        for (int i = 0; i < local_rows; i++) {
            double factor = local_matrix[i * cols + pivot_col];
            for (int j = 0; j < cols; j++) {
                local_matrix[i * cols + j] -= factor * pivoted_row[j];
                local_matrix[i * cols + j] = try_to_convert_to_positive_zero(local_matrix[i * cols + j]);
            }
        }
    }
    else {
        for (int i = 0; i < local_rows; i++) {
            if (i != pivot_local_row) {
                double factor = local_matrix[i * cols + pivot_col];
                for (int j = 0; j < cols; j++) {
                    local_matrix[i * cols + j] -= factor * pivoted_row[j];
                    local_matrix[i * cols + j] = try_to_convert_to_positive_zero(local_matrix[i * cols + j]);
                }
            }
        }
    }
}

void pivot_using_simplex_method(double* local_matrix, int local_rows, int cols, int pivot_row, int pivot_col, int rank, int* displs, double* target_function_row) {
    int start_global_row_index = get_global_row(0, cols, rank, displs);
    int end_global_row_index = get_global_row(local_rows - 1, cols, rank, displs);
    int pivot_local_row = -1;
    double* pivoted_row = (double*)malloc(cols * sizeof(double));

    if (is_local_matrix_contains_pivot_row(start_global_row_index, end_global_row_index, pivot_row)) {
        pivot_local_row = get_local_row(pivot_row, cols, rank, displs);
        double pivot_value = local_matrix[pivot_local_row * cols + pivot_col];

        for (int j = 0; j < cols; j++) {
            double value = local_matrix[pivot_local_row * cols + j] / pivot_value;
            local_matrix[pivot_local_row * cols + j] = try_to_convert_to_positive_zero(value);
            pivoted_row[j] = local_matrix[pivot_local_row * cols + j];
        }

    }
    int pivoted_row_rank = find_rank_by_global_row_index(pivot_row, cols, displs);
    MPI_Bcast(&pivot_local_row, 1, MPI_INT, pivoted_row_rank, MPI_COMM_WORLD);
    MPI_Bcast(pivoted_row, cols, MPI_DOUBLE, pivoted_row_rank, MPI_COMM_WORLD);

    if (!is_local_matrix_contains_pivot_row(start_global_row_index, end_global_row_index, pivot_row)) {
        for (int i = 0; i < local_rows; i++) {
            double factor = local_matrix[i * cols + pivot_col];
            for (int j = 0; j < cols; j++) {
                local_matrix[i * cols + j] -= factor * pivoted_row[j];
                local_matrix[i * cols + j] = try_to_convert_to_positive_zero(local_matrix[i * cols + j]);
            }
        }
    }
    else {
        for (int i = 0; i < local_rows; i++) {
            if (i != pivot_local_row) {
                double factor = local_matrix[i * cols + pivot_col];
                for (int j = 0; j < cols; j++) {
                    local_matrix[i * cols + j] -= factor * pivoted_row[j];
                    local_matrix[i * cols + j] = try_to_convert_to_positive_zero(local_matrix[i * cols + j]);
                }
            }
        }
    }

    double target_function_row_factor = target_function_row[pivot_col];
    for (int j = 0; j < cols; j++) {
        target_function_row[j] -= target_function_row_factor * pivoted_row[j];
        target_function_row[j] = try_to_convert_to_positive_zero(target_function_row[j]);
    }

}


int find_gomory_row_to_cut(double* tableau, int rows, int cols) {
    int row_to_cut = -1;
    double max_fractional_part = 0.0;
    int concerned_rows = rows - 1;

    for (int i = 0; i < rows - 1; i++) {
        double value = tableau[i * cols + cols - 1];
        double fractional_part = value - floor(value);

        if (fractional_part > max_fractional_part) {
            max_fractional_part = fractional_part;
            row_to_cut = i;
        }

    }

    return row_to_cut;
}

int find_gomory_column_to_add(double* tableau, int rows, int cols) {
    int gomory_column_to_add = -1;
    double min_value = INFINITY;

    for (int j = 0; j < cols - 2; j++) {
        double goromy_row_value = tableau[(rows - 2) * cols + j];
        if (goromy_row_value != 0) {
            double last_row_gomory_row_rapport = tableau[(rows - 1) * cols + j] / goromy_row_value;
            if (last_row_gomory_row_rapport <= min_value) {
                gomory_column_to_add = j;
                min_value = last_row_gomory_row_rapport;
            }
        }
    }
    return gomory_column_to_add;
}

double* add_gomory_cut(double* tableau, int old_rows, int old_cols, int row_to_cut, int is_first_time) {

    int new_rows = old_rows + 1;
    int new_cols = old_cols + 1;
    int gomory_row = old_rows - 1;
    int gomory_col = old_cols - 1;

    double* result = (double*)malloc(new_rows * new_cols * sizeof(double));

    int i, j;
    int old_col;
    double value, fractionalPart;
    int old_row;

    omp_set_num_threads(4);

#pragma omp parallel for collapse(2) shared(result, tableau, new_rows, new_cols, old_rows, old_cols, row_to_cut, is_first_time, gomory_row, gomory_col) private(i, j, old_col, value, fractionalPart, old_row)
    for (i = 0; i < new_rows; i++) { 
        for (j = 0; j < new_cols; j++) { 
            old_col = j == old_cols ? j - 1 : j;

            // Строка Гомори
            if (i == gomory_row) {
                if (j == gomory_col) {
                    // Пересечения с столбцом Гомори
                    result[i * new_cols + j] = 1;
                }
                else {
                    value = tableau[row_to_cut * old_cols + old_col];
                    fractionalPart = value - floor(value);
                    result[i * new_cols + j] = fractionalPart != 0 ? -1 * fractionalPart : fractionalPart;
                }
            }
            // Последняя строка
            else if (i == old_rows) {
                if (j == gomory_col) {
                    // Пересечения с столбцом Гомори
                    result[i * new_cols + j] = 0;
                }
                else {
                    old_row = i - 1;

                    if (is_first_time) {
                        result[i * new_cols + j] = tableau[old_row * old_cols + old_col] != 0 ? -1 * tableau[old_row * old_cols + old_col] :
                            tableau[old_row * old_cols + old_col];
                    }
                    else {
                        result[i * new_cols + j] = tableau[old_row * old_cols + old_col];
                    }
                }
            }
            // Other rows
            else {
                if (j == gomory_col) {
                    result[i * new_cols + j] = 0;
                }
                else {
                    result[i * new_cols + j] = tableau[i * old_cols + old_col];
                }
            }
        }
    }

    return result;
}

void scatter_global_matrix(
    int rank,
    int size,
    int rows,
    int cols,
    double* global_matrix,
    int** send_recv_counts,
    int** displs,
    int* local_rows,
    double** local_matrix,
    double** tableau_data,
    double*** tableau
) {

    *send_recv_counts = (int*)malloc(size * sizeof(int));
    *displs = (int*)malloc(size * sizeof(int));

    if (rank == MASTER) {
        *tableau = (double**)malloc(rows * sizeof(double*));
        *tableau_data = (double*)malloc(rows * cols * sizeof(double));
        for (int i = 0; i < rows; i++) {
            (*tableau)[i] = &(*tableau_data)[i * cols];
            for (int j = 0; j < cols; j++) {
                (*tableau)[i][j] = global_matrix[i * cols + j];
            }
        }
    }

    int base_rows = rows / size;
    int extra_rows = rows % size;
    int offset = 0;
    for (int i = 0; i < size; i++) {
        (*send_recv_counts)[i] = (i < extra_rows ? base_rows + 1 : base_rows) * cols;
        (*displs)[i] = offset;
        offset += (*send_recv_counts)[i];
    }

    *local_rows = (*send_recv_counts)[rank] / cols;

    *local_matrix = (double*)malloc((*local_rows) * cols * sizeof(double));

    MPI_Scatterv(
        *tableau_data,
        *send_recv_counts,
        *displs,
        MPI_DOUBLE,
        *local_matrix,
        (*send_recv_counts)[rank],
        MPI_DOUBLE,
        MASTER,
        MPI_COMM_WORLD);
}


void apply_gomory_cuts(double* global_matrix, int rows, int cols, double* basics, int rank, int size) {
    int keep_apply_gomory_cut;
    int is_first_time = 1;

    if (rank == MASTER) {
        keep_apply_gomory_cut = exist_real_value(global_matrix, rows, cols, basics);
        if (keep_apply_gomory_cut) {
            write_simple_text("Apply Gomory\n", rank);
        }
    }

    MPI_Bcast(&keep_apply_gomory_cut, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    int row_to_cut = -1;
    int gomory_col = -1;
    while (keep_apply_gomory_cut) {
        if (rank == MASTER) {
            row_to_cut = find_gomory_row_to_cut(global_matrix, rows, cols);
        }

        MPI_Bcast(&row_to_cut, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

        if (row_to_cut == -1) {
            if (rank == MASTER) {
                write_simple_text("All solutions are integers.\n", rank);
            }
            break;
        }

        if (rank == MASTER) {
            //Добавляем дополнительное ограничение Гомори только в основном процессе
            global_matrix = add_gomory_cut(global_matrix, rows, cols, row_to_cut, is_first_time);
        }

        rows++;
        cols++;

        if (rank == MASTER) {
            gomory_col = find_gomory_column_to_add(global_matrix, rows, cols);
        }
        is_first_time = 0;
        MPI_Bcast(&gomory_col, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
        int gomory_row = rows - 2;
        basics = extend_basics(basics, rows - 2, gomory_col);

        int* send_recv_counts = NULL;
        int* displs = NULL;
        int local_rows = 0;
        double* local_matrix = NULL;
        double* tableau_data = NULL;
        double** tableau = NULL;

        //Разделяем основной матрицы между процессами после добавления дополнительного 
        // ограничения Гомори в основной процессе
        scatter_global_matrix(rank, size, rows, cols, global_matrix, &send_recv_counts, &displs, &local_rows, &local_matrix, &tableau_data, &tableau);

        //Локально применим Pivot
        pivot(local_matrix, local_rows, cols, gomory_row, gomory_col, rank, displs);

        //Собираем результат обратно в основной процесс
        MPI_Gatherv(
            local_matrix,
            send_recv_counts[rank],
            MPI_DOUBLE,
            global_matrix,
            send_recv_counts,
            displs,
            MPI_DOUBLE,
            MASTER,
            MPI_COMM_WORLD);

        if (rank == MASTER) {
            keep_apply_gomory_cut = exist_real_value(global_matrix, rows, cols, basics);
            if (!keep_apply_gomory_cut) {
                write_simple_text("All solutions are integers.\n", rank);
            }
        }
        char title[50];
        snprintf(title, sizeof(title), "Added gomory constraint for row %d ", row_to_cut);
        write_matrix_ordered(local_matrix, local_rows, cols, rank, size, title);
        MPI_Bcast(&keep_apply_gomory_cut, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

        //Пишем в файл только окончательное решение 
        if (!keep_apply_gomory_cut) {
            write_solution(local_matrix, basics, local_rows, rows, cols, rank, size, displs);
        } 
    }
}

void simplex_method(double* local_matrix, int local_rows, int rows, int cols, int rank, int size, int* displs,
    int* send_recv_counts, double* tableau_data, double** tableau, double* target_function_row) {
    double* basics = (double*)malloc((rows - 1) * sizeof(double));
    init_basic(basics, (rows - 1));
    int pivot_col = -1;

    while (1) {
        //Находим свободный столбец во всех процессах - лучше на одни плюс Broadcast
        pivot_col = find_global_pivot_col(target_function_row, cols);
      
        if (pivot_col == -1) {
            //Нашли оптимальное решение с не целочисленными 
            //Отправляем в основной процесс 
            MPI_Gatherv(
                local_matrix,
                send_recv_counts[rank],
                MPI_DOUBLE,
                tableau_data,
                send_recv_counts,
                displs,
                MPI_DOUBLE,
                MASTER,
                MPI_COMM_WORLD);

            write_matrix_ordered(local_matrix, local_rows, cols, rank, size, "Optimal solution");
            write_solution(local_matrix, basics, local_rows, rows, cols, rank, size, displs);

            //Применение метода Гомори
            apply_gomory_cuts(tableau_data, rows, cols, basics, rank, size);

            break;
        }

        int pivot_row = -1;
        pivot_row = find_global_pivot_row(
            local_matrix,
            local_rows,
            cols,
            pivot_col,
            rank,
            size,
            displs);
        if (pivot_row == -1) {
            write_simple_text("Bounded - not solution\n", rank);
            break;
        }
        pivot_using_simplex_method(local_matrix, local_rows, cols, pivot_row, pivot_col, rank, displs, target_function_row);
        //Базисное решение
        basics[pivot_row] = pivot_col;

    }
}


int main(int argc, char* argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    FILE* file = fopen(OUTPUT_FILE_NAME, "a");
    if (file == NULL) {
        fprintf(stderr, "Error opening file.\n");
        return -1;
    }

    double start_time = MPI_Wtime();
    clear_file();

    int rows = ALL_ROWS, cols = ALL_COLS;
    double* flat_mat_a = NULL;
    double** tableau = NULL;
    double* tableau_data = NULL;
    int* send_recv_counts = allocate_array(size);
    int* displs = allocate_array(size);
    double* target_function_row = (double*)malloc(cols * sizeof(double));
   
    if (size > rows) {
        if (rank == MASTER) {
            write_simple_text("The number of processes must be <= the rows of initial matrix\n", rank);
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    if (rank == MASTER) {

        tableau = (double**)malloc(rows * sizeof(double*));
        tableau_data = (double*)malloc(rows * cols * sizeof(double));
        for (int i = 0; i < rows; i++) {
            tableau[i] = &tableau_data[i * cols];
        }


        FILE* file = fopen(DATA_FILE_NAME, "r");
        if (file == NULL) {
            perror("Error opening file");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (fscanf_s(file, "%lf", &tableau[i][j]) != 1) {
                    fprintf(stderr, "Error reading data at [%d][%d]\n", i, j);
                    fclose(file);
                    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                }
            }
        }

        for (int j = 0; j < cols; j++) {
            target_function_row[j] = tableau[rows - 1][j];
        }

        fclose(file);

    }

    //Разделение основной матрицы между всеми процессами 
    int base_rows = rows / size;
    int extra_rows = rows % size;
    int offset = 0;

    for (int i = 0; i < size; i++) {
        send_recv_counts[i] = (i < extra_rows ? base_rows + 1 : base_rows) * cols;
        displs[i] = offset;
        offset += send_recv_counts[i];
    }

    int local_rows = send_recv_counts[rank] / cols;

    //Локальная матрица - матрица каждого процесса
    double* local_matrix = (double*)malloc(local_rows * cols * sizeof(double));

    MPI_Scatterv(
        tableau_data,
        send_recv_counts,
        displs,
        MPI_DOUBLE,
        local_matrix,
        send_recv_counts[rank],
        MPI_DOUBLE,
        MASTER,
        MPI_COMM_WORLD);

    //Отправляем копию строки ЦФ всем процессам
    MPI_Bcast(target_function_row, cols, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    //Локально пишем матрицу в файл
    write_matrix_ordered(local_matrix, local_rows, cols, rank, size, "Initial matrix");

    //Локально применим симплекс-метод
    simplex_method(local_matrix, local_rows, rows, cols, rank, size, displs, send_recv_counts, tableau_data, tableau, target_function_row);

    if (rank == MASTER) {
        free(tableau);
        free(tableau_data);
    }

    free(local_matrix);
    free(target_function_row);

    double end_time = MPI_Wtime();
    if (rank == MASTER) {
        char time_text[50];
        snprintf(time_text, sizeof(time_text), "Execution Time: %f seconds\n", end_time - start_time);

        write_simple_text(time_text, rank);
    }

    MPI_Finalize();
    return 0;

}
