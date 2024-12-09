#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define EPSILON 1e-6

typedef struct {
    double** tableau; // The tableau
    int numRows;
    int numCols;
    double* basics; // The basic values
} Simplex;

// Prototypes
void InitSimplex(Simplex* simplex, double** tableau, int numRows, int numCols);
void FreeSimplex(Simplex* simplex);
void Solve(Simplex* simplex);
int FindPivotColumn(Simplex* simplex);
int FindPivotRow(Simplex* simplex, int pivotCol);
void Pivot(Simplex* simplex, int pivotRow, int pivotCol);
void PrintTableau(Simplex* simplex);
int IsInteger(double value);
int ExistNonIntegralSolution(Simplex* simplex);
double GetFractionPartOf(double value);
void AddGomoryConstraint(Simplex* simplex, double* gomoryConstraint);
double* CreateGomoryConstraint(Simplex* simplex, int constraintRowIndex);
int FindXColumnIndexToAddGomoryConstraint(Simplex* simplex, double* solution);
double* ExtendBasics(Simplex* simplex, int gomoryPivotCol);

// Function principal
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        int numRows = 3;
        int numCols = 5;

        // Initial tableau
        double tableauData[3][5] = {
            {-1, 3, 1, 0, 6},
            {7, 1, 0, 1, 35},
            {-7, -9, 0, 0, 0}
        };

        // Allocate memory for tableau
        double** tableau = (double**)malloc(numRows * sizeof(double*));
        for (int i = 0; i < numRows; i++) {
            tableau[i] = (double*)malloc(numCols * sizeof(double));
            for (int j = 0; j < numCols; j++) {
                tableau[i][j] = tableauData[i][j];
            }
        }

        // Initialize simplex structure
        Simplex simplex;
        InitSimplex(&simplex, tableau, numRows, numCols);

        // Solve the problem
        Solve(&simplex);

        // Free memory
        FreeSimplex(&simplex);
    }

    MPI_Finalize();
    return 0;
}

// Initialize the simplex structure
void InitSimplex(Simplex* simplex, double** tableau, int numRows, int numCols) {
    simplex->numRows = numRows;
    simplex->numCols = numCols;

    simplex->tableau = tableau;

    int numBasics = (numCols - 1) / 2;
    simplex->basics = (double*)malloc(numBasics * sizeof(double));
    for (int i = 0; i < numBasics; i++) {
        simplex->basics[i] = -1;
    }
}

// Free resources used by the simplex
void FreeSimplex(Simplex* simplex) {
    for (int i = 0; i < simplex->numRows; i++) {
        free(simplex->tableau[i]);
    }
    free(simplex->tableau);
    free(simplex->basics);
}

// Solve the linear program
void Solve(Simplex* simplex) {
    while (1) {
        int pivotCol = FindPivotColumn(simplex);
        if (pivotCol == -1) {
            printf("Optimal solution found!\n");
            PrintTableau(simplex);
            return;
        }

        int pivotRow = FindPivotRow(simplex, pivotCol);
        if (pivotRow == -1) {
            printf("The problem is unbounded.\n");
            return;
        }

        Pivot(simplex, pivotRow, pivotCol);
        PrintTableau(simplex);
    }
}

// Find the pivot column
int FindPivotColumn(Simplex* simplex) {
    int pivotCol = -1;
    double mostNegative = 0;

    for (int col = 0; col < simplex->numCols - 1; col++) {
        if (simplex->tableau[simplex->numRows - 1][col] < mostNegative) {
            mostNegative = simplex->tableau[simplex->numRows - 1][col];
            pivotCol = col;
        }
    }
    return pivotCol;
}

// Find the pivot row
int FindPivotRow(Simplex* simplex, int pivotCol) {
    int pivotRow = -1;
    double minRatio = INFINITY;

    for (int row = 0; row < simplex->numRows - 1; row++) {
        if (simplex->tableau[row][pivotCol] > 0) {
            double ratio = simplex->tableau[row][simplex->numCols - 1] / simplex->tableau[row][pivotCol];
            if (ratio < minRatio) {
                minRatio = ratio;
                pivotRow = row;
            }
        }
    }
    return pivotRow;
}

// Perform the pivot operation
void Pivot(Simplex* simplex, int pivotRow, int pivotCol) {
    double pivotValue = simplex->tableau[pivotRow][pivotCol];

    // Normalize the pivot row
    for (int col = 0; col < simplex->numCols; col++) {
        if (simplex->tableau[pivotRow][col] != 0) {
            simplex->tableau[pivotRow][col] /= pivotValue;
        }
    }

    // Eliminate other rows
    for (int row = 0; row < simplex->numRows; row++) {
        if (row != pivotRow) {
            double factor = simplex->tableau[row][pivotCol];
            for (int col = 0; col < simplex->numCols; col++) {
                simplex->tableau[row][col] -= factor * simplex->tableau[pivotRow][col];
            }
        }
    }
}

// Print the tableau
void PrintTableau(Simplex* simplex) {
    printf("\nCurrent Tableau:\n");
    for (int row = 0; row < simplex->numRows; row++) {
        for (int col = 0; col < simplex->numCols; col++) {
            printf("%8.2f ", simplex->tableau[row][col]);
        }
        printf("\n");
    }
    printf("\n");
}

// Check if a value is an integer
int IsInteger(double value) {
    return fabs(value - round(value)) < EPSILON;
}
