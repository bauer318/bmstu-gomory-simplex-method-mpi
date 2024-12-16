#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

// Structure for Simplex Method
typedef struct {
    double** tableau;
    int numRows;
    int numCols;
    double* basics;
} SimplexMethod;

int GetBasicElementsNumber(int numCols) {
    return (numCols - 1) / 2;
}

void InitBasic(double* basics, int length) {
    for (int i = 0; i < length; i++) {
        basics[i] = -1;
    }
}

SimplexMethod* CreateSimplexMethod(double** tableau, int numRows, int numCols) {
    SimplexMethod* simplex = (SimplexMethod*)malloc(sizeof(SimplexMethod));
    simplex->tableau = tableau;
    simplex->numRows = numRows;
    simplex->numCols = numCols;
    int basicsLength = GetBasicElementsNumber(numCols);
    simplex->basics = (double*)malloc(basicsLength * sizeof(double));
    InitBasic(simplex->basics, basicsLength);
    return simplex;
}

int FindPivotColumn(SimplexMethod* simplex) {
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

int FindPivotRow(SimplexMethod* simplex, int pivotCol) {
    int pivotRow = -1;
    double minRatio = INFINITY;
    for (int row = 0; row < simplex->numRows - 1; row++) {
        if (simplex->tableau[row][pivotCol] > 0) {
            double ratio = simplex->tableau[row][simplex->numCols - 1] / simplex->tableau[row][pivotCol];
            if (ratio <= minRatio && ratio > 0) {
                minRatio = ratio;
                pivotRow = row;
            }
        }
    }
    return pivotRow;
}

void Pivot(SimplexMethod* simplex, int pivotRow, int pivotCol) {
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

void PrintTableau(SimplexMethod* simplex) {
    printf("\nCurrent Tableau:\n");
    for (int row = 0; row < simplex->numRows; row++) {
        for (int col = 0; col < simplex->numCols; col++) {
            printf("%8.2f ", simplex->tableau[row][col]);
        }
        printf("\n");
    }
    printf("\n");
}

double* PrintSolutionAndGet(SimplexMethod* simplex) {
    printf("Solution:\n");
    int basicsLength = GetBasicElementsNumber(simplex->numCols + 1);
    double* solution = (double*)calloc(basicsLength, sizeof(double));

    for (int i = 0; i < basicsLength; i++) {
        int xCol = (int)simplex->basics[i];
        if (xCol != -1) {
            solution[xCol] = simplex->tableau[i][simplex->numCols - 1];
        }
    }

    for (int i = 0; i < basicsLength; i++) {
        printf("x%d = %.2f\n", i + 1, solution[i]);
    }

    return solution;
}

void Solve(SimplexMethod* simplex) {
    while (1) {
        int pivotCol = FindPivotColumn(simplex);
        if (pivotCol == -1) {
            printf("Optimal solution found!\n");
            double* solution = PrintSolutionAndGet(simplex);
            free(solution);
            return;
        }

        int pivotRow = FindPivotRow(simplex, pivotCol);
        if (pivotRow == -1) {
            printf("The problem is unbounded.\n");
            return;
        }

        Pivot(simplex, pivotRow, pivotCol);
        PrintTableau(simplex);
        simplex->basics[pivotRow] = pivotCol;
    }
}

void AddGomoryCut(SimplexMethod* simplex, int row, int isFirstTime) {
    int newNumRows = simplex->numRows + 1;
    int newNumCols = simplex->numCols + 1;

    // Realloc rows array to hold one more row
    simplex->tableau = (double**)realloc(simplex->tableau, newNumRows * sizeof(double*));
    if (simplex->tableau == NULL) {
        perror("Failed to reallocate rows");
        exit(EXIT_FAILURE);
    }

    // Shift rows down to make room for the new row at the (numRows - 1) position
    for (int i = newNumRows - 1; i > simplex->numRows - 1; i--) {
        simplex->tableau[i] = simplex->tableau[i - 1];
    }

    // Allocate memory for the new row and initialize it
    simplex->tableau[simplex->numRows - 1] = (double*)malloc(newNumCols * sizeof(double));
    if (simplex->tableau[simplex->numRows - 1] == NULL) {
        perror("Failed to allocate new row");
        exit(EXIT_FAILURE);
    }

    for (int j = 0; j < newNumCols; j++) {
        double value = simplex->tableau[row][j];
        double fractionalPart = value - floor(value);

        simplex->tableau[simplex->numRows - 1][j] = fractionalPart != 0 ? -1 * fractionalPart : fractionalPart;
    }

    // Reallocate all rows to add the new column
    for (int i = 0; i < newNumRows; i++) {
        simplex->tableau[i] = (double*)realloc(simplex->tableau[i], newNumCols * sizeof(double));
        if (simplex->tableau[i] == NULL) {
            perror("Failed to reallocate columns");
            exit(EXIT_FAILURE);
        }
    }

    // Shift data in each row to make room for the new column at (numCols - 1) position
    for (int i = 0; i < newNumRows; i++) {
        for (int j = newNumCols - 1; j > simplex->numCols - 1; j--) {
            simplex->tableau[i][j] = simplex->tableau[i][j - 1];
        }
        if (i != simplex->numRows - 1) {
            simplex->tableau[i][simplex->numCols - 1] = 0.0;
        }
        else {
            simplex->tableau[i][simplex->numCols - 1] = 1;
        }

    }

    if (isFirstTime) {
        for (int j = 0; j < newNumCols; j++) {
            double value = simplex->tableau[newNumRows - 1][j];
            simplex->tableau[newNumRows - 1][j] = value != 0 ? -1 * value : value;
        }
    }

    simplex->numRows = newNumRows;
    simplex->numCols = newNumCols;
}

bool IsInteger(double value) {
    return fabs(value - round(value)) < 1e-6;
}

int findRowToCut(SimplexMethod* simplex) {
    int rowToCut = -1;
    double maxFractionalPart = 0.0;

    for (int i = 0; i < simplex->numRows - 1; i++) {
        double value = simplex->tableau[i][simplex->numCols - 1];
        double fractionalPart = value - floor(value);

        if (fractionalPart > maxFractionalPart) {
            maxFractionalPart = fractionalPart;
            rowToCut = i;
        }
    }

    return rowToCut;
}

int ExistRealValue(SimplexMethod* simplex) {
    int existRealValue = 0;
    for (int i = 0; i < simplex->numRows - 1; i++) {
        if (!IsInteger(simplex->tableau[i][simplex->numCols - 1])) {
            existRealValue = 1;
            break;
        }
    }
    return existRealValue;
}

int FindGomoryColumnToAdd(SimplexMethod* simplex) {
    int gomoryColumnToAdd = -1;
    double minValue = INFINITY;

    for (int j = 0; j < simplex->numCols - 2; j++) {
        double goromyRowValue = simplex->tableau[simplex->numRows - 2][j];
        if (goromyRowValue != 0) {
            double lastRowGomoryRowRapport = simplex->tableau[simplex->numRows - 1][j] / goromyRowValue;
            if (lastRowGomoryRowRapport <= minValue) {
                gomoryColumnToAdd = j;
                minValue = lastRowGomoryRowRapport;
            }
        }
    }
    return gomoryColumnToAdd;
}

void ExtendBasics(SimplexMethod* simplex, int gomoryCol) {

    simplex->basics = (double*)realloc(simplex->basics, (simplex->numRows-1) * sizeof(double));
    if (simplex->basics == NULL) {
        perror("Failed to reallocate memory for basics");
        exit(EXIT_FAILURE);
    }

    simplex->basics[simplex->numRows - 2] = gomoryCol;
}

void ApplyGomoryCuts(SimplexMethod* simplex) {
    Solve(simplex); 
    int existRealValue = ExistRealValue(simplex);
    int isFirstTime = 1;
    while (existRealValue) {
       
        int rowToCut = findRowToCut(simplex);
        if (rowToCut == -1) {
            printf("All solutions are integers.\n");
            break;
        }

        printf("Adding Gomory cut for row %d\n", rowToCut);
        AddGomoryCut(simplex, rowToCut, isFirstTime);
        PrintTableau(simplex);

        int gomoryRow = simplex->numRows - 2;
        int gomoryColumn = FindGomoryColumnToAdd(simplex);

        ExtendBasics(simplex, gomoryColumn);
       
        Pivot(simplex, gomoryRow, gomoryColumn);
     
        PrintTableau(simplex);
        double* solution = PrintSolutionAndGet(simplex);

        existRealValue = ExistRealValue(simplex);
        isFirstTime = 0;
    }
}

void FreeSimplex(SimplexMethod* simplex) {
    for (int i = 0; i < simplex->numRows; i++) {
        free(simplex->tableau[i]);
    }
    free(simplex->tableau);
    free(simplex->basics);
    free(simplex);
}

int main() {
    int numRows = 3;
    int numCols = 5;

    // Example tableau: Maximize z = 7x1 + 9x2
    double tableauData[3][5] = {
        { -1, 3, 1, 0, 6 },
        { 7, 1, 0, 1, 35 },
        { -7, -9, 0, 0, 0 }
    };

    double** tableau = (double**)malloc(numRows * sizeof(double*));
    for (int i = 0; i < numRows; i++) {
        tableau[i] = (double*)malloc(numCols * sizeof(double));
        for (int j = 0; j < numCols; j++) {
            tableau[i][j] = tableauData[i][j];
        }
    }

    SimplexMethod* simplex = CreateSimplexMethod(tableau, numRows, numCols);
    ApplyGomoryCuts(simplex);
    FreeSimplex(simplex);

    return 0;
}
