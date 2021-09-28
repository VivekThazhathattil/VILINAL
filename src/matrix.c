#include <stdio.h>
#include "../include/matrix.h"

matrix_t* zeros(const uint m, uint n){
    if( m <= 0 || n <= 0)
        return NULL;
    /* m pointers for m rows, n double slots for each of n elements in a row */ 
    matrix_t* mat = (matrix_t*) calloc(m, sizeof(matrix_t));
    if(mat == NULL){
        printf("calloc failure; Exiting");
        exit(EXIT_FAIL);
    }

    mat->m = m;
    mat->n = n;
    mat->M = (double**) calloc(m, sizeof(double*));

    for(uint i = 0; i < m; ++i){
        mat->M[i] = (double *) calloc(n, sizeof(double));
        if(mat->M[i] == NULL){
            printf("calloc failure; Exiting");
            freePartiallyFilledMatrix(mat, i);
            exit(EXIT_FAIL);
        }
    }

    return mat;
}

/* collapse the two funcitons below into a single one */
void freePartiallyFilledMatrix(matrix_t *mat, uint failedIdx){
    for(uint i = 0; i < failedIdx; ++i){
        if(mat->M[i] != NULL)
            free(mat->M[i]);
    }
    free(mat->M);
    free(mat);
}

void destroyMatrix(matrix_t *mat){
    for(uint i = 0; i < mat->m; ++i){
        if(mat->M[i] != NULL)
            free(mat->M[i]);
    }
    free(mat->M);
    free(mat);
}
/*---------*/

matrix_t* transpose(matrix_t *mat){
    matrix_t *newMat;
    newMat = zeros(mat->n, mat->m);
    for(uint i = 0; i < mat->m; ++i){
        for(uint j = 0; j < mat->n; ++j){
            newMat->M[j][i] = mat->M[i][j];
        }
    }
    return newMat;
}

matrix_t* makeMatrixFrom2DArray(const double* mat, const uint m, const uint n){
    matrix_t *createdMatrix = zeros(m,n);
    for(uint i = 0; i < m; ++i){
        for(uint j = 0; j < n; ++j){
            createdMatrix->M[i][j] = *((mat + i*n) + j);
        }
    }
    return createdMatrix;
}

matrix_t* columnVector(const matrix_t* inputMatrix, const uint columnNumber){
    if(columnNumber > inputMatrix->n){
        printf("Column number (%d) exceeds the number of columns (%d) of input matrix. Exiting...\n", columnNumber, inputMatrix->n);
        exit(0);
    }
    else if (columnNumber < 1){
        printf("Column number (%d) cannot be less than 1. Exiting...\n", columnNumber);
        exit(0);
    }
    matrix_t *solutionMatrix = zeros(inputMatrix->m, 1);
    for(uint i = 0; i < inputMatrix->m; ++i){
        solutionMatrix->M[i][0] = inputMatrix->M[i][columnNumber-1];
    }
    return solutionMatrix;
}

matrix_t* rowVector(const matrix_t* inputMatrix, const uint rowNumber){
    if(rowNumber > inputMatrix->m){
        printf("Row number (%d) exceeds the number of rows (%d) of input matrix. Exiting...\n", rowNumber, inputMatrix->n);
        exit(0);
    }
    else if (rowNumber < 1){
        printf("Row number (%d) cannot be less than 1. Exiting...\n", rowNumber);
        exit(0);
    }
    matrix_t *solutionMatrix = zeros(1, inputMatrix->n);
    for(uint j = 0; j < inputMatrix->n; ++j){
        solutionMatrix->M[0][j] = inputMatrix->M[rowNumber - 1][j];
    }
    return solutionMatrix;
}

matrix_t *product(const matrix_t *A, const matrix_t *B){
    if(A->n != B->m){
        printf("Cannot multiply the given matrices. A(%d,%d) & B(%d,%d).\n", A->m, A->n, B->m, B->n);
        return NULL;
    }
    matrix_t *prod = zeros(A->m, B->n);
    for(uint i = 0; i < A->m; ++i){
        for(uint j = 0; j < B->n; ++j){
            prod->M[i][j] = 0;
            for(uint k = 0; k < A->n; ++k){
                prod->M[i][j] += A->M[i][k] * B->M[k][j];
            }
        }
    }
    return prod;
}