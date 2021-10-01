#include "../include/matrix.h"
#include <stdio.h>

matrix_t *zeros(const uint m, uint n) {
  if (m <= 0 || n <= 0)
    return NULL;
  /* m pointers for m rows, n double slots for each of n elements in a row */
  matrix_t *mat = (matrix_t *)calloc(m, sizeof(matrix_t));
  if (mat == NULL) {
    printf("calloc failure; Exiting");
    exit(EXIT_FAIL);
  }

  mat->m = m;
  mat->n = n;
  mat->M = (double **)calloc(m, sizeof(double *));

  for (uint i = 0; i < m; ++i) {
    mat->M[i] = (double *)calloc(n, sizeof(double));
    if (mat->M[i] == NULL) {
      printf("calloc failure; Exiting");
      freePartiallyFilledMatrix(mat, i);
      exit(EXIT_FAIL);
    }
  }

  return mat;
}

/* collapse the two funcitons below into a single one */
void freePartiallyFilledMatrix(matrix_t *mat, uint failedIdx) {
  for (uint i = 0; i < failedIdx; ++i) {
    if (mat->M[i] != NULL)
      free(mat->M[i]);
  }
  free(mat->M);
  free(mat);
}

void destroyMatrix(matrix_t *mat) {
  for (uint i = 0; i < mat->m; ++i) {
    if (mat->M[i] != NULL)
      free(mat->M[i]);
  }
  free(mat->M);
  free(mat);
}
/*---------*/

matrix_t *transpose(matrix_t *mat) {
  matrix_t *newMat;
  newMat = zeros(mat->n, mat->m);
  for (uint i = 0; i < mat->m; ++i) {
    for (uint j = 0; j < mat->n; ++j) {
      newMat->M[j][i] = mat->M[i][j];
    }
  }
  return newMat;
}

matrix_t *makeMatrixFrom2DArray(const double *mat, const uint m, const uint n) {
  matrix_t *createdMatrix = zeros(m, n);
  for (uint i = 0; i < m; ++i) {
    for (uint j = 0; j < n; ++j) {
      createdMatrix->M[i][j] = *((mat + i * n) + j);
    }
  }
  return createdMatrix;
}

matrix_t *columnVector(const matrix_t *inputMatrix, const uint columnNumber) {
  if (columnNumber > inputMatrix->n) {
    printf("Column number (%d) exceeds the number of columns (%d) of input "
           "matrix. Exiting...\n",
           columnNumber, inputMatrix->n);
    exit(0);
  } else if (columnNumber < 1) {
    printf("Column number (%d) cannot be less than 1. Exiting...\n",
           columnNumber);
    exit(0);
  }
  matrix_t *solutionMatrix = zeros(inputMatrix->m, 1);
  for (uint i = 0; i < inputMatrix->m; ++i) {
    solutionMatrix->M[i][0] = inputMatrix->M[i][columnNumber - 1];
  }
  return solutionMatrix;
}

matrix_t *rowVector(const matrix_t *inputMatrix, const uint rowNumber) {
  if (rowNumber > inputMatrix->m) {
    printf("Row number (%d) exceeds the number of rows (%d) of input matrix. "
           "Exiting...\n",
           rowNumber, inputMatrix->n);
    exit(0);
  } else if (rowNumber < 1) {
    printf("Row number (%d) cannot be less than 1. Exiting...\n", rowNumber);
    exit(0);
  }
  matrix_t *solutionMatrix = zeros(1, inputMatrix->n);
  for (uint j = 0; j < inputMatrix->n; ++j) {
    solutionMatrix->M[0][j] = inputMatrix->M[rowNumber - 1][j];
  }
  return solutionMatrix;
}

matrix_t *product(const matrix_t *A, const matrix_t *B) {
  if (A->n != B->m) {
    printf("Cannot multiply the given matrices. A(%d,%d) & B(%d,%d).\n", A->m,
           A->n, B->m, B->n);
    return NULL;
  }
  return partial_product(A, B, A->m, B->n);
}

matrix_t *partial_product(const matrix_t *A, const matrix_t *B,
                          const uint partialI, const uint partialJ) {
  matrix_t *prod = zeros(partialI, partialJ);
  for (uint i = 0; i < partialI; ++i) {
    for (uint j = 0; j < partialJ; ++j) {
      prod->M[i][j] = 0;
      for (uint k = 0; k < A->n; ++k) {
        prod->M[i][j] += A->M[i][k] * B->M[k][j];
      }
    }
  }
  return prod;
}

double two_norm(matrix_t *A, uint colIdx) {
  if (colIdx >= A->n) {
    printf("Error encountered in norm(). Column number cannot exceed number of "
           "columns in the matrix. Exiting...\n");
    exit(1);
  }
  double normSquared = 0;
  for (uint i = 0; i < A->m; ++i) {
    normSquared += pow(A->M[i][colIdx], 2);
  }
  if (normSquared == 0) {
    printf("Unexpected: norm equals 0; Exiting...\n");
    exit(1);
  }
  return sqrt(normSquared);
}

void resetToZero(matrix_t *A) {
  setElementsToOneValue(A, 0.0);
  return;
}

void setElementsToOneValue(matrix_t *A, double val) {
  for (uint i = 0; i < A->m; ++i) {
    for (uint j = 0; j < A->n; ++j) {
      A->M[i][j] = val;
    }
  }
  return;
}

void scalarMultiplyMatrix(matrix_t *q, double scalar, int col) {
  if (q->n == 1) {
    col = 0;
  }
  if (col != -1) {
    for (uint i = 0; i < q->m; ++i) {
      q->M[i][col] = q->M[i][col] * scalar;
    }
  } else {
    for (uint i = 0; i < q->m; ++i) {
      for (uint j = 0; j < q->n; ++j) {
        q->M[i][j] = q->M[i][j] * scalar;
      }
    }
  }
  return;
}

matrix_t *create_random(const uint m, const uint n) {
  matrix_t *mat = zeros(m, n);
  for (uint i = 0; i < m; ++i) {
    for (uint j = 0; j < n; ++j) {
      mat->M[i][j] = (10.0 * rand()) / RAND_MAX;
    }
  }
  return mat;
}

matrix_t *pseudoInverse(matrix_t *A) {
  /* formula: pseudoInv(A) = (A^T * A)^(-1) * A^T */
  matrix_t *mat, *AT, *prod, *inv;
  AT = transpose(A);
  prod = product(AT, A);
  inv = inverse(prod);
  mat = product(inv, AT);

  destroyMatrix(AT);
  destroyMatrix(prod);
  destroyMatrix(inv);
  return mat;
}

double determinant(const matrix_t *A) {
  uint rowSize = A->m;
  uint columnSize = A->n;

  if (rowSize != columnSize) {
    printf("Error: determinant(): Dimension mismatch. Operation Not Permitted. "
           "Exiting... \n");
    exit(1);
  } else if (rowSize == 1)
    return (A->M[0][0]);
  else if (rowSize == 2)
    return (A->M[0][0] * A->M[1][1] - A->M[1][0] * A->M[0][1]);
  else {
    matrix_t *minor = zeros(rowSize - 1, columnSize - 1);
    uint rowMinor, columnMinor;
    uint firstRowColumnIndex;
    double sum = 0;
    register uint row, column;
    // exclude first row and current column
    for (firstRowColumnIndex = 0; firstRowColumnIndex < rowSize;
         firstRowColumnIndex++) {
      rowMinor = 0;
      for (row = 1; row < rowSize; row++) {
        columnMinor = 0;
        for (column = 0; column < columnSize; column++) {
          if (column == firstRowColumnIndex)
            continue;
          else
            minor->M[rowMinor][columnMinor] = A->M[row][column];
          columnMinor++;
        }
        rowMinor++;
      }
      if (firstRowColumnIndex % 2 == 0)
        sum += A->M[0][firstRowColumnIndex] * determinant(minor);
      else
        sum -= A->M[0][firstRowColumnIndex] * determinant(minor);
    }
    destroyMatrix(minor);
    return sum;
  }
}

matrix_t *inverse(const matrix_t *A) {
  if (A->m != A->n) {
    printf(
        "Error: inverse(): Given matrix isn't a square matrix. Exiting....\n");
    exit(1);
  }

  if (A->m < 1) {
    printf("Error: Matrix size cannot be less than 1. Exiting... \n");
    exit(1);
  }

  double det = determinant(A);
  if (det == 0) {
    printf("Exception: inverse(): Zero determinant. Matrix is singular. "
           "Exiting....\n");
    exit(1);
  }

  matrix_t *coeffs, *coefMat;
  coefMat = zeros(A->m, A->n);
  if (A->m == 1) {
    coefMat->M[0][0] = (1.0 / pow(det, 2)) * A->M[0][0];
    return coefMat;
  }
  coeffs = zeros(A->m - 1, A->n - 1);

  uint m, n;
  for (uint i = 0; i < A->m; ++i) {
    for (uint j = 0; j < A->n; ++j) {
      m = 0;
      n = 0;
      for (uint ii = 0; ii < A->m; ++ii) {
        if (ii != i) {
          for (uint jj = 0; jj < A->n; ++jj) {
            if (jj != j) {
              coeffs->M[m][n] = A->M[ii][jj];
              ++n;
            }
          }
          ++m;
          n = 0;
        }
      }
      coefMat->M[i][j] = pow(-1, i + j) * determinant(coeffs);
    }
  }
  matrix_t *trans_coef = transpose(coefMat);
  scalarMultiplyMatrix(trans_coef, 1.0 / det, -1);
  destroyMatrix(coefMat);
  destroyMatrix(coeffs);
  return trans_coef;
}

matrix_t *linearCombination(matrix_t *A, matrix_t *B, double alpha,
                            double beta) {
  if (A->m != B->m || A->n != B->n) {
    printf("Matrix order mismatch. Exiting...\n");
    exit(1);
  } else if (A->m < 1 || A->n < 1) {
    printf("Matrix order cannot be less than 1. Exiting... \n");
    exit(1);
  }
  matrix_t *mat = zeros(A->m, A->n);
  for (uint i = 0; i < A->m; ++i) {
    for (uint j = 0; j < A->n; ++j) {
      mat->M[i][j] = alpha * A->M[i][j] + beta * B->M[i][j];
    }
  }
  return mat;
}

matrix_t *subtract(matrix_t *A, matrix_t *B) {
  return linearCombination(A, B, 1, -1);
}

matrix_t *add(matrix_t *A, matrix_t *B) {
  return linearCombination(A, B, 1, 1);
}

matrix_t *createIdentityMatrix(const uint m){
  if(m < 1){
    printf("Error: createIdentityMatrix(): Matrix size cannot be less than 1. Exiting..");
    exit(1);
  }
  matrix_t *mat = zeros(m, m);
  for(uint i = 0; i < m; ++i){
    mat->M[i][i] = 1;
  }
  return mat;
}

matrix_t *copyMatrix(const matrix_t *mat){
  if(mat->m < 1 || mat->n < 1){
    printf("Error: copyMatrix(): Matrix size cannot be less than 1. Exiting..");
    exit(1);
  }
  matrix_t *newMat = zeros(mat->m, mat->n);
  for(uint i = 0; i < mat->m; ++i){
    for(uint j = 0; j < mat->n; ++j){
      newMat->M[i][j] = mat->M[i][j];
    }
  }
  return newMat;
}