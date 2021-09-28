#include "../include/matrix.h"
#include "../include/qr.h"
#include "../include/utils.h"
#include <assert.h>
#include <stdio.h>

void printMatrix(matrix_t *, char *msg);
void matrix_test(void);

int main() {
  matrix_test();
  return 0;
}

void printMatrix(matrix_t *mat, char *msg) {
  printf("%s\n", msg);
  printf("----------------------------------------\n");
  for (uint i = 0; i < mat->m; ++i) {
    for (uint j = 0; j < mat->n; ++j) {
      if (j == mat->n - 1) {
        printf("%lf", mat->M[i][j]);
      } else {
        printf("%lf ", mat->M[i][j]);
      }
    }
    printf("\n");
  }
  printf("========================================\n");
  return;
}

void matrix_test() {
  double a[][5] = {{1, 2, 3, 4, 5}, {6, 7, 8, 9, 10}};
  double b[][3] = {{-1, 35, 1}, {0, 3, 3}, {-3, 0, 2}, {1, 1, 1}, {2, 2, 2}};
  matrix_t *A = makeMatrixFrom2DArray((double *)a, 2, 5);
  matrix_t *B = makeMatrixFrom2DArray((double *)b, 5, 3);
  matrix_t *B_transpose = transpose(B);

  printMatrix(A, "Test #1: makeMatrixFrom2DArray() -> A");

  printMatrix(B, "Test #2a: makeMatrixFrom2DArray() -> B");
  printMatrix(B_transpose, "Test #2b: transpose(B) -> B_transpose");

  matrix_t *row_1_A = rowVector(A, 1);
  printMatrix(row_1_A, "Test #3: rowVector()");

  matrix_t *column_1_A = columnVector(A, 1);
  printMatrix(column_1_A, "Test #4: columnVector()");

  matrix_t *C = product(A, B);
  printMatrix(C, "Test #5: product(A,B) -> C = A x B");

  qr_t *qr = computeQR(B);
  printMatrix(qr->q, "Test #6a: QR factorization computeQR(B) -> q");
  printMatrix(qr->r, "Test #6a: QR factorization computeQR(B) -> r");

  destroyMatrix(A);
  destroyMatrix(B);
  destroyMatrix(B_transpose);
  destroyMatrix(C);
  destroyMatrix(column_1_A);
  destroyMatrix(row_1_A);
  destroyQR(qr);
  return;
}
