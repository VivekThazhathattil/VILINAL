#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../include/matrix.h"
#include "../include/qr.h"
#include "../include/utils.h"

void printMatrix(matrix_t *, char *msg);
void matrix_test(void);

int main() {
  srand(time(0));
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
  printMatrix(qr->r, "Test #6b: QR factorization computeQR(B) -> r");

  matrix_t *D = product(qr->q, qr->r);
  printMatrix(D, "Test #6c: product(q,r): Checking if B = q x r");

  matrix_t *E = create_random(5, 5);
  qr_t *qr2 = computeQR(E);
  printMatrix(qr2->q, "Test #7a: QR factorization computeQR(E) -> q");
  printMatrix(qr2->r, "Test #7b: QR factorization computeQR(E) -> r");
  matrix_t *F = product(qr2->q, qr2->r);
  printMatrix(E, "Test #7c: Original matrix E");
  printMatrix(F, "Test #7d: product(q,r): Checking if E = F = q x r");
  printf("Determinant of E = %lf\n\n", determinant(E));
  matrix_t *E_inv = inverse(E);
  printMatrix(E_inv, "Test #8a: inverse(E) -> E^(-1)");
  matrix_t *G = pseudoInverse(E);
  printMatrix(
      G,
      "Test #8b: pseudoInverse(E) -> should equal (E^T E)^(-1) * E^T = E^(-1)");

  destroyMatrix(A);
  destroyMatrix(B);
  destroyMatrix(B_transpose);
  destroyMatrix(C);
  destroyMatrix(D);
  destroyMatrix(E);
  destroyMatrix(E_inv);
  destroyMatrix(F);
  destroyMatrix(G);
  destroyMatrix(column_1_A);
  destroyMatrix(row_1_A);
  destroyQR(qr);
  destroyQR(qr2);
  return;
}
