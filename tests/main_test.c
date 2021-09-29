#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../include/matrix.h"
#include "../include/matrix_classify.h"
#include "../include/qr.h"
#include "../include/utils.h"

void printMatrix(matrix_t *, char *msg);
void matrix_test(void);
void leastSquares_test(void);
void matrix_classify_test(void);

int main() {
  srand(time(0));
  matrix_test();
  leastSquares_test();
  matrix_classify_test();
  return 0;
}

void leastSquares_test(void) {
  matrix_t *A, *B, *x;
  double a[][2] = {{3, -6}, {4, -8}, {0, 1}};
  double b[][1] = {{-1}, {7}, {2}};

  A = makeMatrixFrom2DArray((double *)a, 3, 2);
  B = makeMatrixFrom2DArray((double *)b, 3, 1);
  x = leastSquares(A, B);

  printMatrix(A, "Test: leastSquares() -> A: Ax = B");
  printMatrix(B, "Test: leastSquares() -> B: Ax = B");
  printMatrix(x, "Test: leastSquares() -> Result: x: Ax = B");

  destroyMatrix(A);
  destroyMatrix(B);
  destroyMatrix(x);

  return;
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

void matrix_test(void) {
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
  matrix_t *q_trans = transpose(qr->q);
  matrix_t *q_prod1 = product(qr->q, q_trans);
  matrix_t *q_prod2 = product(q_trans, qr->q);
  printMatrix(q_prod1, "Test #6d: product(q, q^T): Checking for Identity "
                       "matrix (will satisfy only for symmetric cases");
  printMatrix(q_prod2,
              "Test #6e: product(q^T, q): Checking for Identity matrix");

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
  destroyMatrix(q_trans);
  destroyMatrix(q_prod1);
  destroyMatrix(q_prod2);
  return;
}

void matrix_classify_test(void) {
  double i[][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  double d[][3] = {{3, 0, 0}, {0, 11, 0}, {0, 0, 13}};
  double s[][4] = {{5, 0, 0, 1}, {0, 1, 0, 2}, {0, 0, 13, 5}};
  double z[][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  matrix_t *I, *D, *S, *Z;

  I = makeMatrixFrom2DArray((double *)i, 3, 3);
  D = makeMatrixFrom2DArray((double *)d, 3, 3);
  Z = makeMatrixFrom2DArray((double *)z, 3, 3);
  S = makeMatrixFrom2DArray((double *)s, 3, 4);

  assert(isSquareMatrix(I));
  assert(isSquareMatrix(D));
  assert(isSquareMatrix(Z));
  assert(!isSquareMatrix(S));

  assert(isSymmetricMatrix(I));
  assert(isSymmetricMatrix(D));
  assert(isSymmetricMatrix(Z));
  assert(!isSymmetricMatrix(S));

  assert(isDiagonalMatrix(I));
  assert(isDiagonalMatrix(D));
  assert(isDiagonalMatrix(Z));
  assert(!isDiagonalMatrix(S));

  assert(!isZeroMatrix(I));
  assert(!isZeroMatrix(D));
  assert(!isZeroMatrix(S));
  assert(isZeroMatrix(Z));

  assert(isIdentityMatrix(I));
  assert(!isIdentityMatrix(D));
  assert(!isIdentityMatrix(S));
  assert(!isIdentityMatrix(Z));

  printf("----------------------------------------\n");
  printf("Matrix classification checks successful!\n");
  printf("----------------------------------------\n");

  destroyMatrix(I);
  destroyMatrix(D);
  destroyMatrix(S);
  destroyMatrix(Z);

  return;
}