#include "../include/qr.h"
#include <stdio.h>

static void updateQR(matrix_t *, matrix_t *, matrix_t *, uint);
static void sumQRColumn(matrix_t *, matrix_t *, matrix_t *, matrix_t *, uint);
static double dotEm(matrix_t *, uint, matrix_t *, uint);

qr_t* initializeQR(const uint m, const uint n) {
  qr_t *qr = (qr_t *)malloc(sizeof(qr_t));
  qr->q = zeros(m, n);
  qr->r = zeros(n, n);
  return qr;
}

matrix_t **arrayOfColumnsOfMatrix(matrix_t *mat) {
  matrix_t **cols = (matrix_t **)malloc(sizeof(matrix_t) * mat->n);
  for (uint j = 0; j < mat->n; ++j) {
    cols[j] = columnVector(mat, j);
  }
  return cols;
}

qr_t *computeQR(matrix_t *A) {
  if (A->m < 1 || A->n < 1) {
    printf("width/height of matrix cannot be less than 1.\n");
    return NULL;
  }
  qr_t *qr = initializeQR(A->m, A->n);

  //	/* Remove this and handle it in main algo as 2D ops */
  //	matrix_t **colVectorArray = arrayOfColumnsOfMatrix(A);
  //	matrix_t **q = arrayOfColumnsOfMatrix(qr->q);
  //	matrix_t **r = arrayOfColumnsOfMatrix(qr->r);

  for (uint j = 0; j < A->n; ++j) {
    updateQR(qr->q, qr->r, A, j);
  }

  //	destroyColumnVectorArray(colVectorArray, A->n);
  return qr;
}

void destroyColumnVectorArray(matrix_t **arr, uint arrSize) {
  for (uint i = 0; i < arrSize; ++i) {
    destroyMatrix(arr[i]);
  }
  free(arr);
  return;
}

void destroyQR(qr_t *qr) {
  destroyMatrix(qr->q);
  destroyMatrix(qr->r);
  free(qr);
  return;
}

static void updateQR(matrix_t *q, matrix_t *r, matrix_t *a, uint col) {
  if (q->m != a->m || q->n != a->n || r->m != q->n || r->n != q->n ||
      col >= a->n) {
    printf("Order mismatch in updateQR(). Exiting...\n");
    exit(1);
  }
  /* calculate q_i~ = a_j - q_j x R_ji */
  matrix_t *tempSum = zeros(a->m, 1);
  sumQRColumn(tempSum, a, q, r, col);
  for (uint i = 0; i < a->m; ++i) {
    q->M[i][col] = a->M[i][col] - tempSum->M[i][0];
  }
  r->M[col][col] = two_norm(q, col);
  scalarMultiplyMatrix(q, 1 / r->M[col][col], col);
  destroyMatrix(tempSum);
  return;
}

static void sumQRColumn(matrix_t *temp, matrix_t *a, matrix_t *q, matrix_t *r,
                        uint col) {
  resetToZero(temp);
  for (uint i = 0; i < temp->m; ++i) {
    for (uint j = 0; j < col; ++j) {
      r->M[j][col] = dotEm(q, j, a, col);
      temp->M[i][0] += q->M[i][j] * r->M[j][col];
    }
  }
  return;
}

static double dotEm(matrix_t *q, uint qIdx, matrix_t *a, uint aIdx) {
  if (q->m != a->m) {
    printf("Error: Matrix order mismatch in dotEm(). Exiting...\n");
    exit(1);
  }
  double dotSum = 0;
  for (uint i = 0; i < q->m; ++i) {
    dotSum += q->M[i][qIdx] * a->M[i][aIdx];
  }
  return dotSum;
}