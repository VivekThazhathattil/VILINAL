#include "../include/matrix_classify.h"

uint isSquareMatrix(matrix_t *mat) { return (mat->m == mat->n) ? 1 : 0; }

uint isDiagonalMatrix(matrix_t *mat) {
  if (!isSquareMatrix(mat)) {
    return 0;
  }
  for (uint i = 0; i < mat->m; ++i) {
    for (uint j = 0; j < mat->n; ++j) {
      if (i != j) {
        if (mat->M[i][j] != 0) {
          return 0;
        }
      }
    }
  }
  return 1;
}

uint isSymmetricMatrix(matrix_t *mat) {
  if (!isSquareMatrix(mat)) {
    return 0;
  }
  for (uint i = 0; i < mat->m; ++i) {
    for (uint j = 0; j < mat->n; ++j) {
      if (mat->M[i][j] != mat->M[j][i])
        return 0;
    }
  }
  return 1;
}

uint isZeroMatrix(matrix_t *mat) {
  for (uint i = 0; i < mat->m; ++i) {
    for (uint j = 0; j < mat->n; ++j) {
      if (mat->M[i][j] != 0) {
        return 0;
      }
    }
  }
  return 1;
}

uint isIdentityMatrix(matrix_t *mat) {
  if (!isSquareMatrix(mat)) {
    return 0;
  }

  for (uint i = 0; i < mat->m; ++i) {
    for (uint j = 0; j < mat->n; ++j) {
      if (i == j) {
        if (mat->M[i][j] != 1) {
          return 0;
        }
      } else {
        if (mat->M[i][j] != 0) {
          return 0;
        }
      }
    }
  }
  return 1;
}