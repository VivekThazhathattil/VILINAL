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

uint isOrthogonalMatrix(matrix_t *mat) {
  uint res = 0;
  if (!isSquareMatrix(mat)) {
    return res;
  }

  matrix_t *transMat, *prodLeft, *prodRight, *invMat;
  transMat = transpose(mat);
  prodLeft = product(transMat, mat);
  prodRight = product(mat, transMat);
  invMat = NULL;

  if (isIdentityMatrix(prodLeft) && isIdentityMatrix(prodRight)) {
    res = 1;
  }

  if (determinant(mat) == 0) {
    res = 0;
  }

  else {
    invMat = inverse(mat);
    if (!compareMatrices(invMat, transMat)) {
      res = 0;
    }
  }

  destroyMatrix(transMat);
  destroyMatrix(prodLeft);
  destroyMatrix(prodRight);
  if (invMat != NULL) {
    destroyMatrix(invMat);
  }
  return res;
}

static uint triangularHelper(matrix_t *mat, uint checkDirection) {
  for (uint i = 0; i < mat->m; ++i) {
    for (uint j = 0; j < mat->n; ++j) {
      if (checkDirection) {
        if (i < j && mat->M[i][j] != 0) {
          return 0;
        }
      } else {
        if (i > j && mat->M[i][j] != 0) {
          return 0;
        }
      }
    }
  }
  return 1;
}

uint isUpperTriangular(matrix_t *mat) { return triangularHelper(mat, 0); }

uint isLowerTriangular(matrix_t *mat) { return triangularHelper(mat, 1); }

uint isTriangular(matrix_t *mat) {
  if (isUpperTriangular(mat) || isLowerTriangular(mat)) {
    return 1;
  }
  return 0;
}

uint compareMatrices(matrix_t *mat1, matrix_t *mat2) {
  if (mat1->m != mat2->m || mat1->n != mat2->n) {
    return 0;
  }
  for (uint i = 0; i < mat1->m; ++i) {
    for (uint j = 0; j < mat2->n; ++j) {
      if (mat1->M[i][j] != mat2->M[i][j]) {
        return 0;
      }
    }
  }
  return 1;
}
