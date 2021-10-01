#ifndef __MATRIX_CLASSIFY_H__
#define __MATRIX_CLASSIFY_H__

#include "../include/matrix.h"
#include "../include/utils.h"

uint isSquareMatrix(matrix_t *);
uint isDiagonalMatrix(matrix_t *);
uint isSymmetricMatrix(matrix_t *);
uint isZeroMatrix(matrix_t *);
uint isIdentityMatrix(matrix_t *);
uint compareMatrices(matrix_t *, matrix_t *);

#endif