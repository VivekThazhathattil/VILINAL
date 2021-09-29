#ifndef __MATRIX_H_
#define __MATRIX_H_
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../include/mem.h" // use ony for debug
#include "../include/utils.h"

typedef struct {
  uint m, n;
  double **M;
} matrix_t;

matrix_t *zeros(const uint, const uint);
void destroyMatrix(matrix_t *);
void freePartiallyFilledMatrix(matrix_t *, uint);
matrix_t *transpose(matrix_t *);
matrix_t *makeMatrixFrom2DArray(const double *, const uint, const uint);

matrix_t *columnVector(const matrix_t *inputMatrix, const uint columnNo);
matrix_t *rowVector(const matrix_t *inputMatrix, const uint rowNo);

matrix_t *product(const matrix_t *, const matrix_t *);
matrix_t *partial_product(const matrix_t *, const matrix_t *, const uint,
                          const uint);
matrix_t *add(matrix_t *, matrix_t *);
matrix_t *subtract(matrix_t *, matrix_t *);
matrix_t *linearCombination(matrix_t *, matrix_t *, double, double);

double two_norm(matrix_t *, uint);

void resetToZero(matrix_t *);
void setElementsToOneValue(matrix_t *, double);
void scalarMultiplyMatrix(matrix_t *, double, int);
matrix_t *create_random(const uint, const uint);
double determinant(const matrix_t *);
matrix_t *inverse(const matrix_t *);
matrix_t *pseudoInverse(matrix_t *);
#endif