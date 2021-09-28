#ifndef __MATRIX_H_
#define __MATRIX_H_
    #include <stdio.h>
    #include <math.h>
    #include <stdlib.h>
    #include "../include/utils.h"

    typedef struct{
        uint m, n;
        double **M;
    } matrix_t;

    matrix_t* zeros(const uint, const uint);
    void destroyMatrix(matrix_t* );
    void freePartiallyFilledMatrix(matrix_t*, uint);
    matrix_t* transpose(matrix_t*);
    matrix_t* makeMatrixFrom2DArray(const double*, const uint, const uint);

    matrix_t* columnVector(const matrix_t* inputMatrix, const uint columnNo);
    matrix_t* rowVector(const matrix_t* inputMatrix, const uint rowNo);

    matrix_t *product(const matrix_t *, const matrix_t *);
#endif
