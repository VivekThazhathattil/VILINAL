#ifndef __LU_H__
#define __LU_H__

#include "../include/utils.h"
#include "../include/matrix.h"

typedef struct{
    matrix_t *l, *u;
} lu_t;

lu_t* initializeLU(const matrix_t *mat);
uint checkLUMatrixRestrictions(const matrix_t *);
lu_t* computeLU(const matrix_t *);
void destroyLU(lu_t *);

#endif
