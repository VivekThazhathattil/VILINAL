#include "../include/lu.h"

lu_t *initializeLU(const matrix_t *mat) {
  lu_t *lu = (lu_t *)malloc(sizeof(lu_t));
  lu->l = createIdentityMatrix(mat->m);
  lu->u = copyMatrix(mat);
  return lu;
}

uint checkLUMatrixRestrictions(const matrix_t *mat) {
  if (mat->m < 1) {
    printf(
        "Error: computeLU(): width/height of matrix cannot be less than 1.\n");
    return 0;
  } else if (mat->m != mat->n) {
    printf("Error: computeLU(): square matrix expected.\n");
    return 0;
  }
  return 1;
}

// TODO: This doesn't work if any of the diagonanl entries are 0. Use partial
// pivoting instead
lu_t *computeLU(const matrix_t *mat) {
  if (!checkLUMatrixRestrictions(mat)) {
    return NULL;
  }
  lu_t *lu = initializeLU(mat);
  if (mat->m == 1) {
    return lu;
  }
  double pivot = 0;
  double scaleFactor = 1;
  for (uint i = 0; i < lu->u->m; ++i) {
    pivot = lu->u->M[i][i];
    // uint iTemp = i;
    // if(pivot == 0){
    //  while(i+1 < lu->u->m && pivot == 0){

    //  }
    //}
    for (uint ii = i + 1; ii < lu->u->m; ++ii) {
      scaleFactor = lu->u->M[ii][i] / pivot;
      for (uint jj = i; jj < lu->u->n; ++jj) {
        lu->u->M[ii][jj] -= scaleFactor * lu->u->M[i][jj];
        lu->l->M[ii][jj] -= scaleFactor * lu->l->M[i][jj];
      }
    }
  }
  return lu;
}

void destroyLU(lu_t *lu) {
  destroyMatrix(lu->l);
  destroyMatrix(lu->u);
  free(lu);
  return;
}