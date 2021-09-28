#include  "../include/qr.h"
#include <stdio.h>

static void updateQR(matrix_t *, matrix_t *, matrix_t*, uint);

void initializeQR(qr_t *qr, const uint m, const uint n){
   qr = (qr_t *)malloc(sizeof(qr_t));
   qr->q = zeros(m,n);
   qr->r = zeros(n,n);
   return;
}

matrix_t** arrayOfColumnsOfMatrix(matrix_t *mat){
	matrix_t **cols = (matrix_t **)malloc(sizeof(matrix_t) * mat->n);
	for(uint j = 0; j < mat->n; ++j){
		cols[j] = columnVector(mat, j);
	}
	return cols;
}

qr_t* computeQR(matrix_t *A){
	if(A->m < 1 || A->n < 1){
		printf("width/height of matrix cannot be less than 1.\n");
		return NULL;
	}
	qr_t *qr;
	initializeQR(qr, A->m, A->n);
	
//	/* Remove this and handle it in main algo as 2D ops */
//	matrix_t **colVectorArray = arrayOfColumnsOfMatrix(A);
//	matrix_t **q = arrayOfColumnsOfMatrix(qr->q);
//	matrix_t **r = arrayOfColumnsOfMatrix(qr->r);

	for(uint j = 0; j < A->n; ++j){
		updateQR(qr->q, qr->r, A, j);
	}

//	destroyColumnVectorArray(colVectorArray, A->n);
	return qr;
}

void destroyColumnVectorArray(matrix_t** arr, int arrSize){
	for(uint i = 0; i < arrSize; ++i){
		destroyMatrix(arr[i]);
	}
	free(arr);
	return;
}

static void updateQR(matrix_t *q, matrix_t *r, matrix_t *a, uint col){
	if(q->m != a->m || q->n != a->n || r->m != a->n || r->n != a->n || col >= a->n){
		printf("Order mismatch in updateQR(). Exiting...\n");
		exit(1);
	}
	for(uint i = 0; i < a->m; ++i){
		q->M[i][col] = a->M[i][col] - sum(col)
	}
	return;
}
