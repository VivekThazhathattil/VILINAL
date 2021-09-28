#ifndef __QR_H__
#define __QR_H__

#include "../include/matrix.h"

typedef struct{
    matrix_t *q;
    matrix_t *r;
} qr_t;

void initializeQR(qr_t *, const uint m, const uint n);

#endif
