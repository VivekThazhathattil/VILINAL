#include  "../include/qr.h"

void initializeQR(qr_t * qr, const uint m, const uint n){
   qr = (qr_t *)malloc(sizeof(qr_t));
   qr->q = zeros(m,n);
   qr->r = zeros(n,n);
   return;
}