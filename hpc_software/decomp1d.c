#include <stdio.h>

int decomp1d(int n, int p, int myid, int *s, int *e)
{
    int base = n / p;
    int rem  = n % p;

    if (myid < rem) {
        *s = myid * (base + 1);
        *e = *s + base;
    } else {
        *s = myid * base + rem;
        *e = *s + base - 1;
    }

    return 0;
}

