#include "function.h"


void chinese_remainder_theorem(mpz_t p, mpz_t *sievePrimeList, unsigned long int *candidats, int nbSievePrimes) {
    mpz_t N;
    mpz_t current;
    mpz_t Ni;
    mpz_t inv_Ni;

    mpz_inits(N, current, Ni, inv_Ni, NULL);
    mpz_set_ui(N, 1);


    for (int i = 0; i < nbSievePrimes; i++) {
        mpz_mul(N, N, sievePrimeList[i]);
    }
    gmp_printf("Debug :\nN = %Zu\n", N);
    for (int i = 0; i < nbSievePrimes; i++) {
        mpz_divexact(Ni, N, sievePrimeList[i]);
        mpz_invert(inv_Ni, Ni, sievePrimeList[i]);
        mpz_mul_ui(current, Ni, candidats[i]);
        mpz_mul(current, current, inv_Ni);
        mpz_add(p, p, current);
        mpz_mod(p, p, N);
    }
    mpz_mod(p, p, N);


}


void find_k_first_primes(int size, mpz_t *sievePrimeList) {
    mpz_t current;
    int i = 0, j, isPrime;

    mpz_init(current);
    mpz_set_ui(current, 3);

    while (i < size) {
        isPrime = !mpz_divisible_ui_p(current, 2);

        for (j = 0; j < i && isPrime; j++) {
            if (mpz_divisible_p(current, sievePrimeList[j]))isPrime = 0;
        }

        if (isPrime) {
            mpz_init(sievePrimeList[i]);
            mpz_set(sievePrimeList[i], current);
            i++;
        }

        mpz_add_ui(current, current, 1);
    }


}

double correlation_coeff(double *X, double *Y, int size) {
    double mean_X = 0;
    double mean_Y = 0;

    for (int i = 0; i < size; i++) {
        mean_X = mean_X + X[i];
        mean_Y = mean_Y + Y[i];
    }
    mean_X = mean_X / size;
    mean_Y = mean_Y / size;

    double num = 0;
    double denum;
    double sx = 0;
    double sy = 0;
    for (int i = 0; i < size; i++) {
        num = num + (X[i] - mean_X) * (Y[i] - mean_Y);
        sx = sx + pow(X[i] - mean_X, 2);
        sy = sy + pow(Y[i] - mean_Y, 2);
    }
    denum = sqrt(sx * sy);
    if (denum == 0)return 2;
    return num / denum;
}