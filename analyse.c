#include "gmp.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#define A 2
#define B 3


/*
  Calcul les k premiers et les stocke dans sievePrimeList.
*/
void find_k_first_primes(int size, mpz_t *sievePrimeList) {
    mpz_t current;
    int i = 0, j, isPrime;

    mpz_init(current);
    mpz_set_ui(current, 2);

    while (i < size) {
        isPrime = 1;

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
    if (denum == 0)return 0;
    return num / denum;
}

/*
Renvoie l'indice de la valeur absolue maximale de la liste
*/
unsigned long int argmax(double *L, int size) {
    double max = L[1];
    int imax = 1;
    for (int i = 2; i < size; i++) {
        if (fabs(L[i]) > max) {
            max = fabs(L[i]);
            imax = i;
        }
    }
    return imax;
}
void chinese_remainder_theorem(mpz_t p, mpz_t *sievePrimeList, unsigned long int *candidats, int nbSievePrimes) {
    mpz_t N;
    mpz_t current;
    mpz_t Ni;
    mpz_t inv_Ni;

    mpz_inits(N, current, Ni, inv_Ni, NULL);
    mpz_set_ui(N, 1);

    mpz_t MAX;
    mpz_init(MAX);
    mpz_ui_pow_ui(MAX, 2, 33);


    for (int i = 0; i < nbSievePrimes - 1; i++) {
        mpz_mul(N, N, sievePrimeList[i]);
    }
    gmp_printf("Debug : N = %Zu\n", N);
    gmp_printf("Debug : MAX = %Zu\n", MAX);
    for (int i = 0; i < nbSievePrimes - 1; i++) {
        mpz_divexact(Ni, N, sievePrimeList[i]);
        mpz_invert(inv_Ni, Ni, sievePrimeList[i]);
        mpz_mul_ui(current, Ni, candidats[i]);
        mpz_mul(current, current, inv_Ni);
        mpz_add(p, p, current);
        mpz_mod(p, p, MAX);
    }
    mpz_mod(p, p, MAX);


}

int main(int argc, char const *argv[]) {
    if (argc < 4) {
        printf("Usage : ./analyse [filename] [sieveSize] [nbInt] \n");
        exit(1);
    }

    FILE *fptrHamAndNoise = fopen(argv[1], "r");
    int nbSievePrimes = atoi(argv[3]);
    int nb_candidats = atoi(argv[3]);
    double tmp;
    int small_prime;
    mpz_t z_h, z_m, p;
    mpz_inits(z_h, z_m, p, NULL);
    unsigned long int candidats[nbSievePrimes - 1];
    double mesures[nb_candidats][nbSievePrimes];
    mpz_t *sievePrimeList;

    //setup du gen random
    gmp_randstate_t randstate;
    gmp_randinit_default(randstate);
    gmp_randseed_ui(randstate, time(NULL));


    sievePrimeList = malloc(sizeof(mpz_t) * nbSievePrimes);
    find_k_first_primes(nbSievePrimes, sievePrimeList);



    /*
    Récupération des mesures
    */

    for (int i = 0; i < nb_candidats; i++) {
        for (int j = 0; j < nbSievePrimes; j++) {
            fscanf(fptrHamAndNoise, "%lf", &tmp);
            mesures[i][j] = tmp;
        }
    }

    for (int j = 1; j < nbSievePrimes; j++) {
        small_prime = mpz_get_ui(sievePrimeList[j]);
        double score[small_prime];
        int debug = 0;
        for (int h = 1; h < small_prime; h++) {
            //mpz_set_ui(z_h,h);
            double hypothesis[nb_candidats];
            double comparaisons[nb_candidats];
            for (int i = 0; i < nb_candidats; i++) {
                mpz_set_ui(z_h, h);
                mpz_set_ui(z_m, nb_candidats);
                mpz_sub_ui(z_m, z_m, i);
                mpz_sub_ui(z_m, z_m, 1);
                mpz_mul_ui(z_m, z_m, 2);
                mpz_sub(z_h, z_h, z_m);
                mpz_mod(z_h, z_h, sievePrimeList[j]);
                hypothesis[i] = A * mpz_popcount(z_h) + B; //w(h - (n - i - 1)*2 mod sj)
                comparaisons[i] = mesures[i][j];
            }

            score[h] = correlation_coeff(hypothesis, comparaisons, nb_candidats);
            if (score[h] != 0)debug++;
            else {
                for (int k = 0; k < nb_candidats; k++) {
                    printf("Debug : %lf\n", hypothesis[k]);
                    printf("Debug 2: %lf\n", comparaisons[k]);
                }
            }
            if (j == 1)printf("Score : %lf\n", score[h]);
        }
        if (debug == 0) {
            printf("Error score");
            exit(4);
        }
        candidats[j - 1] = argmax(score, small_prime);
        gmp_printf("p congru à %lu  mod %Zd \n", candidats[j - 1], sievePrimeList[j]);

    }

    chinese_remainder_theorem(p, sievePrimeList, candidats, nbSievePrimes);
    gmp_printf("p = %Zd\n", p);

}