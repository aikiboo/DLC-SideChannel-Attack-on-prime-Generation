#include "gmp.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "function.h"
#define A 2
#define B 3


/*
  Calcul les k premiers et les stocke dans sievePrimeList.
*/


/*
Renvoie l'indice de la valeur absolue maximale de la liste
*/
unsigned long int argmax(double *L, int size) {
    double max = L[0];
    int imax = 0;
    for (int i = 1; i < size; i++) {
        if (fabs(L[i]) > max) {
            max = fabs(L[i]);
            imax = i;
        }
    }
    return max == 0 ? 0 : imax;
}

int main(int argc, char const *argv[]) {
    if (argc < 4) {
        printf("Usage : ./analyse [filename] [sieveSize] [nbInt] \n");
        exit(1);
    }

    FILE *fptrHamAndNoise = fopen(argv[1], "r");
    int nbSievePrimes = atoi(argv[2]) - 1;
    int nb_candidats = atoi(argv[3]);
    double tmp;
    int small_prime;
    mpz_t z_h, z_m, p;
    mpz_inits(z_h, z_m, p, NULL);
    unsigned long int candidats[nbSievePrimes];
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

    for (int j = 0; j < nbSievePrimes; j++) {
        small_prime = mpz_get_ui(sievePrimeList[j]);
        double score[small_prime];
        score[0] = 0;
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
                hypothesis[i] = mpz_popcount(z_h); //w(h - (n - i - 1)*2 mod sj)
                comparaisons[i] = (mesures[i][j] - B) / A;
            }

            score[h] = correlation_coeff(hypothesis, comparaisons, nb_candidats);

        }
        if (j == 10) {
            // for(int l = 0;l<small_prime;l++)printf("Debug : %f \n",score[l]);
        }
        candidats[j] = argmax(score, small_prime);
        gmp_printf("p congru à %lu  mod %Zd \n", candidats[j], sievePrimeList[j]);

    }

    chinese_remainder_theorem(p, sievePrimeList, candidats, nbSievePrimes);
    gmp_printf("p = %Zd\n", p);

}