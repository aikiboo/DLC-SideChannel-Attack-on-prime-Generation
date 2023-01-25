#include "gmp.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "common/gauss.c"

#define A 2
#define B 3
#define PRIME 1
#define COMPOSITE 0

/*
Structure pour notre tableau de résidus
*/
typedef struct {
    int has_a_0;
    mpz_t *array;
} MpzResidues;

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

/*
  génère un nombre de k bits
*/
void gen_k_bits_number(mpz_t out, int size, gmp_randstate_t randstate) {
    mpz_t offset;
    mpz_inits(offset, out, NULL);
    mpz_set_ui(offset, 2);
    mpz_pow_ui(offset, offset, size - 1);
    mpz_urandomb(out, randstate, size - 1);
    mpz_add(out, out, offset);
    mpz_clear(offset);
}

/*
Renvoi un nombre impair, en appellant gen_k_bits_number() et, ajoute 1 si la sortie est pair
*/
void gen_k_bits_number_odd(mpz_t out, int size, gmp_randstate_t randstate) {
    gen_k_bits_number(out, size, randstate);
    if (mpz_divisible_ui_p(out, 2)) {
        mpz_add_ui(out, out, 1);
    }
}

/*
Petmet d'initialiser le tableau de résidus à 0
*/
void update_residues_array(MpzResidues *mpzResidues, mpz_t value, int size, mpz_t *sievePrimeList) {
    for (int i = 0; i < size; i++) {
        mpz_mod(mpzResidues->array[i], value, sievePrimeList[i]);
    }
}

void count_residues_equals_0(MpzResidues *mpzResidues, int size) {
    int tmp = 0;
    for (int i = 0; i < size && !tmp; i++) {
        if (mpz_cmp_ui(mpzResidues->array[i], 0) == 0) {
            tmp = 1;
        }
    }

    mpzResidues->has_a_0 = tmp;
}

void writefiles(MpzResidues *mpzResidues, FILE *fptrValue, FILE *fptrHam, FILE *fptrHam2,
                FILE *fptrHamAndNoise, int size) {
    for (int i = 0; i < size; i++) {
        gmp_fprintf(fptrValue, "%Zu\r\n", mpzResidues->array[i]);
        gmp_fprintf(fptrHam, "%u\n", mpz_popcount(mpzResidues->array[i]));
        gmp_fprintf(fptrHam2, "%u\n", A * mpz_popcount(mpzResidues->array[i]) + B);
        gmp_fprintf(fptrHamAndNoise, "%lf\n", (double) (A * mpz_popcount(mpzResidues->array[i]) + B + gauss()));
    }
}

double correlation_coeff(double *X, double *Y, int size) {
    double mean_X = 0;
    double mean_Y = 0;

    for (int i = 0; i < size; i++) {
        printf("%lf\t %lf\n", X[i], Y[i]);
        mean_X = mean_X + X[i];
        mean_Y = mean_Y + Y[i];
    }
    mean_X = mean_X / size;
    mean_Y = mean_Y / size;
    printf("Debug : %lf\t %lf\n", mean_X, mean_Y);
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
        mpz_mod(p, p, N);
    }
    mpz_mod(p, p, N);


}

int main(int argc, char const *argv[]) {
    if (argc < 3) {
        printf("Usage : ./genDots [nbPrimeInSieve] [sizeOfGeneratedRandom] /\n");
        exit(1);
    }

    int size = 0, nbSievePrimes = 0;
    int nb_candidats = 1;
    int small_prime;
    mpz_t tmpRand;
    mpz_t *sievePrimeList;
    MpzResidues *mpzResidues = malloc(sizeof(MpzResidues));

//setup du gen random
    gmp_randstate_t randstate;
    gmp_randinit_default(randstate);
    gmp_randseed_ui(randstate, time(NULL));

//setup des fichiers de sortie
    FILE *fptrValue, *fptrHam, *fptrHam2, *fptrHamAndNoise;
    fptrValue = fopen("output_value.txt", "w");
    fptrHam = fopen("output_hamming.txt", "w");
    fptrHam2 = fopen("output_hamming_with_func.txt", "w");
    fptrHamAndNoise = fopen("output_hamming_with_func_and_noise.txt", "w");

    if (fptrValue == NULL || fptrHam == NULL) {
        printf("Error!");
        exit(1);
    }

    nbSievePrimes = atoi(argv[1]);
    size = atoi(argv[2]);

//setup du sieve
    sievePrimeList = malloc(sizeof(mpz_t) * nbSievePrimes);
    find_k_first_primes(nbSievePrimes, sievePrimeList);

//on gen un random impair
    gen_k_bits_number_odd(tmpRand, size, randstate);

//setup des MpzResidues
    mpzResidues->array = malloc(sizeof(mpz_t) * nbSievePrimes);
    for (int i = 0; i < nbSievePrimes; i++)mpz_init(mpzResidues->array[i]);
    mpzResidues->has_a_0 = 0;

//execution du crible opti
    update_residues_array(mpzResidues, tmpRand, nbSievePrimes, sievePrimeList);
    writefiles(mpzResidues, fptrValue, fptrHam, fptrHam2, fptrHamAndNoise, nbSievePrimes);
    count_residues_equals_0(mpzResidues, nbSievePrimes);

    while (1) {
        while (mpzResidues->has_a_0) {
            for (int i = 1; i < nbSievePrimes; i++) {
                mpz_add_ui(mpzResidues->array[i], mpzResidues->array[i], 2);
                mpz_mod(mpzResidues->array[i], mpzResidues->array[i], sievePrimeList[i]);
            }
            nb_candidats++;
            gmp_printf("%Zu\n", tmpRand);
            mpz_add_ui(tmpRand, tmpRand, 2);
            writefiles(mpzResidues, fptrValue, fptrHam, fptrHam2, fptrHamAndNoise, nbSievePrimes);
            count_residues_equals_0(mpzResidues, nbSievePrimes);
        }
        if (mpz_probab_prime_p(tmpRand, 10))break;
        mpzResidues->has_a_0 = 1;
    }

    printf("Is prime (1=prob,2=100\%) : %d\n", mpz_probab_prime_p(tmpRand, 10));

    fclose(fptrHamAndNoise);
    fptrHamAndNoise = fopen("output_hamming_with_func_and_noise.txt", "r");
    rewind(fptrHamAndNoise);

    double tmp;
    mpz_t z_h;
    mpz_t z_m;
    mpz_init(z_h);
    mpz_init(z_m);
    unsigned long int candidats[nbSievePrimes - 1];
    mpz_t p;
    mpz_init(p);


    double mesures[nb_candidats][nbSievePrimes];

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

    //gmp_printf("nb of 1 %u\n",mpz_popcount(tmpRand));
    gmp_printf("%Zu\n", tmpRand);

    fclose(fptrValue);
    fclose(fptrHam);
    fclose(fptrHam2);
    fclose(fptrHamAndNoise);

    return 0;
}
