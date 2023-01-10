#include "gmp.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "gauss.c"
#include "common/function.h"

#define A 2
#define B 3

/*
Structure pour notre tableau de résidus
*/
typedef struct {
    int has_a_0;
    mpz_t *array;
} MpzResidues;

/*
  génère un nombre de k bits
*/


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
        gmp_fprintf(fptrValue, "%Zu\n", mpzResidues->array[i]);
        gmp_fprintf(fptrHam, "%u\n", mpz_popcount(mpzResidues->array[i]));
        gmp_fprintf(fptrHam2, "%u\n", A * mpz_popcount(mpzResidues->array[i]) + B);
        gmp_fprintf(fptrHamAndNoise, "%lf\n", A * mpz_popcount(mpzResidues->array[i]) + B + gauss());
    }
}

int main(int argc, char const *argv[]) {
    if (argc < 3) {
        printf("Usage : ./genDots [nbPrimeInSieve] [sizeOfGeneratedRandom] /\n");
        exit(1);
    }

    srand(time(NULL));
    int size = 0, nbSievePrimes = 0,nb_candidats=1;
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
            nb_candidats++;
            for (int i = 0; i < nbSievePrimes; i++) {
                mpz_add_ui(mpzResidues->array[i], mpzResidues->array[i], 2);
                mpz_mod(mpzResidues->array[i], mpzResidues->array[i], sievePrimeList[i]);
            }
            mpz_add_ui(tmpRand, tmpRand, 2);
            writefiles(mpzResidues, fptrValue, fptrHam, fptrHam2, fptrHamAndNoise, nbSievePrimes);
            count_residues_equals_0(mpzResidues, nbSievePrimes);
        }
        if (mpz_probab_prime_p(tmpRand, 10))break;
        mpzResidues->has_a_0 = 1;
    }
    printf("Is prime (1=prob,2=100\%) : %d\n", mpz_probab_prime_p(tmpRand, 10));

    //gmp_printf("nb of 1 %u\n",mpz_popcount(tmpRand));
    printf("Nb de candidats : %d\n",nb_candidats);
    gmp_printf("Premier: \t%Zu\n", tmpRand);

    fclose(fptrValue);
    fclose(fptrHam);
    fclose(fptrHam2);
    fclose(fptrHamAndNoise);

    return 0;
}
