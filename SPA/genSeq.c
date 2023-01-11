#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../common/function.h"

void main(int argc, char const *argv[]) {

    if (argc < 3) {
        printf("Usage : ./genSeq [nbPrimeInSieve] [sizeOfGeneratedRandom] \n");
        exit(1);
    }

    int nbSievePrimes = atoi(argv[1]);
    int size = atoi(argv[2]);
    srand(time(NULL));

    //Setup de la sortie
    FILE *output = fopen("output_sequence.txt", "w");

    //setup du gen random
    gmp_randstate_t randstate;
    gmp_randinit_default(randstate);
    gmp_randseed_ui(randstate, time(NULL));

    //setup du sieve
    mpz_t * sievePrimeList = malloc(sizeof(mpz_t) * nbSievePrimes);
    find_k_first_primes(nbSievePrimes, sievePrimeList);
    //Génération du premier entier impair
    mpz_t tmpRand;
    gen_k_bits_number_odd(tmpRand, size, randstate);
    int isPrime;
    while(1){
        isPrime = 0;
        for(int i =0;i<nbSievePrimes;i++){
            fprintf(output,"c");
            if(mpz_divisible_p(tmpRand,sievePrimeList[i])){
                isPrime =1;
                break;
            }
        }
        //si il a un diviseur dans le sieve
        if(isPrime==0){
            fprintf(output,"a");
            if(mpz_probab_prime_p(tmpRand, 10))break;

        }
        fprintf(output,"b");
        mpz_add_ui(tmpRand,tmpRand,2);
    }
    fclose(output);
}