#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../common/function.h"


void genPrimeAndLog(mpz_t prime,int size,int nbSievePrimes,mpz_t * sievePrimeList,gmp_randstate_t randstate,FILE* outputFile){
    gen_k_bits_number_odd(prime, size, randstate);
    int isPrime;
    while(1){
        isPrime = 0;
        for(int i =0;i<nbSievePrimes;i++){
            fprintf(outputFile,"c");
            if(mpz_divisible_p(prime,sievePrimeList[i])){
                isPrime = 1;
                break;
            }
        }
        //S'il y a un diviseur dans le sieve
        if(isPrime==0){
            fprintf(outputFile,"a");
            if(mpz_probab_prime_p(prime, 10))break;
        }
        fprintf(outputFile,"b");
        mpz_add_ui(prime,prime,2);
    }
}


int main(int argc, char const *argv[]) {

    if (argc < 3) {
        printf("Usage : ./genSeq [nbPrimeInSieve] [sizeOfGeneratedRandom] \n");
        exit(1);
    }

    int nbSievePrimes = atoi(argv[1]);
    int size = atoi(argv[2]);
    srand(time(NULL));
    mpz_t p,q,N,tmpRand;
    mpz_inits(p,q,N,tmpRand,NULL);

    //Setup de la sortie
    FILE *output_p = fopen("output_sequence_p.txt", "w");
    FILE *output_q = fopen("output_sequence_q.txt", "w");
    FILE *module = fopen("module_rsa.txt","w");

    //Setup du gen random
    gmp_randstate_t randstate;
    gmp_randinit_default(randstate);
    gmp_randseed_ui(randstate, time(NULL));

    //Setup du sieve
    mpz_t * sievePrimeList = malloc(sizeof(mpz_t) * nbSievePrimes);
    find_k_first_primes(nbSievePrimes, sievePrimeList);

    genPrimeAndLog(p,size,nbSievePrimes,sievePrimeList,randstate,output_p);
    gmp_printf("Premier p : \t%Zu\n", p);
    genPrimeAndLog(q,size,nbSievePrimes,sievePrimeList,randstate,output_q);
    gmp_printf("Premier q : \t%Zu\n", q);
    mpz_clear(tmpRand);
       



    //Calcul du module
    mpz_mul(N,p,q);

    fclose(output_p);
    fclose(output_q);

    gmp_fprintf(module,"%Zd",N);
    mpz_clears(p,q,N,NULL);
    //on clear le sieve
    for(int x =0;x<nbSievePrimes;x++){
        mpz_clear(sievePrimeList[x]);
    }
    free(sievePrimeList);
    gmp_randclear(randstate);
    fclose(module);
    
}
