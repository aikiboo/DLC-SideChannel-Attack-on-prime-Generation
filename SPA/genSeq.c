#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../common/function.h"

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

    for (int j = 0; j < 2; j++) {
        //Génération du premier entier impair
        gen_k_bits_number_odd(tmpRand, size, randstate);
        int isPrime;
        while(1){
            isPrime = 0;
            for(int i =0;i<nbSievePrimes;i++){
                if (j == 0) fprintf(output_p,"c");
                else fprintf(output_q,"c");
                
                if(mpz_divisible_p(tmpRand,sievePrimeList[i])){
                    isPrime =1;
                    break;
                }
            }
            //S'il y a un diviseur dans le sieve
            if(isPrime==0){
                if (j == 0) fprintf(output_p,"a");
                else fprintf(output_q,"a");
                if(mpz_probab_prime_p(tmpRand, 10))break;

            }
            if (j == 0) fprintf(output_p,"b");
            else fprintf(output_q,"b");
            mpz_add_ui(tmpRand,tmpRand,2);
        }

        printf("Is prime (1=prob,2=100\%) : %d\n", mpz_probab_prime_p(tmpRand, 10));
        if (j == 0) {
            mpz_set(p,tmpRand);
            gmp_printf("Premier p : \t%Zu\n", p);
        }
        else {
            mpz_set(q,tmpRand);
            gmp_printf("Premier q : \t%Zu\n", q);
        }
        mpz_clear(tmpRand);
       

    }

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
