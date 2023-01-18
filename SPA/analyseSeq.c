#include "gmp.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "../common/function.h"

int main(int argc, char const *argv[]) {
     if (argc < 3) {
        printf("Usage : ./analyse [filename] [sieveSize] \n");
        exit(1);
     }

    FILE *output = fopen(argv[1], "r");
    int nbSievePrimes = atoi(argv[2]);
    mpz_t *sievePrimeList;
    mpz_t p;
    int nb_candidats = 1;
    int nb_divisors;
    int taillef;
    char tmp;
    int ind_divisor_prime = -1;

    sievePrimeList = malloc(sizeof(mpz_t) * nbSievePrimes);
    find_k_first_primes(nbSievePrimes, sievePrimeList);


    fseek(output, 0, SEEK_END);
    taillef = ftell(output);
    fseek(output,0,SEEK_SET);


    /*Récupération du nombre de candidats*/
    for(int i = 0; i<taillef; i++) {
    	fscanf(output, "%c", &tmp);
    	if (tmp == 'b'){
    		nb_candidats++;
    	}
    }

    printf("nb_candidats = %d\n", nb_candidats);

    /* Récupération des diviseurs */
    mpz_t *divisors;
    divisors = malloc(sizeof(mpz_t) * nb_candidats);
    mpz_t *congruences;
    congruences = malloc(sizeof(mpz_t) * nb_candidats);
    mpz_t congru;

    fseek(output,0,SEEK_SET);
    int j = 1;
    mpz_inits(divisors[0], congruences[0],NULL);
    mpz_set_ui(divisors[0],2);
    mpz_set_ui(congruences[0],1);


    while(j < nb_candidats){
    	fscanf(output, "%c", &tmp);
    	if (tmp == 'c') {
    		ind_divisor_prime++;
    	}
    	if (tmp == 'b') {
	    	mpz_inits(divisors[j],congruences[j], congru, NULL);
	    	mpz_set(divisors[j],sievePrimeList[ind_divisor_prime]);
	    	mpz_set_ui(congru, nb_candidats);
	    	mpz_sub_ui(congru, congru, j);
	    	mpz_mul_ui(congru, congru, 2);
	    	mpz_mod(congruences[j], congru ,divisors[j]);
	    	j++;
	    	ind_divisor_prime = -1;
    		}
    	if (tmp == 'a') {
    		fscanf(output, "%c", &tmp);
    		ind_divisor_prime = -1;
    		j++;			
    	}
    }

    /*Suppression des doublons et des zéros*/
    nb_divisors = nb_candidats;
    for (int i = 0; i < nb_divisors; i++) {
    	for (int j = i + 1; j < nb_divisors;) {
    		if (mpz_cmp(divisors[j],divisors[i]) == 0) {
    			for (int k = j; k < nb_divisors; k++) {
    				mpz_set(divisors[k],divisors[k+1]);
    				mpz_set(congruences[k], congruences[k+1]);
    			}
    			nb_divisors--;
    		} else {
    			j++;
    		}
    	}
    }

    int izero = 0;
    while (mpz_cmp_ui(divisors[izero],0) != 0) {
    	izero++;
    }
    for(int i = izero; i<nb_divisors-1; i++) {
    	mpz_set(divisors[i],divisors[i+1]);
    	mpz_set(congruences[i],congruences[i+1]);
    }
    nb_divisors--;
    mpz_init(p);
    chinese_remainder_theorem_spa(p, divisors, congruences, nb_divisors);
 	  gmp_printf("p = %Zd\n", p);
}
