
#include "gmp.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "../common/function.h"
typedef struct {
    int size;
    unsigned int *array;
} ArgResults;


/*
Renvoie les indices de la valeur absolue maximale de la liste
*/
unsigned int argmax(double *L, double **H, int size, double norm_comp, int nb_candidats) {
    double** possible_hypothesis = malloc(sizeof(double *)*size);
    double max = L[0];
    ArgResults *results= malloc(sizeof(ArgResults));
    results->size = 1;
    results->array = malloc(sizeof(unsigned int)*size);
    results->array[0] = 0;
    for (int i = 1; i < size; i++) {
        if (fabs(L[i]) > max) {
            max = fabs(L[i]);
            results->array[0] = i;
            results->size = 1;
            possible_hypothesis[0] = malloc(sizeof(double)*nb_candidats);
            for (int k = 0; k < nb_candidats; k++) {
                possible_hypothesis[0][k]= H[i][k];
            }
        }
        else if(fabs(L[i]) == max){
            results->array[results->size] = i;
            possible_hypothesis[results->size] = malloc(sizeof(double)*nb_candidats);
            for (int k = 0; k < nb_candidats; k++) {
                possible_hypothesis[results->size][k]= H[i][k];
            }
            results->size++;
        }
    }
    double diff_norms[results->size];
    for (int j = 0; j<results->size ; j++) {
        diff_norms[j] = fabs(euclidean_norm(possible_hypothesis[j], nb_candidats) - norm_comp);
    }
    int toReturn = results->array[argmin(diff_norms, results->size)];
    free(possible_hypothesis);
    free(results->array);
    free(results);
    return toReturn;
}


double getB(double** mesures,int sizeX,int sizeY){
    double min,total=0;
    int div = sizeX-1;
    for(int i = 0;i<sizeX-1;i++){
        min = 100;
        for(int j = 0;j<sizeY;j++){
            if(mesures[i][j]<min)min = mesures[i][j];
        }
        total+=min;
    }
    return total/div;
}

double getA(double** mesures,int sizeX,int sizeY,double b){
    int minI;
    double min;
    double total=0;
    int div = sizeX-1;
    for(int i = 0;i<sizeX-1;i++){
        min = 100;
        for(int j = 0;j<sizeY;j++){
            if(mesures[i][j]<min){
                min = mesures[i][j];
                minI = j;
            }
        }
        if(min >b) {
            div--;
            continue;
        }
        total += (mesures[i+1][minI]-b);
        //printf("%f\n",(mesures[i+1][minI]-min));
    }

    return total/(div);
}



int main(int argc, char const *argv[]) {
    if (argc < 4) {
        printf("Usage : ./analyse [filename] [sieveSize] [nbInt] \n");
        exit(1);
    }


    int nbSievePrimes = atoi(argv[2]);
    int nb_candidats = atoi(argv[3]);
    if(nb_candidats ==1){
        printf("Trop peu de candidats. Echec\n");
        exit(1);
    }
    FILE *fptrHamAndNoise = fopen(argv[1], "r");
    double norm_comp,tmp;
    int small_prime;
    mpz_t z_h, z_m, p;
    mpz_inits(z_h, z_m, p, NULL);
    unsigned int candidats[nbSievePrimes];
    double** mesures = malloc(sizeof(double*)*nb_candidats);
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
        mesures[i] = malloc(sizeof(double)*nbSievePrimes);
        for (int j = 0; j < nbSievePrimes; j++) {
            fscanf(fptrHamAndNoise, "%lf", &tmp);
            mesures[i][j] = tmp;
            //printf("%f\t",tmp);
        }
        //printf("\n");
    }

    double b = getB(mesures,nb_candidats,nbSievePrimes);
    printf("B = %f\n", b);
    double a = getA(mesures, nb_candidats, nbSievePrimes, b);
    printf("A = %f\n", a);
    int tmp2;
    double *hypothesis,*comparaisons,*score;
    for (int j = 0; j < nbSievePrimes; j++) {
        small_prime = mpz_get_ui(sievePrimeList[j]);
        score =  malloc(sizeof(double)*small_prime);
        double **H;
        H = malloc(sizeof(double*)*small_prime);
        score[0] = 0;
        for (int h = 1; h < small_prime; h++) {
            hypothesis = malloc(sizeof(double)*nb_candidats);
            comparaisons = malloc(sizeof(double)*nb_candidats);
            for (int i = 0; (i < nb_candidats); i++) {
                //m(h − (n − i − 1)τ mod sj )
                tmp2 = h - (nb_candidats-i-1)*2 ;
                mpz_set_si(z_h, tmp2);
                mpz_mod_ui(z_h, z_h, small_prime);
                hypothesis[i] = mpz_popcount(z_h)*a+b; //w(h - (n - i - 1)*2 mod sj)
                comparaisons[i] = mesures[i][j];

            }
           
            score[h] = correlation_coeff(hypothesis, comparaisons, nb_candidats);
            H[h] = malloc(sizeof(double)*nb_candidats);
            for (int k = 0; k<nb_candidats; k++){
                H[h][k] = hypothesis[k];
            }
            norm_comp = euclidean_norm(comparaisons, nb_candidats);
            free(hypothesis);
            free(comparaisons);
        }
        candidats[j] = argmax(score, H, small_prime, norm_comp, nb_candidats);
        gmp_printf("Congru à %d mod %Zd \n",candidats[j],sievePrimeList[j]);
        free(score);
       


    }
    chinese_remainder_theorem_cpa(p, sievePrimeList, candidats, nbSievePrimes);
    gmp_printf("p = %Zd\n", p);
    for (int i = 0; i < nb_candidats; i++)
        free(mesures[i]);
    free(mesures);

    //on clear le sieve
    for(int x =0;x<nbSievePrimes;x++){
        mpz_clear(sievePrimeList[x]);
    }
    free(sievePrimeList);
    mpz_clears(z_h, z_m, p, NULL);
    gmp_randclear(randstate);
}
