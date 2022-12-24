#include "gmp.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "function.h"
#define A 2
#define B 3
 //1353489595387284226541759837701112606791527571
typedef struct {
    int size;
    unsigned int array[10];
} ArgResults;


/*
Renvoie l'indice de la valeur absolue maximale de la liste
*/
ArgResults* argmax(double *L, int size) {
    double max = L[0];
    ArgResults *results= malloc(sizeof(ArgResults));
    results->size = 1;
    results->array[0] = 0;
    for (int i = 1; i < size; i++) {
        if (fabs(L[i]) > max) {
            max = fabs(L[i]);
            results->array[0] = i;
            results->size = 1;
        }
        else if(fabs(L[i]) == max){
            results->array[results->size++] = i;
        }
    }
    return results;
}

void processArgResults(mpz_t p, mpz_t *sievePrimeList, ArgResults**candidats, int nbSievePrimes,int index,unsigned int *toProcessArray)
{
    if(index == nbSievePrimes){
        chinese_remainder_theorem(p, sievePrimeList, toProcessArray, nbSievePrimes);
        gmp_printf("p possible :%Zu\n",p);
    }
    else{
        ArgResults* tmp = candidats[index];
        //printf("Debug size = %d, %p, %d\n",tmp->size,tmp,index);

        for(int i =0;i<tmp->size;i++){
            toProcessArray[index] = tmp->array[i];
     //       gmp_printf("p congru à %lu  mod %Zd \n", toProcessArray[index], sievePrimeList[index]);
            processArgResults(p,sievePrimeList,candidats,nbSievePrimes,index+1,toProcessArray);
        }
    }
}


int main(int argc, char const *argv[]) {
    if (argc < 4) {
        printf("Usage : ./analyse [filename] [sieveSize] [nbInt] \n");
        exit(1);
    }

    FILE *fptrHamAndNoise = fopen(argv[1], "r");
    int nbSievePrimes = atoi(argv[2]);
    int nb_candidats = atoi(argv[3]);
    double tmp;
    int small_prime;
    mpz_t z_h, z_m, p;
    mpz_inits(z_h, z_m, p, NULL);
    ArgResults* candidats[nbSievePrimes];
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
            //printf("%f\t",tmp);
        }
        //printf("\n");
    }
    int tmp2,totalP=1;

    for (int j = 0; j < nbSievePrimes; j++) {
        small_prime = mpz_get_ui(sievePrimeList[j]);
        double score[small_prime];
        score[0] = 0;
        for (int h = 1; h < small_prime; h++) {
            //mpz_set_ui(z_h,h);
            double hypothesis[nb_candidats];
            double comparaisons[nb_candidats];
            for (int i = 0; (i < nb_candidats); i++) {
                //m(h − (n − i − 1)τ mod sj )
                tmp2 = h - (nb_candidats-i-1)*2 ;
                mpz_set_si(z_h, tmp2);
                mpz_mod_ui(z_h, z_h, small_prime);
                hypothesis[i] = mpz_popcount(z_h)*A+B; //w(h - (n - i - 1)*2 mod sj)
                comparaisons[i] = mesures[i][j];
            }

            score[h] = correlation_coeff(hypothesis, comparaisons, nb_candidats);

        }
        candidats[j] = argmax(score, small_prime);
        gmp_printf("Congru à %Zd : %d possibilitées\t",sievePrimeList[j],candidats[j]->size);
        for(int l = 0;l<candidats[j]->size;l++){
            printf("%u: %f\t",candidats[j]->array[l],score[candidats[j]->array[l]]);
        }
        printf("\n");
        totalP*=candidats[j]->size;

    }
    printf("Total : %d possibilitées\n",totalP);
    unsigned int tab[nbSievePrimes];
    processArgResults(p,sievePrimeList,candidats,nbSievePrimes,0,tab);
    gmp_printf("p = %Zd\n", p);

}