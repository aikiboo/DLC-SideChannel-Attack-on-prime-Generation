#include "gmp.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "function.h"
#define A 2
#define B 3
typedef struct {
    int size;
    unsigned int array[10];
} ArgResults;


/*
Renvoie les indices de la valeur absolue maximale de la liste
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
/**
 * Utilise une liste de liste de candidats possibles afin de calculer tous les résultats possibles.
 * @param p
 * @param sievePrimeList
 * @param candidats
 * @param nbSievePrimes
 * @param index
 * @param toProcessArray
 */
void processArgResults(mpz_t p, mpz_t *sievePrimeList, ArgResults**candidats, int nbSievePrimes,int index,unsigned int *toProcessArray)
{
    if(index == nbSievePrimes){
        chinese_remainder_theorem(p, sievePrimeList, toProcessArray, nbSievePrimes);
        //if(mpz_sizeinbase(p,2) == 256)
        gmp_printf("p possible :%Zu\n",p);
        //gmp_printf("size : %zu\n",mpz_sizeinbase(p,2));
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

double getB(double** mesures,int sizeX,int sizeY){
    double min,total=0;
    double tmpMean;
    int div = sizeX-1;
    for(int i = 0;i<sizeX-1;i++){
        min = 100;
        for(int j = 0;j<sizeY;j++){
            if(mesures[i][j]<min)min = mesures[i][j];
        }
        if(i>5 && (min<tmpMean*0.33 || min>tmpMean*3)) {
            div --;
            continue;
        };

        tmpMean = total/(i+1-(sizeX-1-div));
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

    FILE *fptrHamAndNoise = fopen(argv[1], "r");
    int nbSievePrimes = atoi(argv[2]);
    int nb_candidats = atoi(argv[3]);
    double tmp;
    int small_prime;
    mpz_t z_h, z_m, p;
    mpz_inits(z_h, z_m, p, NULL);
    ArgResults* candidats[nbSievePrimes];
    double** mesures;
    mesures = malloc(sizeof(double*)*nb_candidats);
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
    printf("A = %f\n", getA(mesures,nb_candidats,nbSievePrimes,b));
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