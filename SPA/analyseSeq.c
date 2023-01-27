#include "gmp.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "../common/function.h"

/*Renvoie le nombre de candidats produits par l'algorithme de génération de premier*/
int find_nb_candidats(FILE *output) {
    int sizef;
    int nb_candidats = 1;
    char tmp;

    fseek(output, 0, SEEK_END);
    sizef = ftell(output);
    fseek(output, 0, SEEK_SET);

    for (int i = 0; i < sizef; i++) {
        fscanf(output, "%c", &tmp);
        if (tmp == 'b') {
            nb_candidats++;
        }
    }
    return nb_candidats;
}

/*Renvoie la liste des petits premiers qui provoquent changement de candidats (car ils divisent le candidat actuel)*/
int find_divisors(mpz_t *divisors, mpz_t *congruences, mpz_t *divisorForEachCandidate, mpz_t *sievePrimeList,
                  int nb_candidats, FILE *output) {

    int nbDivisors;
    int ind_divisor_prime = -1;
    char tmp;

    mpz_t congru;

    fseek(output, 0, SEEK_SET);
    int j = 0;

    while (j < nb_candidats) {
        fscanf(output, "%c", &tmp);
        if (tmp == 'c') {
            ind_divisor_prime++;
        }
        if (tmp == 'b') {
            mpz_init(congru);
            mpz_set(divisors[j], sievePrimeList[ind_divisor_prime]);
            mpz_set(divisorForEachCandidate[j], sievePrimeList[ind_divisor_prime]);
            mpz_set_ui(congru, nb_candidats);
            mpz_sub_ui(congru, congru, j + 1);
            mpz_mul_ui(congru, congru, 2);
            mpz_mod(congruences[j], congru, divisors[j]);
            j++;
            ind_divisor_prime = -1;
            mpz_clear(congru);
        }
        if (tmp == 'a') {
            fscanf(output, "%c", &tmp);
            ind_divisor_prime = -1;
            j++;
        }
    }

    //Suppression des doublons
    nbDivisors = nb_candidats;
    for (int i = 0; i < nbDivisors; i++) {
        for (int j = i + 1; j < nbDivisors;) {
            if (mpz_cmp(divisors[j], divisors[i]) == 0) {
                for (int k = j; k < nbDivisors - 1; k++) {
                    mpz_set(divisors[k], divisors[k + 1]);
                    mpz_set(congruences[k], congruences[k + 1]);
                }
                nbDivisors--;
            } else {
                j++;
            }
        }
    }

    //Suppression des zéros
    int izero = 0;
    while (mpz_cmp_ui(divisors[izero], 0) != 0) {
        izero++;
    }
    for (int i = izero; i < nbDivisors - 1; i++) {
        mpz_set(divisors[i], divisors[i + 1]);
        mpz_set(congruences[i], congruences[i + 1]);
    }
    //nbDivisors--;
    mpz_set_ui(divisors[nbDivisors - 1], 2);
    mpz_set_ui(congruences[nbDivisors - 1], 1);
    return nbDivisors;
}

/*Renvoie la liste des petits premiers qui ne sont pas dans la liste des diviseurs*/
void find_not_divisors(mpz_t *not_divisors, mpz_t *sievePrimeList, mpz_t *divisors, int nbSievePrimes, int nbDivisors) {
    int k = 0;
    int n = 0;
    int m;
    while (n < nbSievePrimes) {
        m = 0;
        while (m < nbDivisors && mpz_cmp(sievePrimeList[n], divisors[m]) != 0) {
            m++;
        }
        if (m == nbDivisors) {
            mpz_init(not_divisors[k]);
            mpz_set(not_divisors[k], sievePrimeList[n]);
            k++;
        }
        n++;
    }
}

/*Renvoie la liste des éléments auquel le premier cherché n'est pas congru modulo un élément de la liste des non diviseurs*/
int not_congruent(mpz_t *notCongruent, mpz_t *divisorForEachCandidate, unsigned long int nd, int nb_candidats) {
    mpz_init(notCongruent[0]);
    int sizeNotCongruent = 1;
    for (int i = 1; i < nb_candidats; i++) {
        if (mpz_cmp_ui(divisorForEachCandidate[i], nd) > 0 || mpz_cmp_ui(divisorForEachCandidate[i], 0) == 0) {
            mpz_init(notCongruent[sizeNotCongruent]);
            mpz_set_ui(notCongruent[sizeNotCongruent], 2 * (nb_candidats - i));
            mpz_mod_ui(notCongruent[sizeNotCongruent], notCongruent[sizeNotCongruent], nd);
            sizeNotCongruent++;
        }
    }
    return sizeNotCongruent;
}

/*Renvoie la liste des élement auxquels le premier cherché est congru modulo un élément de la liste des non diviseurs*/
void find_hypothesis(mpz_t *hyp, mpz_t *NotCongruent, int sizeNotCongruent, unsigned long int nd) {
    int k = 0;
    int m;
    unsigned int n = 0;

    while (n < nd) {
        m = 0;
        while (m < sizeNotCongruent && mpz_cmp_ui(NotCongruent[m], n) != 0) {
            m++;
        }
        if (m == sizeNotCongruent) {
            mpz_init(hyp[k]);
            mpz_set_ui(hyp[k], n);
            k++;
        }
        n++;
    }
}

void
find_missing_bits(mpz_t p2, int nbMissingBits, mpz_t *not_divisors, mpz_t *divisorForEachCandidate, mpz_t N, mpz_t p,
                  mpz_t s, int nbNotDivisors, int nb_candidats) {
    //Recherche d'un non diviseur qui nous apporte assez d'infos
    int c = nbNotDivisors - 1;
    while (nbMissingBits < mpz_sizeinbase(not_divisors[c], 2) && c > 0) {
        c--;
    }

    unsigned long int nd;
    mpz_t s2;

    while (mpz_divisible_p(N, p2) == 0 && c <= nbNotDivisors) {

        mpz_t *D = malloc(sizeof(mpz_t) * 2), *C = malloc(sizeof(mpz_t) * 2);
        mpz_t gcd;
        mpz_init(gcd);

        mpz_inits(D[0], D[1], C[0], C[1], p2, NULL);

        mpz_set(D[0], s);
        mpz_set(C[0], p);

        mpz_gcd(gcd, s, not_divisors[c]);

        if (mpz_cmp_ui(gcd, 1) == 0) {
            nd = mpz_get_ui(not_divisors[c]);

            //Calcul de ce à quoi p n'est pas congru modulo ce non diviseurs
            mpz_t *notCongruent;
            notCongruent = malloc(sizeof(mpz_t) * nd);
            int sizeNotCongruent = not_congruent(notCongruent, divisorForEachCandidate, nd, nb_candidats);

            //Hypothèses sur p modulo le non diviseur
            mpz_t *hyp;
            int sizeHyp = nd - sizeNotCongruent;
            hyp = malloc(sizeof(mpz_t) * sizeHyp);
            find_hypothesis(hyp, notCongruent, sizeNotCongruent, nd);

            int h = 0;
            mpz_set_ui(D[1], nd);
            mpz_init(s2);
            while (mpz_divisible_p(N, p2) == 0 && h < sizeHyp) {
                mpz_set(C[1], hyp[h]);
                mpz_init(p2);
                chinese_remainder_theorem_spa(p2, s2, D, C, 2);
                h++;
            }
            if (mpz_divisible_p(N, p2) != 0) {
                return;
            }
        }
        c++;
    }
}

void clearMPZArray(int size, mpz_t *array){
    for(int i =0;i<size;i++)mpz_clear(array[i]);
    free(array);

}

void clearMemSieve(int nbSievePrimes,mpz_t *sievePrimeList){
        //on clear le sieve
        for(int x =0;x<nbSievePrimes;x++){
            mpz_clear(sievePrimeList[x]);
        }
        free(sievePrimeList);
}
void clearMemDivisors(int nb_candidatsP,mpz_t *divisors,mpz_t *congruences,mpz_t *divisorForEachCandidats){
    for (int i = 0; i < nb_candidatsP; i++) {
        mpz_clears(divisors[i], congruences[i], divisorForEachCandidats[i], NULL);
    }
    free(divisors);
    free(congruences);
    free(divisorForEachCandidats);
    }

int main(int argc, char const *argv[]) {
    if (argc < 5) {
        printf("Usage : ./analyse [filename output p] [filename output q] [filename module] [sieveSize] \n");
        exit(1);
    }

    FILE *module = fopen(argv[3], "r");
    int nbSievePrimes = atoi(argv[4]);

    //Setup du sieve
    mpz_t *sievePrimeList;
    sievePrimeList = malloc(sizeof(mpz_t) * nbSievePrimes);
    find_k_first_primes(nbSievePrimes, sievePrimeList);

    //Récupération du module
    mpz_t N;
    int sizeN;
    mpz_init(N);
    gmp_fscanf(module, "%Zd\n", N);
    sizeN = mpz_sizeinbase(N, 2);
    fclose(module);

    printf("on tente de retrouver p\n");
    FILE *outputP = fopen(argv[1], "r");

    int nb_candidatsP = find_nb_candidats(outputP);
    printf("nb_candidats = %d\n", nb_candidatsP);

    //Récupération des diviseurs
    mpz_t *divisorsP;
    divisorsP = malloc(sizeof(mpz_t) * nb_candidatsP);
    mpz_t *congruencesP;
    congruencesP = malloc(sizeof(mpz_t) * nb_candidatsP);
    mpz_t *divisorForEachCandidateP;
    divisorForEachCandidateP = malloc(sizeof(mpz_t) * nb_candidatsP);

    for (int i = 0; i < nb_candidatsP; i++) {
        mpz_inits(divisorsP[i], congruencesP[i], divisorForEachCandidateP[i], NULL);
    }

    int nbDivisorsP = find_divisors(divisorsP, congruencesP, divisorForEachCandidateP, sievePrimeList, nb_candidatsP,
                                    outputP);
    fclose(outputP);

    mpz_t ap, sp;
    mpz_inits(ap, sp, NULL);
    //On cherche P avec le CRT
    chinese_remainder_theorem_spa(ap, sp, divisorsP, congruencesP, nbDivisorsP);
    //Si on a trouvé P (N divisible par P et P != 1)
    if (mpz_divisible_p(N, ap) != 0 && mpz_cmp_ui(ap, 1) != 0) { //On a trouvé p
        mpz_t q;
        mpz_init(q);
        mpz_divexact(q, N, ap);
        gmp_printf("p = %Zd\n", ap);
        gmp_printf("q = %Zd\n", q);
        mpz_clear(q);
    }
    else {
        int nb_candidatsQ, nbDivisorsQ;
        mpz_t aq, inv_ap, sq, bq, *divisorsQ, *congruencesQ, *divisorForEachCandidateQ;

        mpz_inits(aq, inv_ap, sq, bq, NULL);

        mpz_invert(inv_ap, ap, sp);
        mpz_mul(aq, inv_ap, N);
        mpz_mod(aq, aq, sp);

        printf("On tente de retrouver q\n");
        FILE *outputQ = fopen(argv[2], "r");

        nb_candidatsQ = find_nb_candidats(outputQ);
        printf("nb_candidats = %d\n", nb_candidatsQ);

        //Récupération des diviseurs
        divisorsQ = malloc(sizeof(mpz_t) * nb_candidatsQ), congruencesQ = malloc(
                sizeof(mpz_t) * nb_candidatsQ), divisorForEachCandidateQ = malloc(sizeof(mpz_t) * nb_candidatsQ);

        for (int i = 0; i < nb_candidatsQ; i++) {
            mpz_inits(divisorsQ[i], congruencesQ[i], divisorForEachCandidateQ[i], NULL);
        }

        nbDivisorsQ = find_divisors(divisorsQ, congruencesQ, divisorForEachCandidateQ, sievePrimeList, nb_candidatsQ,
                                    outputQ);
        fclose(outputQ);

        chinese_remainder_theorem_spa(bq, sq, divisorsQ, congruencesQ, nbDivisorsQ);

        if (mpz_divisible_p(N, bq) != 0 && mpz_cmp_ui(bq, 1) != 0) { //On a trouvé q
            mpz_t p;
            mpz_init(p);
            mpz_divexact(p, N, bq);
            gmp_printf("p = %Zd\n", p);
            gmp_printf("q = %Zd\n", bq);
            mpz_clear(p);
        }
        else {   //On cherche à retrouver p mod ppcm(sp,sq) ou q mod ppcm(sp,sq)

            mpz_t bp, inv_bq, p, q, s, lcmS, a, pgcd, *S = malloc(sizeof(mpz_t) * 2), *P = malloc(
                    sizeof(mpz_t) * 2), *Q = malloc(sizeof(mpz_t) * 2);

            mpz_inits(bp, inv_bq, S[0], S[1], P[0], P[1], Q[0], Q[1], p, q, s, lcmS, a, pgcd, NULL);

            mpz_invert(inv_bq, bq, sq);
            mpz_mul(bp, inv_bq, N);
            mpz_mod(bp, bp, sq);
            mpz_lcm(lcmS, sp, sq);

            mpz_set(P[0], ap);
            mpz_set(P[1], bp);
            mpz_set(Q[0], aq);
            mpz_set(Q[1], bq);
            mpz_set(S[0], sp);
            mpz_divexact(S[1], lcmS, sp);

            mpz_mod(P[0], P[0], S[0]);
            mpz_mod(P[1], P[1], S[1]);

            chinese_remainder_theorem_spa(p, s, S, P, 2);
            gmp_printf("size s = %d\n", mpz_sizeinbase(s, 2));

            if (mpz_divisible_p(N, p) != 0 && mpz_cmp_ui(p, 1) != 0) {
                mpz_divexact(q, N, p);
                gmp_printf("p = %Zd\n", p);
                gmp_printf("q = %Zd\n", q);
            }
            else {
                mpz_mod(Q[0], Q[0], S[0]);
                mpz_mod(Q[1], Q[1], S[1]);

                chinese_remainder_theorem_spa(q, s, S, Q, 2);

                if (mpz_divisible_p(N, q) != 0 && mpz_cmp_ui(q, 1) != 0) {
                    mpz_divexact(p, N, q);
                    gmp_printf("p = %Zd\n", p);
                    gmp_printf("q = %Zd\n", q);
                }
                else {

                    //Calcul du nombre de bits manquants
                    int sizeP = mpz_sizeinbase(p, 2);
                    float nbMissingBitsP = floor(sizeN / 2) - sizeP;
                    printf("nombre de bits manquants pour p : %f\n", nbMissingBitsP);

                    //Calcul des petits premiers non diviseurs
                    int nbNotDivisorsP = nbSievePrimes - nbDivisorsP;
                    mpz_t *not_divisorsP = malloc(sizeof(mpz_t) * nbNotDivisorsP);
                    find_not_divisors(not_divisorsP, sievePrimeList, divisorsP, nbSievePrimes, nbDivisorsP);

                    int bP;
                    bP = mpz_sizeinbase(not_divisorsP[nbNotDivisorsP - 1], 2);
                    if (nbMissingBitsP > bP) {

                        printf("trop de bits manquants pour p, on tente de retrouver q\n");

                        //Calcul du nombre de bits manquants
                        int sizeQ = mpz_sizeinbase(q, 2);
                        float nbMissingBitsQ = floor(sizeN / 2) - sizeQ;
                        printf("nombre de bits manquants pour q : %f\n", nbMissingBitsQ);

                        //Calcul des petits premiers non diviseurs
                        int nbNotDivisorsQ = nbSievePrimes - nbDivisorsQ;
                        mpz_t *not_divisorsQ = malloc(sizeof(mpz_t) * nbNotDivisorsQ);
                        find_not_divisors(not_divisorsQ, sievePrimeList, divisorsQ, nbSievePrimes, nbDivisorsQ);

                        int bQ;
                        bQ = mpz_sizeinbase(not_divisorsQ[nbNotDivisorsQ - 1], 2);

                        if (nbMissingBitsQ > bQ) {
                            printf("échec de l'attaque\n");
                            exit(1);
                        }
                        else{
                            printf("on tente de retrouver les bits manquants de q\n");

                            mpz_t q2;
                            mpz_init(q2);

                            find_missing_bits(q2, nbMissingBitsQ, not_divisorsQ, divisorForEachCandidateQ, N, q, s,
                                              nbNotDivisorsQ, nb_candidatsQ);

                            if (mpz_divisible_p(N, q2) != 0) {
                                mpz_divexact(p, N, q2);
                                gmp_printf("p = %Zd\n", p);
                                gmp_printf("q = %Zd\n", q2);
                                exit(1);
                            }
                            mpz_clear(q2);
                            printf("échec de l'attaque");
                        }
                        clearMPZArray(nbNotDivisorsQ,not_divisorsQ);

                    } else
                    {
                        printf("on tente de retrouver les bits manquants de p\n");

                        mpz_t p2;
                        mpz_init(p2);

                        find_missing_bits(p2, nbMissingBitsP, not_divisorsP, divisorForEachCandidateP, N, p, s,
                                          nbNotDivisorsP, nb_candidatsP);

                        if (mpz_divisible_p(N, p2) != 0) {
                            mpz_divexact(q, N, p2);
                            gmp_printf("p = %Zd\n", p2);
                            gmp_printf("q = %Zd\n", q);
                            exit(1);
                        }
                        mpz_clear(p2);
                        printf("échec de l'attaque\n");
                    }

                    clearMPZArray(nbNotDivisorsP,not_divisorsP);
                }
            }

            mpz_clears(bp, inv_bq, S[0], S[1], P[0], P[1], Q[0], Q[1], p, q, s, lcmS, a, pgcd, NULL);
            free(S);
            free(P);
            free(Q);
        }
        for (int i = 0; i < nb_candidatsQ; i++) {
            mpz_clears(divisorsQ[i], congruencesQ[i], divisorForEachCandidateQ[i], NULL);
        }
        free(divisorsQ);
        free(congruencesQ);
        free(divisorForEachCandidateQ);
        mpz_clears(aq, inv_ap, sq, bq, NULL);
    }


    clearMemSieve(nbSievePrimes,sievePrimeList);
    clearMemDivisors(nb_candidatsP,divisorsP,congruencesP,divisorForEachCandidateP);
    mpz_clears(ap, sp,N, NULL);

}
