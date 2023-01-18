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
   fseek(output,0,SEEK_SET);

   for(int i = 0; i<sizef; i++) {
      fscanf(output, "%c", &tmp);
      if (tmp == 'b'){
         nb_candidats++;
      }
   }
   printf("nb_candidats = %d\n", nb_candidats); 
   return nb_candidats;
}

/*Renvoie la liste des petits premiers qui provoquent changement de candidats (car ils divisent le candidat actuel)*/
int find_divisors(mpz_t *divisors, mpz_t *congruences, mpz_t *divisorForEachCandidate, mpz_t *sievePrimeList, int nb_candidats, FILE *output) {
   int nbDivisors;
   int ind_divisor_prime = -1;
   char tmp;

   mpz_t congru;
    
   fseek(output,0,SEEK_SET);
   int j = 1;
   mpz_set_ui(divisors[0],2);
   mpz_set_ui(congruences[0],1);


   while(j < nb_candidats){
      fscanf(output, "%c", &tmp);
      if (tmp == 'c') {
         ind_divisor_prime++;
      }
      if (tmp == 'b') {
         mpz_init(congru);
         mpz_set(divisors[j],sievePrimeList[ind_divisor_prime]);
         mpz_set(divisorForEachCandidate[j],sievePrimeList[ind_divisor_prime]);
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
    
   //Suppression des doublons
    nbDivisors = nb_candidats;
    for (int i = 0; i < nbDivisors; i++) {
       for (int j = i + 1; j < nbDivisors;) {
          if (mpz_cmp(divisors[j],divisors[i]) == 0) {
             for (int k = j; k < nbDivisors-1; k++) {
                mpz_set(divisors[k],divisors[k+1]);
                mpz_set(congruences[k], congruences[k+1]);
             }
          nbDivisors--;
          } else {
             j++;
          }
       }
    }

   //Suppression des zéros
   int izero = 0;
   while (mpz_cmp_ui(divisors[izero],0) != 0) {
      izero++;
   }
   for(int i = izero; i<nbDivisors-1; i++) {
      mpz_set(divisors[i],divisors[i+1]);
      mpz_set(congruences[i],congruences[i+1]);
   }
   nbDivisors--;
   return nbDivisors;
}


/*Renvoie la liste des petits premiers qui ne sont pas dans la liste des diviseurs*/
void find_not_divisors(mpz_t *not_divisors, mpz_t *sievePrimeList, mpz_t *divisors, int nbSievePrimes, int nbDivisors) {
   int k = 0;
   int n = 0;
   int m;
   while(n < nbSievePrimes) {
      m = 0;
      while(m < nbDivisors && mpz_cmp(sievePrimeList[n],divisors[m]) != 0) {
         m++;
      }
      if (m == nbDivisors) {
         mpz_init(not_divisors[k]);
         mpz_set(not_divisors[k],sievePrimeList[n]);
         k++;
      }
      n++;
   }
}

/*Renvoie la liste des éléments auquel le premier cherché n'est pas congru modulo un élément de la liste des non diviseurs*/
int not_congruent(mpz_t *notCongruent, mpz_t *divisorForEachCandidate, unsigned long int nd, int nb_candidats) {
   mpz_init(notCongruent[0]);
   int sizeNotCongruent = 1;
   for(int i = 1; i<nb_candidats; i++) {
      if(mpz_cmp_ui(divisorForEachCandidate[i], nd) > 0 || mpz_cmp_ui(divisorForEachCandidate[i],0) == 0) {
         mpz_init(notCongruent[sizeNotCongruent]);
         mpz_set_ui(notCongruent[sizeNotCongruent], 2*(nb_candidats-i));
         mpz_mod_ui(notCongruent[sizeNotCongruent],notCongruent[sizeNotCongruent],nd);
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
 
   while(n < nd) {
      m = 0;
      while(m < sizeNotCongruent && mpz_cmp_ui(NotCongruent[m],n) != 0) {
         m++;
      }
      if(m == sizeNotCongruent) {
         mpz_init(hyp[k]);
         mpz_set_ui(hyp[k],n);
         k++;
      }
      n++;
   }
}

int main(int argc, char const *argv[]) {
   if (argc < 5) {
      printf("Usage : ./analyse [filename output p] [filename output q] [filename module] [sieveSize] \n");
      exit(1);
   }

   FILE *module = fopen(argv[3], "r");
   int nbSievePrimes = atoi(argv[4]);

   //Setup du sieve
   mpz_t *sievePrimeList = malloc(sizeof(mpz_t) * nbSievePrimes);
   find_k_first_primes(nbSievePrimes, sievePrimeList);

   //Récupération du module
   mpz_t N;
   int sizeN;
   mpz_init(N);
   gmp_fscanf(module, "%Zd\n", N);
   sizeN = mpz_sizeinbase(N,2);
   fclose(module);

   printf("on tente de retrouver p\n");
   FILE *outputP = fopen(argv[1], "r");
   mpz_t ap;
    
   int nb_candidatsP = find_nb_candidats(outputP);

   //Récupération des diviseurs
   mpz_t *divisorsP;
   divisorsP = malloc(sizeof(mpz_t) * nb_candidatsP);
   mpz_t *congruencesP;
   congruencesP = malloc(sizeof(mpz_t) * nb_candidatsP);
   mpz_t *divisorForEachCandidateP;
   divisorForEachCandidateP = malloc(sizeof(mpz_t) * nb_candidatsP);
   
   for (int i = 0; i < nb_candidatsP; i++){
      mpz_inits(divisorsP[i], congruencesP[i], divisorForEachCandidateP[i],NULL);
   } 

    
   int nbDivisorsP = find_divisors(divisorsP, congruencesP, divisorForEachCandidateP, sievePrimeList, nb_candidatsP, outputP);
   fclose(outputP);
   
   mpz_init(ap);
   mpz_t sp;
   mpz_init(sp);

   chinese_remainder_theorem_spa(ap, sp, divisorsP, congruencesP, nbDivisorsP);

   if (mpz_divisible_p(N,ap) != 0 && mpz_cmp_ui(ap,1) != 0) { //On a trouvé p
      mpz_t q;
      mpz_init(q);
      
      mpz_divexact(q,N,ap);
      gmp_printf("p = %Zd\n", ap);
      gmp_printf("q = %Zd\n", q);     
   }
   else {
      
      mpz_t aq, inv_ap, bq, sq;
      
      mpz_inits(aq, inv_ap, NULL, bq, sq);
      
      mpz_invert(inv_ap, ap, sp);
      mpz_mul(aq,inv_ap,N);
      mpz_mod(aq,aq,sp);


      printf("On tente de retrouver q\n");
      FILE *outputQ = fopen(argv[2], "r");
          
      int nb_candidatsQ = find_nb_candidats(outputQ);

      //Récupération des diviseurs
      mpz_t *divisorsQ = malloc(sizeof(mpz_t) * nb_candidatsQ),
            *congruencesQ = malloc(sizeof(mpz_t) * nb_candidatsQ),
            *divisorForEachCandidateQ = malloc(sizeof(mpz_t) * nb_candidatsQ);

      for (int i = 0; i < nb_candidatsQ; i++){
         mpz_inits(divisorsQ[i], congruencesQ[i], divisorForEachCandidateQ[i],NULL);
      } 
  
      int nbDivisorsQ = find_divisors(divisorsQ, congruencesQ, divisorForEachCandidateQ, sievePrimeList, nb_candidatsQ, outputQ);
      
      fclose(outputQ);

      chinese_remainder_theorem_spa(bq, sq, divisorsQ, congruencesQ, nbDivisorsQ);

      if (mpz_divisible_p(N,bq) != 0 && mpz_cmp_ui(bq,1) != 0) { //On a trouvé q
         mpz_t p;
         mpz_init(p);
         
         mpz_divexact(p,N,bq);
         gmp_printf("p = %Zd\n", p);
         gmp_printf("q = %Zd\n", bq);

      }
      else {

         //On cherche à retrouver p mod ppcm(sp,sq) ou q mod ppcm(sp,sq)
         int i = 0;
         mpz_t p, q, s, lcmS, a, pgcd, bp, inv_bq,
               *S = malloc(sizeof(mpz_t)*2),
               *P = malloc(sizeof(mpz_t)*2),
               *Q = malloc(sizeof(mpz_t)*2);
         
         mpz_inits(bp, inv_bq, S[0], S[1], P[0], P[1], Q[0], Q[1], p, q, s, lcmS, NULL, a, pgcd);

         mpz_invert(inv_bq, bq, sq);
         mpz_mul(bp,inv_bq,N);
         mpz_mod(bp,bp,sq);  
         
         mpz_lcm(lcmS,sp,sq);
         
         mpz_set_ui(S[1],1);
         mpz_set(P[0],ap);
         mpz_set(P[1],bp);
         mpz_set(Q[0],aq);    
         mpz_set(Q[1],bq);
         mpz_set(S[0], lcmS);
         mpz_set(a,S[0]);
         mpz_set_ui(S[1],2);
         
         mpz_gcd(pgcd,a,S[1]);
         
         while(mpz_cmp_ui(pgcd,1) != 0 && i < nbDivisorsP) {
            mpz_divexact(a,S[0],S[1]);
            mpz_set(S[1],divisorsP[i]);
            mpz_gcd(pgcd,a,S[1]);
            i++;
         }
               
         mpz_set(S[0],a);
           
         chinese_remainder_theorem_spa(p, s, S, P, 2);
   
         if (mpz_divisible_p(N,p) != 0) {
            mpz_divexact(q,N,p);
            gmp_printf("p = %Zd\n", p);
            gmp_printf("q = %Zd\n", q);
         }
         else {
            chinese_remainder_theorem_spa(q, s, S, Q, 2);
            if (mpz_divisible_p(N,q) != 0) {
               mpz_divexact(p,N,q);
               gmp_printf("p = %Zd\n", p);
               gmp_printf("q = %Zd\n", q);
            }
            else {
               
               //Calcul du nombre de bits manquants
               int sizeP = mpz_sizeinbase(p,2);
               float nbMissingBitsP = sizeN - ceil(sizeN/2) - sizeP;
               printf("missing bits = %f\n", nbMissingBitsP);

               //Calcul des petits premiers non diviseurs
               int nbNotDivisorsP = nbSievePrimes - nbDivisorsP;
               mpz_t *not_divisorsP = malloc(sizeof(mpz_t) * nbNotDivisorsP);
               find_not_divisors(not_divisorsP, sievePrimeList, divisorsP, nbSievePrimes, nbDivisorsP);

               int bP;
               bP = mpz_sizeinbase(not_divisorsP[nbNotDivisorsP-1],2);
                     
               if (nbMissingBitsP > bP) {
                  
                  printf("trop de bits manquants pour p, on tente de retrouver q\n");
                  
                  //Calcul du nombre de bits manquants
                  int sizeQ = mpz_sizeinbase(q,2);
                  float nbMissingBitsQ = sizeN - ceil(sizeN/2) - sizeQ;
                  printf("missing bits = %f\n", nbMissingBitsQ);

                  //Calcul des petits premiers non diviseurs
                  int nbNotDivisorsQ = nbSievePrimes - nbDivisorsQ;
                  mpz_t *not_divisorsQ = malloc(sizeof(mpz_t) * nbNotDivisorsQ);

                  find_not_divisors(not_divisorsQ, sievePrimeList, divisorsQ, nbSievePrimes, nbDivisorsQ);
                  
                  int bQ;
                  bQ = mpz_sizeinbase(not_divisorsQ[nbNotDivisorsQ-1],2);
               
                  
                  if(nbMissingBitsQ > bQ) {
                     printf("échec de l'attaque\n");
                     exit(1);
                  }
                  
                  else {
                     printf("on tente de retrouver les bits manquants de q\n");
                     
                     //Recherche d'un non diviseur qui nous apporte assez d'infos
                     int c = nbNotDivisorsQ-1;
                     while(nbMissingBitsQ > mpz_sizeinbase(not_divisorsQ[c],2)-1) {
                        c--;
                     }
                           
                     unsigned long int nd;
                     mpz_t s2;

                     while(mpz_divisible_p(N,q) == 0 && c <= nbNotDivisorsP) {
                        mpz_init(pgcd);
                        mpz_gcd(pgcd, s,not_divisorsQ[c]);
                        mpz_t *D;
                        D = malloc(sizeof(mpz_t) * 2);
                        mpz_inits(D[0], D[1], NULL);
                        mpz_t *C;
                        C = malloc(sizeof(mpz_t) * 2);
                        mpz_inits(C[0], C[1], NULL);
                        mpz_set(D[0],s);
                        mpz_set(C[0],q);
                        if (mpz_cmp_ui(pgcd,1) == 0) {
                           nd = mpz_get_ui(not_divisorsQ[c]);

                           //Calcul de ce à quoi p n'est pas congru modulo ce non diviseurs
                           mpz_t *notCongruentQ;
                           notCongruentQ = malloc(sizeof(mpz_t) * nd);
                           int sizeNotCongruentQ = not_congruent(notCongruentQ, divisorForEachCandidateQ, nd, nb_candidatsQ);

                           //Hypothèses sur p modulo le non diviseur
                           mpz_t *hypQ;
                           int sizeHypQ = nd - sizeNotCongruentQ;
                           hypQ = malloc(sizeof(mpz_t) * sizeHypQ);
                           find_hypothesis(hypQ, notCongruentQ, sizeNotCongruentQ, nd);
                           
                           int a = 0;
                           mpz_set_ui(D[1],nd);
                           mpz_init(s2);
                           while(mpz_divisible_p(N,q) == 0 && a<sizeHypQ) {
                              mpz_set(C[1],hypQ[a]);
                              mpz_init(q);
                              chinese_remainder_theorem_spa(q, s2, D, C, 2);
                              a++;
                           }
                           if (mpz_divisible_p(N,q) != 0){
                              mpz_divexact(p,N,q);
                              gmp_printf("p = %Zd\n", p);
                              gmp_printf("q = %Zd\n", q);
                              exit(1);
                           }
                        } 
                        c++;
                     }
                     printf("échec de l'attaque");
                  }
               }
               else {
                  printf("on tente de retrouver les bits manquants de p\n");
                  
                  //Recherche d'un non diviseur qui nous apporte assez d'infos
                  int c = nbNotDivisorsP-1;
                  while(nbMissingBitsP > mpz_sizeinbase(not_divisorsP[c],2)-1) {
                     c--;
                  }

                  unsigned long int nd;
                  mpz_t s2;

                  while(mpz_divisible_p(N,p) == 0 && c <= nbNotDivisorsP) {
                     mpz_init(pgcd);
                     mpz_gcd(pgcd,s,not_divisorsP[c]);
                     mpz_t *D;
                     D = malloc(sizeof(mpz_t) * 2);
                     mpz_inits(D[0], D[1], NULL);
                     mpz_t *C;
                     C = malloc(sizeof(mpz_t) * 2);
                     mpz_inits(C[0], C[1], NULL);
                     mpz_set(D[0],s);
                     mpz_set(C[0],p);
                     if (mpz_cmp_ui(pgcd,1) == 0) {
                        nd = mpz_get_ui(not_divisorsP[c]);

                        //Calcul de ce à quoi p n'est pas congru modulo ce non diviseurs
                        mpz_t *notCongruentP  = malloc(sizeof(mpz_t) * nd);
                        int sizeNotCongruentP = not_congruent(notCongruentP, divisorForEachCandidateP, nd, nb_candidatsP);

                        //Hypothèses sur p modulo le non diviseur
                        mpz_t *hypP;
                        int sizeHypP = nd - sizeNotCongruentP;
                        
                        hypP = malloc(sizeof(mpz_t) * sizeHypP);
                        find_hypothesis(hypP, notCongruentP, sizeNotCongruentP, nd);
                        
                        int a = 0;
                        mpz_set_ui(D[1],nd);
                        mpz_init(s2);
                        while(mpz_divisible_p(N,p) == 0 && a<sizeHypP) {
                           mpz_set(C[1],hypP[a]);
                           mpz_init(p);
                           chinese_remainder_theorem_spa(p, s2, D, C, 2);
                           a++;
                        }
                        
                        if (mpz_divisible_p(N,p) != 0){
                           mpz_divexact(q,N,p);
                           gmp_printf("p = %Zd\n", p);
                           gmp_printf("q = %Zd\n", q);
                           exit(1);
                        }
                     } 
                     c++;
                  }
                  printf("échec de l'attaque\n");
               }       
            }
         }          
      }
   }
}
