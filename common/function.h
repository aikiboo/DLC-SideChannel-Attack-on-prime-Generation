#include "gmp.h"
#include "math.h"

#ifndef DLC_SIDECHANNEL_ATTACK_ON_PRIME_GENERATION_FUNCTION_H
#define DLC_SIDECHANNEL_ATTACK_ON_PRIME_GENERATION_FUNCTION_H

void chinese_remainder_theorem_cpa(mpz_t p, mpz_t *sievePrimeList, unsigned int *candidats, int nbSievePrimes);

void chinese_remainder_theorem_spa(mpz_t p, mpz_t *divisors, mpz_t *congruences, int nb_divisors);

void find_k_first_primes(int size, mpz_t *sievePrimeList);

double correlation_coeff(double *X, double *Y, int size);

double euclidean_norm(double *vector, int size);

unsigned int argmin(double *L, int size);

/*
  génère un nombre de k bits
*/
void gen_k_bits_number(mpz_t out, int size, gmp_randstate_t randstate);

/*
Renvoi un nombre impair, en appellant gen_k_bits_number() et, ajoute 1 si la sortie est pair
*/
void gen_k_bits_number_odd(mpz_t out, int size, gmp_randstate_t randstate);

#endif //DLC_SIDECHANNEL_ATTACK_ON_PRIME_GENERATION_FUNCTION_H
