#include "gmp.h"
#include "math.h"

#ifndef DLC_SIDECHANNEL_ATTACK_ON_PRIME_GENERATION_FUNCTION_H
#define DLC_SIDECHANNEL_ATTACK_ON_PRIME_GENERATION_FUNCTION_H

void chinese_remainder_theorem(mpz_t p, mpz_t *sievePrimeList, unsigned int *candidats, int nbSievePrimes);

void find_k_first_primes(int size, mpz_t *sievePrimeList);

double correlation_coeff(double *X, double *Y, int size);

double euclidean_norm(double *vector, int size);

unsigned int argmin(double *L, int size);

#endif //DLC_SIDECHANNEL_ATTACK_ON_PRIME_GENERATION_FUNCTION_H
