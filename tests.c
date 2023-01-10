#include "common/function.h"
#include "gmp.h"
#include "stdio.h"
#include "stdlib.h"

int main(int argc, char const *argv[]) {

    printf("Chinese test\n");

    mpz_t p;
    mpz_t primes[3];
    unsigned int candidats[3] = {3, 4, 5};

    mpz_inits(p, primes[0], primes[1], primes[2], NULL);
    mpz_set_ui(primes[0], 17);
    mpz_set_ui(primes[1], 11);
    mpz_set_ui(primes[2], 6);
    chinese_remainder_theorem(p, primes, candidats, 3);

    gmp_printf("%Zu\n", p);

    printf("Sieve gen x = 10 \n");
    mpz_t *primes2 = malloc(sizeof(mpz_t) * 50);
    find_k_first_primes(50, primes2);
    for (int i = 0; i < 50; i++)
        gmp_printf("%Zu \t", primes2[i]);
    printf("\n");


    double glaces[3] = {3, 6, 9},
            test3[3] = {3, 3, 3};;
    double heat[3] = {70, 75, 80};
    printf("Correlation de 1  : %f\n ", correlation_coeff(glaces, heat, 3));
    printf("Correlation de 1  : %f\n ", correlation_coeff(test3, heat, 3));
    return 0;
}