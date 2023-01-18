void chinese_remainder_theorem_cpa(mpz_t p, mpz_t *sievePrimeList, unsigned int *candidats, int nbSievePrimes) {
    mpz_t N;
    mpz_t current;
    mpz_t Ni;
    mpz_t inv_Ni, tmp2;

    mpz_inits(N, current, Ni, inv_Ni, tmp2, NULL);
    mpz_set_ui(N, 2);
    mpz_set_ui(tmp2, 2);

    for (int i = 0; i < nbSievePrimes; i++) {
        mpz_mul(N, N, sievePrimeList[i]);
    }
    mpz_divexact(Ni, N, tmp2);
    mpz_invert(inv_Ni, Ni, tmp2);
    mpz_mul_ui(current, Ni, 1);
    mpz_mul(current, current, inv_Ni);
    mpz_add(p, p, current);
    mpz_mod(p, p, N);
    //gmp_printf("Debug :\nN = %Zu\n", N);
    for (int i = 0; i < nbSievePrimes; i++) {
        mpz_divexact(Ni, N, sievePrimeList[i]);
        mpz_invert(inv_Ni, Ni, sievePrimeList[i]);
        mpz_mul_ui(current, Ni, candidats[i]);
        mpz_mul(current, current, inv_Ni);
        mpz_add(p, p, current);
        mpz_mod(p, p, N);
    }
    mpz_mod(p, p, N);
}

void chinese_remainder_theorem_spa(mpz_t p, mpz_t N, mpz_t *divisors, mpz_t *congruences, int nb_divisors) {
    mpz_t current;
    mpz_t Ni;
    mpz_t inv_Ni;

    mpz_inits(N, current, Ni, inv_Ni, NULL);
    mpz_set_ui(N, 1);


    for (int i = 0; i < nb_divisors; i++) {
        mpz_mul(N, N, divisors[i]);
    }
    //gmp_printf("Debug :\nN = %Zu\n", N);
    for (int i = 0; i < nb_divisors; i++) {
        mpz_divexact(Ni, N, divisors[i]);
        mpz_invert(inv_Ni, Ni, divisors[i]);
        mpz_mul(current, Ni, congruences[i]);
        mpz_mul(current, current, inv_Ni);
        mpz_add(p, p, current);
        mpz_mod(p, p, N);
    }
    mpz_mod(p, p, N);
}


void find_k_first_primes(int size, mpz_t *sievePrimeList) {
    mpz_t current;
    int i = 0, j, isPrime;

    mpz_init(current);
    mpz_set_ui(current, 3);

    while (i < size) {
        isPrime = !mpz_divisible_ui_p(current, 2);

        for (j = 0; j < i && isPrime; j++) {
            if (mpz_divisible_p(current, sievePrimeList[j]))isPrime = 0;
        }

        if (isPrime) {
            mpz_init(sievePrimeList[i]);
            mpz_set(sievePrimeList[i], current);
            i++;
        }

        mpz_add_ui(current, current, 1);
    }


}

double correlation_coeff(double *X, double *Y, int size) {
    double mean_X = 0;
    double mean_Y = 0;

    for (int i = 0; i < size; i++) {
        mean_X = mean_X + X[i];
        mean_Y = mean_Y + Y[i];
    }
    mean_X = mean_X / size;
    mean_Y = mean_Y / size;

    double num = 0;
    double denum;
    double sx = 0;
    double sy = 0;
    for (int i = 0; i < size; i++) {
        num = num + (X[i] - mean_X) * (Y[i] - mean_Y);
        sx = sx + pow(X[i] - mean_X, 2);
        sy = sy + pow(Y[i] - mean_Y, 2);
    }
    denum = sqrt(sx * sy);
    if (denum == 0)return 0;
    return num / denum;
}

double euclidean_norm(double *vector, int size) {
    double norm = 0;
    for (int i = 0 ; i < size ; i++) {
        norm = norm + pow(vector[i],2);
    }
    return sqrt(norm);
}

unsigned int argmin(double *L, int size) {
   double min = L[0];
   int ind_min = 0;
   for (int i = 0; i< size ; i++) {
        if( L[i] < min) {
            min = L[i];
            ind_min = i;
        }
   }
   return ind_min;
}

void gen_k_bits_number(mpz_t out, int size, gmp_randstate_t randstate) {
    mpz_t offset;
    mpz_inits(offset, out, NULL);
    mpz_set_ui(offset, 2);
    mpz_pow_ui(offset, offset, size - 1);
    mpz_urandomb(out, randstate, size - 1);
    mpz_add(out, out, offset);
    mpz_clear(offset);
}

void gen_k_bits_number_odd(mpz_t out, int size, gmp_randstate_t randstate) {
    gen_k_bits_number(out, size, randstate);
    if (mpz_divisible_ui_p(out, 2)) {
        mpz_add_ui(out, out, 1);
    }
}
