#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void main(int argc, char const *argv[]) {

    if (argc < 4) {
        printf("Usage : ./analyse [filename] [sieveSize] [nbInt] \n");
        exit(1);
    }

    FILE *fptrHamAndNoise = fopen(argv[1], "r");
    int nbSievePrimes = atoi(argv[2]);
    int nb_candidats = atoi(argv[3]);
    srand(time(NULL));
    FILE *output;
    output = fopen("output_sequence.txt", "w");

}