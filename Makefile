CC = gcc

CFLAGS = -Wall -Wextra

CLIBS = -lgmp -lm

genDots: genDots.c
	$(CC) $(CFLAGS) -o genDots genDots.c $(CLIBS)

analyse: analyse.c
	$(CC) $(CFLAGS) -o analyse function.c analyse.c  $(CLIBS)

all: genDots analyse



clean:
	@rm -f *.o
