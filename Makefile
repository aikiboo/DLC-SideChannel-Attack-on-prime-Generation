CC = gcc

CFLAGS = -Wall -Wextra

CLIBS = -lgmp -lm

genDots: genDots.c
	$(CC) $(CFLAGS) -o genDots common/function.c genDots.c $(CLIBS)

genSeq: SPA/genSeq.c
	$(CC) $(CFLAGS) -o genSeq common/function.c SPA/genSeq.c $(CLIBS)

analyse: analyse.c
	$(CC) $(CFLAGS) -o analyse common/function.c analyse.c  $(CLIBS)

all: genDots analyse

test: tests.c
	$(CC) $(CFLAGS) -o tests common/function.c tests.c  $(CLIBS)

clean:
	@rm -f *.o
