CC = gcc

CFLAGS = -Wall -Wextra

CLIBS = -lgmp -lm

genDots: CPA/genDots.c
	$(CC) $(CFLAGS) -o genDots common/function.c CPA/genDots.c $(CLIBS)

genSeq: SPA/genSeq.c
	$(CC) $(CFLAGS) -o genSeq common/function.c SPA/genSeq.c $(CLIBS)

analyse: CPA/analyse.c
	$(CC) $(CFLAGS) -o analyse common/function.c CPA/analyse.c  $(CLIBS)

analyseSeq: SPA/analyseSeq.c
	$(CC) $(CFLAGS) -o analyseSeq common/function.c SPA/analyseSeq.c  $(CLIBS)

all: genDots analyse genSeq analyseSeq


clean:
	@rm -f *.o
