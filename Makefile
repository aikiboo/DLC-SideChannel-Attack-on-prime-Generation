CC = gcc

NAME = genDots

SRCS = genDots.c

OBJS = $(SRCS:.c=.o)

CFLAGS = -Wall -Wextra

CLIBS = -lgmp -lm

$(NAME): $(OBJS)
	$(CC) $(CFLAGS) -o $(NAME) $(OBJS) $(CLIBS) -lm

all: $(NAME)

%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	@rm -f $(OBJS)
