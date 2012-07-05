FLAGS= -lm -Wall
CC=gcc

td: main.c
	$(CC) $^ -o td $(FLAGS)
histo: makehisto.c
	$(CC) $^ -o makehisto $(FLAGS)
all: td histo
clean:
	rm -rf td makehisto
