FLAGS=-lm -Wall
CC=gcc

all: td histo
td: main.c
	$(CC) $^ -o td $(FLAGS)
histo: makehisto.c
	$(CC) $^ -o makehisto $(FLAGS)
clean:
	rm -rf td makehisto
