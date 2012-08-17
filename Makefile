FLAGS= -lm -Wall -g
CC=gcc
CPP=g++

td: main.c
	$(CPP) $^ -o td $(FLAGS)
histo: makehisto.c
	$(CPP) $^ -o makehisto $(FLAGS)
all: td histo
clean:
	rm -rf td makehisto
