CFLAGS=-Wall -Wextra -O3
LDFLAGS= -lgsl -lgslcblas -lconfig++
CC = g++


main: simulate.cc
	${CC} ${CFLAGS} -o simulate simulate.cc ${LDFLAGS}
