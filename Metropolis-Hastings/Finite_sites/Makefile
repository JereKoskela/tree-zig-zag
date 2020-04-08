CFLAGS=-Wall -Wextra -O3
LDFLAGS= -lgsl -lgslcblas -lconfig++
CC = g++


main: simulate.cc 
	${CC} ${CFLAGS} -o simulate simulate.cc ${LDFLAGS}

gen: gen_data.cc 
	${CC} ${CFLAGS} -o gen_data gen_data.cc ${LDFLAGS}
