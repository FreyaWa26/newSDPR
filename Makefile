
CC = g++

CFLAGS = -g3 -O3 -DHAVE_INLINE -march=native -I/home/gz222/software/gsl/include -std=c++11 

all: SDPR_admix score

SDPR_admix: parse_gen.o mcmc.o regress.o
	${CC} ${CFLAGS} parse_gen.o mcmc.o regress.o -L/home/gz222/software/gsl/lib/ -lgsl -lgslcblas -o SDPR_admix

regress.o: parse_gen.cpp parse_gen.h regress.cpp regress.h
	${CC} ${CFLAGS} -c regress.cpp

parse_gen.o: parse_gen.cpp parse_gen.h
	${CC} ${CFLAGS} -c parse_gen.cpp

mcmc.o: parse_gen.cpp parse_gen.h mcmc.cpp mcmc.h regress.cpp regress.h
	${CC} ${CFLAGS} -c mcmc.cpp

score: score.o
	${CC} ${CFLAGS} score.o -o score

score.o: score.cpp score.h
	${CC} ${CFLAGS} -c score.cpp

clean:
	rm -f *.o
