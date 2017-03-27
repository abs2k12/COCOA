# Copyright by Wojtek Pych, CAMK PAN
# 2014.03.03 version 0.2
#
NAME = sim2obs
#
# the C compiler
CC = gcc
# options for the C compiler
COPT = -O3 -Wall -pedantic
# libraries
LIBS = -lgsl -lgslcblas -lm

${NAME} : ${NAME}.o pfitshead.o swap.o pfitsin.o pfitsio.o
	${CC} -o ${NAME} $^ ${LIBS}

%.o : ${STDH} %.c
	${CC} ${COPT} -c $^
