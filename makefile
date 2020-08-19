#!/bin/bash

BIN = evolve
OBJ = main.o conductivity.o core.o mantle.o fd1d_heat_steady.o
FLAGS = -ffast-math -O3 -lyaml-cpp -std=c++11

REV = $(shell git rev-parse --verify HEAD)
CC = g++
INC = /usr/local/include
LIB = /usr/local/lib

$(BIN): revision.h $(OBJ) *.h
	@rm -rf revision.h
	@echo "string revision = \"$(REV)\";" > revision.h
	$(CC) -o $(BIN) $(OBJ) $(FLAGS)
	@rm -rf revision.h

revision.h:
	@rm -rf revision.h
	@echo "string revision = \"$(REV)\";" > revision.h

$(OBJ): *.h

.cpp.o:
	$(CC) -c $*.cpp $(FLAGS)

clean:
	rm -rf *.txt
	rm -rf *.log
	rm -rf *.o
	rm -rf a.out
