#!/bin/bash

BIN = evolve
SRC = src/main.cpp src/conductivity.cpp src/core.cpp src/mantle.cpp src/fd1d_heat_steady.cpp
OBJ = $(SRC:.cpp=.o)
FLAGS = -ffast-math -O3 -lyaml-cpp -std=c++11

REV = $(shell git rev-parse --verify HEAD)
CC = g++
INC = /usr/local/include
LIB = /usr/local/lib

$(BIN): revision.h $(OBJ)
	@rm -rf revision.h
	@echo "string revision = \"$(REV)\";" > revision.h
	$(CC) -o $(BIN) $(OBJ) $(FLAGS)
	@rm -rf revision.h

.cpp.o:
	$(CC) -c $*.cpp $(FLAGS) -o $@


clean:
	rm -rf *.txt
	rm -rf *.log
	rm -rf src/*.o
	rm -rf a.out
