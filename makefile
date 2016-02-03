#!/bin/bash

REV = $(shell git rev-parse --verify HEAD)
CC = g++
INC = /usr/local/include
LIB = /usr/local/lib

a.out: conductivity.h conductivity.cpp main.h main.cpp
	@rm -rf revision.h
	@echo "string revision = \"$(REV)\";" > revision.h
	$(CC) -O3 -ffast-math -I$(INC) -L$(LIB) -lyaml-cpp conductivity.cpp main.cpp
	@rm -rf revision.h

clean:
	rm -rf *.txt
	rm -rf *.log
	rm -rf a.out
