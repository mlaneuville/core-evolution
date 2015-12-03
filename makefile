#!/bin/bash

REV = $(shell git rev-parse --verify HEAD)
CC = g++

a.out: conductivity.h conductivity.cpp main.h main.cpp
	@rm -rf revision.h
	@echo "string revision = \"$(REV)\";" > revision.h
	$(CC) -O3 -ffast-math conductivity.cpp main.cpp
	@rm -rf revision.h

clean:
	rm -rf *.txt
	rm -rf *.log
	rm -rf a.out
