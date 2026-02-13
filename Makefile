PROGRAM=cleanfasta2


CC=g++
CFLAGS=-Wall -O3 -std=c++17

$PROGRAM: source/main.cpp source/gff.cpp source/sequenceHandler.cpp source/common.cpp
	$(CC) $(CFLAGS) -o $(PROGRAM) source/main.cpp source/gff.cpp source/sequenceHandler.cpp source/common.cpp

