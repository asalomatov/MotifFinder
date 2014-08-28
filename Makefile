CC=g++
#CFLAGS=-c -g -Wall -std=c++0x 
CFLAGS=-c -O3 -s -Wall -std=c++0x 

all: motifFinder.exe

motifFinder.exe: aux.o BioSeq.o mf.o
	$(CC) -o motifFinder.exe aux.o BioSeq.o mf.o

mf.o: mf.cpp
	$(CC) $(CFLAGS) mf.cpp
BioSeq.o: BioSeq.cpp
	$(CC) $(CFLAGS) BioSeq.cpp
aux.o: aux.cpp
	$(CC) $(CFLAGS) aux.cpp

clean:
	rm -rf *.o *.exe

