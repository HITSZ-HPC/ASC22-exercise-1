CC=g++
CPPFLAGS= -Iinclude --std=c++11

.PHONY: default

default: SLIC main

run: default
	./main

SLIC: src/SLIC.cpp include/SLIC.h
	$(CC) $(CPPFLAGS) -c src/SLIC.cpp -o src/SLIC.o

main: main.cpp
	$(CC) $(CPPFLAGS) -c main.cpp -o main.o
	$(CC) $(CPPFLAGS) main.o src/SLIC.o -o main

clean:
	rm *.o
	rm src/*.o
	rm main