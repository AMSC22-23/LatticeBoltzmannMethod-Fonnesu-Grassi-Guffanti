CC = g++
CPPFLAGS = -std=c++20 -Wall -Wextra

.PHONY=run
run: main
	./main

main: main.o
	$(CC) $^ -o $@
main.o: main.cpp
	$(CC) $(CPPFLAGS) -c $^ -o $@


