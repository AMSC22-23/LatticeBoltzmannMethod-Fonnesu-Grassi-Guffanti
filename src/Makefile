CC = g++
CPPFLAGS = -std=c++20 -Wall -Wextra


all: main
	./main
.PHONY: run

main: main.o
	$(CC) $^ -o $@
main.o: main.cpp
	$(CC) $(CPPFLAGS) -c $^ -o $@


clean:
	rm *.o
.PHONY: clean