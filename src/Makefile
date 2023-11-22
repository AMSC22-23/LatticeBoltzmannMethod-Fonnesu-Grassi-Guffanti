CC = g++
CPPFLAGS = -std=c++20 -Wall -Wextra
SRCS = main.cpp velocitySet.cpp
OBJS = $(SRCS:.cpp=.o)

all: main
	./main
.PHONY: run

main: $(OBJS)
	$(CC) -o $@ $(OBJS)

%.o: %.cpp
	$(CC) $(CPPFLAGS) -c $< -o $@


clean:
	rm *.o
.PHONY: clean