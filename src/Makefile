CC          := g++

CXXFLAGS        := -O3 -Wall -std=c++11

../bin/hw2: hw2.cpp
	$(CC) $(CXXFLAGS) -pthread -o $@ $? -fopenmp
clean:
	rm -f ../bin/hw2