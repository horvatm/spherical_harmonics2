CXX = g++

DEPS = $(wildcard include/*.h)
INCLUDE = -I../include -I.    
LDFLAGS = 
CXXFLAGS = -O3 -Wall -fopenmp -std=c++17

.PHONY: all clean

all: test_spher_harm.e test_dspher_harm.e

%.e : %.cpp ${DEPS}
	${CXX} ${CXXFLAGS} ${INCLUDE} $< -o $@ ${LDFLAGS}

clean:
	rm -rf *.e *.dat 
