CXX = g++
CXXFLAGS = -Wall -ansi -pedantic -O3 -g --std=c++0x
LDFLAGS = -Wl,-rpath,/usr/local/lib
OBJECTS = main.o nwdp.o
BINPATH = ../bin
#DPALIGNPATH = /media/data/phd/projects/DotAligner/src

all: $(BINPATH)/DotAligner Makefile

$(BINPATH)/DotAligner: Makefile $(OBJECTS) 
		$(CXX) -o $(BINPATH)/DotAligner $(LDFLAGS) $(OBJECTS)

main.o: main.cpp
		$(CXX) $(CXXFLAGS) -c main.cpp
nwdp.o: nwdp.cpp
		$(CXX) $(CXXFLAGS) -c nwdp.cpp

clean: 
		rm -f -v *.o $(BINPATH)/DotAligner

.PHONY: clean
