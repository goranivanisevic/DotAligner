CXX = g++
#CXXFLAGS = -Wall -ansi -pedantic -O3 -g --std=c++0x
CXXFLAGS = -Wall -ansi -pedantic -O0 -g --std=c++0x #debug mode
LDFLAGS = -Wl,-rpath,/usr/local/lib
OBJECTS = main_v2.o nwdp_v2.o
BINPATH = ../bin
#DPALIGNPATH = /media/data/phd/projects/DotAligner/src

all: $(BINPATH)/DotAligner2 Makefile

$(BINPATH)/DotAligner2: Makefile $(OBJECTS) 
		$(CXX) -o $(BINPATH)/DotAligner2 $(LDFLAGS) $(OBJECTS)

main_v2.o: main_v2.cpp
		$(CXX) $(CXXFLAGS) -c main_v2.cpp
nwdp_v2.o: nwdp_v2.cpp
		$(CXX) $(CXXFLAGS) -c nwdp_v2.cpp

clean: 
		rm -f -v *.o $(BINPATH)/DotAligner2

.PHONY: clean
