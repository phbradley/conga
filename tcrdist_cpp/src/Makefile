# the compiler
CC = g++

## compiler flags:

# debugging
#CCFLAGS  = -g -std=c++11 -Wall

# production
CCFLAGS  = -O3 -std=c++11 -Wall

# 
INCLUDES = -I ../include/

# recompile if any .hh files changed
HHS = *.hh

all: ../bin/find_neighbors ../bin/calc_distributions ../bin/count_matches_single_chain ../bin/count_matches_paired


../bin/find_neighbors:  find_neighbors.cc  $(HHS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o ../bin/find_neighbors find_neighbors.cc

../bin/calc_distributions:  calc_distributions.cc  $(HHS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o ../bin/calc_distributions calc_distributions.cc

../bin/count_matches_single_chain:  count_matches_single_chain.cc  $(HHS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o ../bin/count_matches_single_chain count_matches_single_chain.cc

../bin/count_matches_paired:  count_matches_paired.cc  $(HHS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o ../bin/count_matches_paired count_matches_paired.cc


clean:
	-rm ./bin/*
