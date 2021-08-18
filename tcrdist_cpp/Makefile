# the compiler
CC = g++

## compiler flags:

# debugging
#CCFLAGS  = -g -std=c++11 -Wall

# production
CCFLAGS  = -O3 -std=c++11 -Wall

# EDIT THIS to point to the location where you installed boost
# this should be the main directory containing boost/ doc/ libs/ etc as subdirectories
# this is only used by the (unused, untested) find_cdr3_motifs code right now
# so you don't have to worry about it
BOOSTDIR = /home/pbradley/include/boost_1_67_0

# 
INCLUDES = -I ./include/

# recompile if any .hh files changed
HHS = ./src/*.hh

all: ./bin/find_neighbors ./bin/find_neighbors_single_chain ./bin/calc_distributions ./bin/count_matches_single_chain ./bin/count_matches_paired ./bin/find_paired_matches ./bin/compute_single_chain_distance_matrix

# this is still under development:
motifs: ./bin/find_cdr3_motifs


./bin/find_neighbors:  ./src/find_neighbors.cc  $(HHS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o ./bin/find_neighbors ./src/find_neighbors.cc

./bin/find_neighbors_single_chain:  ./src/find_neighbors_single_chain.cc  $(HHS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o ./bin/find_neighbors_single_chain ./src/find_neighbors_single_chain.cc

./bin/calc_distributions:  ./src/calc_distributions.cc  $(HHS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o ./bin/calc_distributions ./src/calc_distributions.cc

./bin/count_matches_single_chain:  ./src/count_matches_single_chain.cc  $(HHS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o ./bin/count_matches_single_chain ./src/count_matches_single_chain.cc

./bin/count_matches_paired:  ./src/count_matches_paired.cc  $(HHS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o ./bin/count_matches_paired ./src/count_matches_paired.cc

./bin/find_paired_matches:  ./src/find_paired_matches.cc  $(HHS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o ./bin/find_paired_matches ./src/find_paired_matches.cc

./bin/compute_single_chain_distance_matrix:  ./src/compute_single_chain_distance_matrix.cc  $(HHS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o ./bin/compute_single_chain_distance_matrix  ./src/compute_single_chain_distance_matrix.cc

./bin/find_cdr3_motifs:  ./src/find_cdr3_motifs.cc  $(HHS)
	$(CC) $(CCFLAGS) $(INCLUDES) -I $(BOOSTDIR) -o ./bin/find_cdr3_motifs ./src/find_cdr3_motifs.cc



clean:
	-rm ./bin/*
