SHELL    = /bin/bash
CXX      = g++
CXXFLAGS = -O3 -g -std=c++11
src_dir  = ../src
src      = $(wildcard $(src_dir)/*.cpp)
obj      = $(src:.cpp=.o)
exe      = $(wildcard *.cpp)
targets  = test prob_2_part_a prob_2_part_b prob_2_part_c


# make reference : http://nuclear.mutantstargoat.com/articles/make/
# @ is a built-in make variable containing the target of each rule
# ^ is a built-in variable containing all the dependencies of each rule
# < is a built-in variable containing only the first element of the dependencies


all: $(exe) $(targets)

test: test.cpp $(obj)
	$(CXX) $(CXXFLAGS) -I $(src_dir) -o $@ $^

prob_2_part_a: prob_2_part_a.cpp $(obj)
	$(CXX) $(CXXFLAGS) -I $(src_dir) -o $@ $^

prob_2_part_b: prob_2_part_b.cpp $(obj)
	$(CXX) $(CXXFLAGS) -I $(src_dir) -o $@ $^

prob_2_part_c: prob_2_part_c.cpp $(obj)
	$(CXX) $(CXXFLAGS) -I $(src_dir) -o $@ $^


.PHONY: clean
clean:
	rm -f $(targets) $(obj)


