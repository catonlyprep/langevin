CXX ?= g++
CXXFLAGS = -O3 -fPIC -Wall -std=c++11 -march=native
LFLAGS = -lgsl -lgslcblas -lhdf5_cpp -lhdf5

#####################################
# No need to modify below this line #
#####################################

# Build targets
TARGETEXE = langevin
TARGETLIB = liblangevin.so

# Source & object files
SRCS = src/langevin.cpp
OBJS = $(SRCS:.cpp=.o)

.PHONY: all
all: $(TARGETEXE) $(TARGETLIB)

$(TARGETEXE): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LFLAGS)
	
$(TARGETLIB): $(OBJS)
	$(CXX) $(CXXFLAGS) -shared $^ -o $@ $(LFLAGS)

.PHONY: clean

clean:
	rm *~ *.o  $(TARGETEXE) $(TARGETLIB)
