APP_NAME=t0

CXX=g++ -m64 -std=c++11
CXXFLAGS= -I. -O3  -Wall -pedantic -fopenmp -Wno-unknown-pragmas

OBJS =  t0.o util/colorConv.o

%.o: %.cpp
	$(CXX) $< $(CXXFLAGS) -c -o $@

# Compilation Commands
$(APP_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

clean:
	/bin/rm -rf *~ *.o $(APP_NAME) *.class