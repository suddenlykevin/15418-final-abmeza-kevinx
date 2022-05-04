APP_NAME=t0

CXX=g++
CXXFLAGS= -pedantic -lm -I. -Wall -fopenmp -Wno-unknown-pragmas

SRC =  t0.cpp pixImage.cpp util/colorConv.cpp util/superPixel.cpp util/CycleTimer.h

# Compilation Commands
$(APP_NAME): $(SRC)
	$(CXX) $(SRC) $(CXXFLAGS) -o $@

clean:
	/bin/rm -rf *~ *.o $(APP_NAME) *.class