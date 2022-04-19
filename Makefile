APP_NAME=t0

CXX=g++
CXXFLAGS= -Wall -pedantic -lm

SRC =  t0.cpp util/colorConv.cpp util/superPixel.cpp

# Compilation Commands
$(APP_NAME): $(SRC)
	$(CXX) $(SRC) $(CXXFLAGS) -o $@

clean:
	/bin/rm -rf *~ *.o $(APP_NAME) *.class