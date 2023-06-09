CC=g++
CFLAGS=-Wall -O3 -Ofast
LIBS=-lm

SRC=solver.cpp problem.cpp vector_utils.cpp
OBJ=$(SRC:%.cpp=%.o)

%.o: %.cpp %.hpp
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@

exec: main.cpp $(OBJ)
	$(CC) $(CFLAGS) $(INCDIR) $(LDFLAGS) $^ $(LIBS) -o $@

clean:
	rm -f *.o
	rm -f exec
