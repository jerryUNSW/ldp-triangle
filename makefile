CC = g++
CPPFLAGS = -I. -I/usr/local/include  
LDFLAGS = -g -fopenmp -L/usr/local/lib  
LDLIBS = -lgsl -lgslcblas -lm -lgomp  # Link against GSL libraries and libgomp

DEPS = graph.h utility.h counting.h
OBJ = graph.o main.o utility.o counting.o

%.o: %.cpp $(DEPS)
	$(CC) -std=c++1y -c -O3 -o $@ $< $(CPPFLAGS) $(LDFLAGS)

counting: $(OBJ)
	$(CC) -std=c++1y -O3 -pthread -o $@ $^ $(CPPFLAGS) $(LDFLAGS) $(LDLIBS)

clean:
	-rm -f counting *.o

