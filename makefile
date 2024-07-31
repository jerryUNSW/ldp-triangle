# CC=g++
# CPPFLAGS=-I.
# LDFLAGS=-g -fopenmp 
# DEPS = bigraph.h utility.h abcore.h
# OBJ = bigraph.o main.o utility.o abcore.o

# %.o: %.cpp $(DEPS)
# 	$(CC) -std=c++1y -c -O3 -o $@ $< $(CPPFLAGS) $(LDFLAGS)  # Include LDFLAGS here

# abcore: $(OBJ)
# 	$(CC) -std=c++1y -O3 -pthread -o $@ $^ $(CPPFLAGS) $(LDFLAGS) -lgomp  # Include LDFLAGS and link against libgomp

# clean:
# 	-rm -f abcore *.o

CC = g++
CPPFLAGS = -I. -I/usr/local/include  
LDFLAGS = -g -fopenmp -L/usr/local/lib  
LDLIBS = -lgsl -lgslcblas -lm -lgomp  # Link against GSL libraries and libgomp

DEPS = bigraph.h utility.h abcore.h
OBJ = bigraph.o main.o utility.o abcore.o

%.o: %.cpp $(DEPS)
	$(CC) -std=c++1y -c -O3 -o $@ $< $(CPPFLAGS) $(LDFLAGS)

abcore: $(OBJ)
	$(CC) -std=c++1y -O3 -pthread -o $@ $^ $(CPPFLAGS) $(LDFLAGS) $(LDLIBS)

clean:
	-rm -f abcore *.o

