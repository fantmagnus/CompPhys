
CC = gcc
CFLAGS = -O3 -Wall -std=c99
LIBS = -lm -lgsl -lgslcblas

HEADERS = func.h initfcc.h alpotential.h
OBJECTS = func.o initfcc.o alpotential.o MD_main.o
PROGRAM = MD

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

