
CC = gcc
CFLAGS = -O3 -Wall -std=c99
LIBS = -lm -lgsl -lgslcblas

HEADERS = 
OBJECTS = DMC_main.o
PROGRAM = run

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

test:
	make
	./run

