
CC = gcc
CFLAGS = -O3 -Wall -std=c99
LIBS = -lm -lgsl -lgslcblas

HEADERS = 
OBJECTS = main1.o
PROGRAM = run

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

