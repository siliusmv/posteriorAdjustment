c-files = $(wildcard *.c)
header-files = $(wildcard *.h)
compiled-files = $(c-files:.c=.so)

include_flags = $(shell gsl-config --cflags)
lib_flags = $(shell gsl-config --libs)
cflags = $(include_flags) $(lib_flags)

%.so: %.c
	gcc-11 $(cflags) -fopenmp -Wall -fpic -O -c -o tmp.o $<
	gcc-11 $(cflags) -fopenmp -shared -o $@ tmp.o
	rm tmp.o

%.o: %.c $(header-files)
	gcc-11 $(cflags) -fopenmp -Wall -fpic -O -o $@ $<

