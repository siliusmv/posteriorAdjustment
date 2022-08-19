GCC = gcc-12
c-files = $(wildcard *.c)
header-files = $(wildcard *.h)
compiled-files = $(c-files:.c=.so)
GSL_CONFIG = gsl/bin/gsl-config

gsl_flags = $(shell $(GSL_CONFIG) --cflags) $(shell $(GSL_CONFIG) --libs)

%.so: %.c $(header-files)
	$(GCC) $(gsl_flags) -Wall -fpic -O -c -o tmp.o $<
	$(GCC) $(gsl_flags) -shared -o $@ tmp.o
	rm tmp.o

%.o: %.c $(header-files)
	$(GCC) $(gsl_flags) -Wall -fpic -O -o $@ $<

