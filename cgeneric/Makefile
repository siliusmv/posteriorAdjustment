GCC = gcc-12
c-files = $(wildcard *.c)
header-files = $(wildcard *.h)
compiled-files = $(c-files:.c=.so)

gsl_flags = $(shell gsl-config --cflags) $(shell gsl-config --libs)

%.so: %.c $(header-files)
	$(GCC) $(gsl_flags) -Wall -fpic -O -c -o tmp.o $<
	$(GCC) $(gsl_flags) -shared -o $@ tmp.o
	rm tmp.o

%.o: %.c $(header-files)
	$(GCC) $(gsl_flags) -Wall -fpic -O -o $@ $<

