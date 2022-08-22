GCC = gcc-12
c-files = $(wildcard *.c)
header-files = $(wildcard *.h)

.PHONY: test

%.so: %.c $(header-files)
	$(GCC) -Wall -fpic -O -c -o tmp.o $<
	$(GCC) -shared -o $@ tmp.o
	rm tmp.o

%.o: %.c %.h
	$(GCC) -Wall -fpic -O -c -o $@ $<

test: test.c smat-operations.o spde-precision.o
	$(GCC) test.c smat-operations.o spde-precision.o -o test.o
