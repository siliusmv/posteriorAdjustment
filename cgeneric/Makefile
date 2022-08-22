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

test: test1.c test2.c smat-operations.o spde-precision.o test-funcs.o
	$(GCC) test1.c smat-operations.o spde-precision.o test-funcs.o -o test1.o
	$(GCC) test2.c smat-operations.o spde-precision.o test-funcs.o -o test2.o
