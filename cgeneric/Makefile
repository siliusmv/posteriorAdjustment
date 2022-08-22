GCC = gcc-12

.PHONY: test clean all

all: a.so b.so test

a.so: a.c cgeneric.h
	$(GCC) -Wall -fpic -O -c -o a.o a.c
	$(GCC) -shared -o a.so a.o

b.so: b.c cgeneric.h smat-operations.c smat-operations.h spde-precision.c spde-precision.h
	$(GCC) -Wall -fpic -O -c -o b.o b.c
	$(GCC) -shared -o b.so b.o

test: test1.c test2.c smat-operations.o spde-precision.o test-funcs.o cgeneric.h
	$(GCC) test1.c smat-operations.o spde-precision.o test-funcs.o -o test1.o
	$(GCC) test2.c smat-operations.o spde-precision.o test-funcs.o -o test2.o

%.o: %.c %.h
	$(GCC) -Wall -fpic -O -c -o $@ $<

clean:
	rm *.o
	rm *.so
