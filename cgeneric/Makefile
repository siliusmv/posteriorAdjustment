GCC = gcc-12

.PHONY: test clean all

all: a.so b.so test

a.so: a.c cgeneric.h
	$(GCC) -Wall -fpic -O -c -o a.o -lm a.c
	$(GCC) -shared -o a.so -lm a.o

b.so: b.c cgeneric.h smat-operations.c smat-operations.h spde-precision.c spde-precision.h
	$(GCC) -Wall -fpic -O -c -o b.o -lm b.c
	$(GCC) -shared -o b.so -lm b.o

test: test1.c test2.c smat-operations.o spde-precision.o test-funcs.o cgeneric.h
	$(GCC) -o test1.o -lm test1.c smat-operations.o spde-precision.o test-funcs.o
	$(GCC) -o test2.o -lm test2.c smat-operations.o spde-precision.o test-funcs.o

%.o: %.c %.h
	$(GCC) -Wall -fpic -O -c -o $@ -lm $<

clean:
	rm *.o
	rm *.so
