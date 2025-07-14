# Makefile for 
# A. Ustun 11.09.2023
#
# For a successful compilation
# find directory name by 
# $(dirname $(find /usr -name gmt.h))
# and replace it with default path of GMTINC variable
#

SHELL := /bin/bash
ICC     = icc
GCC     = gcc
MCC     = icc

CDEBUG = -g
CFLAGS = -Wall
LDLIBS = -lm
LDFLAGS= 
ICFLAGS= -O2 -xCORE-AVX512 -qopt-zmm-usage=high -no-multibyte-chars
IMFLAGS= -O2 -xCORE-AVX2
GCXFLAG= -O2 -march=native -ffast-math
OMFLAGS= -qopenmp -parallel
GOFLAGS= -fopenmp
GCFLAGS= -O3 -march=native -ffast-math
GMTINC = $(shell gmt-config --cflags)
GMTLIB = $(shell gmt-config --libs)


BINDIR = /usr/local/bin
MANDIR = /usr/local/share/man
OPTREP = -g -qopt-report=5 -qopt-report-phase=all

MACFLAG= -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk

ICCOBJ = tcc.oi tc-cylinder.oi
GCCOBJ = tcc.og tc-cylinder.og
MACOBJ = tcc.om tc-cylinder.om

all: tc-cylinder-icc tc-cylinder-gcc

tc-cylinder-icc: $(ICCOBJ)
	$(ICC) $^ -o $@ $(ICFLAGS) $(OMFLAGS) $(LDLIBS) $(GMTLIB)

tc-cylinder-gcc: $(GCCOBJ)
	$(GCC) $^ -o $@ $(GCFLAGS) $(GOFLAGS) $(LDLIBS) $(GMTLIB)

tc-cylinder-mac: $(MACOBJ)
	$(MCC) $^ -o $@ $(IMFLAGS) $(OMFLAGS) $(LDLIBS) $(GMTLIB) $(MACFLAG) 

.SUFFIXES: .c .oi .og .om

.c.oi:
	$(ICC) -c $(GMTINC) $(OMFLAGS) $(CFLAGS) $(ICFLAGS) $(OPTREP) -o $@ $< 

.c.og:
	$(GCC) -c $(GMTINC) $(GCFLAGS) $(GOFLAGS) $(CFLAGS) $(LDLIBS) -o $@ $<

.c.om:
	$(MCC) -c $(GMTINC) $(OMFLAGS) $(CFLAGS) $(IMFLAGS) $(OPTREP) $(MACFLAG) -o $@ $< 



install: tc-cylinder-icc
	cp -u tc-cylinder-icc $(BINDIR)/tc-cylinder
	cp -u ../man/tc-cylinder.1 $(MANDIR)/man1 

clean:
	rm -rf *.oi *.og *.om tc-cylinder

run-icctest: tc-cylinder-icc
	cat ../TestData/testpois.lfH | ./tc-cylinder-icc \
		-G../TestData/testDEM15s.grd -I../TestData/testDEM1s.grd \
		-Z56 -p3 -r0.5s > testout.txt 


run-icc: tc-cylinder-icc
	time head -n1 ../TestData/testpois.lfH | ./tc-cylinder-icc \
		-G../TestData/testDEM15s.grd -I../TestData/testDEM1s.grd \
		-Z56 -p3 -r0.5s -V | (head; tail)  

run-gcc: tc-cylinder-gcc
	time head -n1 ../TestData/testpois.lfH | ./tc-cylinder-gcc \
		-G../TestData/testDEM15s.grd -I../TestData/testDEM1s.grd \
		-Z56 -p3 -r0.5s | (head; tail)

run-mac: tc-cylinder-mac
	time head ../TestData/testpois.lfH | ./tc-cylinder-mac \
		-G../TestData/testDEM15s.grd -I../TestData/testDEM1s.grd \
		-Z56 -p3 -r0.5s
run:
	time head -n1 ../TestData/testpois.lfH | ./tc-cylinder-icc \
		-G../TestData/testDEM15s.grd -I../TestData/testDEM1s.grd \
		-Z56 -p3 -r0.5s

speed-icc:
	for run in {1..5}; do time head -n50 ../TestData/testpois.lfH | \
		./tc-cylinder-icc -G../TestData/testDEM15s.grd -I../TestData/testDEM1s.grd \
		-Z56 -p3 -r0.5s; done 2> speed-icc-050.test 1> /dev/null
	for run in {1..5}; do time head -n100 ../TestData/testpois.lfH | \
		./tc-cylinder-icc -G../TestData/testDEM15s.grd -I../TestData/testDEM1s.grd \
		-Z56 -p3 -r0.5s; done 2> speed-icc-100.test 1> /dev/null
	paste speed-icc*.test

speed-gcc:
	for run in {1..5}; do time head -n50 ../TestData/testpois.lfH | \
		./tc-cylinder-gcc -G../TestData/testDEM15s.grd -I../TestData/testDEM1s.grd \
		-Z56 -p3 -r0.5s; done 2> speed-gcc-050.test 1> /dev/null
	for run in {1..5}; do time head -n100 ../TestData/testpois.lfH | \
		./tc-cylinder-gcc -G../TestData/testDEM15s.grd -I../TestData/testDEM1s.grd \
		-Z56 -p3 -r0.5s; done 2> speed-gcc-100.test 1> /dev/null
	paste speed-gcc*.test

speed-mac:
	for run in {1..5}; do time head -n50 ../TestData/testpois.lfH | \
		./tc-cylinder-mac -G../TestData/testDEM15s.grd -I../TestData/testDEM1s.grd \
		-Z56 -p3 -r0.5s; done 2> speed-mac-050.test 1> /dev/null
	for run in {1..5}; do time head -n100 ../TestData/testpois.lfH | \
		./tc-cylinder-mac -G../TestData/testDEM15s.grd -I../TestData/testDEM1s.grd \
		-Z56 -p3 -r0.5s; done 2> speed-mac-100.test 1> /dev/null
	paste speed-mac.test

