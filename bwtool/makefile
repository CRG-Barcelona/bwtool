include ../local.mk
include ${KENT}/inc/common.mk

L += -lm -lz
MYLIBDIR = ${KENT}/lib/${MACHTYPE}
MYLIBS = ../lib/andycrg.a ${MYLIBDIR}/jkweb.a
HG_INC=-I${KENT}/inc -I../include

A = bwtool
O = aggregate.o bwtool.o chromgraph.o distrib.o fill.o find.o lift.o matrix.o paste.o rand.o remove.o sax.o shift.o summarize.o

bwtool: ${O} ${MYLIBS}
	${CC} ${COPT} -o ${DESTDIR}${BINDIR}/${A}${EXE} $O ${MYLIBS} $L

compile:: ${O}
	${CC} ${COPT} -o ${A}${EXE} ${O} ${MYLIBS} $L

clean::
	rm -f ${A}${EXE} ${O}
