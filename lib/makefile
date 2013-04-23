include ../local.mk
include ${KENT}/inc/common.mk

KENTLIBDIR = ${KENT}/lib/${MACHTYPE}
HG_INC=-I${KENT}/inc -I${GSL}/include -I${SAMTOOL} -I../include

O = bigs.o cluster.o extrema.o fastq_stuff.o favs.o metaBig.o random_coord.o sax.o

andycrg.a: $(O)
	ar rcus andycrg.a $(O)

clean:
	rm -f ${O} andycrg.a
