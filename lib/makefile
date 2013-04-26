include ../local.mk
include ${KENT}/inc/common.mk

HG_INC=-I${KENT}/inc -I../include
ifeq (${USE_BAM},1)
	L+=${SAMDIR}/libbam.a
	HG_INC+=-I${SAMDIR}
endif

O = bigs.o cluster.o extrema.o fastq_stuff.o favs.o metaBig.o metaBigBam.o mt19937-64.o random_coord.o sax.o stuff.o

andycrg.a: $(O)
	ar rcus andycrg.a $(O)

clean:
	rm -f ${O} andycrg.a
