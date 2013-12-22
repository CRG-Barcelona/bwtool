check_PROGRAMS = bwmake wigmake
bwmake_CFLAGS = -I$(top_srcdir)/libjkweb 
bwmake_LDADD = $(top_srcdir)/libjkweb/libjkweb.a
bwmake_SOURCES = bwmake.c
wigmake_CFLAGS = -I$(top_srcdir)/libjkweb -I$(top_srcdir)/lib
wigmake_LDADD = $(top_srcdir)/lib/libbeato.a $(top_srcdir)/libjkweb/libjkweb.a
wigmake_SOURCES = wigmake.c
TESTS = \
	scripts/aggregate_main.wig_agg1.bed.1.sh \
	scripts/aggregate_main.wig_agg1.bed.2.sh \
	scripts/aggregate_main.wig_agg1.bed.3.sh \
	scripts/aggregate_2_and_2.sh \
	scripts/chromgraph_main_every_5.sh \
	scripts/extract_main_agg1.sh \
	scripts/find_main.bw_morethan4.sh \
	scripts/shift_main.bw_+3.sh \
	scripts/shift_main.bw_-2.sh \
	scripts/fill_main.bw_zero.sh \
	scripts/lift_to_new_main.sh \
	scripts/lift_main_inversion.1.sh \
	scripts/lift_main_inversion.2.sh \
	scripts/matrix_two_wigs.1.sh \
	scripts/matrix_two_wigs.2.sh \
	scripts/matrix_two_wigs.3.sh \
	scripts/matrix_two_wigs.4.sh \
	scripts/paste_main.bw_second.bw.1.sh \
	scripts/paste_main.bw_second.bw.2.sh \
	scripts/remove_main.bw_agg1.bed.sh \
	scripts/summary_main_every3.sh \
	scripts/summary_main_every10_fillzero_wsum.sh
