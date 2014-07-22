bin_PROGRAMS = hgaprec
hgaprec_SOURCES = env.hh ratings.hh ratings.cc matrix.hh main.cc log.cc log.hh hgaprec.cc hgaprec.hh gpbase.hh
hgaprec_LDADD = -larpack -llapack -lblas
#
# AUTOCONF 2.63 has a bug!
# temp fix: uncomment to use these options
#
#if HAVE_NMFLIB
#hgaprec_LDADD += ./nmflib/lib/libnnmf.a
#endif
#if DEBUG
#AM_CFLAGS = -g  -O0
#AM_CXXFLAGS = -g -O0
#endif
#
