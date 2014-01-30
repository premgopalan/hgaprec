bin_PROGRAMS = hgaprec
hgaprec_SOURCES = env.hh ratings.hh ratings.cc matrix.hh main.cc log.cc log.hh hgaprec.cc hgaprec.hh gpbase.hh 
hgaprec_LDADD = ./nmflib/lib/libnnmf.a -larpack -llapack -lblas
if DEBUG
AM_CFLAGS = -g  -O0
AM_CXXFLAGS = -g -O0
endif

