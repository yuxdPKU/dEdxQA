##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libdedxqa.la

#
libdedxqa_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib64 \
  -L$(OFFLINE_MAIN)/lib

libdedxqa_la_LIBADD = \
  -lActsCore \
  -lActsPluginTGeo \
  -lActsExamplesCommon \
  -lActsExamplesDetectorTGeo \
  -lActsExamplesFramework \
  -lSubsysReco \
  -ltrackbase_historic_io \
  -ltrack_io \
  -ltrack_reco \
  -lglobalvertex_io \
  -lg4detectors_io \
  -lphparameter

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -DHomogeneousField

pkginclude_HEADERS = \
  dEdx.h

libdedxqa_la_SOURCES = \
  dEdx.cc

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_dedxqa

BUILT_SOURCES  = testexternals.cc

testexternals_dedxqa_SOURCES = testexternals.cc
testexternals_dedxqa_LDADD = libdedxqa.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f $(BUILT_SOURCES)
