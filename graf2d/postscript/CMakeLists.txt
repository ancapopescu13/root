############################################################################
# CMakeLists.txt file for building ROOT graf2d/postcript package
# @author Pere Mato, CERN
############################################################################
set(libname PostScript)

ROOT_USE_PACKAGE(graf2d/graf)
ROOT_USE_PACKAGE(graf2d/mathtext)

#ROOT_GENERATE_DICTIONARY(G__PostScript *.h MODULE Postscript LINKDEF LinkDef.h OPTIONS "-writeEmptyRootPCM")
ROOT_GENERATE_DICTIONARY(G__${libname} *.h MODULE ${libname} LINKDEF LinkDef.h OPTIONS "-writeEmptyRootPCM")

ROOT_LINKER_LIBRARY(${libname} *.cxx G__${libname}.cxx LIBRARIES Vcshared DEPENDENCIES Graf)
ROOT_INSTALL_HEADERS()
