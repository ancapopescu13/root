############################################################################
# CMakeLists.txt file for building ROOT hist/spectrumpainter package
############################################################################

set(libname SpectrumPainter)

ROOT_USE_PACKAGE(hist/hist)


ROOT_GENERATE_DICTIONARY(G__${libname} *.h MODULE ${libname} LINKDEF LinkDef.h OPTIONS "-writeEmptyRootPCM")

#ROOT_LINKER_LIBRARY(${libname} *.cxx G__SpectrumPainter.cxx LIBRARIES Vcshared DEPENDENCIES Graf Hist )

ROOT_LINKER_LIBRARY(${libname} *.cxx G__${libname}.cxx LIBRARIES Vcshared DEPENDENCIES Graf Hist)
ROOT_INSTALL_HEADERS()
