############################################################################
# CMakeLists.txt file for building ROOT montecarlo/vmc package
# @author Pere Mato, CERN
############################################################################
set(libname VMC)

ROOT_USE_PACKAGE(montecarlo/eg)
ROOT_USE_PACKAGE(geom/geom)

#ROOT_STANDARD_LIBRARY_PACKAGE(VMC DEPENDENCIES EG Geom MathCore)

ROOT_GENERATE_DICTIONARY(G__${libname} *.h MODULE ${libname} LINKDEF LinkDef.h OPTIONS "-writeEmptyRootPCM")

ROOT_LINKER_LIBRARY(${libname} *.cxx G__${libname}.cxx LIBRARIES Vcshared DEPENDENCIES EG Geom MathCore)
ROOT_INSTALL_HEADERS()
