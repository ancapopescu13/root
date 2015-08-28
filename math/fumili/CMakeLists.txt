############################################################################
# CMakeLists.txt file for building ROOT math/fumili package
############################################################################
set(libname Fumili)

ROOT_USE_PACKAGE(math/matrix)
ROOT_USE_PACKAGE(math/mathcore)
ROOT_USE_PACKAGE(hist/hist)

ROOT_GENERATE_DICTIONARY(G__${libname} *.h MODULE ${libname} LINKDEF LinkDef.h OPTIONS "-writeEmptyRootPCM")

#ROOT_STANDARD_LIBRARY_PACKAGE(Fumili LIBRARIES Vcshared ${_srcs} DEPENDENCIES Graf Hist MathCore)

ROOT_LINKER_LIBRARY(${libname} *.cxx G__${libname}.cxx LIBRARIES Vcshared DEPENDENCIES Graf Hist MathCore)
ROOT_INSTALL_HEADERS()
