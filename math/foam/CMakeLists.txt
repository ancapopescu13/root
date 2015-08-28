############################################################################
# CMakeLists.txt file for building ROOT math/foam package
############################################################################
set(libname Foam)

ROOT_USE_PACKAGE(hist/hist)

ROOT_GENERATE_DICTIONARY(G__${libname} *.h MODULE ${libname} LINKDEF LinkDef.h OPTIONS "-writeEmptyRootPCM")

ROOT_LINKER_LIBRARY(${libname} *.cxx G__${libname}.cxx LIBRARIES Vcshared ${_srcs} DEPENDENCIES Hist MathCore)
ROOT_INSTALL_HEADERS()
