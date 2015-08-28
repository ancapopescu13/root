############################################################################
# CMakeLists.txt file for building ROOT math/physics package
############################################################################
set(libname MLP)

ROOT_USE_PACKAGE(tree/tree)
ROOT_USE_PACKAGE(tree/treeplayer)
ROOT_USE_PACKAGE(hist/hist)

#ROOT_STANDARD_LIBRARY_PACKAGE(MLP DEPENDENCIES Hist Matrix Tree Graf Gpad TreePlayer MathCore)

ROOT_GENERATE_DICTIONARY(G__${libname} *.h MODULE ${libname} LINKDEF LinkDef.h OPTIONS "-writeEmptyRootPCM")

ROOT_LINKER_LIBRARY(${libname} *.cxx G__${libname}.cxx LIBRARIES Vcshared DEPENDENCIES Hist Matrix Tree Graf Gpad TreePlayer MathCore)
ROOT_INSTALL_HEADERS()

