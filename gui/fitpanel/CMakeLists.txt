############################################################################
# CMakeLists.txt file for building ROOT gui/fitpanel package
############################################################################
set(libname FitPanel)

ROOT_USE_PACKAGE(gui/gui)
ROOT_USE_PACKAGE(tree/treeplayer)

#ROOT_STANDARD_LIBRARY_PACKAGE(FitPanel DEPENDENCIES Gui Gpad Hist Graf Tree MathCore DICTIONARY_OPTIONS "-writeEmptyRootPCM")

ROOT_GENERATE_DICTIONARY(G__${libname} *.h MODULE ${libname} LINKDEF LinkDef.h OPTIONS "-writeEmptyRootPCM")

ROOT_LINKER_LIBRARY(${libname} *.cxx G__${libname}.cxx LIBRARIES Vcshared DEPENDENCIES Gui Gpad Hist Graf Tree MathCore)
ROOT_INSTALL_HEADERS()
