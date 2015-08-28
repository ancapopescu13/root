############################################################################
# CMakeLists.txt file for building ROOT gui/ged package
############################################################################
set(libname Ged)

ROOT_USE_PACKAGE(gui/gui)
ROOT_USE_PACKAGE(tree/treeplayer)

ROOT_GENERATE_DICTIONARY(G__${libname} *.h MODULE ${libname} LINKDEF LinkDef.h OPTIONS "-writeEmptyRootPCM")

#ROOT_STANDARD_LIBRARY_PACKAGE(Ged DEPENDENCIES Hist Gpad Graf Gui Tree TreePlayer)

ROOT_LINKER_LIBRARY(${libname} *.cxx G__${libname}.cxx LIBRARIES Vcshared DEPENDENCIES Hist Gpad Graf Gui Tree TreePlayer)
ROOT_INSTALL_HEADERS()
