############################################################################
# CMakeLists.txt file for building ROOT gui/sessionviewer package
############################################################################
set(libname SessionViewer)

ROOT_USE_PACKAGE(gui/gui)
ROOT_USE_PACKAGE(tree/treeplayer)
ROOT_USE_PACKAGE(proof/proof)

#ROOT_STANDARD_LIBRARY_PACKAGE(SessionViewer DEPENDENCIES Proof Gui Hist Gpad Graf Tree MathCore)

ROOT_GENERATE_DICTIONARY(G__${libname} T*.h MODULE ${libname} LINKDEF LinkDef.h)

ROOT_LINKER_LIBRARY(${libname} *.cxx G__${libname}.cxx LIBRARIES Vcshared DEPENDENCIES Proof Gui Hist Gpad Graf Tree MathCore)
ROOT_INSTALL_HEADERS()

