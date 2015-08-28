############################################################################
# CMakeLists.txt file for building ROOT geom/gdml package
############################################################################
set(libname Gdml)

ROOT_USE_PACKAGE(geom/geom)
ROOT_USE_PACKAGE(io/xml)

ROOT_GENERATE_DICTIONARY(G__${libname} *.h MODULE ${libname} LINKDEF LinkDef.h OPTIONS "-writeEmptyRootPCM")

#ROOT_LINKER_LIBRARY(Gdml *.cxx G__Gdml.cxx LIBRARIES Vcshared ${_srcs} DEPENDENCIES Geom XMLIO Hist RIO)

ROOT_LINKER_LIBRARY(${libname} *.cxx G__${libname}.cxx LIBRARIES Vcshared DEPENDENCIES Geom XMLIO Hist RIO)
ROOT_INSTALL_HEADERS()
