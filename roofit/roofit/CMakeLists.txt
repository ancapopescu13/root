############################################################################
# CMakeLists.txt file for building ROOT rootfit package
# @author Pere Mato, CERN
############################################################################
set(libname RooFit)

ROOT_GENERATE_DICTIONARY(G__RooFit *.h MODULE RooFit LINKDEF LinkDef1.h OPTIONS "-writeEmptyRootPCM")

ROOT_LINKER_LIBRARY(RooFit  *.cxx G__RooFit.cxx LIBRARIES Vcshared Core 
                           DEPENDENCIES RooFitCore Tree RIO Matrix MathCore ${ROOT_MATHMORE_LIBRARY} )
ROOT_INSTALL_HEADERS()
