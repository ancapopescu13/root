############################################################################
# CMakeLists.txt file for building ROOT html package
############################################################################

ROOT_USE_PACKAGE(graf2d/gpad)
ROOT_USE_PACKAGE(core/thread)

ROOT_STANDARD_LIBRARY_PACKAGE(Html LIBRARIES Vcshared ${_srcs} DEPENDENCIES Graf Thread)

