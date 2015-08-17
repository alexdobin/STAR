# CMake configuration for getopt library
project("getopt" C)

# add __GNU_LIBRARY__ and HAVE_STRING_H preprocessors for getopt  
if (WIN32)
add_definitions(-D__GNU_LIBRARY__ -DHAVE_STRING_H=1)
endif()

include_directories("${CMAKE_SOURCE_DIR}/common")

file (GLOB GET_OPT_SOURCE_LIST [a-zA-Z0-9] *.c *.h)
add_library(getopt OBJECT ${GET_OPT_SOURCE_LIST})