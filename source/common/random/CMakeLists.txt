# CMake configuration for common/random folder.
include_directories(${Boost_INCLUDE_DIRS})

include_directories(${BOOST_INCLUDE_DIR})
file (GLOB RANDOM_SOURCE_LIST [a-zA-Z0-9] *.cpp *.h)
add_library(random OBJECT ${RANDOM_SOURCE_LIST})