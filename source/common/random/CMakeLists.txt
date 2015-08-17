# CMake configuration for common/random folder.
#TODO : Correct this path by checking if Boost installed and using BOOST_INCLUDE_DIR variable
set(BOOST_INCLUDE_DIR "${CMAKE_BINARY_DIR}/bootstrap/boost/include/boost-1_56")

include_directories(${BOOST_INCLUDE_DIR})
file (GLOB RANDOM_SOURCE_LIST [a-zA-Z0-9] *.cpp *.h)
add_library(random OBJECT ${RANDOM_SOURCE_LIST})