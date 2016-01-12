# Program Name
set(THIS_PROGRAM  "STAR_ReadFile")

include_directories (${ZLIB_DIR})
#for unistd.h 
include_directories(${CMAKE_SOURCE_DIR}/common)

add_executable (${THIS_PROGRAM} STAR_ReadFile.cpp)
target_link_libraries (${THIS_PROGRAM} ${ZLIB_LIBRARY})
add_dependencies(${THIS_PROGRAM} ${THIS_PROJECT}_opt) 

set_target_properties(${THIS_PROGRAM} PROPERTIES  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


