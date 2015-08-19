################################################################################
# Set the compilation flages to report all warnings and enable debugging 
# in the generated objects and executable files
cmake_minimum_required(VERSION 2.8.0)

set(THIS_PROGRAM "STAR")
set(CMAKE_C_FLAGS "-Wall -g")
set(THIS_MODULE_DIR "${CMAKE_SOURCE_DIR}/cmake")
set(THIS_REDIST_DIR "${CMAKE_SOURCE_DIR}/../redist")

project(${THIS_PROGRAM} CXX)

# redist includes
#include("${THIS_MODULE_DIR}/boost.cmake")

add_subdirectory(opt) 
add_subdirectory(common)
add_subdirectory(STAR_ReadFile)


include_directories (${Boost_INCLUDE_DIRS})
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/common)
include_directories (${ZLIB_DIR})
include_directories (${HTSLIB_DIR})
include_directories (${CMAKE_CURRENT_BINARY_DIR}/opt)

file (GLOB SOURCE_LIST *.cpp *.c)
foreach (SOURCE_FILE ${SOURCE_LIST})
    get_filename_component(SOURCE_NAME ${SOURCE_FILE} NAME_WE)
	if (WIN32)
		SET_SOURCE_FILES_PROPERTIES( ${SOURCE_FILE} PROPERTIES LANGUAGE CXX )
	endif()
endforeach ()

## Define the target application executable 
##
# exclude SharedMemory.cpp from build in WIN32, new file SharedMemorySegment added for Windows.
if (WIN32)
	list(REMOVE_ITEM SOURCE_LIST ${CMAKE_CURRENT_SOURCE_DIR}/SharedMemory.cpp)
	list(REMOVE_ITEM SOURCE_LIST ${CMAKE_CURRENT_SOURCE_DIR}/SharedMemory.h)
endif ()

set (LIBS_DEPENDS ${HTSLIB_LIBRARY} ${ZLIB_LIBRARY} ${Boost_LIBRARIES})	  
# add winsock2 lib for Windows.
if (WIN32)
	set (LIBS_DEPENDS ${LIBS_DEPENDS} ws2_32)	  
endif ()

add_executable(${THIS_PROGRAM} ${SOURCE_LIST} $<TARGET_OBJECTS:getopt> $<TARGET_OBJECTS:random>)

target_link_libraries (${THIS_PROGRAM} ${LIBS_DEPENDS})
