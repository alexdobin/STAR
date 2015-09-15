################################################################################
# Set the compilation flages to report all warnings and enable debugging 
# in the generated objects and executable files
cmake_minimum_required(VERSION 2.8.0)

set(THIS_PROJECT "STAR")
set(THIS_PROGRAM "STAR")
set(THIS_MODULE_DIR "${CMAKE_SOURCE_DIR}/cmake")
set(THIS_REDIST_DIR "${CMAKE_SOURCE_DIR}/../redist")

if (NOT WIN32)
	set(CMAKE_CXX_COMPILER "g++")
	set(CMAKE_CXX_FLAGS "-pipe -std=c++11 -Wall -Wextra -fopenmp")
endif()

#project(${THIS_PROJECT} CXX)

# redist includes
include("${THIS_MODULE_DIR}/boost.cmake")

add_subdirectory(opt)

if (WIN32) 
	add_subdirectory(common)
	add_subdirectory(STAR_ReadFile)
	include_directories (${CMAKE_CURRENT_SOURCE_DIR}/common)
endif()


include_directories (${Boost_INCLUDE_DIRS})
include_directories (${CMAKE_CURRENT_SOURCE_DIR})

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
# exclude SharedMemory.cpp, we have SharedMemorySegment.cpp that implement Boost Shared Memory.
list(REMOVE_ITEM SOURCE_LIST ${CMAKE_CURRENT_SOURCE_DIR}/SharedMemory.cpp)

# exclude IfstreamReadIn.cpp if NOT WIN32, this file implement pipe and ifstream wrapper over that in Windows.  
if (NOT WIN32)
	list(REMOVE_ITEM SOURCE_LIST ${CMAKE_CURRENT_SOURCE_DIR}/IfstreamReadIn.cpp)
	list(REMOVE_ITEM SOURCE_LIST ${CMAKE_CURRENT_SOURCE_DIR}/DirFunctions.cpp)
endif ()

set (LIBS_DEPENDS ${HTSLIB_LIBRARY} ${ZLIB_LIBRARY} ${Boost_LIBRARIES})	  
# add winsock2 lib for Windows.
if (WIN32)
	set (LIBS_DEPENDS ${LIBS_DEPENDS} ws2_32)	  
endif ()
if (WIN32)
	add_executable(${THIS_PROGRAM} ${SOURCE_LIST} $<TARGET_OBJECTS:getopt> $<TARGET_OBJECTS:random>)
else()
	add_executable(${THIS_PROGRAM} ${SOURCE_LIST})
endif()
if (WIN32)
	target_link_libraries (${THIS_PROGRAM} ${LIBS_DEPENDS})
else () 
	target_link_libraries (${THIS_PROGRAM} ${LIBS_DEPENDS} rt)
endif() 
add_dependencies(${THIS_PROGRAM} ${THIS_PROGRAM}_opt) 
