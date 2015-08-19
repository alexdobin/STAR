################################################################################
# Set the compilation flages to report all warnings and enable debugging 
# in the generated objects and executable files
cmake_minimum_required(VERSION 2.8.0)

set(THIS_PROGRAM "STAR")
set(CMAKE_C_FLAGS "-Wall -g")
set(HTSLIB_PREFIX "htslib-1.2.1")
set(THIS_MODULE_DIR "${CMAKE_SOURCE_DIR}/cmake")
set(THIS_REDIST_DIR "${CMAKE_SOURCE_DIR}/../redist")

project(${THIS_PROGRAM} CXX)

# redist includes
#include("${THIS_MODULE_DIR}/boost.cmake")

#add_subdirectory(${HTSLIB_PREFIX})
add_subdirectory(common)
add_subdirectory(STAR_ReadFile)
add_subdirectory(opt) 

include_directories (${Boost_INCLUDE_DIRS})
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/common)

# required support for gzip compression
# find_package(ZLIB)
# #static_find_library(ZLIB zlib.h z)
# if    (ZLIB_FOUND)
    # include_directories(${ZLIB_INCLUDE_DIRS})
	# if (WIN32)
		# set  (THIS_ADDITIONAL_LIB ${THIS_ADDITIONAL_LIB} ${ZLIB_LIBRARIES})
	# else ()
		# set  (THIS_ADDITIONAL_LIB ${THIS_ADDITIONAL_LIB} z)
	# endif()
	# message(STATUS "zlib found")
# else  ()
    # set(TMP_MSG "zlib library not found")
    # if (WIN32)
        # set (TMP_MSG "${TMP_MSG}. On win32 this can be installed as part of GnuWin32: http://gnuwin32.sourceforge.net/downlinks/zlib.php")
    # endif ()
    # #message(FATAL_ERROR "${TMP_MSG}")
	# set (ZLIB_INCLUDE_DIRS "C:/Program Files (x86)/zlib-1.2.5-bin-x64/zlib/include")
	# set (ZLIB_LIBRARIES "C:/Program Files (x86)/zlib-1.2.5-bin-x64/zlib/lib")
	
	# include_directories(${ZLIB_INCLUDE_DIRS})
	# if (WIN32)
		# set  (THIS_ADDITIONAL_LIB ${THIS_ADDITIONAL_LIB} ${ZLIB_LIBRARIES})
	# else ()
		# set  (THIS_ADDITIONAL_LIB ${THIS_ADDITIONAL_LIB} z)
	# endif()
	
# endif ()

# The following File GLOB command will search for all the files in the current directory
# that matches the specified expressions (*.cpp  and *.h) and assigning 
# the resulting files to the variable SOURCE_LIST 

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

set (LIBS_REQUIRED ${HTSLIB_LIBRARY} ${THIS_ADDITIONAL_LIB} ${Boost_LIBRARIES})	  
# add winsock2 lib for Windows.
if (WIN32)
	set (LIBS_REQUIRED ${LIBS_REQUIRED} ws2_32)	  
endif ()

add_executable(${THIS_PROGRAM} ${SOURCE_LIST} $<TARGET_OBJECTS:getopt> $<TARGET_OBJECTS:random>)

target_link_libraries (${THIS_PROGRAM} ${LIBS_REQUIRED})
