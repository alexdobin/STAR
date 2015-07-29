################################################################################
# Set the compilation flages to report all warnings and enable debugging 
# in the generated objects and executable files
cmake_minimum_required(VERSION 2.8.0)

set(THIS_PROJECT_NAME "STAR")
set(CMAKE_C_FLAGS "-Wall -g")
set(HTSLIB_PREFIX "htslib-1.2.1")
set(THIS_MODULE_DIR "${CMAKE_SOURCE_DIR}/cmake")
set(THIS_REDIST_DIR "${CMAKE_SOURCE_DIR}/../redist")

# redist includes
include("${THIS_MODULE_DIR}/boost.cmake")

add_subdirectory(${HTSLIB_PREFIX})
add_subdirectory(common)

include_directories (${Boost_INCLUDE_DIRS})
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/common)

# required support for gzip compression
find_package(ZLIB)
#static_find_library(ZLIB zlib.h z)
if    (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
	if (WIN32)
		set  (THIS_ADDITIONAL_LIB ${THIS_ADDITIONAL_LIB} ${ZLIB_LIBRARIES})
	else ()
		set  (THIS_ADDITIONAL_LIB ${THIS_ADDITIONAL_LIB} z)
	endif()
	message(STATUS "zlib found")
else  ()
    set(TMP_MSG "zlib library not found")
    if (WIN32)
        set (TMP_MSG "${TMP_MSG}. On win32 this can be installed as part of GnuWin32: http://gnuwin32.sourceforge.net/downlinks/zlib.php")
    endif ()
    message(FATAL_ERROR "${TMP_MSG}")
endif ()

# The following File GLOB command will search for all the files in the current directory
# that matches the specified expressions (*.cpp  and *.h) and assigning 
# the resulting files to the variable SOURCE_LIST 

file (GLOB SOURCE_LIST *.cpp *.h)

## Define the target application executable 
##
# exclude SharedMemory.cpp from build in WIN32, new file SharedMemorySegment added for Windows.
if (WIN32)
list(REMOVE_ITEM SOURCE_LIST ${CMAKE_CURRENT_SOURCE_DIR}/SharedMemory.cpp)
endif ()

add_executable(${THIS_PROJECT_NAME} ${SOURCE_LIST})
set (LIBS_REQUIRED ${THIS_PROJECT_NAME}_${HTSLIB_PREFIX} ${THIS_PROJECT_NAME}_common ${THIS_ADDITIONAL_LIB} ${Boost_LIBRARIES})	  

# add winsock2 lib for Windows.
if (WIN32)
set (LIBS_REQUIRED ${LIBS_REQUIRED} ws2_32)	  
endif ()

target_link_libraries (${THIS_PROJECT_NAME} ${LIBS_REQUIRED})

