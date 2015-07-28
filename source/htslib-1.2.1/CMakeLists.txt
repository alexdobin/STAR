#
# Manta
# Copyright (c) 2013-2015 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/sequencing/licenses/>
#

###############################################################################
##
## Configuration file for the htslib-1.2.1 subfolder
##
##
################################################################################

include_directories (BEFORE SYSTEM ${THIS_CXX_BEFORE_SYSTEM_INCLUDES})
include_directories (${THIS_CXX_ALL_INCLUDES})
include_directories (${CMAKE_CURRENT_BINARY_DIR})
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories (${THIS_CXX_CONFIG_H_DIR})
include_directories (${CMAKE_SOURCE_DIR})
include_directories (${CMAKE_SOURCE_DIR}/common)

get_filename_component(CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
message (STATUS "Adding c++ library subdirectory: ${CURRENT_DIR_NAME}")

##
## Some generators (VS) require all targets to be unique across the project.
## Therefore, a unique prefix is needed to create the target names which are
## shared across libraries
##

string(REGEX REPLACE ${CMAKE_SOURCE_DIR}/c[+][+]/ "" TMP1 ${CMAKE_CURRENT_SOURCE_DIR}/)
string(REGEX REPLACE "/" "_" THIS_UNIQUE_PREFIX ${TMP1})

##
## build the library
##

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

# The version number.
set(PACKAGE_VERSION  "1.2.1")
  
FILE(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/version.h
"\#define HTS_VERSION \"@PACKAGE_VERSION@\"\n"
)

file(GLOB_RECURSE THIS_LIBRARY_SOURCES *.cpp *.c)
# TODO : change this based on compile flag to include iRODS library 
list(REMOVE_ITEM THIS_LIBRARY_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/hfile_irods.c)
foreach (SOURCE_FILE ${THIS_LIBRARY_SOURCES})
    get_filename_component(SOURCE_NAME ${SOURCE_FILE} NAME_WE)
	if (WIN32)
		SET_SOURCE_FILES_PROPERTIES( ${SOURCE_FILE} PROPERTIES LANGUAGE CXX )
	endif()
    if (${SOURCE_NAME}_COMPILE_FLAGS)
		set_source_files_properties(${SOURCE_FILE} PROPERTIES COMPILE_FLAGS ${${SOURCE_NAME}_COMPILE_FLAGS})
    endif ()
endforeach ()

# add cram source  
#add_subdirectory(cram) 

# we don't need to add headers to the library for the build to work, but we
# add headers on win32 to make the Visual Studio solution easier to work with:
if (WIN32)
    file(GLOB THIS_LIBRARY_HEADERS *.hh)
    set (THIS_LIBRARY_SOURCES ${THIS_LIBRARY_SOURCES} ${THIS_LIBRARY_HEADERS})
endif ()

if (THIS_LIBRARY_SOURCES)
    set (LIB_TARGET_NAME "${THIS_PROJECT_NAME}_${CURRENT_DIR_NAME}")
    add_library     (${LIB_TARGET_NAME} STATIC ${THIS_LIBRARY_SOURCES})
    #add_dependencies(${LIB_TARGET_NAME} ${THIS_OPT})

    # make the target project use folders when applying cmake IDE generators like Visual Studio
    file(RELATIVE_PATH THIS_RELATIVE_LIBDIR "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    set_property(TARGET ${LIB_TARGET_NAME} PROPERTY FOLDER "${THIS_RELATIVE_LIBDIR}")
endif()

##
## build the unit tests if a "test" subdirectory is found:
##
if (IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/test")
    message (STATUS "Adding c++ test subdirectory:    ${CURRENT_DIR_NAME}/test")
    add_subdirectory (test)
endif ()

