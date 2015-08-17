
# Program Name
set(THIS_PROGRAM  "STAR_ReadFile")
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# required support for gzip compression
find_package(ZLIB)
#static_find_library(ZLIB zlib.h z)
if  (ZLIB_FOUND)
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

file (GLOB THIS_PROGRAM_SOURCE_LIST [a-zA-Z0-9] *.c *.cpp)
# we don't need to add headers to the library for the build to work, but we
# add headers on win32 to make the Visual Studio solution easier to work with:
if (WIN32)
    file(GLOB THIS_PROGRAM_HEADERS *.h)
    set (THIS_PROGRAM_SOURCE_LIST ${THIS_PROGRAM_SOURCE_LIST} ${THIS_PROGRAM_HEADERS})
endif ()

if (THIS_PROGRAM_SOURCE_LIST)
    add_executable (${THIS_PROGRAM} ${THIS_PROGRAM_SOURCE_LIST})
	target_link_libraries (${THIS_PROGRAM} ${THIS_ADDITIONAL_LIB})
endif()


