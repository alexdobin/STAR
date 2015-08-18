
message(STATUS "Creating external tools in subdirectory opt")

# convenience macro to set in both current and parent scope:
macro(superset symbol value)
    set(${symbol} "${value}")
    set(${symbol} "${value}" PARENT_SCOPE)
endmacro()

# TODO : Get this based on build config 
set (VS_BUILD_CONFIG "Release")
set (THIS_REDIST_DIR ${CMAKE_SOURCE_DIR}/redist)

#
# zlib
#
set(ZLIB_PREFIX "zlib-1.2.8")
superset(ZLIB_DIR "${CMAKE_CURRENT_BINARY_DIR}/${ZLIB_PREFIX}")
if (WIN32)
	superset(ZLIB_LIBRARY "${ZLIB_DIR}/${VS_BUILD_CONFIG}/zlib.lib")
else ()
	# TODO : Check this library name
	superset(ZLIB_LIBRARY "${ZLIB_DIR}/z.a")
endif ()

add_custom_command(
    OUTPUT ${ZLIB_DIR}
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${ZLIB_DIR}"
    COMMAND ${CMAKE_COMMAND} -E tar xjf "${THIS_REDIST_DIR}/${ZLIB_PREFIX}.tar.gz"
	COMMENT "Unpacking zlib package")

set (ZLIB_FINAL_TASK ${ZLIB_DIR})

if (NOT WIN32)
    add_custom_command(
        OUTPUT ${ZLIB_LIBRARY}
        COMMAND $(MAKE) >| zlib.build.log
        WORKING_DIRECTORY ${ZLIB_DIR}
		DEPENDS ${ZLIB_DIR}
        COMMENT "Building zlib package")
		
	set (ZLIB_FINAL_TASK ${ZLIB_LIBRARY})
endif ()

if (WIN32)
	add_custom_command(
        OUTPUT ${ZLIB_LIBRARY}
        COMMAND cmake -G "Visual Studio 12 2013" -A "win32" "${ZLIB_DIR}"
		COMMAND msbuild /p:Configuration=${VS_BUILD_CONFIG} /p:Platform=Win32 ${SAMTOOLS_DIR}/zlib.sln
        DEPENDS ${ZLIB_DIR}
        WORKING_DIRECTORY ${ZLIB_DIR}
        COMMENT "Building zlib package")
		
	set (ZLIB_FINAL_TASK ${ZLIB_LIBRARY})
endif()

add_custom_target(${THIS_PROJECT_NAME}_zlib DEPENDS "${ZLIB_FINAL_TASK}")

#
# htslib
#
set(HTSLIB_PREFIX "htslib-1.2.1")
superset(HTSLIB_DIR "${CMAKE_CURRENT_BINARY_DIR}/${HTSLIB_PREFIX}")
if (WIN32)
	superset(HTSLIB_LIBRARY "${HTSLIB_DIR}/${VS_BUILD_CONFIG}/libhts.lib")
else ()
	superset(HTSLIB_LIBRARY "${HTSLIB_DIR}/libhts.a")
endif ()

add_custom_command(
    OUTPUT ${HTSLIB_DIR}
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${HTSLIB_DIR}"
    COMMAND ${CMAKE_COMMAND} -E tar xjf "${THIS_REDIST_DIR}/${HTSLIB_PREFIX}.tar.bz2"
	COMMENT "Unpacking htslib library")

set (HTSLIB_FINAL_TASK ${HTSLIB_DIR})

if (NOT WIN32)
    add_custom_command(
        OUTPUT ${HTSLIB_LIBRARY}
        COMMAND $(MAKE) -C "${HTSLIB_DIR}" lib-static tabix bgzip 2>| htslib.build.log
        DEPENDS ${HTSLIB_DIR}
		DEPENDS ${ZLIB_LIBRARY}
        COMMENT "Building htslib library")
		
    set (HTSLIB_FINAL_TASK ${HTSLIB_LIBRARY})
endif ()

if (WIN32)
	add_custom_command(
        OUTPUT ${HTSLIB_LIBRARY}
		WORKING_DIRECTORY ${HTSLIB_DIR}
        COMMAND cmake -G "Visual Studio 12 2013" -A "win32" "${HTSLIB_DIR}"
		COMMAND msbuild /p:Configuration=${VS_BUILD_CONFIG} /p:Platform=Win32 ${HTSLIB_DIR}/htslib.sln
        DEPENDS ${HTSLIB_DIR}
		DEPENDS ${ZLIB_LIBRARY}
		COMMENT "Building htslib library")
		
	set (HTSLIB_FINAL_TASK ${HTSLIB_LIBRARY})
endif()

add_custom_target(${THIS_PROJECT_NAME}_htslib DEPENDS "${HTSLIB_FINAL_TASK}")
