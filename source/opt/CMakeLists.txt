message(STATUS "Creating external tools in subdirectory opt")
# convenience macro to set in both current and parent scope:
macro(superset symbol value)
    set(${symbol} "${value}")
    set(${symbol} "${value}" PARENT_SCOPE)
endmacro()

if (NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE Release) 
endif ()

find_package(ZLIB)
#static_find_library(ZLIB zlib.h z)
if  (ZLIB_FOUND)
	include_directories(${ZLIB_INCLUDE_DIRS})
	if (WIN32)
		set  (THIS_ADDITIONAL_LIB ${THIS_ADDITIONAL_LIB} ${ZLIB_LIBRARIES})
	else ()
		# TODO : get zlib library name in Linux
		set  (ZLIB_LIBRARIES z)
	endif()
	message(STATUS "zlib found")
else ()
	message(STATUS "zlib NOT found, building from source")
	
	set(ZLIB_PREFIX "zlib-1.2.8")
	superset(ZLIB_DIR "${CMAKE_CURRENT_BINARY_DIR}/${ZLIB_PREFIX}")
	if (WIN32)
		if (CMAKE_BUILD_TYPE STREQUAL "Debug")
			superset(ZLIB_LIBRARY "${ZLIB_DIR}/${CMAKE_BUILD_TYPE}/zlibstaticd.lib")
		else()
			superset(ZLIB_LIBRARY "${ZLIB_DIR}/${CMAKE_BUILD_TYPE}/zlibstatic.lib")
		endif ()
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
			COMMENT "Building zlib package")
			
		set (ZLIB_FINAL_TASK ${ZLIB_LIBRARY})
	endif ()

	if (WIN32)
		add_custom_command(
			OUTPUT ${ZLIB_LIBRARY}
			COMMAND cmake -G "Visual Studio 12 2013" -A "x64" "${ZLIB_DIR}"
			COMMAND msbuild /p:Configuration=${CMAKE_BUILD_TYPE} /p:Platform=x64 ${ZLIB_DIR}/zlib.sln
			WORKING_DIRECTORY ${ZLIB_DIR}
			DEPENDS ${ZLIB_DIR}
			COMMENT "Building zlib package")
			
		set (ZLIB_FINAL_TASK ${ZLIB_LIBRARY})
	endif()
endif ()

#
# htslib
#
set(HTSLIB_PREFIX "htslib-1.2.1")
superset(HTSLIB_DIR "${CMAKE_CURRENT_BINARY_DIR}/${HTSLIB_PREFIX}")
if (WIN32)
	superset(HTSLIB_LIBRARY "${HTSLIB_DIR}/${CMAKE_BUILD_TYPE}/libhts.lib")
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
        COMMENT "Building htslib library")
		
    set (HTSLIB_FINAL_TASK ${HTSLIB_LIBRARY})
endif ()

if (WIN32)
	add_custom_command(
        OUTPUT ${HTSLIB_LIBRARY}
		WORKING_DIRECTORY ${HTSLIB_DIR}
        COMMAND cmake -G "Visual Studio 12 2013" -A "x64" "${HTSLIB_DIR}" -DBOOST_ROOT:STRING=${Boost_INCLUDE_DIRS} -DZLIB_DIR:STRING=${ZLIB_DIR}
		COMMAND msbuild /p:Configuration=${CMAKE_BUILD_TYPE} /p:Platform=x64 ${HTSLIB_DIR}/Project.sln
		DEPENDS ${ZLIB_LIBRARY}
        DEPENDS ${HTSLIB_DIR}
		COMMENT "Building htslib library")
		
	set (HTSLIB_FINAL_TASK ${HTSLIB_LIBRARY})
endif()

add_custom_target(${THIS_PROJECT}_opt DEPENDS "${HTSLIB_FINAL_TASK}")
