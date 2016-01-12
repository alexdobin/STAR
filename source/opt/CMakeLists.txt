message(STATUS "Creating external tools in subdirectory opt")
# convenience macro to set in both current and parent scope:
macro(superset symbol value)
    set(${symbol} "${value}")
    set(${symbol} "${value}" PARENT_SCOPE)
endmacro()

if (NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE Release) 
endif ()

set (THIS_PROJECT_NAME "STAR")

#
# zlib
#
set(ZLIB_PREFIX "zlib-1.2.8")
superset(ZLIB_DIR "${CMAKE_CURRENT_BINARY_DIR}/${ZLIB_PREFIX}")
if (WIN32)
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		superset (ZLIB_LIBRARY ${ZLIB_DIR}/${CMAKE_BUILD_TYPE}/zlibstaticd.lib)
	else()
		superset (ZLIB_LIBRARY ${ZLIB_DIR}/${CMAKE_BUILD_TYPE}/zlibstatic.lib)
	endif()
else ()
    superset(ZLIB_LIBRARY "${ZLIB_DIR}/libz.a")
endif ()

add_custom_command(
	OUTPUT ${ZLIB_DIR}
	COMMAND ${CMAKE_COMMAND} -E remove_directory "${ZLIB_DIR}"
	COMMAND ${CMAKE_COMMAND} -E tar xjf "${THIS_REDIST_DIR}/${ZLIB_PREFIX}.tar.bz2"
	COMMENT "Unpacking zlib package")

if (NOT WIN32)
    set (ASM_MAKE_OPTION LOC=-DASMV OBJA=match.o)
    if (THIS_ARCH MATCHES "^.*86$")
        set(ZLIB_ARCH "")
        set(ASM_FILE_COPY cp contrib/asm686/match.S match.S)
    elseif (THIS_ARCH MATCHES "^x86_64$")
        set(ZLIB_ARCH "--64")
        set(ASM_FILE_COPY cp contrib/amd64/amd64-match.S match.S)
    else ()
        set (ZLIB_ARCH "")
        set (ASM_FILE_COPY "")
        set (ASM_MAKE_OPTION "")
    endif ()

    add_custom_command(
        OUTPUT ${ZLIB_LIBRARY}
        COMMAND CC=${CMAKE_C_COMPILER} ${ZLIB_DIR}/configure --prefix="${ZLIB_DIR}" --static ${ZLIB_ARCH} >zlib.config.log
        COMMAND ${ASM_FILE_COPY}
        COMMAND $(MAKE) ${ASM_MAKE_OPTION} >zlib.build.log 2>zlib.build.error.log
        WORKING_DIRECTORY ${ZLIB_DIR}
        DEPENDS ${ZLIB_DIR}
        COMMENT "Building zlib package")
else ()
    # Sometimes the platform (x64,win32...) is specified as part of the target name
    # this captures the cases where it is specified separately:
    set (ZLIB_PLATFORM "")
    if (CMAKE_GENERATOR_PLATFORM)
        set (ZLIB_PLATFORM -A "${CMAKE_GENERATOR_PLATFORM}")
    endif ()
    add_custom_command(
        OUTPUT ${ZLIB_LIBRARY}
        COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" ${ZLIB_PLATFORM} -DCMAKE_C_COMPILER="${CMAKE_C_COMPILER}" ${ZLIB_DIR} >zlib.config.log
        COMMAND ${CMAKE_MAKE_PROGRAM} "${ZLIB_DIR}/zlib.sln" /p:Configuration=${CMAKE_BUILD_TYPE} >zlib.build.log
        WORKING_DIRECTORY ${ZLIB_DIR}
        DEPENDS ${ZLIB_DIR}
        COMMENT "Building zlib package")
endif ()

set(THIS_ZLIB "${THIS_PROJECT_NAME}_zlib")
add_custom_target(${THIS_ZLIB} DEPENDS "${ZLIB_LIBRARY}")

#
# htslib
#
set(HTSLIB_PREFIX "htslib")
superset(HTSLIB_DIR "${CMAKE_CURRENT_BINARY_DIR}/${HTSLIB_PREFIX}")
if (WIN32)
	superset(HTSLIB_LIBRARY "${HTSLIB_DIR}/${CMAKE_BUILD_TYPE}/hts.lib")
else ()
	superset(HTSLIB_LIBRARY "${HTSLIB_DIR}/libhts.a")
endif ()

add_custom_command(
    OUTPUT ${HTSLIB_DIR}
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${HTSLIB_DIR}"
    COMMAND ${CMAKE_COMMAND} -E tar xjf "${THIS_REDIST_DIR}/${HTSLIB_PREFIX}.tar.bz2"
    DEPENDS ${THIS_ZLIB}
    COMMENT "Unpacking htslib library")

set (HTSLIB_FINAL_TASK ${HTSLIB_DIR})

if (NOT WIN32)
    # TODO : note that htslib ./configure CFLAGS="custom" will overwrite the standard "-g -O2 -Wall" CFLAGS, so
    # we need to restore at least the -O2 for reasonable library performance:
    add_custom_command(
	OUTPUT ${HTSLIB_LIBRARY}
	COMMAND ${CMAKE_COMMAND} ${HTSLIB_DIR} -DZLIB_DIR:STRING=${ZLIB_DIR}
	COMMAND $(MAKE) >htslib.build.log 2>htslib.build.error.log
	WORKING_DIRECTORY ${HTSLIB_DIR}
	DEPENDS ${HTSLIB_DIR}
	COMMENT "Building htslib library")
    set (HTSLIB_FINAL_TASK ${HTSLIB_LIBRARY})
else ()
	add_custom_command(
        OUTPUT ${HTSLIB_LIBRARY}
		WORKING_DIRECTORY ${HTSLIB_DIR}
        COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" -A "x64" "${HTSLIB_DIR}" -DBOOST_ROOT:STRING=${Boost_INCLUDE_DIRS} -DZLIB_DIR:STRING=${ZLIB_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
		COMMAND msbuild /p:Configuration=${CMAKE_BUILD_TYPE} /p:Platform=x64 ${HTSLIB_DIR}/Project.sln
		DEPENDS ${ZLIB_LIBRARY}
        DEPENDS ${HTSLIB_DIR}
		COMMENT "Building htslib library")
		
	set (HTSLIB_FINAL_TASK ${HTSLIB_LIBRARY})
endif()

add_custom_target(${THIS_PROJECT}_opt DEPENDS "${HTSLIB_FINAL_TASK}")
