# Program Name
set(THIS_PROGRAM  "STAR_ReadFile")

include_directories (${ZLIB_DIR})

add_executable (${THIS_PROGRAM} STAR_ReadFile.cpp)
add_dependencies(${THIS_PROGRAM} ${ZLIB_LIBRARY})
target_link_libraries (${THIS_PROGRAM} ${ZLIB_LIBRARY})


