# tetgen package

include(ExternalProject)

set(TETGEN_SOURCE_DIR ${MAIN_SOURCE_DIR}/tetgen)
set(TETGEN_BINARY_DIR ${MAIN_BINARY_DIR}/tetgen)
set(TETGEN_INSTALL_PREFIX ${MAIN_INSTALL_PREFIX}/tetgen)

message(STATUS "Install tetgen release 1.0.0")

ExternalProject_Add(tetgen
    GIT_REPOSITORY https://github.com/fvicini/tetgen.git
    GIT_TAG 0669eb14f527bf0654c347607a5fce24c3492fee # release 1.0.0
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
    SOURCE_DIR ${TETGEN_SOURCE_DIR}
    BINARY_DIR ${TETGEN_BINARY_DIR}
    CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${TETGEN_INSTALL_PREFIX}
    )
