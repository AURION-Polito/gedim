cmake_minimum_required(VERSION 3.9)

project(GeDiM LANGUAGES CXX C VERSION 1.0.0)

# Set C++ Standard
###########################################################################################
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create library variables
###############################################################################################
set(GEDIM_LIBRARY_SOURCES ${PROJECT_NAME}_sources)
set(GEDIM_LIBRARY_HEADERS ${PROJECT_NAME}_headers)
set(GEDIM_LIBRARY_PUBLIC_HEADERS ${PROJECT_NAME}_public_headers)
set(GEDIM_LIBRARY_LINKED_LIBRARIES ${PROJECT_NAME}_linked_libraries)
set(GEDIM_LIBRARY_INCLUDE ${PROJECT_NAME}_include)
set(GEDIM_LIBRARY_ADDITIONAL_INCLUDE ${PROJECT_NAME}_additional_include)
set(GEDIM_LIBRARY_COMPILE_DEFINITIONS ${PROJECT_NAME}_compile_definitions)
set(GEDIM_INSTALL_DIR ${PROJECT_NAME})
set(GEDIM_INSTALL_BINARY_DIR bin/)
set(GEDIM_INSTALL_INCLUDE_DIR include/)
set(GEDIM_INSTALL_LIB_DIR lib/)
set(GEDIM_INSTALL_ARCHIVE_DIR lib/)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()
Message(STATUS "Gedim Build configuration: " ${CMAKE_BUILD_TYPE})

# IMPOSE WARNINGS ON DEBUG
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic-errors")

# IMPOSE CXX FLAGS FOR WINDOWS
if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif (WIN32)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${GEDIM_INSTALL_DIR} CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
Message(STATUS "Gedim Library will be installed in: " ${CMAKE_INSTALL_PREFIX})

set(USE_MPI OFF CACHE BOOL "Enable MPI support")
set(ENABLE_METIS OFF CACHE BOOL "Enable Metis extension")
set(ENABLE_VTK OFF CACHE BOOL "Enable VTK extension")
set(ENABLE_MKL OFF CACHE BOOL "Enable MKL extension")
set(ENABLE_VORO OFF CACHE BOOL "Enable VORO extension")
set(ENABLE_SUITESPARSE OFF CACHE BOOL "Enable SuiteSparse extension")
set(ENABLE_PETSC OFF CACHE BOOL "Enable PETSc extension")
set(VERBOSE 3 CACHE STRING "Select verbosity level (0, 1, 2 3 or 4)")
set_property(CACHE VERBOSE PROPERTY STRINGS 0 1 2 3 4)
if(NOT ${VERBOSE} MATCHES 0 AND NOT ${VERBOSE} MATCHES 1 AND NOT ${VERBOSE} MATCHES 2 AND NOT ${VERBOSE} MATCHES 3 AND NOT ${VERBOSE} MATCHES 4)
    message(SEND_ERROR "VERBOSE variable must be either 0 (None), 1 (Error), 2 (Warning), 3 (Info) or 4 (Debug)")
endif()
set(LOGGING 3 CACHE STRING "Select logging level (0, 1, 2 or 3)")
set_property(CACHE LOGGING PROPERTY STRINGS 0 1 2 3)
if(NOT ${LOGGING} MATCHES 0 AND NOT ${LOGGING} MATCHES 1 AND NOT ${LOGGING} MATCHES 2 AND NOT ${LOGGING} MATCHES 3)
    message(SEND_ERROR "LOGGING variable must be either 0 (None), 1 (Console), 2 (LogFile) or 3 (Console and LogFile)")
endif()

# Add dependencies
###########################################################################################
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

## MPI
if(${USE_MPI})
    find_package(MPI REQUIRED)
    list(APPEND ${GEDIM_LIBRARY_LINKED_LIBRARIES} PUBLIC MPI::MPI_CXX)
endif()

## Eigen3
find_package(Eigen3 CONFIG REQUIRED)
list(APPEND ${GEDIM_LIBRARY_LINKED_LIBRARIES} PUBLIC Eigen3::Eigen)

## Blas
find_package(BLAS REQUIRED)
list(APPEND ${GEDIM_LIBRARY_LINKED_LIBRARIES} PUBLIC BLAS::BLAS)

## Lapack
find_package(LAPACK REQUIRED)
list(APPEND ${GEDIM_LIBRARY_LINKED_LIBRARIES} PUBLIC LAPACK::LAPACK)

## MKL
if(${ENABLE_MKL})
    set(MKL_INTERFACE lp64) # Eigen requires lp64 interface
    find_package(MKL CONFIG REQUIRED)
    list(APPEND ${GEDIM_LIBRARY_LINKED_LIBRARIES} MKL::MKL)
endif()

## triangle
find_package(triangle CONFIG REQUIRED)
list(APPEND ${GEDIM_LIBRARY_LINKED_LIBRARIES} PUBLIC triangle::triangle)

## tetgen
find_package(tetgen CONFIG REQUIRED)
list(APPEND ${GEDIM_LIBRARY_LINKED_LIBRARIES} PUBLIC tetgen::tetgen)

## VORO
if(${ENABLE_VORO})
    find_package(VORO CONFIG REQUIRED)
    list(APPEND ${GEDIM_LIBRARY_LINKED_LIBRARIES} PUBLIC VORO::voro++)
endif()

## VTK
if(${ENABLE_VTK})
    find_package(VTK CONFIG REQUIRED)
    list(APPEND ${GEDIM_LIBRARY_LINKED_LIBRARIES} PUBLIC ${VTK_LIBRARIES})
endif()

## METIS
if(${ENABLE_METIS})
    find_package(metis CONFIG REQUIRED)
    list(APPEND ${GEDIM_LIBRARY_LINKED_LIBRARIES} METIS::metis)
endif()

## SuiteSparse
if(${ENABLE_SUITESPARSE})
    find_package(SuiteSparse CONFIG REQUIRED)
    message("Found SuiteSparse: ${SUITESPARSE_LIBRARIES}")
    list(APPEND ${GEDIM_LIBRARY_LINKED_LIBRARIES}  ${SUITESPARSE_LIBRARIES})
    list(APPEND ${GEDIM_LIBRARY_INCLUDE}  ${SUITESPARSE_INCLUDE_DIRS})
endif()

if(${ENABLE_PETSC})
    # PkgConfig
    find_package(PkgConfig)

    # PETSc
    if (PKG_CONFIG_FOUND)
        pkg_check_modules(PETSC PETSc)
    endif()

    if (PETSC_FOUND)
        message("Found PETSc: ${PETSC_LINK_LIBRARIES}")
        list(APPEND COMPILE_OPTIONS ${PETSC_CFLAGS})
        list(APPEND ${GEDIM_LIBRARY_INCLUDE}  ${PETSC_INCLUDE_DIRS})
        list(APPEND ${GEDIM_LIBRARY_LINKED_LIBRARIES}  ${PETSC_LINK_LIBRARIES})
        set(LINK_FLAGS "${LINK_FLAGS} ${PETSC_LDFLAGS}")
        set(CMAKE_REQUIRED_FLAGS ${PETSC_CFLAGS})
        set(CMAKE_REQUIRED_INCLUDES "${PETSC_INCLUDE_DIRS}")
    endif()
endif()

# Insert Sources
###############################################################################################
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)
list(APPEND ${GEDIM_LIBRARY_SOURCES} ${gedim_src_sources})
list(APPEND ${GEDIM_LIBRARY_HEADERS} ${gedim_src_headers})
list(APPEND ${GEDIM_LIBRARY_INCLUDE} ${gedim_src_includes})

# Create public headers
###############################################################################################
set(GEDIM_LIBRARY_ADDITIONAL_INCLUDE_PATH ${CMAKE_CURRENT_BINARY_DIR}/${GEDIM_LIBRARY_ADDITIONAL_INCLUDE})

foreach(header ${${GEDIM_LIBRARY_HEADERS}})
    get_filename_component(headerDirectory ${header} DIRECTORY)
    get_filename_component(headerFileName ${header} NAME_WE)
    get_filename_component(headerExtension ${header} EXT)

    if(${headerExtension} MATCHES ".hpp.in")
        set(headerExtension ".hpp")
        set(headerDirectory ${GEDIM_LIBRARY_ADDITIONAL_INCLUDE_PATH})

        configure_file(${header} ${headerDirectory}/${headerFileName}${headerExtension} @ONLY)
    elseif(NOT ${headerExtension} MATCHES ".hpp" AND NOT ${headerExtension} MATCHES ".h")
        message(SEND_ERROR "Header extension not recognized for file ${headerFileName}${headerExtension}")
    endif()

    list(APPEND ${GEDIM_LIBRARY_PUBLIC_HEADERS} ${headerDirectory}/${headerFileName}${headerExtension})
endforeach()

list(APPEND ${GEDIM_LIBRARY_INCLUDE} ${GEDIM_LIBRARY_ADDITIONAL_INCLUDE_PATH})

# Add variables for clang anc cppcheck
###############################################################################################
list(APPEND CLANGFORMAT_FILES ${${GEDIM_LIBRARY_SOURCES}})
list(APPEND CLANGFORMAT_FILES ${${GEDIM_LIBRARY_HEADERS}})
set(CLANGFORMAT_FILES ${CLANGFORMAT_FILES} PARENT_SCOPE)

list(APPEND CPPCHECK_FILES ${${GEDIM_LIBRARY_SOURCES}})
list(APPEND CPPCHECK_FILES ${${GEDIM_LIBRARY_HEADERS}})
set(CPPCHECK_FILES ${CPPCHECK_FILES} PARENT_SCOPE)

# Create library
###############################################################################################
add_library(${PROJECT_NAME} STATIC ${${GEDIM_LIBRARY_SOURCES}} ${${GEDIM_LIBRARY_HEADERS}})
add_library(GeDiM::GeDiM ALIAS ${PROJECT_NAME})

target_link_libraries(${PROJECT_NAME} ${${GEDIM_LIBRARY_LINKED_LIBRARIES}})

target_include_directories(${PROJECT_NAME} PRIVATE ${${GEDIM_LIBRARY_INCLUDE}})
target_include_directories(${PROJECT_NAME} PUBLIC $<INSTALL_INTERFACE:${GEDIM_INSTALL_INCLUDE_DIR}>)
target_compile_definitions(${PROJECT_NAME} PUBLIC ${${GEDIM_LIBRARY_COMPILE_DEFINITIONS}})

target_compile_options(${PROJECT_NAME} PUBLIC -fPIC)

# Create install target
###############################################################################################
set_target_properties(${PROJECT_NAME} PROPERTIES 
    PUBLIC_HEADER "${${GEDIM_LIBRARY_PUBLIC_HEADERS}}"
    EXPORT_NAME GeDiM)

install(
    TARGETS ${PROJECT_NAME}
    EXPORT GeDiMTargets
    RUNTIME DESTINATION ${GEDIM_INSTALL_BINARY_DIR} # runtime executable
    LIBRARY DESTINATION ${GEDIM_INSTALL_LIB_DIR} # dynamic libraries
    ARCHIVE DESTINATION ${GEDIM_INSTALL_ARCHIVE_DIR} # static libraries
    PUBLIC_HEADER DESTINATION ${GEDIM_INSTALL_INCLUDE_DIR} # headers
)


include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/GeDiMConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/GeDiMConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_PREFIX}/${GEDIM_INSTALL_LIB_DIR}/cmake/${GEDIM_INSTALL_DIR}
    PATH_VARS GEDIM_INSTALL_INCLUDE_DIR
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/GeDiMConfigVersion.cmake
    VERSION ${GeDiM_VERSION}
    COMPATIBILITY AnyNewerVersion
)

export(
    TARGETS ${PROJECT_NAME}
    NAMESPACE GeDiM::
    FILE ${CMAKE_CURRENT_BINARY_DIR}/GeDiMTargets.cmake
)

install(
    EXPORT GeDiMTargets
    NAMESPACE GeDiM::
    DESTINATION ${CMAKE_INSTALL_PREFIX}/${GEDIM_INSTALL_LIB_DIR}/cmake/${GEDIM_INSTALL_DIR}
)

install(
    FILES
    ${CMAKE_CURRENT_BINARY_DIR}/GeDiMConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/GeDiMConfigVersion.cmake
    DESTINATION
    ${CMAKE_INSTALL_PREFIX}/${GEDIM_INSTALL_LIB_DIR}/cmake/${GEDIM_INSTALL_DIR}
)
