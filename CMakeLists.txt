cmake_minimum_required(VERSION 3.3)

project(Gedim LANGUAGES CXX C VERSION 0.4)

# CACHE VARIABLES DEFINITION
set(GEDIM_EXTERNALS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external CACHE STRING "Directory where GeDiM external components were built (absolute path)")
set(USE_MPI OFF CACHE BOOL "Enable MPI support")
set(MPI_COMPILER_C mpicc CACHE STRING "C compiler used when USE_MPI is ON")
set(MPI_COMPILER_CXX mpic++ CACHE STRING "C++ compiler used when USE_MPI is ON")

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

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()
Message(STATUS "Build configuration: " ${CMAKE_BUILD_TYPE})

 set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/library/ CACHE STRING "Install path prefix, prepended onto install directories." FORCE)
Message(STATUS "Library will be installed in: " ${CMAKE_INSTALL_PREFIX})

# IMPOSE CXX STANDARD
set(CMAKE_CXX_STANDARD 11)

# ADD EXTERNAL DIRECTORY
add_subdirectory(${GEDIM_EXTERNALS_DIR})

# ADD SOURCES
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)

# ADD TESTS
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/test)

# ADD Eigen
#######################################
list(APPEND ${PROJECT_NAME}_include ${GEDIM_EXTERNALS_DIR}/eigen)

# ADD Google Testing Environment
####################################### Start
add_definitions(-DGTEST_LANGUAGE_CXX11)

find_package(Threads REQUIRED)
include(GoogleTest)

if ($ENV{GOOGLETEST_DIR})
    set(GOOGLETEST_DIR $ENV{GOOGLETEST_DIR})
else ()
    if (NOT ${GEDIM_EXTERNALS_DIR}/googletest STREQUAL "")
        message(VERBOSE "Using googletest src dir specified at Qt Creator wizard")
    endif ()
    set(GOOGLETEST_DIR ${GEDIM_EXTERNALS_DIR}/googletest)
endif ()

if (EXISTS ${GOOGLETEST_DIR})
    set(GTestSrc ${GOOGLETEST_DIR}/googletest)
    set(GMockSrc ${GOOGLETEST_DIR}/googlemock)
elseif (UNIX AND EXISTS /usr/src/gtest)
    set(GTestSrc /usr/src/gtest)
    message(WARNING "Using gtest from system")
    if (EXISTS /usr/src/gmock)
        set(GMockSrc /usr/src/gmock)
    endif ()
else ()
    message(FATAL_ERROR "No googletest src dir found - set GOOGLETEST_DIR to enable!")
endif ()

set(GTestFiles ${GTestSrc}/src/gtest-all.cc)
set(GTestIncludes ${GTestSrc} ${GTestSrc}/include)
if (NOT ${GMockSrc} STREQUAL "")
    list(APPEND GTestFiles ${GMockSrc}/src/gmock-all.cc)
    list(APPEND GTestIncludes ${GMockSrc} ${GMockSrc}/include)
endif ()

list(APPEND ${PROJECT_NAME}_unittest_include ${GTestIncludes})
####################################### End

# ADD external libs
set(${PROJECT_NAME}_libs ${${PROJECT_NAME}_libs} "test")

# ADD code
list(APPEND ${PROJECT_NAME}_include ${src_include})
list(APPEND ${PROJECT_NAME}_unittest_include ${test_include})

# CONFIGURE PUBLIC_HEADERS FOR LIBRARY INSTALL
####################################### Start
set(${PROJECT_NAME}_additional_include ${CMAKE_CURRENT_BINARY_DIR}/additional_include)

foreach(header ${src_headers})
	get_filename_component(headerDirectory ${header} DIRECTORY)
	get_filename_component(headerFileName ${header} NAME_WE)
	get_filename_component(headerExtension ${header} EXT)
	  
  if(${headerExtension} MATCHES ".hpp.in")
    set(headerExtension ".hpp")
    set(headerDirectory ${${PROJECT_NAME}_additional_include})
    
    configure_file(${header} ${headerDirectory}/${headerFileName}${headerExtension})
  elseif(NOT ${headerExtension} MATCHES ".hpp" AND NOT ${headerExtension} MATCHES ".h")
    message(SEND_ERROR "Header extension not recognized for file ${headerFileName}${headerExtension}")
  endif()
    
  list(APPEND ${PROJECT_NAME}_public_headers ${headerDirectory}/${headerFileName}${headerExtension})
endforeach()

list(APPEND ${PROJECT_NAME}_include ${${PROJECT_NAME}_additional_include})
####################################### End

# Add additional files

# CREATE LIBRARY
add_library(${PROJECT_NAME} STATIC ${src_sources} ${src_headers})

target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${${PROJECT_NAME}_include})
if(${${PROJECT_NAME}_libs})
	add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_libs})
	target_link_libraries(${PROJECT_NAME} ${${PROJECT_NAME}_libs})
endif()
target_compile_options(${PROJECT_NAME} PUBLIC -fPIC )
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${${PROJECT_NAME}_public_headers}")

# ADD INSTALL DIRECTIVE
install(TARGETS ${PROJECT_NAME}
		RUNTIME DESTINATION bin # runtime executable
        LIBRARY DESTINATION lib # dynamic libraries
        ARCHIVE DESTINATION lib/static # static libraries
        PUBLIC_HEADER DESTINATION include # headers
)

# CREATE TESTS EXECUTABLE
add_executable(${PROJECT_NAME}_unitTests EXCLUDE_FROM_ALL mainUnitTests.cpp ${src_sources} ${src_headers} ${test_sources} ${test_headers} ${GTestFiles})

target_include_directories(${PROJECT_NAME}_unitTests SYSTEM PUBLIC ${${PROJECT_NAME}_include})
target_include_directories(${PROJECT_NAME}_unitTests SYSTEM PRIVATE ${${PROJECT_NAME}_unittest_include})
if(${${PROJECT_NAME}_libs})
	add_dependencies(${PROJECT_NAME}_unitTests ${${PROJECT_NAME}_libs})
endif()
target_link_libraries(${PROJECT_NAME}_unitTests Threads::Threads)

# Set Unit Tests
gtest_discover_tests(${PROJECT_NAME}_unitTests)

# the parent scope (if any) will have the ${${PROJECT_NAME}_include} containing the directories to be included for using GeDiM
get_directory_property(parentDir PARENT_DIRECTORY)
string(LENGTH "${parentDir}" hasParent)
if(${hasParent})
  set(${PROJECT_NAME}_include ${${PROJECT_NAME}_include} PARENT_SCOPE)
  set(${PROJECT_NAME}_libs ${${PROJECT_NAME}_libs} PARENT_SCOPE)
  set(${PROJECT_NAME}_unittest_include ${${PROJECT_NAME}_unittest_include} PARENT_SCOPE)
endif()
