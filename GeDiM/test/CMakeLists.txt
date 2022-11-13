cmake_minimum_required(VERSION 3.9)

project(GedimUnitTest LANGUAGES CXX C)

# Create library test variables
###############################################################################################
set(GEDIM_LIBRARY_TEST_SOURCES ${PROJECT_NAME}_sources)
set(GEDIM_LIBRARY_TEST_HEADERS ${PROJECT_NAME}_headers)
set(GEDIM_LIBRARY_TEST_PUBLIC_HEADERS ${PROJECT_NAME}_public_headers)
set(GEDIM_LIBRARY_TEST_INCLUDE ${PROJECT_NAME}_include)
set(GEDIM_LIBRARY_TEST_LINKED_LIBRARIES ${PROJECT_NAME}_linked_libraries)

# Add dependencies
###############################################################################################
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

## Threads
find_package(Threads REQUIRED)
list(APPEND ${GEDIM_LIBRARY_TEST_LINKED_LIBRARIES} PRIVATE Threads::Threads)

## GTest
if (CMAKE_VERSION VERSION_LESS 3.20)
    find_package(GTest REQUIRED)
    find_package(GMock REQUIRED)

    list(APPEND ${GEDIM_LIBRARY_TEST_LINKED_LIBRARIES} PRIVATE GTest::GTest)
    list(APPEND ${GEDIM_LIBRARY_TEST_LINKED_LIBRARIES} PRIVATE GTest::Main)
    list(APPEND ${GEDIM_LIBRARY_TEST_LINKED_LIBRARIES} PRIVATE ${GMOCK_BOTH_LIBRARIES})
else()
    find_package(GTest REQUIRED)

    list(APPEND ${GEDIM_LIBRARY_TEST_LINKED_LIBRARIES} PRIVATE GTest::gtest)
    list(APPEND ${GEDIM_LIBRARY_TEST_LINKED_LIBRARIES} PRIVATE GTest::gtest_main)
    list(APPEND ${GEDIM_LIBRARY_TEST_LINKED_LIBRARIES} PRIVATE GTest::gmock)
    list(APPEND ${GEDIM_LIBRARY_TEST_LINKED_LIBRARIES} PRIVATE GTest::gmock_main)
endif()


## GeDiM
if (NOT TARGET GeDiM::GeDiM)
	find_package(GeDiM CONFIG REQUIRED)
	list(APPEND ${GEDIM_LIBRARY_TEST_LINKED_LIBRARIES} PRIVATE GeDiM::GeDiM)
else()
	get_target_property(GEDIM_LIBRARY_HEADERS GeDiM::GeDiM PUBLIC_HEADER)
	get_target_property(GEDIM_LIBRARY_NAME GeDiM::GeDiM EXPORT_NAME)
	get_target_property(GEDIM_LIBRARY_LINKED_LIBRARIES GeDiM::GeDiM LINK_LIBRARIES)
	get_target_property(GEDIM_LIBRARY_INCLUDE_DIRECTORIES GeDiM::GeDiM INCLUDE_DIRECTORIES)
	get_target_property(GEDIM_LIBRARY_SOURCES GeDiM::GeDiM SOURCES)

	list(APPEND ${GEDIM_LIBRARY_TEST_LINKED_LIBRARIES} PRIVATE ${GEDIM_LIBRARY_NAME} ${GEDIM_LIBRARY_LINKED_LIBRARIES})
	list(APPEND ${GEDIM_LIBRARY_TEST_HEADERS} ${GEDIM_LIBRARY_HEADERS})
	list(APPEND ${GEDIM_LIBRARY_TEST_INCLUDE} ${GEDIM_LIBRARY_INCLUDE_DIRECTORIES})
	list(APPEND ${GEDIM_LIBRARY_TEST_LINKED_LIBRARIES} PRIVATE GeDiM::GeDiM)
endif()

## Eigen3
find_package(Eigen3 CONFIG REQUIRED)

## Blas
find_package(BLAS REQUIRED)

## Lapack
find_package(LAPACK REQUIRED)

## VTK
find_package(VTK CONFIG REQUIRED)

## triangle
find_package(triangle CONFIG REQUIRED)

## tetgen
find_package(tetgen CONFIG REQUIRED)

# Insert Sources
###############################################################################################
add_subdirectory(Geometry)
add_subdirectory(Algebra)
add_subdirectory(Mesh)
add_subdirectory(IO)
add_subdirectory(Quadrature)

list(APPEND ${GEDIM_LIBRARY_TEST_PUBLIC_HEADERS} ${gedim_test_src_headers})
list(APPEND ${GEDIM_LIBRARY_TEST_SOURCES} ${gedim_test_src_sources})
list(APPEND ${GEDIM_LIBRARY_TEST_HEADERS} ${gedim_test_src_headers})
list(APPEND ${GEDIM_LIBRARY_TEST_INCLUDE} ${gedim_test_src_includes})

# Create library test executable
###############################################################################################
add_executable(${PROJECT_NAME}
    main.cpp
    ${${GEDIM_LIBRARY_TEST_SOURCES}}
    ${${GEDIM_LIBRARY_TEST_HEADERS}})

target_link_libraries(${PROJECT_NAME} ${${GEDIM_LIBRARY_TEST_LINKED_LIBRARIES}})
target_include_directories(${PROJECT_NAME} PRIVATE ${${GEDIM_LIBRARY_TEST_INCLUDE}})

set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${${GEDIM_LIBRARY_TEST_PUBLIC_HEADERS}}")

gtest_discover_tests(${PROJECT_NAME})
