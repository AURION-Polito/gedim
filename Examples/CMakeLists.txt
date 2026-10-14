cmake_minimum_required(VERSION 3.12)
project(GeDiM_Examples LANGUAGES CXX C VERSION 1.0)

add_subdirectory(Excavate_Domain_3D)

set(CLANGFORMAT_FILES ${CLANGFORMAT_FILES} PARENT_SCOPE)
