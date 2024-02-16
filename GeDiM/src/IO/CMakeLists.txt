list(APPEND gedim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/IOUtilities.cpp)
list(APPEND gedim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/FileTextReader.cpp)
list(APPEND gedim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/StringsUtilities.cpp)
list(APPEND gedim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/Configurations.cpp)
list(APPEND gedim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/VTKUtilities.cpp)
list(APPEND gedim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/UCDUtilities.cpp)

list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/IOUtilities.hpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/FileTextReader.hpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/IFileTextReader.hpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/StringsUtilities.hpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/Configurations.hpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/ConfigurationPropertySupportedTypes.hpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/VTKUtilities.hpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/UCDUtilities.hpp)

list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/TimeUtilities.hpp)

list(APPEND gedim_src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(gedim_src_sources ${gedim_src_sources} PARENT_SCOPE)
set(gedim_src_headers ${gedim_src_headers} PARENT_SCOPE)
set(gedim_src_includes ${gedim_src_includes} PARENT_SCOPE)
