file(GLOB_RECURSE ALL_SOURCE_FILES *.cpp *.h)

set(ENABLE_CLANGFORMAT OFF CACHE BOOL "Enable clangformat")

if(${ENABLE_CLANGFORMAT})
	add_custom_target(
		clangformat
		ALL
		COMMAND clang-format
		-style=GNU
		-i
		${CLANGFORMAT_FILES}
	)
endif()
