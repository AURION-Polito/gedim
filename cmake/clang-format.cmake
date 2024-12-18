set(ENABLE_CLANGFORMAT OFF CACHE BOOL "Enable clangformat")

set(CLANGFORMAT_STYLE "{BasedOnStyle: Microsoft, AlignAfterOpenBracket: Align, BinPackArguments: 'false', BinPackParameters: 'false'}")

if(${ENABLE_CLANGFORMAT})

	list(FILTER CLANGFORMAT_FILES INCLUDE REGEX ".cpp|.hpp")

	add_custom_target(
		clangformat
		ALL
		COMMAND clang-format
		-style=${CLANGFORMAT_STYLE}
		-i
		${CLANGFORMAT_FILES}
	)
endif()
