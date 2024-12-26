set(ENABLE_CLANGFORMAT OFF CACHE BOOL "Enable clangformat")

set(CLANGFORMAT_STYLE "{BasedOnStyle: Microsoft, AllowAllArgumentsOnNextLine: 'false', AllowAllParametersOfDeclarationOnNextLine: 'false', AlwaysBreakAfterDefinitionReturnType: None, AlwaysBreakAfterReturnType: None, BinPackArguments: 'false', BinPackParameters: 'false', ColumnLimit: '120', PenaltyBreakAssignment: '5', PenaltyBreakBeforeFirstCallParameter: '20', PenaltyBreakComment: '0', PenaltyBreakFirstLessLess: '1', PenaltyBreakString: '5', PenaltyExcessCharacter: '1', PenaltyReturnTypeOnItsOwnLine: '150'}")

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
