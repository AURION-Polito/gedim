set(ENABLE_PREDEND_LICENSE OFF CACHE BOOL "Enable prepend_license")

if(${ENABLE_PREDEND_LICENSE})

	list(FILTER CLANGFORMAT_FILES INCLUDE REGEX ".cpp|.hpp")
		
	foreach(FILE ${CLANGFORMAT_FILES})
		list(APPEND COMMANDS COMMAND bash -c "echo \"Processing ${FILE}\"")
		list(APPEND COMMANDS COMMAND bash -c "${CMAKE_CURRENT_LIST_DIR}/prepend_license.sh ${FILE}")
	endforeach()

	add_custom_target(
		prepend_license
		ALL
		${COMMANDS}
	)
endif()
