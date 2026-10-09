.PHONY: create distclean clean run

main_folder = not-set
external_folder = not-set
mesh_folder = not-set
initial_conditions_folder = not-set

ifeq ($(main_folder), not-set)
	MAKEFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
	MAKEFILE_DIR  := $(dir $(MAKEFILE_PATH))
	main_folder = ${MAKEFILE_DIR}..
endif

clean:
	@docker container prune
	@docker rmi $$(docker images --filter "dangling=true" -q --no-trunc)

create:
	@docker build --progress=plain -f Dockerfile --target build_windows . -t build_windows:1.0.0

run:
	@echo "Shared folder: ${main_folder}"
	@docker run --rm -it -v ${main_folder}:/main_folder -p 8080:8080 build_windows:1.0.0 bash

distclean:
	@docker container prune
	@docker rmi build_windows:1.0.0
	@docker rmi $$(docker images --filter "dangling=true" -q --no-trunc)
	@docker builder prune
	@docker buildx prune -f
