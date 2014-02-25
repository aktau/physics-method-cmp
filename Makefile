SOURCE_DATA=data.tmp

determine_lua = \
		lua=; \
		for x in luajit lua lua52; do \
		  if type "$${x%% *}" >/dev/null 2>/dev/null; then lua=$$x; break; fi; \
		done; \
		if [ -z "$$lua" ]; then echo 1>&2 "Unable to find a lua binary"; exit 2; fi

default: comparison.html

$(SOURCE_DATA): integrators.lua
	@$(determine_lua); \
	$$lua integrators.lua > "$@"

comparison.html: $(SOURCE_DATA)
	./plot.sh "$<"
