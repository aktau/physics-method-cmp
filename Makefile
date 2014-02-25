SOURCE_DATA=data.tmp
SOURCE_META=meta.plt

determine_lua = \
		lua=; \
		for x in luajit lua lua52; do \
		  if type "$${x%% *}" >/dev/null 2>/dev/null; then lua=$$x; break; fi; \
		done; \
		if [ -z "$$lua" ]; then echo 1>&2 "Unable to find a lua binary"; exit 2; fi

default: comparison.html

$(SOURCE_DATA): integrators.lua
	@echo "rebuilding data..."
	@$(determine_lua); \
	$$lua integrators.lua data 2>&1 > "$@"

$(SOURCE_META): integrators.lua
	@echo "rebuilding meta..."
	@$(determine_lua); \
	$$lua integrators.lua meta 2>&1 > "$@"

comparison.html: $(SOURCE_DATA) $(SOURCE_META)
	./plot.sh $^

open: comparison.html
	open "$<"

clean:
	rm -f $(SOURCE_DATA) || true
	rm -f $(SOURCE_META) || true
	rm -f comparison.html || true

.PHONY: clean open
