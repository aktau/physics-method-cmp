SOURCE_DATA=data.tmp

default: comparison.html

$(SOURCE_DATA): integrators.lua
	./integrators.lua > "$@"

comparison.html: $(SOURCE_DATA)
	./int "$<"

