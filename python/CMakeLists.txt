
message("In python/CMakeLists.txt")

# Find external packages
find_package( ROOT COMPONENTS Core Matrix MathCore Physics)

# Find Python
find_package(Python3 COMPONENTS Interpreter Development)

message("Found Python? ${Python3_FOUND}")
message("Python version: ${Python3_VERSION}")

# Set variables
set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
set(SETUP_PY "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
set(PY_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/build")

# Substitute ${VAR}/@VAR@ into the input file
configure_file(${SETUP_PY_IN} ${SETUP_PY})

# Add build commands
add_custom_command(OUTPUT ${PY_OUTPUT}
	COMMAND ${Python3_EXECUTABLE} ${SETUP_PY} build_ext --build-lib ${PY_OUTPUT}
	COMMAND ${CMAKE_COMMAND} -E touch ${PY_OUTPUT}
)

add_custom_target(event_shapes ALL DEPENDS ${PY_OUTPUT})

install(CODE "execute_process(
		COMMAND ${Python3_EXECUTABLE} -m pip install --use-feature=in-tree-build ${CMAKE_CURRENT_BINARY_DIR}
	)"
)