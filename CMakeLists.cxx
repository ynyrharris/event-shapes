
# Find external packages
find_package( ROOT COMPONENTS Core Matrix MathCore Physics)

ROOT_GENERATE_DICTIONARY( G__EventShapesLib
	EventShapes/EventShapes.h
	LINKDEF Root/LinkDef.h
)

add_library( EventShapesLib SHARED
	Root/EventShapes.cxx EventShapes/EventShapes.h
	G__EventShapesLib.cxx
)

target_include_directories( EventShapesLib PUBLIC
	${CMAKE_SOURCE_DIR}/EventShapes
	${ROOT_INCLUDE_DIRS}
)

target_link_libraries( EventShapesLib
	${ROOT_LIBRARIES}
)

install(TARGETS EventShapesLib DESTINATION ${PROJECT_BINARY_DIR}/lib)
install(FILES EventShapes/EventShapes.h DESTINATION ${PROJECT_BINARY_DIR}/include)

# Write shell script for exporting environment variables
file(WRITE ${CMAKE_BINARY_DIR}/setup_env.sh
"#!$ENV{SHELL}
export LD_LIBRARY_PATH=\"${CMAKE_BINARY_DIR}/EventShapes:\${LD_LIBRARY_PATH}\"
")
