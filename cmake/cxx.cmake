
# Find external packages
find_package( ROOT COMPONENTS Core Matrix MathCore Physics)
find_package( Eigen3 3.3 REQUIRED NO_MODULE )

#ROOT_GENERATE_DICTIONARY( G__EventShapesLib
#	EventShapes/EventShapes.h
#	LINKDEF Root/LinkDef.h
#)

add_library( EventShapesLib SHARED
	Root/EventShapes.cxx
	Root/Thrust.cxx
	EventShapes/EventShapes.h
#	G__EventShapesLib.cxx
)

target_include_directories( EventShapesLib PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}
	${ROOT_INCLUDE_DIRS}
)

target_link_libraries( EventShapesLib
	${ROOT_LIBRARIES}
	Eigen3::Eigen
)

ROOT_GENERATE_DICTIONARY(G__EventShapesLib
	EventShapes/EventShapes.h
	MODULE EventShapesLib
	LINKDEF Root/LinkDef.h
)

# Add executable for testing the C++
add_executable( estest ${CMAKE_CURRENT_SOURCE_DIR}/tests/main.cxx )
target_include_directories( estest PUBLIC "${PROJECT_BINARY_DIR}" )
target_link_libraries( estest PUBLIC EventShapesLib )

install(TARGETS EventShapesLib DESTINATION ${PROJECT_BINARY_DIR}/lib)
install(FILES EventShapes/EventShapes.h DESTINATION ${PROJECT_BINARY_DIR}/include)

# Write shell script for exporting environment variables
file(WRITE ${CMAKE_BINARY_DIR}/setup_env.sh
"#!$ENV{SHELL}
export LD_LIBRARY_PATH=\"${CMAKE_BINARY_DIR}/EventShapes:\${LD_LIBRARY_PATH}\"
")
