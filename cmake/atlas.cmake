
# The name of the package
atlas_subdir( EventShapes )

# Find external packages
find_package( ROOT COMPONENTS Core Matrix MathCore Physics )

# Generate a CINT dictionary source file
atlas_add_root_dictionary( EventShapesLib EventShapesDictSource
	ROOT_HEADERS EventShapes/*.h Root/LinkDef.h
)

atlas_add_library( EventShapesLib
	EventShapes/*.h Root/*.cxx ${EventShapesDictSource}
	PUBLIC_HEADERS EventShapes
	INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
	LINK_LIBRARIES ${ROOT_LIBRARIES}
)