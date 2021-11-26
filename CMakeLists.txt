cmake_minimum_required( VERSION 3.17 )

project(EVENTSHAPES VERSION 2.0)

# Set up build options
option(PYTHON_SETUP  "Build package using Cython" YES)
option(CXX_SETUP     "Build package as a standalone C++ library" NO)

# Find external packages required in every build scheme
find_package( ROOT COMPONENTS Core Matrix MathCore Physics)
find_package( Eigen3 3.3 REQUIRED NO_MODULE )

# Select the build scheme based on user's configuration
if (PYTHON_SETUP AND (NOT CXX_SETUP))
	include(cmake/python.cmake)
elseif (CXX_SETUP)
	include(cmake/cxx.cmake)
endif()
