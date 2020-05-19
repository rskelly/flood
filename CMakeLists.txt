cmake_minimum_required(VERSION 3.10)

project (flood)
add_subdirectory (libgeo)

enable_language(C)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_CXX_FLAGS "-Wall -D_GLIBCXX_PARALLEL")
set (CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set (CMAKE_CXX_FLAGS_RELEASE "-O3")
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set (CMAKE_INSTALL_RPATH "/usr/local/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set (CMAKE_CXX_COMPILER_LAUNCHER ccache)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wno-unknown-pragmas -fopenmp")
set (FS_LIB stdc++fs) # For filesystem features in c++17.

find_package (Eigen3 REQUIRED)
find_package (GDAL 2 REQUIRED)
find_package (Qt5Widgets REQUIRED)


# Configure directories###################################################################################

include_directories(BEFORE SYSTEM 
	include
	libgeo/include
	ann/include 
	${EIGEN3_INCLUDE_DIR} 
	${GDAL_INCLUDE_DIR}
)
	 
link_directories(./build/lib ./build/bin 
	${GDAL_LIBRARYDIR} 
)

# Build libraries ####################################################################################################

#add_library(floodplugin SHARED src/flood_qgis.cpp ${MOCSRC} ${UISRC})
#target_link_libraries(floodplugin ${QT_LIBRARIES} qgis_app qgis_core)

# Build programs #####################################################################################################

add_executable (flood src/flood_app.cpp src/flood.cpp)
target_link_libraries (flood ${GDAL_LIBRARY} geogrid geoutil)

# Install ###########################################################################################################

install (TARGETS flood DESTINATION bin)

