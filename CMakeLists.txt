cmake_minimum_required(VERSION 3.1)

enable_language(C)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)

if(APPLE)
	# This does not use GCC.
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++ -Wno-unknown-pragmas -Wall -D_GLIBCXX_PARALLEL")
else()
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -fopenmp -Wno-unknown-pragmas -Wall -D_GLIBCXX_PARALLEL")
endif()

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set (CMAKE_INSTALL_RPATH "/usr/local/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set (CMAKE_CXX_COMPILER_LAUNCHER ccache)

set (CMAKE_AUTOUIC ON)
set (CMAKE_AUTOMOC ON)
set (CMAKE_INCLUDE_CURRENT_DIR ON)
set (CMAKE_AUTOMOC_OPTIONS "-Iinclude -Isrc")

set (THREADS_PREFER_PTHREAD_FLAG ON)

find_package (GDAL 2.2 REQUIRED)
find_package (Qt5Widgets REQUIRED)
find_library (GEOS_LIBRARY geos PATHS /opt/lib /usr/local/lib /usr/lib NO_CMAKE_SYSTEM_PATH)

project (flood)
add_subdirectory (libgeo)

# Configure directories###################################################################################

include_directories(BEFORE SYSTEM ./include ./libgeo/include ./libgeo/ann/include /opt/include /usr/include/gdal /usr/include/eigen3) 
link_directories(. /opt/lib /usr/local/lib /usr/lib)

# Build libraries ####################################################################################################

#add_library(floodplugin SHARED src/flood_qgis.cpp ${MOCSRC} ${UISRC})
#target_link_libraries(floodplugin ${QT_LIBRARIES} qgis_app qgis_core)

# Build programs #####################################################################################################

add_executable (flood src/flood.cpp)
target_link_libraries (flood ${GDAL_LIBRARY} geoann georaster geoutil)

# Install ###########################################################################################################

install (TARGETS flood DESTINATION bin)

