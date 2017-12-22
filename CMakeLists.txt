cmake_minimum_required(VERSION 3.1)

enable_language(C)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -D_GLIBCXX_PARALLEL")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")

#set (Boost_DEBUG OFF)
#set (Boost_USE_STATIC_LIBS OFF)
#set (Boost_USE_STATIC_RUNTIME OFF)

set (CMAKE_AUTOUIC ON)
set (CMAKE_AUTOMOC ON)
set (CMAKE_INCLUDE_CURRENT_DIR ON)
set (CMAKE_AUTOMOC_OPTIONS "-Iinclude -Isrc")

project (flood)
add_subdirectory (libgeo)

find_package (GDAL 2.2 REQUIRED)
find_package (OpenMP REQUIRED) # Implies pthread?
find_package (Qt4 REQUIRED QtCore QtGui QtXml)

include (${QT_USE_FILE})

qt4_wrap_ui(UISRC src/flood_qgis.ui)
qt4_wrap_cpp(MOCSRC include/flood_qgis.hpp)

#set (GEOS_LIBRARY /usr/lib/x86_64-linux-gnu/libgeos.a)

if(OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
	set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
endif()

# Configure directories
include_directories (BEFORE SYSTEM ${GDAL_INCLUDE_DIR})
include_directories (BEFORE SYSTEM include)
include_directories (BEFORE SYSTEM /usr/include/qgis)
include_directories (BEFORE SYSTEM libgeo/include)
include_directories (BEFORE SYSTEM libgeo/ann/include)
include_directories (BEFORE SYSTEM /usr/include/eigen3)
#link_directories (. ${Boost_LIBRARY_DIRS} ${GDAL_LIBRARY})

# Build libraries ####################################################################################################

add_library(floodplugin SHARED src/flood_qgis.cpp ${MOCSRC} ${UISRC})
target_link_libraries(floodplugin ${QT_LIBRARIES} qgis_app qgis_core)

# Build programs #####################################################################################################

add_executable (flood src/flood.cpp)
target_link_libraries (flood ${GDAL_LIBRARY} geoann georaster geoutil)

# Install ###########################################################################################################

install (TARGETS flood DESTINATION bin)
