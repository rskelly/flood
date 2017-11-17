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

project (flood)
add_subdirectory (../libgeo libgeo)

find_package (GDAL 2.1 REQUIRED)
find_package (OpenMP REQUIRED) # Implies pthread?

set (GEOS_LIBRARY /usr/lib/x86_64-linux-gnu/libgeos.a)

if(OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
	set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
endif()

# Configure directories
link_directories (. ${Boost_LIBRARY_DIRS} ${GDAL_LIBRARY})
include_directories (BEFORE SYSTEM ../libgeo/include)
include_directories (BEFORE SYSTEM ${GDAL_INCLUDE_DIR})
include_directories (BEFORE SYSTEM /usr/include/eigen3)
include_directories (BEFORE SYSTEM /opt/include)
link_directories (BEFORE SYSTEM /opt/lib)


# Build programs #####################################################################################################

add_executable (flood src/flood.cpp)
target_link_libraries (flood ${GDAL_LIBRARY} ${GEOS_LIBRARY} raster util)

# Install ###########################################################################################################

install (TARGETS flood DESTINATION bin)
