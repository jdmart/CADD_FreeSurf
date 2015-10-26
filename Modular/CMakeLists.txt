cmake_minimum_required(VERSION 2.6)

project(CADD)
enable_language(Fortran)
enable_language(CXX OPTIONAL)

# list of source files
FILE(GLOB SRC_FILES *.f)
FILE(GLOB INC_FILES *.par)
INCLUDE_DIRECTORIES(../Disl 
	${CMAKE_BINARY_DIR}/Disl)

# Create Library for linking into main executable
add_library(mdlro SHARED ${SRC_FILES} ${INC_FILES})

#add_library(mdlra STATIC ${SRC_FILES} ${INC_FILES})


# set the linker flags that are not libraries

