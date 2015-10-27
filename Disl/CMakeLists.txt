cmake_minimum_required(VERSION 2.6)

project(CADD_Disl)
enable_language(Fortran)
enable_language(CXX OPTIONAL)

# list of source files
FILE(GLOB SRC_FILES *.f)
FILE(GLOB INC_FILES *.par)

### Checking to see if compile flags are inherited

# Create Library for linking into main executable
add_library(dislo SHARED ${SRC_FILES} ${INC_FILES})

#add_library(disla STATIC ${SRC_FILES} ${INC_FILES})


# set the linker flags that are not libraries
if(BLAS_LINKER_FLAGS)
  set_target_properties(qsdd PROPERTIES LINK_FLAGS ${BLAS_LINKER_FLAGS})
endif()

#if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
#  message(FATAL_ERROR "Sorry qsdd does not compile with gfortran. Please choose an other compiler by setting the variable CMAKE_Fortran_COMPILER")
#endif()
