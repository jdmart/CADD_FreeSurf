cmake_minimum_required(VERSION 2.6)

project(CADD)
enable_language(Fortran)
enable_language(CXX OPTIONAL)

# list of source files
INCLUDE_DIRECTORIES(Disl 
  Modular 
  ${CMAKE_BINARY_DIR}/Disl 
  ${CMAKE_BINARY_DIR}/Modular )

FILE(GLOB SRC_FILES *.f)
FILE(GLOB INC_FILES *.par)

SUBDIRS(Disl 
  Modular)

#===============================================================================
#Macro to add some compilation flags
macro(add_ff_flags flag)

  set(type ${ARGN})
  if(type)
    set(_var CMAKE_Fortran_FLAGS_${type})
  else()
    set(_var CMAKE_Fortran_FLAGS)
  endif()
  if(NOT ${_var} MATCHES "${flag} ")
    set(${_var} "${flag} ${${_var}}"
      CACHE STRING "Flags used by the compiler during all build types." FORCE)
  endif()
endmacro()

macro(remove_ff_flags flag)
  set(type ${ARGN})
  if(type)
    set(_var CMAKE_Fortran_FLAGS_${type})
  else()
    set(_var CMAKE_Fortran_FLAGS)
  endif()
  string(REPLACE "${flag} " "" ${_var} "${${_var}}")
  set(${_var} "${${_var}}"
    CACHE STRING "Flags used by the compiler during all build types." FORCE)
endmacro()
#===============================================================================

#  message("adding flags")
include_directories(${CMAKE_BINARY_DIR})

# to add ff flags for a given build type
#add_ff_flags(-fastsse RELEASE)
#add_ff_flags(-fast RELEASE)
#12 add_ff_flags(-g-traceback)
#12 add_ff_flags(-checkbounds -g)
#1 add_ff_flags(-mcmodel=medium)
#1 add_ff_flags(-shared-intel)
#1 add_ff_flags(-traceback)
#12 add_ff_flags(-heap-arrays)

#1 add_ff_flags(-fPIC)
#1 add_ff_flags(-w)

# you can also define a profile that you can select with CMAKE_BUILD_TYPE
set(CMAKE_Fortran_FLAGS_SUNF90 "-fast -xmodel=medium -xtarget=barcelona -m64" CACHE STRING "Flags used by the compiler during sunf90 builds")
set(CMAKE_EXE_LINKER_FLAGS_SUNF90 "" CACHE STRING "Flags used by the linker during sunf90 builds")
mark_as_advanced(CMAKE_Fortran_FLAGS_SUNF90)
mark_as_advanced(CMAKE_EXE_LINKER_FLAGS_SUNF90)

# Add executable 
add_executable(CADD ${SRC_FILES} ${INC_FILES})
target_link_libraries(CADD dislo mdlro)

#if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
#  message(FATAL_ERROR "Sorry CADD does not compile with gfortran. Please choose an other compiler by setting the variable CMAKE_Fortran_COMPILER")
#endif()
