#################################################################
# HEADER
#################################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

PROJECT(original_f90_code Fortran)
MESSAGE(STATUS "PROJECT_NAME: ${PROJECT_NAME}")

SET(${PROJECT_NAME}_VERSION 0.0.1)
SET(${PROJECT_NAME}_SOVERSION 1)

#################################################################
# PROJECT PATHS
#################################################################

SET(ROOT_PATH ${CMAKE_SOURCE_DIR})
#SET(SRC_PATH ${ROOT_PATH}/Sources)
#SET(LIB_PATH ${SRC_PATH}/Lib)
#SET(DRIVERS_PATH ${SRC_PATH}/Drivers)
SET(CMAKE_MODULE_PATH "${ROOT_PATH}/CMake/Modules") 

#################################################################
# BUILD PATHS
#################################################################

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
INCLUDE_DIRECTORIES(${CMAKE_Fortran_MODULE_DIRECTORY})

#################################################################
# CONFIGURATION TYPES & BUILD MODE & TESTING
#################################################################
SET(CMAKE_CONFIGURATION_TYPES DEBUG RELEASE) # Specifies the available build types
IF(NOT CMAKE_BUILD_TYPE)                     # If build type was not specified, set DEBUG as default
  SET(CMAKE_BUILD_TYPE DEBUG CACHE STRING
      "Choose the type of build, options are: NONE DEBUG RELEASE"
      FORCE)

  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS NONE DEBUG RELEASE)
ENDIF(NOT CMAKE_BUILD_TYPE)

#################################################################
# FFLAGS depending on the compiler and the build type
#################################################################

GET_FILENAME_COMPONENT(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
MESSAGE(STATUS "COMPILER INFO: ${CMAKE_Fortran_COMPILER_ID} - ${Fortran_COMPILER_NAME}")

IF (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" OR Fortran_COMPILER_NAME MATCHES "gfortran*")
  # gfortran 
  SET(FORTRAN_FLAGS "-fdefault-real-8 -ffree-line-length-0 -cpp -Wimplicit-interface ${EXTRA_FLAGS} ")
  SET(CMAKE_Fortran_FLAGS "${FORTRAN_FLAGS}" CACHE STRING "" FORCE)
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3" CACHE STRING "" FORCE)
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-g -fbacktrace -fbounds-check -fprofile-arcs -ftest-coverage -Wimplicit-interface" CACHE STRING "" FORCE)
ELSEIF (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" OR Fortran_COMPILER_NAME MATCHES "ifort*")
  # ifort
  SET(FORTRAN_FLAGS "-r8 -fpp -W1 ${EXTRA_FLAGS} ")
  SET(CMAKE_Fortran_FLAGS "${FORTRAN_FLAGS} " CACHE STRING "" FORCE)
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3" CACHE STRING "" FORCE)
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -traceback -g -debug all -check all -ftrapuv -warn nointerfaces" CACHE STRING "" FORCE)
ELSEIF (${CMAKE_Fortran_COMPILER_ID} STREQUAL "XL" OR Fortran_COMPILER_NAME MATCHES "xlf*")
  # xlf
  SET(FORTRAN_FLAGS "-q64 -qrealsize=8 -qsuffix=f=f90:cpp=f90 ${EXTRA_FLAGS}")
  SET(CMAKE_Fortran_FLAGS "${FORTRAN_FLAGS}" CACHE STRING "" FORCE)
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -qstrict" CACHE STRING "" FORCE)
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -qfullpath -qkeepparm" CACHE STRING "" FORCE)
ELSE ()
  MESSAGE("Unknown compiler. No compiler-dependent flags could be used, just standard -O2 and -O0 -g for release/debug, resp.")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
ENDIF ()

message (STATUS "CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
message (STATUS "CMAKE_Fortran_FLAGS: " ${CMAKE_Fortran_FLAGS})
message (STATUS "CMAKE_Fortran_FLAGS_RELEASE: " ${CMAKE_Fortran_FLAGS_RELEASE})
message (STATUS "CMAKE_Fortran_FLAGS_DEBUG: " ${CMAKE_Fortran_FLAGS_DEBUG})

#################################################################
# ADD SUBDIRECTORIES AND DRIVER
#################################################################

#ADD_SUBDIRECTORY(${LIB_PATH})
#ADD_SUBDIRECTORY(${DRIVERS_PATH})

FILE(GLOB_RECURSE DRIVER_SRC *.f90 *.F90)
#MESSAGE(STATUS ${DRIVER_SRC})
SET(EXE_NAME example1)
ADD_EXECUTABLE(${EXE_NAME} ${DRIVER_SRC})
