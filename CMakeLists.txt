cmake_minimum_required(VERSION 3.9)
project(EigBench)



set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake-modules)

################################################################
### Get git version number                                   ###
### Generates a header gitversion/gitversion.h               ###
### Include it using #include <gitversion.h>                 ###
### Gives a namespace GIT:: with several git version numbers.###
################################################################
include(cmake-modules/gitversion.cmake)


##################################################################
### Define variables for directory paths accessible at runtime ###
##################################################################
include(cmake-modules/directory.cmake)


##################################################################
### Add a debugging type message, only active in Debug-mode   ###
##################################################################
include(cmake-modules/debug_msg.cmake)




##################################################################
### Get operating system modes into a variable            ###
##################################################################
set(GET_OS_INFO_CMD lsb_release -a)
if(${CMAKE_HOST_APPLE})
    set(GET_OS_INFO_CMD "sw_vers")
endif()
execute_process(COMMAND ${GET_OS_INFO_CMD}
        OUTPUT_VARIABLE OS_PROPERTIES
        OUTPUT_STRIP_TRAILING_WHITESPACE)
message("========= DETECTED OS =========")
message("${OS_PROPERTIES}")
message("===============================")


###################
### Setup flags ###
###################
message("Compiling with: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
set(COMMON_OPTIONS -g -Wall -march=haswell)
set(DEBUG_OPTIONS       -Wpedantic -Wextra -O0)
set(RELEASE_OPTIONS     -DNDEBUG -O3)

##########################
### OS-dependent flags ###
##########################
if(${CMAKE_HOST_APPLE})
    list(APPEND COMMON_OPTIONS -mmacosx-version-min=10.13)
else()
    list(APPEND COMMON_OPTIONS -mavx ) #-mfma makes things slower?
endif()

################################
### Compiler-dependent flags ###
################################
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" )
    list(APPEND COMMON_OPTIONS -Wno-implicit-fallthrough -Wno-deprecated-declarations -Wno-ignored-attributes -Wno-int-in-bool-context -Wno-maybe-uninitialized -Wno-enum-compare)
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    list(APPEND COMMON_OPTIONS -stdlib=libstdc++ -Wno-deprecated-declarations -Wno-ignored-attributes -Wno-invalid-partial-specialization -Wno-missing-braces -Wno-overloaded-virtual -Wno-uninitialized)
endif()


################################
### Print current CXX flags  ###
################################
message(STATUS "COMMON CXX flags:  ${COMMON_OPTIONS} ${CMAKE_CXX_FLAGS}")
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "RELEASE CXX flags: ${RELEASE_OPTIONS}")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "DEBUG CXX flags: ${DEBUG_OPTIONS}")
endif()





###########################################
###  Add all source files recursively   ###
###########################################
file(GLOB_RECURSE SOURCES "source/*.cpp")
file(GLOB_RECURSE HEADERS "source/*.h")
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
target_include_directories(${PROJECT_NAME} PUBLIC source)

#Force cmake to find .a library suffixes
set(STATIC_BUILD ON)
set(CUSTOM_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
set(CMAKE_FIND_LIBRARY_SUFFIXES ${CUSTOM_SUFFIX} ${CMAKE_FIND_LIBRARY_SUFFIXES})
#target_link_libraries  (${PROJECT_NAME} PRIVATE -static)                                             ### Static linkage
#target_link_libraries  (${PROJECT_NAME} PRIVATE -static-libgcc -static-libstdc++)                                             ### Static linkage




###########################################
###  Apply RELEASE/DEBUG compile flags  ###
###########################################
target_compile_options (${PROJECT_NAME} PRIVATE ${COMMON_OPTIONS})                                   ### Common options
target_compile_options (${PROJECT_NAME} PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")               ### Debug build options
target_compile_options (${PROJECT_NAME} PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")           ### Release build options

###########################################
###  Apply RELEASE/DEBUG compile flags  ###
###########################################
target_compile_options (${PROJECT_NAME} PRIVATE ${COMMON_OPTIONS})                                   ### Common options
target_compile_options (${PROJECT_NAME} PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")               ### Debug build options
target_compile_options (${PROJECT_NAME} PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")           ### Release build options
target_link_libraries  (${PROJECT_NAME} PRIVATE -flto)                                               ### Link time optimization flag
target_link_libraries  (${PROJECT_NAME} PRIVATE -lstdc++fs)                                          ### For experimental::filesystem
target_link_libraries  (${PROJECT_NAME} PRIVATE -lpthread )



##############################################
###  Optional OpenMP support               ###
###  Note that Clang has some              ###
###  trouble with static openmp            ###
###  Also, libiomp5 might help for         ###
###  shared linking.                       ###
##############################################
#find_package(OpenMP)
#if (OpenMP_FOUND)
#    #    target_compile_options(${PROJECT_NAME} PRIVATE -fopenmp=libiomp5)
#    if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" )
#        message(STATUS "OpenMP ON" )
#        target_compile_options    (${PROJECT_NAME} PRIVATE ${OpenMP_CXX_FLAGS})
#        target_link_libraries     (${PROJECT_NAME} PRIVATE ${OpenMP_CXX_LIBRARIES})
#        target_compile_definitions(${PROJECT_NAME} PRIVATE -DOpenMP_AVAILABLE)
#    endif()
#    #    get_cmake_property(_variableNames VARIABLES)
#    #    foreach (_variableName ${_variableNames})
#    #        if("${_variableName}" MATCHES "openmp" OR "${_variableName}" MATCHES "OpenMP")
#    #            message(STATUS "${_variableName}=${${_variableName}}")
#    #        endif()
#    #    endforeach()
#else()
#    message(STATUS "OpenMP OFF" )
#endif()


###############################
###  Require c++17 support  ###
###############################
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
set_target_properties(${PROJECT_NAME}
        PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF)





###################
####  Link MPI  ###
###################
find_package(MPI REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE ${MPI_LIBRARIES})
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
target_include_directories(${PROJECT_NAME} PRIVATE ${MPI_INCLUDE_PATH})



###########################################
###  Find or fetch all dependencies   ###
###########################################
target_link_libraries(${PROJECT_NAME} PRIVATE -lpthread)

set(INSTALL_DIRECTORY ${PROJECT_SOURCE_DIR}/libs)
include(cmake-modules/FindGFortran.cmake)
#include(cmake-modules/Find_dont_install_INTELMKL.cmake)     # MKL - Intel's Math Kernel Library, use the BLAS implementation in Eigen and Arpack. Includes lapack.
include(cmake-modules/Fetch_OpenBLAS.cmake)                 # If MKL is not found openblas will be used instead. Includes lapack.
include(cmake-modules/FindLAPACKE.cmake)
include(cmake-modules/Fetch_Eigen3.cmake)                   # Eigen library
include(cmake-modules/Fetch_Flens.cmake)                    # Flens library
#include(cmake-modules/Fetch_Elemental.cmake)                # Elemental library
#include(cmake-modules/Fetch_HDF5.cmake)                     # HDF5 binary file format for output data_struct.
include(cmake-modules/Fetch_SuperLU.cmake)                  #
#include(cmake-modules/Fetch_SuiteSparse.cmake)              #
include(cmake-modules/Fetch_Arpack-ng.cmake)                # Iterative Eigenvalue solver for a few eigenvalues/eigenvectors using Arnoldi method.
include(cmake-modules/Fetch_Arpackpp.cmake)                 # C++ frontend for arpack-ng
include(cmake-modules/Fetch_GSL.cmake)                      # GSL for computing some numerical integrals.
include(cmake-modules/Fetch_Armadillo.cmake)                # Armadillo library
#include(cmake-modules/Fetch_PETSc.cmake)                    # PETSc library




######################################
###  Require experimental headers  ###
######################################
# Use cmake_policy CMP0067 to let try_compile honor current CMAKE_CXX_STANDARD setting.
# See here
# https://stackoverflow.com/questions/47213356/cmake-using-corrext-c-standard-when-checking-for-header-files
cmake_policy(SET CMP0067 NEW)
include(CheckIncludeFileCXX)
set(CMAKE_REQUIRED_FLAGS -std=c++17  -lstdc++fs -l${PATH_TO_LIBSTDCXXFS})
check_include_file_cxx(experimental/filesystem    has_experimental_filesystem)
check_include_file_cxx(experimental/type_traits   has_experimental_type_traits)
if(NOT has_experimental_filesystem OR NOT has_experimental_type_traits )
    message(FATAL_ERROR "\n\
        Missing one or more experimental headers.\n\
        Consider using a newer compiler (GCC 7 or above, Clang 6 or above),\n\
        or checking the compiler flags.")
endif()






message("=============================================================================")
message("=== SUMMARY OF CMAKE CONFIGURATION:  THE FOLLOWING LIBRARIES WILL BE USED ===")
message("=============================================================================")
message("    INTEL MKL:   ${MKL_LIBRARIES}")
message("    BLAS     :   ${BLAS_LIBRARIES}")
message("    LAPACK   :   ${LAPACK_LIBRARIES}")
message("    LAPACKE  :   ${LAPACKE_LIBRARIES}")
message("    EIGEN3   :   ${EIGEN3_INCLUDE_DIR}")
message("    ARMADILLO:   ${ARMADILLO_INCLUDE_DIRS}")
message("    BLAZE    :   ${BLAZE_INCLUDE_DIRS}")
message("    ELEMENTAL:   ${ELEMENTAL_LIBRARIES}")
message("    HDF5     :   ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_CXX_HL_LIBRARIES}")
message("    ARPACK   :   ${ARPACK_LIBRARIES}")
message("    ARPACKPP :   ${ARPACKPP_INCLUDE_DIR}")
message("    GSL      :   ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES}")
message("    GFORTRAN :   ${GFORTRAN_LIB}")
message("    QUADMATH :   ${QUADMATH_LIB}")
message("=============================================================================")
message("")




###################
### INFORMATION ###
###################
message("")
message("=============================================================================")
message("")
message("      To clean cmake files, run './build.sh -c'")
message("      To clean cmake files downloaded libraries, run './build.sh -l '")
message("      To launch this program,  './run.sh'")
message("")
message("=============================================================================")
message("")

