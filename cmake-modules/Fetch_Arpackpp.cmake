



message(STATUS "SEARCHING FOR ARPACK++ IN SYSTEM...")
find_library(ARPACKPP_LIBRARIES
        NAMES arpackpp arpack++ libarpack2++ libarpack++ libarpackpp
        PATH_SUFFIXES lib lib32 lib64
        )
find_path(ARPACKPP_INCLUDE_DIR
        NAMES ardsnsym.h
        PATHS /usr/include/arpack++ /usr/include /usr/local/
        )

if (NOT ARPACKPP_LIBRARIES OR NOT ARPACKPP_INCLUDE_DIR)
    # Try finding arpack as module library
    set(ARPACKPP_LIBRARIES "")
    message(STATUS "SEARCHING FOR ARPACK IN LOADED MODULES")
    find_path(ARPACKPP_INCLUDE_DIR
            NAMES ardsnsym.h arpack++
            PATHS $ENV{ARPACKPP_DIR}/include
            NO_DEFAULT_PATH
            )
endif()


message(STATUS "Note that old versions of Arpack++ (e.g. the default in Ubuntu Trusty 14.04 LTS) may fail to compile, requiring '-fpermissive'.")
if (ARPACKPP_LIBRARIES OR ARPACKPP_INCLUDE_DIR AND NOT "${OS_PROPERTIES}" MATCHES "trusty" )
    message(STATUS "Arpack++ library found in system: ${ARPACKPP_LIBRARIES}")
    message(STATUS "Arpack++ include found in system: ${ARPACKPP_INCLUDE_DIR}")
    add_library(arpack++ INTERFACE)
    set_target_properties(arpack++ PROPERTIES
            INTERFACE_LINK_LIBRARIES "blas;lapack;arpack;${ARPACKPP_LIBRARIES}"
            INTERFACE_INCLUDE_DIRECTORIES "${ARPACKPP_INCLUDE_DIR}")
    target_link_libraries(${PROJECT_NAME} PRIVATE arpack++)

else()
    message(STATUS "Arpack++ will be installed into ${INSTALL_DIRECTORY}/arpackpp on first build.")
    #####################################################################
    ### Prepare lists with generator expressions, replacing all semicolons.
    ### Otherwise, passing raw lists results  in only the first element
    ### of the list to be passed.
    ####################################################################
    string (REPLACE ";" "$<SEMICOLON>" BLAS_LIBRARIES_GENERATOR     "${BLAS_LIBRARIES}")
    string (REPLACE ";" "$<SEMICOLON>" LAPACK_LIBRARIES_GENERATOR   "${LAPACK_LIBRARIES}")
    string (REPLACE ";" "$<SEMICOLON>" FC_LDLAGS_GENERATOR          "${FC_LDLAGS}")
    ####################################################################
    include(ExternalProject)
    ExternalProject_Add(library_ARPACK++
            GIT_REPOSITORY      https://github.com/m-reuter/arpackpp.git
            GIT_TAG             master
            PREFIX              "${INSTALL_DIRECTORY}/arpack++"

            UPDATE_COMMAND ""


            CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
            -DEXAMPLES=ON
            -DCMAKE_BUILD_TYPE=Release
            -DBLAS_LIB=${BLAS_LIBRARIES_GENERATOR}
            -DOPENBLAS_LIB=${BLAS_LIBRARIES_GENERATOR}
            -DARPACK_LIB=${ARPACK_LIBRARIES}
            -DSUPERLU=ON
            -DSUPERLU_LIB=${SUPERLU_LIBRARIES}
            -DSUPERLU_INC=${SUPERLU_INCLUDE_DIRS}
            -DUMFPACK=OFF
            DEPENDS blas lapack arpack superlu gfortran
#            -DCMAKE_FIND_ROOT_PATH=<INSTALL_DIR>;${BLAS_LIBRARIES_GENERATOR};${ARPACK_LIBRARIES};${SUPERLU_LIBRARIES}
#            -DCMAKE_C_FLAGS=-w -m64 -fPIC
#            -DCMAKE_Fortran_FLAGS=-w -m64 -fPIC
#            -DEXTRA_LDLAGS=${FC_LDLAGS_GENERATOR}

            #            -DMPI=OFF
#            -DINTERFACE64=OFF
#            -DBUILD_SHARED_LIBS=${ARPACK_SHARED}

            #            UPDATE_COMMAND ""
#            TEST_COMMAND ""
#            INSTALL_COMMAND ""
#            CONFIGURE_COMMAND ""
#            BUILD_COMMAND
#            ${CMAKE_COMMAND} -E make_directory <INSTALL_DIR>/include && find <INSTALL_DIR>/include -maxdepth 1 -type l -delete &&
#            ${CMAKE_COMMAND} -E create_symlink <SOURCE_DIR>/include <INSTALL_DIR>/include/arpack++
            )

    ExternalProject_Get_Property(library_ARPACK++ INSTALL_DIR)
    add_library(arpack++ INTERFACE)
    set(ARPACKPP_INCLUDE_DIR ${INSTALL_DIR}/include)
    set_target_properties(arpack++ PROPERTIES
            INTERFACE_LINK_LIBRARIES "arpack;blas;lapack;superlu;gfortran"
            INTERFACE_INCLUDE_DIRECTORIES ${ARPACKPP_INCLUDE_DIR}
            )
    add_dependencies(arpack++ library_ARPACK++ blas lapack)
    target_link_libraries(${PROJECT_NAME} PRIVATE arpack++)
    target_include_directories(${PROJECT_NAME} PRIVATE ${ARPACKPP_INCLUDE_DIR})
endif()



