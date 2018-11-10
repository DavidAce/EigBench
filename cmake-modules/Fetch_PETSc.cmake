  message(STATUS "PETSc will be installed into ${INSTALL_DIRECTORY}/petsc on first build.")

  #####################################################################
  ### Prepare lists with generator expressions, replacing all semicolons.
  ### Otherwise, passing raw lists results  in only the first element
  ### of the list to be passed.
  ####################################################################
  string (REPLACE ";" "$<SEMICOLON>" BLAS_LIBRARIES_GENERATOR     "${BLAS_LIBRARIES}")
  string (REPLACE ";" "$<SEMICOLON>" LAPACK_LIBRARIES_GENERATOR   "${LAPACK_LIBRARIES}")
  string (REPLACE ";" "$<SEMICOLON>" FC_LDLAGS_GENERATOR          "${FC_LDLAGS}")
  ####################################################################
  if(${STATIC_BUILD})
      set(PETSC_SHARED OFF)
  else()
      set(PETSC_SHARED ON)
  endif()

  include(ExternalProject)
  ExternalProject_Add(library_PETSC
          URL      http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.10.2.tar.gz
          PREFIX              "${INSTALL_DIRECTORY}/petsc"
          UPDATE_COMMAND ""
          CONFIGURE_COMMAND
             cd <SOURCE_DIR> &&
             ./configure
                    --with-mpi=false
                    --prefix=<INSTALL_DIR>
          BUILD_COMMAND
            cd <SOURCE_DIR> &&
            ${CMAKE_MAKE_PROGRAM}
          INSTALL_COMMAND
            cd <SOURCE_DIR> &&
            ${CMAKE_MAKE_PROGRAM}install

#
          DEPENDS blas lapack gfortran
          )
  ExternalProject_Get_Property(library_PETSC INSTALL_DIR)
  set(PETSC_INCLUDE_DIRS ${INSTALL_DIR}/include)
  add_library(petsc STATIC IMPORTED)
#  set_target_properties(arpack
#          PROPERTIES
#          IMPORTED_LOCATION "${INSTALL_DIR}/lib/libarpack${CUSTOM_SUFFIX}"
#          INTERFACE_LINK_LIBRARIES "blas;lapack;gfortran"
#          INTERFACE_LINK_FLAGS            "-lpthread"
#          )
#  #            INTERFACE_LINK_FLAGS      )

  add_dependencies(arpack library_PETSC blas lapack gfortran )
  target_link_libraries(${PROJECT_NAME} PRIVATE petsc)
  target_include_directories(${PROJECT_NAME} PRIVATE ${PETSC_INCLUDE_DIRS})
  #For convenience, define these variables
  get_target_property(PETSC_LIBRARIES petsc IMPORTED_LOCATION)