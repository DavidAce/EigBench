  message(STATUS "SuiteSparse will be installed into ${INSTALL_DIRECTORY}/suitesparse on first build.")

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
      set(SUITESPARSE_SHARED OFF)
  else()
      set(SUITESPARSE_SHARED ON)
  endif()
  set(OLD_PREFIX "/usr/local/lib/")
  include(ExternalProject)
  ExternalProject_Add(library_SUITESPARSE
          URL       http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-5.3.0.tar.gz
#          URL       http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.4.5.tar.gz
          PREFIX              "${INSTALL_DIRECTORY}/suitesparse"
          UPDATE_COMMAND ""

          CONFIGURE_COMMAND
          sed -i "/BLAS = -lopenblas/c\\BLAS = ${BLAS_LIBRARIES_GENERATOR} -lpthread" <SOURCE_DIR>/SuiteSparse_config/SuiteSparse_config.mk
#          sed -i "/INSTALL_LIB = /usr/local/lib/c\\INSTALL_LIB = ${INSTALL_DIRECTORY}/suitesparse/lib" <SOURCE_DIR>/SuiteSparse_config/SuiteSparse_config.mk &&
#          sed -i "/INSTALL_INCLUDE = /usr/local/lib/c\\INSTALL_INCLUDE = ${INSTALL_DIRECTORY}/suitesparse/include" <SOURCE_DIR>/SuiteSparse_config/SuiteSparse_config.mk

          BUILD_IN_SOURCE 1
          BUILD_COMMAND $(MAKE) BLAS=${BLAS_LIBRARIES_GENERATOR} LAPACK=${LAPACK_LIBRARIES_GENERATOR} library
          INSTALL_COMMAND
          $(MAKE) install INSTALL=${INSTALL_DIRECTORY}/suitesparse
          DEPENDS blas lapack gfortran

          )
  ExternalProject_Get_Property(library_SUITESPARSE INSTALL_DIR)
  set(SUITESPARSE_INCLUDE_DIRS ${INSTALL_DIR}/include)

  add_library(suitesparse UNKNOWN IMPORTED)
  set_target_properties(suitesparse
          PROPERTIES
          IMPORTED_LOCATION "${INSTALL_DIR}/lib/libsuitesparse${CUSTOM_SUFFIX}"
          INTERFACE_LINK_LIBRARIES "blas;lapack;gfortran"
          INTERFACE_LINK_FLAGS            "-lpthread"
          )

  add_dependencies(suitesparse library_SUITESPARSE blas lapack gfortran )
  target_link_libraries(${PROJECT_NAME} PRIVATE suitesparse -lpthread -m64)
  target_include_directories(${PROJECT_NAME} PRIVATE ${SUITESPARSE_INCLUDE_DIRS})
  #For convenience, define these variables
  get_target_property(SUITESPARSE_LIBRARIES suitesparse IMPORTED_LOCATION)

