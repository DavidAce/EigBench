  message(STATUS "SuperLU will be installed into ${INSTALL_DIRECTORY}/superlu on first build.")

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
      set(SUPERLU_SHARED OFF)
  else()
      set(SUPERLU_SHARED ON)
  endif()

  include(ExternalProject)
  ExternalProject_Add(library_SUPERLU
          URL       http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_5.2.1.tar.gz
          PREFIX              "${INSTALL_DIRECTORY}/superlu"
          UPDATE_COMMAND ""
          CMAKE_ARGS
          -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
          -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
          -DEXAMPLES=ON
          -DCMAKE_BUILD_TYPE=Release
          DEPENDS blas lapack gfortran
          #-DCMAKE_C_FLAGS=-w -m64 -fPIC
          # -DCMAKE_Fortran_FLAGS=-w -m64 -fPIC
          )
  ExternalProject_Get_Property(library_SUPERLU INSTALL_DIR)
  set(SUPERLU_INCLUDE_DIRS ${INSTALL_DIR}/include)

  add_library(superlu UNKNOWN IMPORTED)
  set_target_properties(superlu
          PROPERTIES
          IMPORTED_LOCATION "${INSTALL_DIR}/lib/libsuperlu${CUSTOM_SUFFIX}"
          INTERFACE_LINK_LIBRARIES "blas;lapack;gfortran"
          INTERFACE_LINK_FLAGS            "-lpthread"
          )

  add_dependencies(superlu library_SUPERLU blas lapack gfortran )
  target_link_libraries(${PROJECT_NAME} PRIVATE superlu -lpthread -m64)
  target_include_directories(${PROJECT_NAME} PRIVATE ${SUPERLU_INCLUDE_DIRS})
  #For convenience, define these variables
  get_target_property(SUPERLU_LIBRARIES superlu IMPORTED_LOCATION)

