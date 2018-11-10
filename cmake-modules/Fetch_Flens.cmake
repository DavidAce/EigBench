


#find_package(Flens 3.3)
if(FLENS_FOUND)
    message(STATUS "EIGEN FOUND IN SYSTEM: ${FLENS_INCLUDE_DIR}")
    add_library(FLENS INTERFACE)
else()
    message(STATUS "Flens will be installed into ${INSTALL_DIRECTORY}/flens on first build.")

    include(ExternalProject)
    ExternalProject_Add(library_FLENS
            GIT_REPOSITORY https://github.com/michael-lehn/FLENS.git
            GIT_TAG public
            GIT_PROGRESS 1
            PREFIX "${INSTALL_DIRECTORY}/flens"
            CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#            UPDATE_COMMAND ""
#            TEST_COMMAND ""
#            INSTALL_COMMAND ""
#            CONFIGURE_COMMAND ""
#            BUILD_COMMAND
#            ${CMAKE_COMMAND} -E make_directory <INSTALL_DIR>/include/flens && find <INSTALL_DIR>/include/flens -maxdepth 1 -type l -delete &&
#            ln -s <SOURCE_DIR>/Eigen/ <SOURCE_DIR>/unsupported/ <INSTALL_DIR>/include/flens/
            )


    ExternalProject_Get_Property(library_FLENS INSTALL_DIR)
    add_library(FLENS INTERFACE)
    set(FLENS_INCLUDE_DIR ${INSTALL_DIR}/include)
    add_dependencies(FLENS library_FLENS)
endif()



set_target_properties(FLENS PROPERTIES
        INTERFACE_INCLUDE_DIRECTORY     "${FLENS_INCLUDE_DIR}"
        INTERFACE_COMPILE_OPTIONS       "${FLENS_COMPILER_FLAGS}"
        )
target_link_libraries(${PROJECT_NAME} PRIVATE FLENS)
# Add SYSTEM flag to suppress warnings
target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE ${FLENS_INCLUDE_DIR})
target_compile_options(${PROJECT_NAME} PRIVATE ${FLENS_COMPILER_FLAGS})
