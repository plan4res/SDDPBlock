# --------------------------------------------------------------------------- #
#    Custom CMake find module for Eigen3                                      #
#                                                                             #
#    This module finds Eigen3 include directories.                            #
#    Use it by invoking find_package() with the form:                         #
#                                                                             #
#        find_package(Eigen3 [version] [EXACT] [REQUIRED])                    #
#                                                                             #
#    The results are stored in the following variables:                       #
#                                                                             #
#        Eigen3_FOUND         - True if headers are found                     #
#        Eigen3_INCLUDE_DIRS  - Include directories                           #
#        Eigen3_VERSION       - Version number                                #
#                                                                             #
#    The search results are saved in these persistent cache entries:          #
#                                                                             #
#        Eigen3_INCLUDE_DIR   - Directory containing headers                  #
#                                                                             #
#    This module can read a search path from the variable:                    #
#                                                                             #
#        Eigen3_ROOT          - Preferred Eigen3 location                     #
#                                                                             #
#    The following IMPORTED target is also defined:                           #
#                                                                             #
#        Eigen3::Eigen                                                        #
#                                                                             #
#    This find module is provided because the user may not have Eigen3        #
#    configuration file installed in the system default directories, which    #
#    happens to be the case for some of our main developers and testers.      #
#                                                                             #
#                              Niccolo' Iardella                              #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# Check if already in cache
if (Eigen3_INCLUDE_DIR)
    set(Eigen3_FOUND TRUE)
else ()

    # ----- Find the headers ------------------------------------------------ #
    # Note that find_path() also creates a cache entry
    find_path(Eigen3_INCLUDE_DIR Eigen/Dense
              PATH_SUFFIXES eigen3
              DOC "Eigen3 include directory.")

    # ----- Parse the version ----------------------------------------------- #
    if (Eigen3_INCLUDE_DIR)
        file(STRINGS
             "${Eigen3_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h"
             _eigen3_version_lines REGEX "#define EIGEN_(WORLD|MAJOR|MINOR)_VERSION")

        string(REGEX REPLACE ".*EIGEN_WORLD_VERSION *\([0-9]*\).*" "\\1" _eigen3_version_major "${_eigen3_version_lines}")
        string(REGEX REPLACE ".*EIGEN_MAJOR_VERSION *\([0-9]*\).*" "\\1" _eigen3_version_minor "${_eigen3_version_lines}")
        string(REGEX REPLACE ".*EIGEN_MINOR_VERSION *\([0-9]*\).*" "\\1" _eigen3_version_patch "${_eigen3_version_lines}")

        set(Eigen3_VERSION "${_eigen3_version_major}.${_eigen3_version_minor}.${_eigen3_version_patch}")
        unset(_eigen3_version_lines)
        unset(_eigen3_version_major)
        unset(_eigen3_version_minor)
        unset(_eigen3_version_patch)
    endif ()

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET, REQUIRED and version-related
    # options passed to find_package(). It also sets <PackageName>_FOUND if
    # REQUIRED_VARS are set.
    # REQUIRED_VARS should be cache entries and not output variables. See:
    # https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
    find_package_handle_standard_args(Eigen3
                                      REQUIRED_VARS Eigen3_INCLUDE_DIR
                                      VERSION_VAR Eigen3_VERSION)
endif ()

# ----- Export the target --------------------------------------------------- #
if (Eigen3_FOUND)
    set(Eigen3_INCLUDE_DIRS "${Eigen3_INCLUDE_DIR}")

    if (NOT TARGET Eigen3::Eigen)
        add_library(Eigen3::Eigen INTERFACE IMPORTED)
        set_target_properties(
                Eigen3::Eigen PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${Eigen3_INCLUDE_DIR}")
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(Eigen3_INCLUDE_DIR Eigen3_VERSION)

# --------------------------------------------------------------------------- #
