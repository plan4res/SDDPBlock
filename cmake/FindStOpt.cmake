# --------------------------------------------------------------------------- #
#    CMake find module for StOpt                                              #
#                                                                             #
#    This module finds StOpt include directories and libraries.               #
#    Use it by invoking find_package() with the form:                         #
#                                                                             #
#        find_package(StOpt [version] [EXACT] [REQUIRED])                     #
#                                                                             #
#    The results are stored in the following variables:                       #
#                                                                             #
#        <LibraryName>_FOUND         - True if headers are found              #
#        <LibraryName>_INCLUDE_DIRS  - Include directories                    #
#        <LibraryName>_LIBRARIES     - Libraries to be linked                 #
#        <LibraryName>_VERSION       - Version number                         #
#                                                                             #
#    The search results are saved in these persistent cache entries:          #
#                                                                             #
#        <LibraryName>_INCLUDE_DIR   - Directory containing headers           #
#        <LibraryName>_LIBRARY       - The found library                      #
#                                                                             #
#    for each <LibraryName> in: StOpt, StOpt_geners.                          #
#                                                                             #
#    This module can read a search path from the variable:                    #
#                                                                             #
#        StOpt_ROOT          - Preferred StOpt location                       #
#                                                                             #
#    The following IMPORTED targets are also defined:                         #
#                                                                             #
#        StOpt::StOpt                                                         #
#        StOpt::geners                                                        #
#                                                                             #
#    This find module is provided because StOpt does not provide              #
#    a CMake configuration file on its own.                                   #
#                                                                             #
#                              Niccolo' Iardella                              #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# ----- Requirements -------------------------------------------------------- #
find_package(BZip2 REQUIRED QUIET)
find_package(ZLIB REQUIRED QUIET)
find_package(Boost REQUIRED COMPONENTS system timer QUIET)

# This will try first with Eigen3 own configuration file,
# then with the find module we provide.
find_package(Eigen3 QUIET CONFIG)
if (NOT Eigen3_FOUND)
    find_package(Eigen3 REQUIRED)
endif ()

# Check if already in cache
if (StOpt_geners_INCLUDE_DIR AND StOpt_geners_LIBRARY AND
    StOpt_INCLUDE_DIR AND StOpt_LIBRARY)
    set(StOpt_FOUND TRUE)
else ()

    # ----- Find the geners library ----------------------------------------- #
    # Note that find_path() creates a cache entry
    find_path(StOpt_geners_INCLUDE_DIR
              NAMES geners/uriUtils.hh
              DOC "geners include directory.")

    # Note that find_library() creates a cache entry
    find_library(StOpt_geners_LIBRARY
                 NAMES geners
                 DOC "geners library.")

    # ----- Find the StOpt library ------------------------------------------ #
    # Note that find_path() creates a cache entry
    find_path(StOpt_INCLUDE_DIR
              NAMES StOpt/sddp/OptimizerSDDPBase.h
              DOC "StOpt include directory.")

    # Note that find_library() creates a cache entry
    find_library(StOpt_LIBRARY
                 NAMES StOpt
                 DOC "StOpt library.")

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET, REQUIRED and version-related
    # options passed to find_package(). It also sets <PackageName>_FOUND if
    # REQUIRED_VARS are set.
    # REQUIRED_VARS should be cache entries and not output variables. See:
    # https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
    find_package_handle_standard_args(
            StOpt
            REQUIRED_VARS
            StOpt_geners_INCLUDE_DIR StOpt_geners_LIBRARY
            StOpt_LIBRARY StOpt_INCLUDE_DIR)
endif ()

# ----- Export the targets -------------------------------------------------- #
if (StOpt_FOUND)
    set(StOpt_geners_INCLUDE_DIRS "${StOpt_geners_INCLUDE_DIR}")
    set(StOpt_geners_LIBRARIES "${StOpt_geners_LIBRARY}")

    if (NOT TARGET StOpt::geners)
        add_library(StOpt::geners UNKNOWN IMPORTED)
        set_target_properties(
                StOpt::geners PROPERTIES
                IMPORTED_LOCATION "${StOpt_geners_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${StOpt_geners_INCLUDE_DIR}")
    endif ()

    set(StOpt_INCLUDE_DIRS "${StOpt_INCLUDE_DIR}")
    set(StOpt_LIBRARIES "${StOpt_LIBRARY}")

    if (NOT TARGET StOpt::StOpt)
        add_library(StOpt::StOpt UNKNOWN IMPORTED)
        set_target_properties(
                StOpt::StOpt PROPERTIES
                IMPORTED_LOCATION "${StOpt_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${StOpt_INCLUDE_DIR}"
                INTERFACE_LINK_LIBRARIES "StOpt::geners;Eigen3::Eigen;BZip2::BZip2;ZLIB::ZLIB;Boost::system;Boost::timer")
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(StOpt_geners_INCLUDE_DIR
                 StOpt_geners_LIBRARY
                 StOpt_INCLUDE_DIR
                 StOpt_LIBRARY)

# --------------------------------------------------------------------------- #
