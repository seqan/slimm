# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/slimm/CMakeLists.txt
#
# CMakeLists.txt file for slimm.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_apps_slimm)
message (STATUS "Configuring apps/slimm")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------
# Search SeqAn and select dependencies.
find_package(ZLIB)
find_package(BZip2)
find_package(SeqAn REQUIRED)

if (NOT CXX11_FOUND)
    message (FATAL_ERROR "  C++11 not found: Not building SLIMM.")
    return ()
endif (NOT CXX11_FOUND)

# Do not build without Zlib.
if (NOT ZLIB_FOUND)
  message (WARNING "WARNING: Zlib not found! SLIMM will not work compressed files including BAM files")
endif (NOT ZLIB_FOUND)

# Warn if BZip2 was not found.
if (NOT BZIP2_FOUND)
  message (WARNING "WARNING: BZip2 not found. Building SLIMM without support for bzipped input and output.")
endif (NOT BZIP2_FOUND)

if(CMAKE_COMPILER_IS_GNUCXX)
  if (491 GREATER _GCC_VERSION)
    message (FATAL_ERROR "Your GCC version is too old. Minimum version is GCC-4.9.1!")
    return ()
  endif(491 GREATER _GCC_VERSION)
endif(CMAKE_COMPILER_IS_GNUCXX)

# ----------------------------------------------------------------------------
# App-Level Configuration
# ----------------------------------------------------------------------------

set (SEQAN_APP_VERSION "0.1.9")

# ----------------------------------------------------------------------------
# Change to approprate compiler flags
# ----------------------------------------------------------------------------
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()


# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")

# Update the list of file names below if you add source files to your application.
add_executable(slimm    slimm.cpp
                        slimm.h
                        timer.h
                        load_mapping_files.h
                        librf.h)


# Add dependencies found by find_package (SeqAn).
target_link_libraries (slimm ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (slimm)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install slimm in ${PREFIX}/bin directory
install (TARGETS slimm
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/slimm for SeqAn release builds.
install (FILES LICENSE README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (slimm)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable slimm in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} slimm CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:slimm")
  set (CPACK_PACKAGE_NAME "slimm")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "slimm - a tool to detect microbes from shotgun metagenomes")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Temesgen H. Dadi <temesgen.dadi@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (slimm "slimm")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:slimm")
