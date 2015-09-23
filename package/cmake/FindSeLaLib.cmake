# - Try to find SeLaLib
# Once done this will define
#
#  SELALIB_FOUND        - system has SeLaLib
#  SELALIB_INCLUDES     - the SeLaLib include directories
#  SELALIB_LIBRARIES    - Link these to use SeLaLib
#
#  Usage:
#  find_package(SeLaLib)
#
SET(TRIAL_PATHS /usr /usr/local /opt/local)

FIND_PATH(SELALIB_INCLUDE_DIRS NAMES selalib.h 
                           HINTS ${TRIAL_PATHS} 
                           PATH_SUFFIXES include/selalib
                           DOC "SeLaLib include path")

FIND_LIBRARY(SELALIB_LIBRARIES NAMES selalib
                               HINTS ${TRIAL_PATHS} 
                               PATH_SUFFIXES lib 
                               DOC "SeLaLib libraries")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SELALIB DEFAULT_MSG
                                  SELALIB_INCLUDE_DIRS SELALIB_LIBRARIES)

