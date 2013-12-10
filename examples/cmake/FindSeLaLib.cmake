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
SET(TRIAL_PATHS $ENV{SELALIB_ROOT}/usr
                /usr/local)

FIND_PATH (SELALIB_ROOT include/selalib.h
           HINTS ENV SELALIB_ROOT
           PATHS /usr/local $ENV{HOME}/local
           DOC "SeLaLib directory")

FIND_PATH(SELALIB_INCLUDES NAMES selalib.h 
                           HINTS ${TRIAL_PATHS} 
                           PATH_SUFFIXES include 
                           DOC "SeLaLib include path")

FIND_LIBRARY(SELALIB_LIBRARIES NAMES selalib 
                               HINTS ${TRIAL_PATHS}   
                               PATH_SUFFIXES lib 
                               DOC "SeLaLib libraries")

IF(SELALIB_LIBRARIES AND SELALIB_INCLUDES)
   SET(SELALIB_FOUND TRUE)
ELSE()
   SET(SELALIB_FOUND FALSE)
ENDIF()
