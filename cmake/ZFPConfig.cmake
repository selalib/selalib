# Find the ZFP compression library,
# http://computation.llnl.gov/projects/floating-point-compression

if(DEFINED ENV{ZFP_ROOT})
  set(ZFP_ROOT
      $ENV{ZFP_ROOT}
      CACHE PATH "ZFP installation location")
endif()

find_path(
  ZFP_INCLUDE_DIR
  NAMES zfp.h
  HINTS ${ZFP_ROOT}
  PATH_SUFFIXES include)

find_library(
  ZFP_LIBRARIES
  NAMES zfp libzfp
  HINTS ${ZFP_ROOT}
  PATH_SUFFIXES lib lib64)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ZFP DEFAULT_MSG ZFP_INCLUDE_DIR ZFP_LIBRARIES)

if(ZFP_FOUND)
  include_directories(${ZFP_INCLUDE_DIR})
  message(STATUS "ZFP_INCLUDE_DIR:${ZFP_INCLUDE_DIR}")
  message(STATUS "ZFP_LIBRARIES:${ZFP_LIBRARIES}")
  add_definitions(-DUSE_ZFP)
endif()
