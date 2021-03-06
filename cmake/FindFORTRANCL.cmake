set(FORTRANCL_ROOT
    "/usr/local"
    CACHE PATH "Root dir for fortrancl library")

find_path(
  FORTRANCL_INCLUDE_DIRS
  NAMES cl.mod
  HINTS /opt/local ${FORTRANCL_ROOT}
  PATH_SUFFIXES include
  DOC "PATH TO cl.mod")

find_library(
  FORTRANCL_LIBRARIES
  NAMES fortrancl
  HINTS /opt/local ${FORTRANCL_ROOT}
  PATH_SUFFIXES lib
  DOC "PATH TO libfortrancl.a")

find_library(OPENCL_LIBRARIES OpenCL DOC "OpenCL lib")

if(FORTRANCL_INCLUDE_DIRS
   AND FORTRANCL_LIBRARIES
   AND OPENCL_LIBRARIES)
  set(FORTRANCL_FOUND YES)
endif()
