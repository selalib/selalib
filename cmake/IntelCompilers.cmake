# compiler flags for intel 8.x, 9.x, 10.x, 11.x check the version and warn users
# if older compilers are being used SSE4.2 option is available for 11.1  and
# higher

# common options for intel compilers enable Interprocedural (IP) Optimizations
set(INTEL_COMPILER 1)
add_definitions(-DADD_ -DINLINE_ALL=inline)
set(ENABLE_OPENMP 1)
# set(INTEL_OPTS "-g -restrict -unroll -ansi-alias -O3 -ip -openmp")
set(INTEL_OPTS "-g -restrict -unroll -O3 -ip -openmp")

# MESSAGE(STATUS "CPU_IDENTITY ${CPU_IDENTITY}")

# grep icpc version and determine what options to be used
exec_program(
  icpc ARGS
  -v
  OUTPUT_VARIABLE ICPC_OUTPUT
  RETURN_VALUE ICPC_INFO)
string(REGEX MATCH "[0-9]+" ICPC_VERSION "${ICPC_OUTPUT}")
string(REGEX REPLACE "[0-9]+\\.([0-9]+)" "\\1" ICPC_MINOR_VERSION
                     "${ICPC_OUTPUT}")

# Use deprecated options prior to 11.1
set(ICC_DEPRECATED_OPTS FALSE)
if(ICPC_VERSION LESS 11)
  set(ICC_DEPRECATED_OPTS TRUE)
else(ICPC_VERSION LESS 11)
  if(ICPC_MINOR_VERSION LESS 1)
    set(ICC_DEPRECATED_OPTS TRUE)
  endif(ICPC_MINOR_VERSION LESS 1)
endif(ICPC_VERSION LESS 11)

# -ipo_obj force generation of real object files (requires -ipo)
# set(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3 -Ob=1 -ansi -ipo
# -ipo_obj -cxxlib-icc") set(CMAKE_CC_FLAGS "-restrict -unroll -fno-alias -O3
# -Ob=1 -ansi -ipo -ipo_obj") set(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias
# -g -ansi") set(CMAKE_CC_FLAGS "-restrict -unroll -fno-alias -g -ansi")
# set(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3 -cxxlib-icc")
# set(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3 -ansi -fno-fnalias
# -ivdep_parallel -Ob=2") set(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3
# -ivdep_parallel -Ob=2") set(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3
# -Ob=2 -qp") set(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3 -Ob=2
# -cxxlib-icc") set(CMAKE_CXX_FLAGS "-restrict -unroll -fno-alias -O3")

# check if -ftz is accepted
set(CMAKE_TRY_ICC_FLAGS "-ftz")
check_c_compiler_flag(${CMAKE_TRY_ICC_FLAGS} INTEL_CC_FLAGS)
if(INTEL_CC_FLAGS)
  set(INTEL_OPTS "${INTEL_OPTS} ${CMAKE_TRY_ICC_FLAGS}")
endif(INTEL_CC_FLAGS)

if(ICC_DEPRECATED_OPTS)
  set(INTEL_OPTS "${INTEL_OPTS} -prefetch")
else(ICC_DEPRECATED_OPTS)
  set(INTEL_OPTS "${INTEL_OPTS} -opt-prefetch")
endif(ICC_DEPRECATED_OPTS)

exec_program(
  grep ARGS
  flags /proc/cpuinfo
  OUTPUT_VARIABLE CPU_FLAGS
  RETURN_VALUE CPUINFO)

set(HAVE_SSE4_2 0)
set(HAVE_SSE4_1 0)
set(HAVE_SSSE3 0)
set(HAVE_SSE3 0)
set(HAVE_SSE2 0)
set(HAVE_SSE 0)

set(SSE_OPT_SET FALSE)
if(CPU_FLAGS MATCHES "sse4_2")
  if(ICC_DEPRECATED_OPTS)
    message(WARNING "SSE4.2 needs version 11.1 and higher.")
  else(ICC_DEPRECATED_OPTS)
    set(HAVE_SSE4_2 1)
    set(SSE_OPT_SET TRUE)
    set(INTEL_OPTS "${INTEL_OPTS} -xSSE4.2")
  endif(ICC_DEPRECATED_OPTS)
endif(CPU_FLAGS MATCHES "sse4_2")

if(CPU_FLAGS MATCHES "sse4_1")
  set(HAVE_SSE4_1 1)
  if(NOT SSE_OPT_SET)
    set(SSE_OPT_SET TRUE)
    if(ICC_DEPRECATED_OPTS)
      set(INTEL_OPTS "${INTEL_OPTS} -xS")
    else(ICC_DEPRECATED_OPTS)
      set(INTEL_OPTS "${INTEL_OPTS} -xSSE4.1")
    endif(ICC_DEPRECATED_OPTS)
  endif(NOT SSE_OPT_SET)
endif(CPU_FLAGS MATCHES "sse4_1")

if(CPU_FLAGS MATCHES "ssse3")
  set(HAVE_SSSE3 1)
  if(NOT SSE_OPT_SET)
    set(SSE_OPT_SET TRUE)
    if(ICC_DEPRECATED_OPTS)
      set(INTEL_OPTS "${INTEL_OPTS} -xT")
    else(ICC_DEPRECATED_OPTS)
      set(INTEL_OPTS "${INTEL_OPTS} -xSSSE3")
    endif(ICC_DEPRECATED_OPTS)
  endif(NOT SSE_OPT_SET)
endif(CPU_FLAGS MATCHES "ssse3")

if(CPU_FLAGS MATCHES "sse3")
  set(HAVE_SSE3 1)
  if(NOT SSE_OPT_SET)
    set(SSE_OPT_SET TRUE)
    if(ICC_DEPRECATED_OPTS)
      set(INTEL_OPTS "${INTEL_OPTS} -xP")
    else(ICC_DEPRECATED_OPTS)
      set(INTEL_OPTS "${INTEL_OPTS} -xSSE3")
    endif(ICC_DEPRECATED_OPTS)
  endif(NOT SSE_OPT_SET)
endif(CPU_FLAGS MATCHES "sse3")

if(CPU_FLAGS MATCHES "sse2")
  set(HAVE_SSE2 1)
  if(NOT SSE_OPT_SET)
    set(SSE_OPT_SET TRUE)
    if(ICC_DEPRECATED_OPTS)
      set(INTEL_OPTS "${INTEL_OPTS} -xN")
    else(ICC_DEPRECATED_OPTS)
      set(INTEL_OPTS "${INTEL_OPTS} -xSSE2")
    endif(ICC_DEPRECATED_OPTS)
  endif(NOT SSE_OPT_SET)
endif(CPU_FLAGS MATCHES "sse2")

# SET(CMAKE_TRY_ICC_FLAGS "-xSSE4.1")
# CHECK_C_COMPILER_FLAG(${CMAKE_TRY_ICC_FLAGS} INTEL_CC_FLAGS)
# IF(INTEL_CC_FLAGS) SET(HAVE_SSE4_2 1) IF(NOT SSE_OPT_SET) SET(INTEL_OPTS
# "${INTEL_OPTS} ${CMAKE_TRY_ICC_FLAGS}") SET(SSE_OPT_SET TRUE) ENDIF(NOT
# SSE_OPT_SET) ENDIF(INTEL_CC_FLAGS)
#
# SET(CMAKE_TRY_ICC_FLAGS "-xSSSE3")
# CHECK_C_COMPILER_FLAG(${CMAKE_TRY_ICC_FLAGS} INTEL_CC_FLAGS)
# IF(INTEL_CC_FLAGS) SET(HAVE_SSSE3 1) IF(NOT SSE_OPT_SET) SET(INTEL_OPTS
# "${INTEL_OPTS} ${CMAKE_TRY_ICC_FLAGS}") SET(SSE_OPT_SET TRUE) ENDIF(NOT
# SSE_OPT_SET) ENDIF(INTEL_CC_FLAGS)

# SET(CMAKE_TRY_ICC_FLAGS "-xSSE3") CHECK_C_COMPILER_FLAG(${CMAKE_TRY_ICC_FLAGS}
# INTEL_CC_FLAGS) IF(INTEL_CC_FLAGS) SET(HAVE_SSE3 1) IF(NOT SSE_OPT_SET)
# SET(INTEL_OPTS "${INTEL_OPTS} ${CMAKE_TRY_ICC_FLAGS}") SET(SSE_OPT_SET TRUE)
# ENDIF(NOT SSE_OPT_SET) ENDIF(INTEL_CC_FLAGS)
#
# SET(CMAKE_TRY_ICC_FLAGS "-xSSE2") CHECK_C_COMPILER_FLAG(${CMAKE_TRY_ICC_FLAGS}
# INTEL_CC_FLAGS) IF(INTEL_CC_FLAGS) SET(HAVE_SSE2 1) IF(NOT SSE_OPT_SET)
# SET(INTEL_OPTS "${INTEL_OPTS} ${CMAKE_TRY_ICC_FLAGS}") SET(SSE_OPT_SET TRUE)
# ENDIF(NOT SSE_OPT_SET) ENDIF(INTEL_CC_FLAGS)
#
# SET(CMAKE_TRY_ICC_FLAGS "-xSSE") CHECK_C_COMPILER_FLAG(${CMAKE_TRY_ICC_FLAGS}
# INTEL_CC_FLAGS) IF(INTEL_CC_FLAGS) SET(HAVE_SSE 1) IF(NOT SSE_OPT_SET)
# SET(INTEL_OPTS "${INTEL_OPTS} ${CMAKE_TRY_ICC_FLAGS}") SET(SSE_OPT_SET TRUE)
# ENDIF(NOT SSE_OPT_SET) ENDIF(INTEL_CC_FLAGS)

# if($ENV{CXX} matches cmpi) set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ccl icpc")
# endif($ENV{CXX} matches cmpi)

# IF(${CMAKE_SYSTEM_PROCESSOR} matches "unknown") IF(CPU_IDENTITY matches "E5")
# set(CMAKE_CXX_FLAGS "${INTEL_OPTS} -xT ") set(CMAKE_C_FLAGS "${INTEL_OPTS}
# -xT") ELSE(CPU_IDENTITY matches "E5") set(CMAKE_CXX_FLAGS "${INTEL_OPTS} -xP
# ") set(CMAKE_C_FLAGS "${INTEL_OPTS} -xP") endif(CPU_IDENTITY matches "E5")
# set(HAVE_SSE 1) set(HAVE_SSE2 1) set(HAVE_SSE3 1) set(HAVE_SSSE3 1)
# #set(USE_PREFETCH 1) #set(PREFETCH_AHEAD 12) endif(${CMAKE_SYSTEM_PROCESSOR}
# matches "unknown")

set(CMAKE_CXX_FLAGS "$ENV{CXX_FLAGS} ${INTEL_OPTS} -Wno-deprecated")
set(CMAKE_C_FLAGS "$ENV{CC_FLAGS} ${INTEL_OPTS} -std=c99 -Wno-deprecated")

# ifc -> ifort set(FORTRAN_LIBS " -lifcore -lifport")
set(F77 ifort)
set(F77OPTFLAGS -fpp2 -O3)
set(F77FLAGS ${F77OPTFLAGS})

# ##############################################################################
# KCC needs to be used to build static libraries
# ##############################################################################
# set(CMAKE_AR xild) set(CMAKE_CXX_CREATE_STATIC_LIBRARY "<CMAKE_AR> -lib -o
# <TARGET> <OBJECTS>")
