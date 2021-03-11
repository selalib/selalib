#
# checking processor name to pick optimization flags for linux and mac os x
#
if(CMAKE_SYSTEM_NAME MATCHES "Linux")
  exec_program(
    grep ARGS
    model /proc/cpuinfo
    OUTPUT_VARIABLE CPUINFO_MODEL
    RETURN_VALUE CPUINFO)

  if(CPUINFO_MODEL MATCHES "E5")
    set(CPU_IDENTITY "E5xxx")
    message("-- CPU Model Name Quad-core Xeon")
  endif(CPUINFO_MODEL MATCHES "E5")

  if(CPUINFO_MODEL MATCHES "E4")
    set(CPU_IDENTITY "E4xxx")
    message("-- CPU Model Name Core 2 Duo")
  endif(CPUINFO_MODEL MATCHES "E4")

  if(CPUINFO_MODEL MATCHES "Itanium")
    set(CPU_IDENTITY "IA64")
    message("-- CPU Model Name Itanium family")
  endif(CPUINFO_MODEL MATCHES "Itanium")

endif(CMAKE_SYSTEM_NAME MATCHES "Linux")

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  exec_program(
    system_profiler ARGS
    -detailLevel basic
    OUTPUT_VARIABLE MACSYSINFO_MODEL
    RETURN_VALUE CPUINFO)
  if(MACSYSINFO_MODEL MATCHES "Core 2 Duo")
    set(CPU_IDENTITY "E4xxx")
    message("-- CPU Model Name Core 2 Duo ")
  endif(MACSYSINFO_MODEL MATCHES "Core 2 Duo")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
