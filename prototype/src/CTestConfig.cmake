FIND_PROGRAM(LSB_RELEASE_COMMAND NAMES lsb_release)
if(LSB_RELEASE_COMMAND)
   EXECUTE_PROCESS(
      COMMAND ${LSB_RELEASE_COMMAND} -is
      COMMAND tr -d '\n'
      OUTPUT_VARIABLE LSB_ID
   )
   EXECUTE_PROCESS(
      COMMAND ${LSB_RELEASE_COMMAND} -rs
      COMMAND tr -d '\n'
      OUTPUT_VARIABLE LSB_VER
   )

   set(LINUX_NAME "${LSB_ID}-${LSB_VER}")
else()
if(EXISTS "/etc/issue")
      set(LINUX_NAME "")
      file(READ "/etc/issue" LINUX_ISSUE)
      # Fedora case
      if(LINUX_ISSUE MATCHES "Fedora")
        string(REGEX MATCH "release ([0-9]+)" FEDORA "${LINUX_ISSUE}")
        set(LINUX_NAME "Fedora")
	set(LINUX_VER "${CMAKE_MATCH_1}")
      endif(LINUX_ISSUE MATCHES "Fedora")
      # CentOS case
      if(LINUX_ISSUE MATCHES "CentOS")
        string(REGEX MATCH "release ([0-9]+\\.[0-9]+)" CENTOS "${LINUX_ISSUE}")
        set(LINUX_NAME "CentOS")        
	set(LINUX_VER "${CMAKE_MATCH_1}")
      endif(LINUX_ISSUE MATCHES "CentOS")
      # Redhat case
      # Red Hat Enterprise Linux Server release 5 (Tikanga)
      if(LINUX_ISSUE MATCHES "Red Hat")
        string(REGEX MATCH "release ([0-9]+\\.*[0-9]*)" REDHAT "${LINUX_ISSUE}")
        set(LINUX_NAME "RedHat")        
	set(LINUX_VER "${CMAKE_MATCH_1}")
      endif(LINUX_ISSUE MATCHES "Red Hat")
   set(LINUX_NAME "${LINUX_NAME}-${LINUX_VER}")
endif()
endif()

MESSAGE(STATUS "LINUX_NAME:${LINUX_NAME}")

SET(CTEST_PROJECT_NAME "Selalib")
SET(CTEST_NIGHTLY_START_TIME "00:00:00 EST")
SET(CTEST_DROP_METHOD "http")
SET(CTEST_DROP_SITE_CDASH TRUE)
SET(CTEST_DROP_SITE "cdash.inria.fr")
SET(CTEST_DROP_LOCATION "/CDash/submit.php?project=Selalib")
SET(CTEST_TRIGGER_SITE "")
FIND_PROGRAM(CTEST_GIT_COMMAND NAMES git)
SET(UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

IF(LINUX_NAME)
SET(BUILDNAME "${Fortran_COMPILER_NAME}-${Fortran_COMPILER_VERSION}-${LINUX_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
ELSE()
SET(BUILDNAME "${Fortran_COMPILER_NAME}-${Fortran_COMPILER_VERSION}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
ENDIF()
MESSAGE(STATUS "BUILDNAME:${BUILDNAME}")
