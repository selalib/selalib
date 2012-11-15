<<<<<<< HEAD:cmake/DartConfig.cmake
SET(CTEST_PROJECT_NAME "Pipol")
SET(NIGHTLY_START_TIME "20:00:00 EST")
SET(DROP_METHOD "http")
SET(DROP_SITE "cdash.inria.fr")
SET(DROP_LOCATION "/CDash/submit.php?project=Pipol")
SET(DROP_SITE_CDASH TRUE)
=======
SITE_NAME(HOSTNAME)
MESSAGE(STATUS "HOSTNAME:${HOSTNAME}")
IF(HOSTNAME MATCHES "inrialpes")
   SET(CTEST_PROJECT_NAME "Pipol")
   SET(NIGHTLY_START_TIME "20:00:00 EST")
   SET(DROP_METHOD "http")
   SET(DROP_SITE "cdash.inria.fr")
   SET(DROP_LOCATION "/CDash/submit.php?project=Pipol")
   SET(DROP_SITE_CDASH TRUE)
ENDIF()
>>>>>>> origin/phdf5-3d:cmake/DartConfig.cmake
