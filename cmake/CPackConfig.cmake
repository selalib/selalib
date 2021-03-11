set(CPACK_BINARY_DEB "ON")
set(CPACK_PACKAGE_VENDOR "IRMA")
set(CPACK_PACKAGE_NAME "selalib")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME}")
set(CPACK_PACKAGE_INSTALL_REGISTRY_KEY
    "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${LIB_TYPE}")
set(CPACK_PACKAGE_VERSION
    "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}"
)

set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/RELEASE.md")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING.txt")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_RELOCATABLE TRUE)

set(CPACK_PACKAGING_INSTALL_PREFIX "/usr/lib/selalib")
set(CPACK_COMPONENTS_ALL_IN_ONE_PACKAGE ON)
set(CPACK_COMPONENTS_ALL libraries headers fortheaders)

set(CPACK_PACKAGE_MAINTAINER "selalib-user@lists.gforge.inria.fr")
set(CPACK_PACKAGE_CONTACT "pierre.navaro@univ-rennes1.fr")

site_name(CPACK_HOSTNAME)

if(CPACK_HOSTNAME MATCHES "irma-suse")
  set(CPACK_GENERATOR "DEB")
  set(CPACK_DEBIAN_PACKAGE_DEPENDS
      "libhdf5-openmpi-dev,libopenmpi-dev,libfftw3-mpi-dev,liblapack-dev,openmpi-bin"
  )
  set(CPACK_PACKAGE_ARCHITECTURE "amd64")
  set(CPACK_PACKAGE_FILE_NAME
      "${CPACK_PACKAGE_NAME}_${CPACK_PACKAGE_VERSION}_${CPACK_PACKAGE_ARCHITECTURE}"
  )
endif()

if(CPACK_HOSTNAME MATCHES "selalib-fedora18")
  set(CPACK_GENERATOR "RPM")
  set(CPACK_PACKAGE_ARCHITECTURE "x86_64")
  set(CPACK_RPM_PACKAGE_REQUIRES
      "gcc-gfortran,hdf5-mpich2-devel,mpich2-devel,fftw-devel")
  set(CPACK_PACKAGE_FILE_NAME
      "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}.fc18.${CPACK_PACKAGE_ARCHITECTURE}"
  )
endif()

set(CPACK_PACKAGE_URL "http://selalib.gforge.inria.fr")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "Modular library for the gyrokinetic simulation model by a semi-Lagrangian method."
)
set(CPACK_PACKAGE_DESCRIPTION
    " As part of french action Fusion, Calvi INRIA project developed
       in collaboration with CEA Cadarache GYSELA simulation code for
       gyrokinetic simulation of plasma turbulence in Tokamaks. Development
       and testing of new numerical methods is generally done on different
       simplified models before its integration into GYSELA. No specification
       is implemented for this aim, which involves several rewriting code
       lines before the final integration, which is cumbersome and inefficient.
       SeLaLib is an API for components of basic numerical implementation also
       in GYSELA, a framework for parallel single streamline built in order
       to improve reliability and simplify the work of development.  ")

include(InstallRequiredSystemLibraries)

include(CPack)
