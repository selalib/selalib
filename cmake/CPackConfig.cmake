OPTION (SELALIB_NO_PACKAGES "CPACK - Disable packaging" OFF)

MARK_AS_ADVANCED (SELALIB_NO_PACKAGES)

#-----------------------------------------------------------------------------
# Set the cpack variables
#-----------------------------------------------------------------------------
IF (NOT SELALIB_NO_PACKAGES)

  SET (CPACK_BINARY_DEB                "ON")
  SET (CPACK_GENERATOR                 "DEB")
  SET (CPACK_PACKAGE_VENDOR            "IRMA")
  SET (CPACK_PACKAGE_NAME              "selalib")
  SET (CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME}")
  SET (CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${LIB_TYPE}")
  SET (CPACK_PACKAGE_VERSION_MAJOR     "1")
  SET (CPACK_PACKAGE_VERSION_MINOR     "0")
  SET (CPACK_PACKAGE_VERSION_PATCH     "0")
  SET (CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

  SET (CPACK_PACKAGE_ARCHITECTURE "amd64")
  SET (CPACK_DEBIAN_PACKAGE_DEPENDS "libhdf5-openmpi-dev,libopenmpi-dev,libfftw3-mpi-dev,liblapack-dev,openmpi-bin")

  SET (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/RELEASE.txt")
  SET (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING.txt")
  SET (CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
  SET (CPACK_PACKAGE_RELOCATABLE TRUE)

  SET (CPACK_PACKAGING_INSTALL_PREFIX "/usr")
  SET (CPACK_COMPONENTS_ALL_IN_ONE_PACKAGE ON)

  SET (CPACK_PACKAGE_MAINTAINER "selalib-user@lists.gforge.inria.fr")
  SET (CPACK_PACKAGE_CONTACT "pierre.navaror@inria.fr")
  SET (CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}_${CPACK_PACKAGE_VERSION}_${CPACK_PACKAGE_ARCHITECTURE}")
    
  SET (CPACK_PACKAGE_URL "http://selalib.gforge.inria.fr")
  SET (CPACK_PACKAGE_SUMMARY "Modular library for the gyrokinetic simulation model by a semi-Lagrangian method.")
  SET (CPACK_PACKAGE_DESCRIPTION 
        " As part of french action Fusion, Calvi INRIA project developed 
          in collaboration with CEA Cadarache GYSELA simulation code for 
          gyrokinetic simulation of plasma turbulence in Tokamaks. Development 
          and testing of new numerical methods is generally done on different 
          simplified models before its integration into GYSELA. No specification 
          is implemented for this aim, which involves several rewriting code 
          lines before the final integration, which is cumbersome and inefficient. 
          SeLaLib is an API for components of basic numerical implementation also 
          in GYSELA, a framework for parallel single streamline built in order 
          to improve reliability and simplify the work of development.  "
    )
  
  INCLUDE(InstallRequiredSystemLibraries)

  INCLUDE (CPack)

  #---------------------------------------------------------------------------
  # Now list the cpack commands
  #---------------------------------------------------------------------------
  CPACK_ADD_COMPONENT (libraries DISPLAY_NAME "SELALIB Libraries"
                       GROUP Runtime)
  CPACK_ADD_COMPONENT (headers   DISPLAY_NAME "SELALIB Headers" 
                       DEPENDS libraries
                       GROUP Development)
  CPACK_ADD_COMPONENT (hdfdocuments DISPLAY_NAME "SELALIB Documents"
                       GROUP Documents)
  CPACK_ADD_COMPONENT (configinstall DISPLAY_NAME "SELALIB CMake files" 
                       DEPENDS libraries
                       GROUP Development)

ENDIF(NOT SELALIB_NO_PACKAGES)
