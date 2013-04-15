# SELALIB_INSTALL_LIB_DIR, SELALIB_INSTALL_INCLUDE_DIR :
#   Customize the 'lib', and 'include' installation directories.
#
OPTION (SELALIB_NO_PACKAGES "CPACK - Disable packaging" OFF)
MARK_AS_ADVANCED (SELALIB_NO_PACKAGES)

#-----------------------------------------------------------------------------
# Set the core names of all the libraries
#-----------------------------------------------------------------------------
SET (SELALIB_LIB_CORENAME             "selalib")
SET (SELALIB_MEMORY_LIB_CORENAME      "sll_memory")
SET (SELALIB_ASSERT_LIB_CORENAME      "sll_assert")
SET (SELALIB_WP_LIB_CORENAME          "sll_working_precision")
SET (SELALIB_CONSTANTS_LIB_CORENAME   "sll_constants")
SET (SELALIB_UTLITIES_LIB_CORENAME    "sll_utilities")

IF (NOT SELALIB_INSTALL_LIB_DIR)
  SET (SELALIB_INSTALL_LIB_DIR lib)
ENDIF (NOT SELALIB_INSTALL_LIB_DIR)

IF (NOT SELALIB_INSTALL_INCLUDE_DIR)
  SET (SELALIB_INSTALL_INCLUDE_DIR include)
ENDIF (NOT SELALIB_INSTALL_INCLUDE_DIR)

SET (SELALIB_PACKAGE "selalib")
SET (SELALIB_PACKAGE_NAME "SELALIB")
SET (SELALIB_PACKAGE_VERSION "1")
SET (SELALIB_PACKAGE_VERSION_MINOR "0")
SET (SELALIB_PACKAGE_VERSION_MAJOR "0")
SET (SELALIB_PACKAGE_URL "http://selalib.gforge.inria.fr")
SET (SELALIB_PACKAGE_BUGREPORT "selalib-dev@gforge.inria.fr")

#-----------------------------------------------------------------------------
# Set the cpack variables
#-----------------------------------------------------------------------------
IF (NOT SELALIB_NO_PACKAGES)

  SET (CPACK_PACKAGE_VENDOR "IRMA")
  SET (CPACK_PACKAGE_NAME "${SELALIB_PACKAGE_NAME}")
  SET (CPACK_PACKAGE_INSTALL_DIRECTORY "${SELALIB_PACKAGE_NAME}")
  SET (CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${SELALIB_PACKAGE_NAME}-${SELALIB_PACKAGE_VERSION}-${LIB_TYPE}")
  SET (CPACK_PACKAGE_VERSION "${SELALIB_PACKAGE_VERSION}")
  SET (CPACK_PACKAGE_VERSION_MAJOR "${SELALIB_PACKAGE_VERSION_MAJOR}")
  SET (CPACK_PACKAGE_VERSION_MINOR "${SELALIB_PACKAGE_VERSION_MINOR}")
  SET (CPACK_PACKAGE_VERSION_PATCH "")
  SET (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/RELEASE.txt")
  SET (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING.txt")
  SET (CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
  SET (CPACK_PACKAGE_RELOCATABLE TRUE)

  SET (CPACK_PACKAGING_INSTALL_PREFIX "/usr")
  SET (CPACK_COMPONENTS_ALL_IN_ONE_PACKAGE ON)

  SET (CPACK_DEBIAN_PACKAGE_SECTION "Libraries")
  SET (CPACK_DEBIAN_PACKAGE_MAINTAINER "${SELALIB_PACKAGE_BUGREPORT}")
    
  SET (CPACK_RPM_COMPONENT_INSTALL ON)
  SET (CPACK_RPM_PACKAGE_RELOCATABLE ON)
  SET (CPACK_RPM_PACKAGE_LICENSE "BSD-style")
  SET (CPACK_RPM_PACKAGE_GROUP "Development/Libraries")
  SET (CPACK_RPM_PACKAGE_URL "${SELALIB_PACKAGE_URL}")
  SET (CPACK_RPM_PACKAGE_SUMMARY "Modular library for the gyrokinetic simulation model by a semi-Lagrangian method.")
  SET (CPACK_RPM_PACKAGE_DESCRIPTION 
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

  SET (CPACK_INSTALL_CMAKE_PROJECTS "${SELALIB_BINARY_DIR};SELALIB;ALL;/")
  
  INCLUDE (CPack)

  #---------------------------------------------------------------------------
  # Now list the cpack commands
  #---------------------------------------------------------------------------
  CPACK_ADD_COMPONENT (libraries 
      DISPLAY_NAME "SELALIB Libraries"
      GROUP Runtime
  )
  CPACK_ADD_COMPONENT (headers 
      DISPLAY_NAME "SELALIB Headers" 
      DEPENDS libraries
      GROUP Development
  )
  CPACK_ADD_COMPONENT (hdfdocuments 
      DISPLAY_NAME "SELALIB Documents"
      GROUP Documents
  )
  CPACK_ADD_COMPONENT (configinstall 
      DISPLAY_NAME "SELALIB CMake files" 
      DEPENDS libraries
      GROUP Development
  )

ENDIF(NOT SELALIB_NO_PACKAGES)
