# SLL_INSTALL_LIB_DIR, SLL_INSTALL_INCLUDE_DIR :
#   Customize the 'lib', and 'include' installation directories.
#
OPTION (SLL_NO_PACKAGES "CPACK - Disable packaging" OFF)
MARK_AS_ADVANCED (SLL_NO_PACKAGES)

#-----------------------------------------------------------------------------
# Set the core names of all the libraries
#-----------------------------------------------------------------------------
SET (SLL_LIB_CORENAME             "selalib")
SET (SLL_MEMORY_LIB_CORENAME      "sll_memory")
SET (SLL_ASSERT_LIB_CORENAME      "sll_assert")
SET (SLL_WP_LIB_CORENAME          "sll_working_precision")
SET (SLL_CONSTANTS_LIB_CORENAME   "sll_constants")
SET (SLL_UTLITIES_LIB_CORENAME    "sll_utilities")

IF (NOT SLL_INSTALL_LIB_DIR)
  SET (SLL_INSTALL_LIB_DIR lib)
ENDIF (NOT SLL_INSTALL_LIB_DIR)

IF (NOT SLL_INSTALL_INCLUDE_DIR)
  SET (SLL_INSTALL_INCLUDE_DIR include)
ENDIF (NOT SLL_INSTALL_INCLUDE_DIR)

SET (SLL_PACKAGE "selalib")
SET (SLL_PACKAGE_NAME "SLL")
SET (SLL_PACKAGE_VERSION "1")
SET (SLL_PACKAGE_VERSION_MINOR "0")
SET (SLL_PACKAGE_VERSION_MAJOR "0")
SET (SLL_PACKAGE_URL "http://selalib.gforge.inria.fr")
SET (SLL_PACKAGE_BUGREPORT "selalib-dev@gforge.inria.fr")

#-----------------------------------------------------------------------------
# Set the cpack variables
#-----------------------------------------------------------------------------
IF (NOT SLL_NO_PACKAGES)

  SET (CPACK_PACKAGE_VENDOR "IRMA")
  SET (CPACK_PACKAGE_NAME "${SLL_PACKAGE_NAME}")
  SET (CPACK_PACKAGE_INSTALL_DIRECTORY "${SLL_PACKAGE_NAME}")
  SET (CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${SLL_PACKAGE_NAME}-${SLL_PACKAGE_VERSION}-${LIB_TYPE}")
  SET (CPACK_PACKAGE_VERSION "${SLL_PACKAGE_VERSION}")
  SET (CPACK_PACKAGE_VERSION_MAJOR "${SLL_PACKAGE_VERSION_MAJOR}")
  SET (CPACK_PACKAGE_VERSION_MINOR "${SLL_PACKAGE_VERSION_MINOR}")
  SET (CPACK_PACKAGE_VERSION_PATCH "")
  SET (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/RELEASE.txt")
  SET (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING.txt")
  SET (CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
  SET (CPACK_PACKAGE_RELOCATABLE TRUE)

  SET (CPACK_PACKAGING_INSTALL_PREFIX "/usr")
  SET (CPACK_COMPONENTS_ALL_IN_ONE_PACKAGE ON)

  SET (CPACK_DEBIAN_PACKAGE_SECTION "Libraries")
  SET (CPACK_DEBIAN_PACKAGE_MAINTAINER "${SLL_PACKAGE_BUGREPORT}")
    
  SET (CPACK_RPM_COMPONENT_INSTALL ON)
  SET (CPACK_RPM_PACKAGE_RELOCATABLE ON)
  SET (CPACK_RPM_PACKAGE_LICENSE "BSD-style")
  SET (CPACK_RPM_PACKAGE_GROUP "Development/Libraries")
  SET (CPACK_RPM_PACKAGE_URL "${SLL_PACKAGE_URL}")
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

  SET (CPACK_INSTALL_CMAKE_PROJECTS "${SLL_BINARY_DIR};SLL;ALL;/")
  
  INCLUDE (CPack)

  #---------------------------------------------------------------------------
  # Now list the cpack commands
  #---------------------------------------------------------------------------
  CPACK_ADD_COMPONENT (libraries 
      DISPLAY_NAME "SLL Libraries"
      GROUP Runtime
  )
  CPACK_ADD_COMPONENT (headers 
      DISPLAY_NAME "SLL Headers" 
      DEPENDS libraries
      GROUP Development
  )
  CPACK_ADD_COMPONENT (hdfdocuments 
      DISPLAY_NAME "SLL Documents"
      GROUP Documents
  )
  CPACK_ADD_COMPONENT (configinstall 
      DISPLAY_NAME "SLL CMake files" 
      DEPENDS libraries
      GROUP Development
  )

ENDIF(NOT SLL_NO_PACKAGES)
