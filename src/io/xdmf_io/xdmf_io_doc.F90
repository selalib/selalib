!> @defgroup xdmf sll_xdmf
!!
!! @authors Yaman Güçlü    - <yaman.guclu@gmail.com>
!! @authors Marco Restelli - <marco.restelli@gmail.com>
!!
!! @brief 
!! Library for writing simulation output data in the form of an XDMF database.
!! 
!! @details
!! Typically, an XDMF database is composed of a sequence of time dependent
!! results. At each snapshot in time, one lightweight XML file is produced,
!! together with a heavier HDF5 file.
!! See <a href= http://www.xdmf.org/index.php/Main_Page target=_blank>
!! here </a>.
!!
!! @note
!! This library is work-in progress. Only the XML module is currently stable.
!!
!! <h4> Header files available </h4>
!!   + No header files.
!! <!--
!!   + *sll_m_xdmf.h*
!! -->
!!
!! <h4> How to use it </h4>
!! <!--
!!  + Include header file:
!!    \code #include "sll_vector_space.h" \endcode
!! -->
!!  + Import Fortran modules:
!!    \code
!!      use sll_m_xml
!!    \endcode
!!  + Add dependency to *CMakeLists.txt*:
!!    \code
!!      target_link_libraries( <my_lib/exec> sll_m_xdmf ... )
!!    \endcode
!!
!! <h4> Examples </h4>
!!   TODO: Add some fortran lines to explain how to use the library
!! \code
!! \endcode
!!

