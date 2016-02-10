! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory 
! in file Doxyfile.in (line 691) if it is excluded. 
! Type 'make doc' in build directory.
! To check the results, open : 
! selalib/doc/build/html/doxygen/html/defgroup pic_particle_typess.html 
! The following lines will be read by doxygen to generate documentation:


!> @defgroup pic_particle_types sll_pic_particle_types 
!> @brief Contains basic particle and particle group types 
!>        for 2d and 4d PIC simulations
!> @author Sever Hirstoaga
!> You can add a contact, do not put your email to prevent spam.
!> @details  The sll_particle_representation.h file contains Macros
!>           for conversion between standard (position,velocity) 
!>           particle type and the type defined in the module
!>           sll_m_particle_representations.
!> <b> Headers file available </b>
!>  - sll_pic_particle_types.h
!>
!> <b> Modules available </b>
!>  List fortran module available
!>  - sll_pic_particle_types
!>
!> <b> How to use it </b>
!> - Header file : \code #include 'sll_pic_particle_types.h' \endcode
!> - Link with   <code>-lsll_%s</code>
!> - Add <code> use sll_pic_particle_types </code>
!>
!> <b> Examples </b>
!> -Add some fortran lines to explain how ti use the library
!> \code
!> call initialize(my_type, arg_1, arg_2)
!> call solve(my_type, your_result)
!> \endcode
!>
