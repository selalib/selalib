! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory
! in file Doxyfile.in (line 691) if it is excluded.
! Type 'make doc' in build directory.
! To check the results, open :
! selalib/doc/build/html/doxygen/html/defgroup characteristicss.html
! The following lines will be read by doxygen to generate documentation:

!> @defgroup characteristics sll_characteristics
!> @brief
!> computes the characteristics for advection
!> @author Michel Mehrenberger
!> You can add a contact, do not put your email to prevent spam.
!> @details
!> computes the characteristics
!> it is then used for the advection
!> in combination with interpolation for the BSL method
!> there exists conservative versions
!> that are used with CSL (more in development)

!> <b> Headers file available </b>
!>  - sll_characteristics.h
!>
!> <b> Modules available </b>
!>  List fortran module available
!>  - sll_characteristics
!>
!> <b> How to use it </b>
!> - Header file : \code #include 'sll_characteristics.h' \endcode
!> - Link with   <code>-lsll_%s</code>
!> - Add <code> use sll_characteristics </code>
!>
!> <b> Examples </b>
!> -Add some fortran lines to explain how ti use the library
!> \code
!> call initialize(my_type, arg_1, arg_2)
!> call solve(my_type, your_result)
!> \endcode
!>
