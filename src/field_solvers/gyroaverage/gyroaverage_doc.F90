! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory
! in file Doxyfile.in (line 691) if it is excluded.
! Type 'make doc' in build directory.
! To check the results, open :
! selalib/doc/build/html/doxygen/html/defgroup gyroaverages.html
! The following lines will be read by doxygen to generate documentation:

!> @defgroup gyroaverage sll_gyroaverage
!> @brief
!> Compute the gyroaverage operator in a polar mesh.
!> @author Selalib team
!> Michel Mehrenberger
!> Christophe Steiner
!> @details
!> Three methods for computing the gyroaverage operator in a 2D polar mesh are implemented ( Hermite, Pade , splines). The description of these methods are given in :
!> Gyroaverage for a polar mesh
!> C. Steiner, M. Mehrenberger, N. Crouseilles, V. Grandgirard, G. Latu, F. Rozar
!> EPJD 2015
!> Two methods for solving the quasi-neutrality equation are implemented ( Pade , circles ). The description of these methods are given in :
!> Résolution numérique de l'opérateur de gyromoyenne, schémas d'advection et couplage. Applications à l'équation de Vlasov.
!> C. Steiner
!> Thesis 2014, Chapter 9
!> <b> Headers file available </b>
!>  - sll_gyroaverage.h
!>
!> <b> Modules available </b>
!>  List fortran module available
!>  - sll_gyroaverage
!>
!> <b> How to use it </b>
!> - Header file : \code #include 'sll_gyroaverage.h' \endcode
!> - Link with   <code>-lsll_%s</code>
!> - Add <code> use sll_gyroaverage </code>
!>
!> <b> Examples </b>
!> -Add some fortran lines to explain how ti use the library
!> \code
!> call initialize(my_type, arg_1, arg_2)
!> call solve(my_type, your_result)
!> \endcode
!>
