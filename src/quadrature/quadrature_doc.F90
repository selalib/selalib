! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory 
! in file Doxyfile.in (line 691) if it is excluded. 
! Type 'make doc' in build directory.
! To check the results, open : 
! selalib/doc/build/html/doxygen/html/defgroup integrations.html 
! The following lines will be read by doxygen to generate documentation:


!!         

!> @defgroup quadrature sll_quadrature 
!> @brief 
!> Implementation of various quadrature rules for integration of a function.
!> @authors Edwin Chacon-Golcher, Laura Mendoza and Pierre Navaro.
!! @details 
!> This module aims at providing a single interface to the process of 
!> integrating a function on a given interval.
!!   - Gauss-Legendre points and weights
!!   - Gauss-Lobatto points and weights
!!
!> <b> How to use it </b>
!> - Link with   <code>-lsll_quadrature</code>
!> - Add <code> use sll_m_gauss_legendre_integration </code>
!> - Add <code> use sll_m_gauss_lobatto_integration </code>
!>
!> <b> Examples </b>
!> \code
!> sll_f_gauss_legendre_points_and_weights(5,-1.0_f64,1.0_f64)
!>  x = gauss_lobatto_points( 10, -1._f64, 1._f64)
!>  w = gauss_lobatto_weights(10)
!> \endcode
!>
