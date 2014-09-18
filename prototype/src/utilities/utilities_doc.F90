! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory 
! in file Doxyfile.in (line 691) if it is excluded. 
! Type 'make doc' in build directory.
! To check the results, open : 
! selalib/prototype/documentation/build/html/doxygen/html/namespaces.html 
! The following lines will be read by doxygen to generate documentation:


!> @namespace sll_utilities 
!> @brief 
!> Some numerical low level utilities
!> @details
!> - Cyclic tridiagonal system solver
!> - Pentadiagonal system solver
!> - Skyline format matrix solver
!> - Display functions for matrix and vector
!> - integer to string function
!> @author Selalib team 
!>
!> <b> Headers file available </b>
!>  - sll_utilities.h
!>
!> <b> Modules available </b>
!>  - sll_utilities
!>  - sll_tridiagonal
!>  - sll_toep_penta_diagonal
!>  - skyline
!>
!> <b> How to use it </b>
!> - Header file : \code #include 'sll_utilities.h' \endcode
!> - Link with   <code>-lsll_utilities</code>
!> - Add <code> use sll_utilities </code>
!>
!> <b> Tridiagonal system solver </b>
!> 
!> To solve systems of the form Ax=b, where A is a tridiagonal matrix, Selalib 
!! offers a native, robust tridiagonal system solver. The present implementation 
!! contains only a serial version.
!! The algorith is based on an LU factorisation of a given matrix,
!! with row pivoting. The tridiagonal matrix must be given as a single array,
!! with a memory layout shown next.
!> \f[ A = \begin{bmatrix}
!! a(2) & a(3) & & & & & a(1)
!! \\ a(4) & a(5) & a(6) & & & &
!! \\ & a(7) & a(8) & a(9) & & &
!! \\ & & \ddots & \ddots & \ddots & &
!! \\ & & & \ddots & \ddots & \ddots &
!! \\ & & & & a(3n-5) & a(3n-4)&a(3n-3)
!! \\ a(3n)& & & & & a(3n-2) & a(3n-1)
!! \end{bmatrix} \f]
!>
!> Usage:
!>
!> To solve a tridiagonal system, first: \n
!> -# Assemble the matrix 'a' as a single array with the layout just 
!>    described above
!> -# Use 'setup_cyclic_tridiag( a, n, cts, ipiv )' to factorize the system
!>    -# In 'setup_cyclic_tridag', a is the array to be factorized, stored 
!>        with the layout shown above. 'n' is essentially the problem size.
!>        cts and ipiv are respectively real and integer arrays of size 7*n 
!>        and n that are needed to return factorization information. ipiv
!>        is the usual 'pivot' array.
!> -# To solve the system, make a call to 
!>         'solve_cyclic_tridiag(cts, ipiv, b, n, x)'
!>    Here, cts and ipiv are the ones returned by setup_cyclic_tridiag. The
!>    function returns the solution to Ax = b, storing the results in 'x'.
!>    In case that an 'in-place' computation is desired, it is acceptable to
!>    make the call like: 
!>          solve_cyclic_tridiag(cts, ipiv, b, n, b)
!>
