!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE: sll_tridiagonal
!
! DESCRIPTION:
!> @file tridiagonal.F90
!> @namespace sll_tridiagonal
!> @author Module Author Name and Affiliation
!> @brief Tridiagonal system solver.
!> @details To solve systems of the form Ax=b, where A is a tridiagonal matrix, Selalib 
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
!
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------ 
module sll_tridiagonal
#include "sll_working_precision.h"
implicit none

  interface solve_cyclic_tridiag
     module procedure solve_cyclic_tridiag_double, solve_cyclic_tridiag_complex
  end interface
contains

! careful with side-effects here
#define SWP(aval,bval) swp=(aval); aval=(bval); bval=swp

  
  !---------------------------------------------------------------------------  
  !
  ! Implementation notes:
  ! **********************
  ! (Adapted) description for the C implementation:
  !
  ! setup_cyclic_tridiag computes the LU factorization of a cylic
  !  tridiagonal matrix specified in a band-diagonal representation.
  !  solve_cyclic_tridiag uses this factorization to solve of the
  !  system to solve the system Ax = b quickly and robustly.
  !
  !  Unambigously, the input tridiagonal system is specified by:
  !
  !       + a(2)    a(3)                                      a(1) +
  !       | a(4)    a(5)    a(6)                                   |
  !   A = |         a(7)    a(8)    a(9)                           |
  !       |                         ...                            |
  !       |                                a(3n-5) a(3n-4) a(3n-3) |
  !       + a(3n)                                  a(3n-2) a(3n-1) +
  ! 
  ! The LU factorization does partial (row) pivoting, so the
  ! factorization requires slightly more memory to hold than standard
  ! non-pivoting tridiagonal (or cyclic tridiagonal) solution methods
  ! The benefit is that this routine can accomodate systems that are
  ! not diagonally dominant.
  !
  ! The output factorization is stored in "cts" and "ipiv" (the pivot
  ! information).  In general, all you need to know about "cts" is that 
  ! it is what you give to solve_cyclic_tridiag. However, for the 
  ! masochistic, the final factorization is stored in seven vectors 
  ! of length n which are packed into the vector cts in the order: 
  ! d u v q r l m. The L and U factors of A are built out of vectors 
  ! and have the structure:
  !
  !       + 1                               +
  !       | l0  1                           |
  !       |     l1  1                       |
  !       |         ...                     |
  !   L = |             ...                 |
  !       |                 ln4 1           |
  !       |                     ln3 1       |
  !       + m0  m1  m2 ...  mn4 mn3 ln2 1   +
  !
  ! and:
  !
  !       + d0  u0  v0                          q0  r0  +
  !       |     d1  u1  v1                      q1  r1  |
  !       |         d2  u2  v2                  q2  r2  |
  !       |             ...                     :   :   |
  !   U = |                 ...                 :   :   |
  !       |                     dn6 un6 vn6     qn6 rn6 |
  !       |                         dn5 un5 vn5 qn5 rn5 |
  !       |                             dn4 un4 vn4 rn4 |
  !       |                                 dn3 un3 vn3 |
  !       |                                     dn2 un2 |
  !       +                                         dn1 +
  !
  ! Such that LU = PA. (The vector p describes the pivoting matrix.
  ! p[i] = j indicates that in step i of the factorization, row i was
  ! swapped with row j.  See Golub & van Loan. Matrix Computations.
  ! Section 3.4.3 for more details.)  For the convenience of other
  ! functions that use the output, the elements of d,l,m,u,v,q,r not
  ! shown explicitly above are set to zero.
  !
  ! Note that if a zero pivot is encountered (i.e. the det(A) is
  ! indistinguishable from zero as far as the computer can tell),
  ! this routine returns an error."
  !
  ! *************************************************************************

!---------------------------------------------------------------------------  
!> @brief Give the factorization of the matrix in argument.
!> @details setup_cyclic_tridiag has been adapted from the C version written
!>          by Kevin Bowers for the Desmond molecular dynamics code.
!>          For the Fortran implementation, we have adjusted the algorithm
!>          such that it is compatible with the 1-based array indexing:
!>
!> @param a the matrix to be factorized
!> @param[in] n the problem size (A is nXn)     
!> @param[out] ipiv an ineteger array of length n on wich pivoting information will be returned
!> @param[out] cts a real array of size 7n where factorization information will be returned 
subroutine setup_cyclic_tridiag( a, n, cts, ipiv )
  intrinsic :: abs
  sll_real64, dimension(:) :: a    ! 3*n size allocation
  sll_int32,  intent(in)   :: n    ! a is nXn
  sll_int32,  intent(out), dimension(1:n)           :: ipiv
  sll_real64, intent(out), dimension(1:7*n), target :: cts  ! 7*n allocation

  ! The following variables represent a scratch space where the local
  ! computations are made.
  sll_real64 :: s11, s12, s13, s14, s15
  sll_real64 :: s21, s22, s23, s24, s25
  sll_real64 :: s31, s32, s33, s34, s35
  sll_real64 :: t1, t2, t3, swp

  sll_real64, pointer :: d(:)
  sll_real64, pointer :: u(:)
  sll_real64, pointer :: v(:)
  sll_real64, pointer :: q(:)
  sll_real64, pointer :: r(:)
  sll_real64, pointer :: l(:)
  sll_real64, pointer :: m(:)
  sll_int32 :: i
  
  s11 = 0._f64
  s12 = 0._f64
  s13 = 0._f64
  s14 = 0._f64
  s15 = 0._f64

  s21 = 0._f64
  s22 = 0._f64
  s23 = 0._f64
  s24 = 0._f64
  s25 = 0._f64

  s31 = 0._f64
  s32 = 0._f64
  s33 = 0._f64
  s34 = 0._f64
  s35 = 0._f64
  
  
  cts(1:7*n) = 0._f64
  d => cts(    1:n  )
  u => cts(  n+1:2*n)
  v => cts(2*n+1:3*n)
  q => cts(3*n+1:4*n)
  r => cts(4*n+1:5*n)
  l => cts(5*n+1:6*n)
  m => cts(6*n+1:7*n)

  ! The following indexing allows us to dereference the array a as if it
  ! were a rank-2 array. Note that a(1,0) = a(1), i.e.: the cyclic term 
  ! in the first equation and that a(n,n+1) = a(n), i.e.: the cyclic term
  ! in the last equation.
#define aa(ii,jj) a(2*((ii)-1)+(jj)+1)
  select case (n)
     case (3) ! we just reproduce the input array on the scratch space
        s11 = aa(1,1); s12 = aa(1,2); s13 = aa(1,0)
        s21 = aa(2,1); s22 = aa(2,2); s23 = aa(2,3)
        s31 = aa(3,4); s32 = aa(3,2); s33 = aa(3,3)
     case (4) 
        s11 = aa(1,1); s12 = aa(1,2); s13 = 0.0;     s14 = aa(1,0)
        s21 = aa(2,1); s22 = aa(2,2); s23 = aa(2,3); s24 = 0.0
        s31 = aa(4,5); s32 = 0.0;     s33 = aa(4,3); s34 = aa(4,4)
     case default
        s11 = aa(1,1);  s12 = aa(1,2);s13 = 0.0;    s14 = 0.0;     s15 = aa(1,0)
        s21 = aa(2,1);  s22 = aa(2,2);s23 = aa(2,3);s24 = 0.0;     s25 = 0.0
        s31 = aa(n,n+1);s32 = 0.0;    s33 = 0.0;    s34 = aa(n,n-1);s35= aa(n,n)
     end select

     ! Factor the matrix with row pivoting
     do i=1,n
        if(i .lt. n-3) then
           ! Matrix thus far (zeros not shown)
           ! 
           !   |  1   2   3   4   5 ... i-1   i   i+1   i+2  ...   n-1    n
           ! --+-------------------------------------------------------------
           ! 1 | d1  u1  v1                                         q1   r1
           ! 2 | l1  d2  u2  v2                                     q2   r2
           ! 3 |     l2  d3  u3  v4                                 q3   r3
           ! : |                    ...                              :    :
           ! : |                        ...                          :    :
           !i-1|                             di1  ui1   vi1         qi   ri1
           ! i |                             li1  s11   s12  s13... s14  s15
           !i+1|                                  s21   s22  s23... s24  s25
           !   |
           !   |                     Untouched rows of A
           !   |
           ! n | m1  m2  m3  m4  m5  ...     mi1  s31   s32  s33 ...s34  s35

           ! Pivot

           t1 = abs( s11 )
           t2 = abs( s21 )
           t3 = abs( s31 )
           if( (t1 < t2) .or. (t1 < t3) ) then
              if( t3 > t2 ) then ! swap rows i and n
                 SWP(s11,s31) 
                 SWP(s12,s32) 
                 SWP(s13,s33)
                 SWP(s14,s34)
                 SWP(s15,s35)
                 ipiv(i) = n
              else ! swap rows i and i+1
                 SWP(s11,s21) 
                 SWP(s12,s22) 
                 SWP(s13,s23)
                 SWP(s14,s24)
                 SWP(s15,s25)
                 ipiv(i) = i+1
              end if
           else ! no swapping necessary
              ipiv(i) = i
           end if
           ! Eliminate the current column of A below the diagonal. The 
           ! column of L will be stored in place of the zeros.
           if( s11 == 0.0 ) print *, 'zero determinant' ! FIX THIS
           s21 = s21/s11
           s22 = s22 - s21*s12
           s23 = s23 - s21*s13
           s24 = s24 - s21*s14
           s25 = s25 - s21*s15
           s31 = s31/s11
           s32 = s32 - s31*s12
           s33 = s33 - s31*s13
           s34 = s34 - s31*s14
           s35 = s35 - s31*s15

           ! Save the column of L an row of U

           d(i) = s11; u(i) = s12; v(i) = s13; q(i) = s14; r(i) = s15
           l(i) = s21
           m(i) = s31
           
           ! Advance the scratch space for the next iteration

           if( i<(n-4) ) then
              s11=s22;        s12=s23;        s13=0.0;        s14=s24; s15=s25
              s21=aa(i+2,i+1);s22=aa(i+2,i+2);s23=aa(i+2,i+3);s24=0.0; s25=0.0
              s31=s32;        s32=s33;        s33=0.0;        s34=s34; s35=s35
           else
              s11=s22;        s12=s23;        s13=s24;        s14=s25;
              s21=aa(i+2,i+1);s22=aa(i+2,i+2);s23=aa(i+2,i+3);s24=0.0
              s31=s32;        s32=s33;        s33=s34;        s34=s35
           end if
        else if( i==(n-3) ) then
           ! Matrix so far (zeros not shown):
           !
           !         | 1   2   3   4   5  ... n-5 n-4 n-3 n-2 n-1  n
           !     ----+---------------------------------------------------
           !     1   | d1  u1  v1                             q1  r1
           !     2   | l1  d2  u2  v2                         q2  r2
           !     3   |     l2  d3  u3  v3                     q3  r3
           !     :   |                ...                     :   :
           !     :   |                    ...                 :   :
           !     n-5 |                        dn5 un5 vn5     qn5 rn5
           !     n-4 |                        ln5 dn4 un4 vn4 qn4 rn4
           !     n-3 |                            ln4 s11 s12 s13 s14
           !     n-2 |                                s21 s22 s23 s24
           !     n-1 |             Untouched row of A
           !     n   | m1  m2  m3  m4  m5 ... mn5 mn4 s31 s32 s33 s34 
           !
           ! Pivot
           t1 = abs( s11 )
           t2 = abs( s21 )
           t3 = abs( s31 )
           if( (t1 < t2) .or. (t1 < t3) ) then
              if( t3 > t2 ) then ! swap rows i and n
                 SWP(s11,s31) 
                 SWP(s12,s32) 
                 SWP(s13,s33)
                 SWP(s14,s34)
                 ipiv(i) = n
              else ! swap rows i and i+1
                 SWP(s11,s21) 
                 SWP(s12,s22) 
                 SWP(s13,s23)
                 SWP(s14,s24)
                 ipiv(i) = i+1
              end if
           else
              ipiv(i) = i
           end if

           ! Eliminate the current column of A below the diagonal. The column
           ! of L will be stored in place of the zeros.

           if( s11==0.0 ) print *, 'Zero determinant' ! FIX: Do something else!
           s21 = s21/s11
           s22 = s22 - s21*s12
           s23 = s23 - s21*s13
           s24 = s24 - s21*s14
           s31 = s31/s11
           s32 = s32 - s31*s12
           s33 = s33 - s31*s13
           s34 = s34 - s31*s14
           
           ! Save the column of L and row of U
           d(i) = s11; u(i) = s12; v(i) = s13; r(i) = s14
           l(i) = s21
           m(i) = s31
           q(i) = 0.0
           
           ! Advance the scratch space for the next iteration
           s11 = s22;        s12 = s23;         s13 = s24
           s21 = aa(i+2,i+1);s22 = aa(i+2,i+2); s23 = aa(i+2,i+3)
           s31 = s32;        s32 = s33;         s33 = s34
        else if( i==n-2 ) then
           ! Matrix so far (zeros not shown):
           !
           !        | 1   2   3   4   5  ... n-5 n-4 n-3 n-2 n-1 n
           !    ----+---------------------------------------------------
           !    1   | d1  u1  v1                             q1  r1
           !    2   | l1  d2  u2  v2                         q2  r2
           !    3   |     l2  d3  u3  v3                     q3  r3
           !    :   |                ...                     :   :
           !    :   |                    ...                 :   :
           !    n-5 |                        dn5 un5 vn5     qn5 rn5
           !    n-4 |                        ln5 dn4 un4 vn4 qn4 rn4
           !    n-3 |                            ln4 dn3 un3 vn3 rn3
           !    n-2 |                                ln3 s11 s12 s13
           !    n-1 |                                    s21 s22 s23
           !    n   | m1  m2  m3  m4  m5 ... mn5 mn4 mn3 s31 s32 s33 
           !
           !       Pivot 
           t1 = abs( s11 )
           t2 = abs( s21 )
           t3 = abs( s31 )
           if( (t1 < t2) .or. (t1 < t3) ) then
              if( t3 > t2 ) then ! swap rows i and n
                 SWP(s11,s31) 
                 SWP(s12,s32) 
                 SWP(s13,s33)
                 ipiv(i) = n
              else ! swap rows i and i+1
                 SWP(s11,s21) 
                 SWP(s12,s22) 
                 SWP(s13,s23)
                 ipiv(i) = i+1
              end if
           else
              ipiv(i) = i
           end if
           ! Eliminate the current column of A below the diagonal. The
           ! column of L will be sotred inplace of the zeros.
           if( s11==0.0 )  print *, 'zero determinant' ! FIX THIS
           s21 = s21/s11
           s22 = s22 - s21*s12
           s23 = s23 - s21*s13
           s31 = s31/s11
           s32 = s32 - s31*s12
           s33 = s33 - s31*s13
           
           ! Save the column of L and row of U
           d(i) = s11; u(i) = s12; v(i) = s13 
           l(i) = s21
           m(i) = s31
           
           q(i) = 0.0
           r(i) = 0.0
           
           ! Advance the scratch space for the next iteration
           s11 = s22;        s12 = s23
           s21 = s32;        s22 = s33
        else if( i==n-1 ) then

           ! Matrix so far (zeros not shown):
           !
           !     | 1   2   3   4   5  ... n-5 n-4 n-3 n-2 n-1 n
           ! ----+---------------------------------------------------
           ! 1   | d1  u1  v1                             q1  r1
           ! 2   | l1  d2  u2  v2                         q2  r2
           ! 3   |     l2  d3  u3  v3                     q3  r3
           ! :   |                ...                     :   :
           ! :   |                    ...                 :   :
           ! n-5 |                        dn5 un5 vn5     qn5 rn5
           ! n-4 |                        ln5 dn4 un4 vn4 qn4 rn4
           ! n-3 |                            ln4 dn3 un3 vn3 rn3
           ! n-2 |                                ln3 dn2 un2 vn2
           ! n-1 |                                    ln2 s11 s12
           ! n   | m1  m2  m3  m4  m5 ... mn5 mn4 mn3 mn2 s21 s22 
           !
           !    Pivot
           t1 = abs( s11 )
           t2 = abs( s21 )
           if( t1 < t2 ) then 
              SWP(s11,s21) 
              SWP(s12,s22) 
              ipiv(i) = i+1
           else ! swap rows i and i+1
              ipiv(i) = i
           end if
           
           ! Eliminate the current column of A below the diagonal.
           ! The column of L will be stored in place of the zeros.
           if( s11==0.0 ) print *, 'Zero determinant' ! FIX THIS
           s21 = s21/s11
           s22 = s22 - s21*s12

           ! Save the column of L and row of U

           d(i) = s11; u(i) = s12
           l(i) = s21

           v(i) = 0.0
           r(i) = 0.0
           q(i) = 0.0
           m(i) = 0.0

           ! Advance the scratch space for the next iteration

           s11 = s22

        else ! i==n
           ! Matrix so far (zeros not shown):
           !
           !       | 1   2   3   4   5  ... n-5 n-4 n-3 n-2 n-1 n
           !   ----+---------------------------------------------------
           !   1   | d1  u1  v1                             q1  r1
           !   2   | l1  d2  u2  v2                         q1  r1
           !   3   |     l2  d3  u3  v3                     q3  r3
           !   :   |                ...                     :   :
           !   :   |                    ...                 :   :
           !   n-5 |                        dn5 un5 vn5     qn5 rn5
           !   n-4 |                        ln5 dn4 un4 vn4 qn4 rn4
           !   n-3 |                            ln4 dn3 un3 vn3 rn3
           !   n-2 |                                ln3 dn2 un2 vn2
           !   n-1 |                                    ln2 dn1 un1
           !   n   | m1  m2  m3  m4  m5 ... mn5 mn4 mn3 mn2 ln1 s11 

           if( s11 == 0.0 ) print *, 'zero determinant' ! FIX THIS
           ipiv(i) = i
           d(i) = s11
           u(i) = 0.0
           v(i) = 0.0
           r(i) = 0.0
           q(i) = 0.0
           l(i) = 0.0
           m(i) = 0.0
        end if
     end do

   end subroutine setup_cyclic_tridiag
#undef aa

   !---------------------------------------------------------------------------  
   !> @author Routine Author Name and Affiliation.
   !> @brief Solves tridiagonal system. 
   !> @details Computes the solution of:
   !>
   !>          <center> <b> A x = b </b> </center>
   !> 
   !>          For a cyclic tridiagonal matrix A. The matrix cts is filled with
   !>          the output of the function setup_cyclic_tridiag.  Note that the
   !>          call:
   !>
   !>          solve_cyclic_tridiag( cts, ipiv, b, n, b )
   !>
   !>          is valid if you want run in-place and overwrite the right hand side
   !>          with the solution.
   !>
   !> @param cts a real array of size 7n where factorization information will be returned
   !> @param[in] ipiv an ineteger array of length n on wich pivoting information will be returned
   !> @param b the second member of the equation
   !> @param n the problem size
   !> @param x the solution vector
   !--------------------------------------------------------------------------- 
 
   subroutine solve_cyclic_tridiag_double( cts, ipiv, b, n, x )
     ! size of the allocations:
     ! x:     n
     ! cts:  7n
     ! ipiv:  n
     ! b:     n

     sll_int32,  intent(in)                 :: n    ! matrix size
     sll_real64, dimension(1:7*n), target   :: cts  ! 7*n size allocation
     sll_int32,  dimension(1:n), intent(in) :: ipiv
     sll_real64, target                     :: b(n)
     sll_real64, target                     :: x(n)  
     sll_real64, pointer, dimension(:)                    :: bptr
     sll_real64, pointer, dimension(:)                    :: xptr  
     sll_real64                             :: swp
     sll_int32                              :: i
     sll_int32                              :: inew
     sll_real64, pointer                    :: d(:)
     sll_real64, pointer                    :: u(:)
     sll_real64, pointer                    :: v(:)
     sll_real64, pointer                    :: q(:)
     sll_real64, pointer                    :: r(:)
     sll_real64, pointer                    :: l(:)
     sll_real64, pointer                    :: m(:)

     d => cts(    1:n  )
     u => cts(  n+1:2*n)
     v => cts(2*n+1:3*n)
     q => cts(3*n+1:4*n)
     r => cts(4*n+1:5*n)
     l => cts(5*n+1:6*n)
     m => cts(6*n+1:7*n)

     !do i=1,n
     !  print *,'duvqrlm=',i,d(i),u(i),v(i),q(i),r(i),l(i),m(i) 
     !enddo     

     bptr =>b(1:n)
     xptr =>x(1:n)
     ! FIX: ADD SOME ERROR CHECKING ON ARGUMENTS
     if( .not. associated(xptr, target=bptr) ) then
        do i=1,n
           x(i) = b(i)
        end do
     end if
     ! 'x' contains now the informatin in 'b', in case that it was given as 
     ! a different array.
     !
     ! Overwrite x with the solution of Ly = Pb
     do i=1,n-1
        SWP(x(i),x(ipiv(i)))
        x(i+1) = x(i+1) - l(i)*x(i)
        x(n)   = x(n) - m(i)*x(i)
     end do

     ! Overwrite x with the solution of Ux = y
     i    = n
     x(i) = x(i)/d(i)
     i    = i-1
     x(i) = (x(i) - u(i)*x(i+1))/d(i)
     !i    = i-1
     inew    = i-1
     do i=inew,1,-1
     !do i=i,1,-1
        x(i) = (x(i)-(u(i)*x(i+1) + v(i)*x(i+2) + &
                      q(i)*x(n-1) + r(i)*x(n) ))/d(i)
     end do
   end subroutine solve_cyclic_tridiag_double

   subroutine solve_cyclic_tridiag_complex( cts, ipiv, b, n, x )
     ! size of the allocations:
     ! x:     n
     ! cts:  7n
     ! ipiv:  n
     ! b:     n

     sll_int32,  intent(in)                 :: n    ! matrix size
     sll_real64, dimension(1:7*n), target   :: cts  ! 7*n size allocation
     sll_int32,  dimension(1:n), intent(in) :: ipiv
     sll_comp64, target                     :: b(n)
     sll_comp64, target                     :: x(n)  
     sll_comp64, pointer, dimension(:)      :: bptr
     sll_comp64, pointer, dimension(:)      :: xptr  
     sll_comp64                             :: swp
     sll_int32                              :: i,inew
     sll_real64, pointer                    :: d(:)
     sll_real64, pointer                    :: u(:)
     sll_real64, pointer                    :: v(:)
     sll_real64, pointer                    :: q(:)
     sll_real64, pointer                    :: r(:)
     sll_real64, pointer                    :: l(:)
     sll_real64, pointer                    :: m(:)

     d => cts(    1:n  )
     u => cts(  n+1:2*n)
     v => cts(2*n+1:3*n)
     q => cts(3*n+1:4*n)
     r => cts(4*n+1:5*n)
     l => cts(5*n+1:6*n)
     m => cts(6*n+1:7*n)

     bptr =>b(1:n)
     xptr =>x(1:n)
     
     ! FIX: ADD SOME ERROR CHECKING ON ARGUMENTS
     if( .not. associated(xptr, target=bptr) ) then
        do i=1,n
           x(i) = b(i)
        end do
     end if
     ! 'x' contains now the informatin in 'b', in case that it was given as 
     ! a different array.
     !
     ! Overwrite x with the solution of Ly = Pb
     do i=1,n-1
        SWP(x(i),x(ipiv(i)))
        x(i+1) = x(i+1) - l(i)*x(i)
        x(n)   = x(n) - m(i)*x(i)
     end do

     ! Overwrite x with the solution of Ux = y
     i    = n
     x(i) = x(i)/d(i)
     i    = i-1
     x(i) = (x(i) - u(i)*x(i+1))/d(i)
     inew    = i-1
     do i=inew,1,-1
        x(i) = (x(i)-(u(i)*x(i+1) + v(i)*x(i+2) + &
                      q(i)*x(n-1) + r(i)*x(n) ))/d(i)
     end do
   end subroutine solve_cyclic_tridiag_complex

   ! @brief Solves the system ax=b
   ! param[in] a Global matrix
   ! param[in] b Second member
   ! param[in] n Problem size (a is nxn)
   ! param[out] x Solution vector
   !SUBROUTINE solve_tridiag(a, b, n ,x)
   ! sll_int32, intent(in)                   :: n
   ! sll_real64, DIMENSION(n), intent(in)    :: a
   ! sll_real64, DIMENSION(n), intent(inout) :: b
   ! sll_real64, DIMENSION(1:7*n)            :: cts 
   ! sll_int32,  DIMENSION(1:n)              :: ipiv
   ! sll_real64, DIMENSION(:)                :: x
   ! CALL setup_cyclic_tridiag( a, n, cts, ipiv )
   ! CALL solve_cyclic_tridiag( cts, ipiv, b, n, x )
   !END SUBROUTINE solve_tridiag

end module sll_tridiagonal
