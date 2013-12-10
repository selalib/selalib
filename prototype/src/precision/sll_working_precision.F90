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
!
! MODULE: sll_working_precision
!
!> @author
!> unknown author
!
! DESCRIPTION: 
!> Define the kind type parameter for intern type data.
!> @brief
!> If you want to select the kind parameter \a n, you need to write \n
!> \code
!> real(kind=n) :: var1
!> real*n       :: var2 
!> \endcode
!> The two entries \a var1 and \a var2 are equivalents.
!> You can also define the constant like this \a 23.455_n.
!>
!> <b> How to use this module: </b> \n
!> ****************************
!>
!> First, call the module \a sll_woring_precision like that
!> \code #include "sll_working_precision.h" \endcode
!> Now, you can use the types:
!> \code
!> sll_int32  :: i        !integer simple precision
!> sll_int64  :: N        !integer double precision
!> sll_real32 :: theta    !real simple precision
!> sll_real64 :: my_pi    !real double precision
!>
!> my_pi = 3.1415926535897932384626433_f64
!> theta = 2.0*real(N,f32)
!> \endcode
!
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------

! From the Fortran Standard (2.4.1.1): "The kind type parameter indicates the 
! decimal exponent range for the integer type (4.4.1), the decimal precision 
! and exponent range for the real and complex types (4.4.2, 4.4.3), and the 
! representation methods for the character and logical types (4.4.4, 4.4.5)."
module sll_working_precision
  implicit none
  intrinsic :: kind, selected_real_kind

  ! The intent is that i32 will hold values up to 2**32-1
  !> i32 is the kind type for 32-bit integers
  integer, parameter :: i32 = kind(0)
  !> i64 is the kind type for 64-bit integers
  integer, parameter :: i64 = kind(2_8**32) !i64=kind(1.0d0) should be specific enough
  !> f32 is the kind type for 32-bit reals (simple precision)
  integer, parameter :: f32 = selected_real_kind(1,37)
  !> f64 is the kind type for 64-bit reals (double precision)
  integer, parameter :: f64 = selected_real_kind(1,99)

end module sll_working_precision
