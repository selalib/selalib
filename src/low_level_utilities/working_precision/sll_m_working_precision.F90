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

! From the Fortran Standard (2.4.1.1): "The kind type parameter indicates the 
! decimal exponent range for the integer type (4.4.1), the decimal precision 
! and exponent range for the real and complex types (4.4.2, 4.4.3), and the 
! representation methods for the character and logical types (4.4.4, 4.4.5)."

!> @ingroup working_precision
!> @brief
!> Module to select the kind parameter. 
!> @details
!> In future developement we consider to put here some ISO_C bindings
!> to call selalib from python.
module sll_m_working_precision

  implicit none
!  intrinsic :: kind, selected_real_kind
  intrinsic :: selected_int_kind, selected_real_kind

  public :: i32, i64, f32, f64
  private

  ! The intent is that i32 will hold values up to 2**32-1
  !> i32 is the kind type for 32-bit integers
!  integer, parameter :: i32 = kind(0)
  integer, parameter :: i32 = selected_int_kind(9)
  !> i64 is the kind type for 64-bit integers
!  integer, parameter :: i64 = kind(2_8**32) !i64=kind(1.0d0) should be specific enough
  integer, parameter :: i64 = selected_int_kind(18)
  !> f32 is the kind type for 32-bit reals (simple precision)
  integer, parameter :: f32 = selected_real_kind(1,37)
  !> f64 is the kind type for 64-bit reals (double precision)
  integer, parameter :: f64 = selected_real_kind(1,99)

end module sll_m_working_precision
