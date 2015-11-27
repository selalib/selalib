
!> @details Abramovic and Stegun, Chapter 25.2

module sll_m_lagrange_fast
#include "sll_working_precision.h"
  implicit none

  ! --- compile-time constants to avoid run-time division
  sll_real64, parameter :: inv_6       = 1./6.
  sll_real64, parameter :: inv_12      = 1./12.
  sll_real64, parameter :: inv_24      = 1./24.
  sll_real64, parameter :: inv_36      = 1./36.
  sll_real64, parameter :: inv_48      = 1./48.
  sll_real64, parameter :: inv_120     = 1./120.
  sll_real64, parameter :: inv_576     = 1./576.
  sll_real64, parameter :: inv_720     = 1./720.
  sll_real64, parameter :: inv_1440    = 1./1440.
  sll_real64, parameter :: inv_5040    = 1./5040.
  sll_real64, parameter :: inv_14400   = 1./14400.
  sll_real64, parameter :: inv_17280   = 1./17280.
  sll_real64, parameter :: inv_30240   = 1./30240.
  sll_real64, parameter :: inv_40320   = 1./40320.
  sll_real64, parameter :: inv_80640   = 1./80640.
  sll_real64, parameter :: inv_362880  = 1./362880.
  sll_real64, parameter :: inv_3628800 = 1./3628800.

contains

  subroutine lagr_3pt_coeff(pp, p)
    implicit none
    sll_real64, intent(out) :: pp(3) !< Lagrange interpolations coefficients
    sll_real64, intent(in) :: p      !< offset in units of grid spacing
    pp(1) = p*(p-1.)*0.5
    pp(2) = 1. - p*p
    pp(3) = p*(p+1.)*0.5
  end subroutine

  ! --- single point 3-pt-lagrange interpolation
  function lagr_3pt(fm1, f0, f1, p)
    implicit none
    sll_real64 :: lagr_3pt !< interpolated value
    sll_real64, intent(in) :: fm1, f0, f1, p !< known function values at point -1, 0, 1 (relative to where we want to interpolate)
    sll_real64 :: pp(3)
    call lagr_3pt_coeff(pp, p)
    lagr_3pt = pp(1) * fm1 &
      + pp(2) * f0  &
      + pp(3) * f1
  end function lagr_3pt

  ! --- vectorizable 3-pt-lagrange interpolation
  subroutine lagr_3pt_vec(fi, fp, p)
    implicit none
    sll_real64, intent(in) :: fi(:), p
    sll_real64, intent(out) :: fp(:)
    sll_real64 :: pp(3)
    sll_int32 :: i, n
    call lagr_3pt_coeff(pp, p)
    n = size(fi)
    do i=2,n-1
      fp(i) = pp(1) * fi(i-1) &
        + pp(2) * fi(i)   &
        + pp(3) * fi(i+1)
    enddo
  end subroutine lagr_3pt_vec



  subroutine lagr_5pt_coeff(pp, p)
    implicit none
    sll_real64, intent(out) :: pp(5)
    sll_real64, intent(in) :: p
    pp(1) = (p*p-1.)*p*(p-2.)*inv_24
    pp(2) = -(p-1.)*p*(p*p-4.)*inv_6
    pp(3) = (p*p-1.)*(p*p-4.)*0.25
    pp(4) = -(p+1.)*p*(p*p-4.)*inv_6
    pp(5) = (p*p-1.)*p*(p+2.)*inv_24
  end subroutine

  ! --- single point 5-pt-lagrange interpolation
  function lagr_5pt(fm2, fm1, f0, f1, f2, p)
    implicit none
    sll_real64 :: lagr_5pt
    sll_real64, intent(in) :: fm2, fm1, f0, f1, f2, p
    sll_real64 :: pp(5)
    call lagr_5pt_coeff(pp, p)
    lagr_5pt = pp(1) * fm2 &
      + pp(2) * fm1 &
      + pp(3) * f0  &
      + pp(4) * f1  &
      + pp(5) * f2
  end function lagr_5pt

  ! --- vectorizable 5-pt-lagrange interpolation
  subroutine lagr_5pt_vec(fi, fp, p)
    implicit none
    sll_real64, intent(in) :: fi(:), p
    sll_real64, intent(out) :: fp(:)
    sll_real64 :: pp(5)
    sll_int32 :: i, n
    call lagr_5pt_coeff(pp, p)
    n = size(fi)
    do i=3, n-2
      fp(i) = pp(1) * fi(i-2) &
        + pp(2) * fi(i-1) &
        + pp(3) * fi(i)   &
        + pp(4) * fi(i+1) &
        + pp(5) * fi(i+2)
    enddo
  end subroutine lagr_5pt_vec



  subroutine lagr_7pt_coeff(pp, p)
    implicit none
    sll_real64, intent(out) :: pp(7)
    sll_real64, intent(in) :: p
    pp(1) = p*(p-3)*(p**2-4)*(p**2-1)*inv_720
    pp(2) = -p*(p-2)*(p**2-9)*(p**2-1)*inv_120
    pp(3) = p*(p-1)*(p**2-9)*(p**2-4)*inv_48
    pp(4) = -(p**2-9)*(p**2-4)*(p**2-1)*inv_36
    pp(5) = (p+1)*p*(p**2-9)*(p**2-4)*inv_48
    pp(6) = -(p+2)*p*(p**2-9)*(p**2-1)*inv_120
    pp(7) = (p+3)*p*(p**2-4)*(p**2-1)*inv_720
  end subroutine

  ! --- single point 7-pt-lagrange interpolation
  function lagr_7pt(fm3, fm2, fm1, f0, f1, f2, f3, p)
    implicit none
    sll_real64 :: lagr_7pt
    sll_real64, intent(in) :: fm3, fm2, fm1, f0, f1, f2, f3, p
    sll_real64 :: pp(7)
    call lagr_7pt_coeff(pp, p)
    lagr_7pt = pp(1) * fm3 &
      + pp(2) * fm2 &
      + pp(3) * fm1 &
      + pp(4) * f0  &
      + pp(5) * f1  &
      + pp(6) * f2  &
      + pp(7) * f3
  end function lagr_7pt

  ! --- vectorizable 7-pt-lagrange interpolation
  subroutine lagr_7pt_vec(fi, fp, p)
    implicit none
    sll_real64, intent(in) :: fi(:), p
    sll_real64, intent(out) :: fp(:)
    sll_real64 :: pp(7)
    sll_int32 :: i, n
    call lagr_7pt_coeff(pp, p)
    n = size(fi)
    do i=4, n-3
      fp(i) = pp(1) * fi(i-3) &
        + pp(2) * fi(i-2) &
        + pp(3) * fi(i-1) &
        + pp(4) * fi(i)   &
        + pp(5) * fi(i+1) &
        + pp(6) * fi(i+2) &
        + pp(7) * fi(i+3)
    enddo
  end subroutine lagr_7pt_vec




  subroutine lagr_9pt_coeff(pp, p)
    implicit none
    sll_real64, intent(out) :: pp(9)
    sll_real64, intent(in) :: p
    pp(1) = p*(p-4)*(p**2-9)*(p**2-4)*(p**2-1)*inv_40320
    pp(2) = -p*(p-3)*(p**2-16)*(p**2-4)*(p**2-1)*inv_5040
    pp(3) = p*(p-2)*(p**2-16)*(p**2-9)*(p**2-1)*inv_1440
    pp(4) = -p*(p-1)*(p**2-16)*(p**2-9)*(p**2-4)*inv_720
    pp(5) = (p**2-16)*(p**2-9)*(p**2-4)*(p**2-1)*inv_576
    pp(6) = -(p+1)*p*(p**2-16)*(p**2-9)*(p**2-4)*inv_720
    pp(7) = (p+2)*p*(p**2-16)*(p**2-9)*(p**2-1)*inv_1440
    pp(8) = -(p+3)*p*(p**2-16)*(p**2-4)*(p**2-1)*inv_5040
    pp(9) = (p+4)*p*(p**2-9)*(p**2-4)*(p**2-1)*inv_40320
  end subroutine

  ! --- single point 9-pt-lagrange interpolation
  function lagr_9pt(fm4, fm3, fm2, fm1, f0, f1, f2, f3, f4, p)
    implicit none
    sll_real64 :: lagr_9pt
    sll_real64, intent(in) :: fm4, fm3, fm2, fm1, f0, f1, f2, f3, f4, p
    sll_real64 :: pp(9)
    call lagr_9pt_coeff(pp, p)
    lagr_9pt = pp(1) * fm4 &
      + pp(2) * fm3 &
      + pp(3) * fm2 &
      + pp(4) * fm1 &
      + pp(5) * f0  &
      + pp(6) * f1  &
      + pp(7) * f2  &
      + pp(8) * f3  &
      + pp(9) * f4
  end function lagr_9pt

  ! --- vectorizable 9-pt-lagrange interpolation
  subroutine lagr_9pt_vec(fi, fp, p)
    implicit none
    sll_real64, intent(in) :: fi(:), p
    sll_real64, intent(out) :: fp(:)
    sll_real64 :: pp(9)
    sll_int32 :: i, n
    call lagr_9pt_coeff(pp, p)
    n = size(fi)
    do i=5, n-4
      fp(i) = pp(1) * fi(i-4) &
        + pp(2) * fi(i-3) &
        + pp(3) * fi(i-2) &
        + pp(4) * fi(i-1) &
        + pp(5) * fi(i)   &
        + pp(6) * fi(i+1) &
        + pp(7) * fi(i+2) &
        + pp(8) * fi(i+3) &
        + pp(9) * fi(i+4)
    enddo
  end subroutine lagr_9pt_vec




  subroutine lagr_11pt_coeff(pp, p)
    implicit none
    sll_real64, intent(out) :: pp(11)
    sll_real64, intent(in) :: p
    ! generated using Maple
    pp(1)  = p*(p-5)*(p**2-16)*(p**2-9)*(p**2-4)*(p**2-1)*inv_3628800
    pp(2)  = -p*(p-4)*(p**2-25)*(p**2-9)*(p**2-4)*(p**2-1)*inv_362880
    pp(3)  = p*(p-3)*(p**2-25)*(p**2-16)*(p**2-4)*(p**2-1)*inv_80640
    pp(4)  = -p*(p-2)*(p**2-25)*(p**2-16)*(p**2-9)*(p**2-1)*inv_30240
    pp(5)  = p*(p-1)*(p**2-25)*(p**2-16)*(p**2-9)*(p**2-4)*inv_17280
    pp(6)  = -(p**2-25)*(p**2-16)*(p**2-9)*(p**2-4)*(p**2-1)*inv_14400
    pp(7)  = (p+1)*p*(p**2-25)*(p**2-16)*(p**2-9)*(p**2-4)*inv_17280
    pp(8)  = -(p+2)*p*(p**2-25)*(p**2-16)*(p**2-9)*(p**2-1)*inv_30240
    pp(9)  = (p+3)*p*(p**2-25)*(p**2-16)*(p**2-4)*(p**2-1)*inv_80640
    pp(10) = -(p+4)*p*(p**2-25)*(p**2-9)*(p**2-4)*(p**2-1)*inv_362880
    pp(11) = (p+5)*p*(p**2-16)*(p**2-9)*(p**2-4)*(p**2-1)*inv_3628800
  end subroutine

  ! --- single point 11-pt-lagrange interpolation
  function lagr_11pt(fm5, fm4, fm3, fm2, fm1, f0, f1, f2, f3, f4, f5, p)
    implicit none
    sll_real64 :: lagr_11pt
    sll_real64, intent(in) :: fm5, fm4, fm3, fm2, fm1, f0, f1, f2, f3, f4, f5, p
    sll_real64 :: pp(11)
    call lagr_11pt_coeff(pp, p)
    lagr_11pt = pp(1) * fm5 &
      + pp(2) * fm4 &
      + pp(3) * fm3 &
      + pp(4) * fm2 &
      + pp(5) * fm1 &
      + pp(6) * f0  &
      + pp(7) * f1  &
      + pp(8) * f2  &
      + pp(9) * f3  &
      + pp(10)* f4  &
      + pp(11)* f5
  end function lagr_11pt

  ! --- vectorizable 11-pt-lagrange interpolation
  subroutine lagr_11pt_vec(fi, fp, p)
    implicit none
    sll_real64, intent(in) :: fi(:), p
    sll_real64, intent(out) :: fp(:)
    sll_real64 :: pp(11)
    sll_int32 :: i, n
    call lagr_11pt_coeff(pp, p)
    n = size(fi)
    do i=6, n-5
      fp(i) = pp(1) * fi(i-5) &
        + pp(2) * fi(i-4) &
        + pp(3) * fi(i-3) &
        + pp(4) * fi(i-2) &
        + pp(5) * fi(i-1) &
        + pp(6) * fi(i)   &
        + pp(7) * fi(i+1) &
        + pp(8) * fi(i+2) &
        + pp(9) * fi(i+3) &
        + pp(10)* fi(i+4) &
        + pp(11)* fi(i+5)
    enddo
  end subroutine lagr_11pt_vec


  ! --- Lagrange interpolation, without boundary conditions ---
  ! fi(:)      input array of length n
  ! fp(:)      output array of length n
  ! p          offset in units of dx
  ! stencil    number of points in fi used for interpolation
  subroutine lagrange(fi, fp, p, stencil)
    implicit none
    sll_real64, intent(in) :: fi(:)
    sll_real64, intent(out) :: fp(:)
    sll_real64, intent(in) :: p
    sll_int32, intent(in) :: stencil
    sll_int32 :: n, i

    n = size(fi)

    select case (stencil)
      case(3)
        i=1
        fp(i) = lagr_3pt(fi(i), fi(i+1), fi(i+2), p-1.)
        call lagr_3pt_vec(fi, fp, p)
        i=n
        fp(i) = lagr_3pt(fi(i-2), fi(i-1), fi(i), p+1.)

      case(5)
        i=1
        fp(i) = lagr_5pt(fi(i), fi(i+1), fi(i+2), fi(i+3), fi(i+4), p-2.)
        i=2
        fp(i) = lagr_5pt(fi(i-1), fi(i), fi(i+1), fi(i+2), fi(i+3), p-1.)
        call lagr_5pt_vec(fi, fp, p)
        i=n-1
        fp(i) = lagr_5pt(fi(i-3), fi(i-2), fi(i-1), fi(i), fi(i+1), p+1.)
        i=n
        fp(i) = lagr_5pt(fi(i-4), fi(i-3), fi(i-2), fi(i-1), fi(i), p+2.)

      case default
        write(*,*) "Lagrange stencil not implemented."
    end select

  end subroutine lagrange


  ! --- Lagrange interpolation, periodic boundary conditions ---
  ! fi(:)      input array of length n
  ! fp(:)      output array of length n
  ! p          offset in units of dx
  ! stencil    number of points in fi used for interpolation
  subroutine lagrange_periodic(fi, fp, p, stencil)
    implicit none
    sll_real64, intent(in) :: fi(:)
    sll_real64, intent(out) :: fp(:)
    sll_real64, intent(in) :: p
    sll_int32, intent(in) :: stencil

    sll_int32 :: n
    n = size(fi)

    select case (stencil)

      case(3)
        fp(1) = lagr_3pt(fi(n), fi(1), fi(2), p)
        call lagr_3pt_vec(fi, fp, p)
        fp(n) = lagr_3pt(fi(n-1), fi(n), fi(1), p)

      case(5)
        fp(1) = lagr_5pt(fi(n-1), fi(n), fi(1), fi(2), fi(3), p)
        fp(2) = lagr_5pt(fi(n), fi(1), fi(2), fi(3), fi(4), p)
        call lagr_5pt_vec(fi, fp, p)
        fp(n-1) = lagr_5pt(fi(n-3), fi(n-2), fi(n-1), fi(n), fi(1), p)
        fp(n) = lagr_5pt(fi(n-2), fi(n-1), fi(n), fi(1), fi(2), p)

      case default
        write(*,*) "Lagrange stencil not implemented."
    end select
  end subroutine lagrange_periodic


  ! --- Lagrange interpolation with halo cell boundaries, ---
  !     where the input array already contains halo cells.
  !     ==> This is what you would typically use for an
  !         MPI decomposition with ghost cells.
  !
  ! fi(:)      input array of length n, including the halos
  ! fp(:)      output array of length n, only the inner part is overwritten
  !            (ie boundaries of half stencil width are untouched)
  ! p          offset in units of dx
  ! stencil    number of points {3,5,7,9,11} in fi used for interpolation
  subroutine lagrange_halo_cells(fi, fp, p, stencil)
    implicit none
    sll_real64, intent(in) :: fi(:)
    sll_real64, intent(out) :: fp(:)
    sll_real64, intent(in) :: p
    sll_int32, intent(in) :: stencil

    select case (stencil)
      case(3)
        call lagr_3pt_vec(fi, fp, p)

      case(5)
        call lagr_5pt_vec(fi, fp, p)

      case(7)
        call lagr_7pt_vec(fi, fp, p)

      case(9)
        call lagr_9pt_vec(fi, fp, p)

      case(11)
        call lagr_11pt_vec(fi, fp, p)

      case default
        write(*,*) "Lagrange stencil not implemented."
    end select

  end subroutine lagrange_halo_cells


end module sll_m_lagrange_fast
