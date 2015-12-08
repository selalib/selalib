!> @ingroup lagrange_interpolation
!> @author Klaus Reuter, MPCDF
!> @contact Katharina Kormann, IPP
!> @brief
!> Module for 1D Lagrange interpolation on a uniform grid (only odd order)
!> @details
!> This is an alternative implementation of the Lagrange interpolation for equidistant grids. The only function implemented is an interpolation for a given displacement (interpolate_array_disp). The purpose of this implementation is to provide a fast alternative that exploits the simplifications in this special case.
!> Note: The implementation is based on the formulas in Abramowitz and Stegun: Handbook of Mathematical Functions, Chapter 25.2
module sll_m_lagrange_fast
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"
  implicit none

  public :: sll_s_interpolate_array_disp_lagrange_fixed_no_bc, &
       sll_s_interpolate_array_disp_lagrange_fixed_periodic, &
       sll_s_interpolate_array_disp_lagrange_fixed_halo_cells

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! --- compile-time constants to avoid run-time division
  sll_real64, parameter :: inv_6       = 1._f64/6._f64
  sll_real64, parameter :: inv_12      = 1._f64/12._f64
  sll_real64, parameter :: inv_24      = 1._f64/24._f64
  sll_real64, parameter :: inv_36      = 1._f64/36._f64
  sll_real64, parameter :: inv_48      = 1._f64/48._f64
  sll_real64, parameter :: inv_120     = 1._f64/120._f64
  sll_real64, parameter :: inv_576     = 1._f64/576._f64
  sll_real64, parameter :: inv_720     = 1._f64/720._f64
  sll_real64, parameter :: inv_1440    = 1._f64/1440._f64
  sll_real64, parameter :: inv_5040    = 1._f64/5040._f64
  sll_real64, parameter :: inv_14400   = 1._f64/14400._f64
  sll_real64, parameter :: inv_17280   = 1._f64/17280._f64
  sll_real64, parameter :: inv_30240   = 1._f64/30240._f64
  sll_real64, parameter :: inv_40320   = 1._f64/40320._f64
  sll_real64, parameter :: inv_80640   = 1._f64/80640._f64
  sll_real64, parameter :: inv_362880  = 1._f64/362880._f64
  sll_real64, parameter :: inv_3628800 = 1._f64/3628800._f64

contains

!------------------------------------------------------------------------------!
  !> Compute coefficients for Lagrange interpolation for normalized displacement \a p
  subroutine lagr_3pt_coeff(pp, p)
    sll_real64, intent(out) :: pp(3) !< Lagrange interpolations coefficients
    sll_real64, intent(in)  :: p     !< displacement in units of grid spacing

    pp(1) = p*(p-1.)*0.5
    pp(2) = 1. - p*p
    pp(3) = p*(p+1.)*0.5
  end subroutine



!------------------------------------------------------------------------------!
  !> single point 3-pt-lagrange interpolation
  function lagr_3pt(fm1, f0, f1, p)
    sll_real64 :: lagr_3pt        !< interpolated value
    sll_real64, intent(in) :: fm1 !< known function values at point -1 (relative to where we want to interpolate)
    sll_real64, intent(in) :: f0  !< known function values at point 0 (relative to where we want to interpolate)
    sll_real64, intent(in) :: f1  !< known function values at point 1 (relative to where we want to interpolate)
    sll_real64, intent(in) :: p   !< displacement in units of grid spacing

    sll_real64 :: pp(3)
    call lagr_3pt_coeff(pp, p)
    lagr_3pt = pp(1) * fm1 &
      + pp(2) * f0  &
      + pp(3) * f1
  end function lagr_3pt


!------------------------------------------------------------------------------!
  !> vectorizable 3-pt-lagrange interpolation
  subroutine lagr_3pt_vec(fi, fp, p)
    sll_real64, intent(in)  :: fi(:) !< known function values 
    sll_real64, intent(out) :: fp(:) !< interpolated function values
    sll_real64, intent(in)  :: p     !< displacement in units of grid spacing

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


!------------------------------------------------------------------------------!
  !> Compute coefficients for Lagrange interpolation for normalized displacement \a p
  subroutine lagr_5pt_coeff(pp, p)
    sll_real64, intent(out) :: pp(5)  !< Lagrange interpolations coefficients
    sll_real64, intent(in)  :: p      !< displacement in units of grid spacing

    pp(1) = (p*p-1.)*p*(p-2.)*inv_24
    pp(2) = -(p-1.)*p*(p*p-4.)*inv_6
    pp(3) = (p*p-1.)*(p*p-4.)*0.25
    pp(4) = -(p+1.)*p*(p*p-4.)*inv_6
    pp(5) = (p*p-1.)*p*(p+2.)*inv_24
  end subroutine
 

!------------------------------------------------------------------------------!
  !> single point 5-pt-lagrange interpolation
  function lagr_5pt(fm2, fm1, f0, f1, f2, p)
    sll_real64 :: lagr_5pt        !< interpolated value
    sll_real64, intent(in) :: fm2 !< known function values at point -2 (relative to where we want to interpolate)
    sll_real64, intent(in) :: fm1 !< known function values at point -1 (relative to where we want to interpolate)
    sll_real64, intent(in) :: f0  !< known function values at point 0 (relative to where we want to interpolate)
    sll_real64, intent(in) :: f1  !< known function values at point 1 (relative to where we want to interpolate))
    sll_real64, intent(in) :: f2  !< known function values at point 2 (relative to where we want to interpolate)
    sll_real64, intent(in) :: p   !< displacement in units of grid spacing

    sll_real64 :: pp(5)

    call lagr_5pt_coeff(pp, p)
    lagr_5pt = pp(1) * fm2 &
      + pp(2) * fm1 &
      + pp(3) * f0  &
      + pp(4) * f1  &
      + pp(5) * f2
  end function lagr_5pt


!------------------------------------------------------------------------------!
  !> vectorizable 5-pt-lagrange interpolation
  subroutine lagr_5pt_vec(fi, fp, p)
    sll_real64, intent(in)  :: fi(:) !< known function values 
    sll_real64, intent(out) :: fp(:) !< interpolated function values
    sll_real64, intent(in)  :: p     !< displacement in units of grid spacing

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



!------------------------------------------------------------------------------!
  !> Compute coefficients for Lagrange interpolation for normalized displacement \a p
  subroutine lagr_7pt_coeff(pp, p)
    sll_real64, intent(out) :: pp(7) !< Lagrange interpolations coefficients
    sll_real64, intent(in) :: p      !< displacement in units of grid spacing

    pp(1) = p*(p-3)*(p**2-4)*(p**2-1)*inv_720
    pp(2) = -p*(p-2)*(p**2-9)*(p**2-1)*inv_120
    pp(3) = p*(p-1)*(p**2-9)*(p**2-4)*inv_48
    pp(4) = -(p**2-9)*(p**2-4)*(p**2-1)*inv_36
    pp(5) = (p+1)*p*(p**2-9)*(p**2-4)*inv_48
    pp(6) = -(p+2)*p*(p**2-9)*(p**2-1)*inv_120
    pp(7) = (p+3)*p*(p**2-4)*(p**2-1)*inv_720
  end subroutine


!------------------------------------------------------------------------------!
  !> single point 7-pt-lagrange interpolation
  function lagr_7pt(fm3, fm2, fm1, f0, f1, f2, f3, p)
    implicit none
    sll_real64 :: lagr_7pt !< interpolated value
    sll_real64, intent(in) :: fm3 !< known function values at point -3 (relative to where we want to interpolate)
    sll_real64, intent(in) :: fm2 !< known function values at point -2 (relative to where we want to interpolate)
    sll_real64, intent(in) :: fm1 !< known function values at point -1 (relative to where we want to interpolate)
    sll_real64, intent(in) :: f0  !< known function values at point 0 (relative to where we want to interpolate)
    sll_real64, intent(in) :: f1  !< known function values at point 1 (relative to where we want to interpolate))
    sll_real64, intent(in) :: f2  !< known function values at point 2 (relative to where we want to interpolate)
    sll_real64, intent(in) :: f3  !< known function values at point 3 (relative to where we want to interpolate)
    sll_real64, intent(in) :: p   !< displacement in units of grid spacing
   
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


!------------------------------------------------------------------------------!
  !> vectorizable 7-pt-lagrange interpolation
  subroutine lagr_7pt_vec(fi, fp, p)
    sll_real64, intent(in)  :: fi(:) !< known function values 
    sll_real64, intent(out) :: fp(:) !< interpolated function values
    sll_real64, intent(in)  :: p     !< displacement in units of grid spacing

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


!------------------------------------------------------------------------------!
  !> Compute coefficients for Lagrange interpolation for normalized displacement \a p
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


!------------------------------------------------------------------------------!
  !> single point 9-pt-lagrange interpolation
  function lagr_9pt(fm4, fm3, fm2, fm1, f0, f1, f2, f3, f4, p)
    implicit none
    sll_real64 :: lagr_9pt !< interpolated value
    sll_real64, intent(in) :: fm4 !< known function values at point -4 (relative to where we want to interpolate)
    sll_real64, intent(in) :: fm3 !< known function values at point -3 (relative to where we want to interpolate)
    sll_real64, intent(in) :: fm2 !< known function values at point -2 (relative to where we want to interpolate)
    sll_real64, intent(in) :: fm1 !< known function values at point -1 (relative to where we want to interpolate)
    sll_real64, intent(in) :: f0  !< known function values at point 0 (relative to where we want to interpolate)
    sll_real64, intent(in) :: f1  !< known function values at point 1 (relative to where we want to interpolate))
    sll_real64, intent(in) :: f2  !< known function values at point 2 (relative to where we want to interpolate)
    sll_real64, intent(in) :: f3  !< known function values at point 3 (relative to where we want to interpolate)
    sll_real64, intent(in) :: f4  !< known function values at point 4 (relative to where we want to interpolate)
    sll_real64, intent(in) :: p   !< displacement in units of grid spacing

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


!------------------------------------------------------------------------------!
  !> vectorizable 9-pt-lagrange interpolation
  subroutine lagr_9pt_vec(fi, fp, p)
    sll_real64, intent(in)  :: fi(:) !< known function values 
    sll_real64, intent(out) :: fp(:) !< interpolated function values
    sll_real64, intent(in)  :: p     !< displacement in units of grid spacing

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


!------------------------------------------------------------------------------!
  !> Compute coefficients for Lagrange interpolation for normalized displacement \a p
  subroutine lagr_11pt_coeff(pp, p)
    sll_real64, intent(out) :: pp(11) !< Lagrange interpolations coefficients
    sll_real64, intent(in)  :: p      !< displacement in units of grid spacing

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


!------------------------------------------------------------------------------!
  !> single point 11-pt-lagrange interpolation
  function lagr_11pt(fm5, fm4, fm3, fm2, fm1, f0, f1, f2, f3, f4, f5, p)
    implicit none
    sll_real64 :: lagr_11pt !< interpolated value
    sll_real64, intent(in) :: fm5 !< known function values at point -5 (relative to where we want to interpolate)
    sll_real64, intent(in) :: fm4 !< known function values at point -4 (relative to where we want to interpolate)
    sll_real64, intent(in) :: fm3 !< known function values at point -3 (relative to where we want to interpolate)
    sll_real64, intent(in) :: fm2 !< known function values at point -2 (relative to where we want to interpolate)
    sll_real64, intent(in) :: fm1 !< known function values at point -1 (relative to where we want to interpolate)
    sll_real64, intent(in) :: f0  !< known function values at point 0 (relative to where we want to interpolate)
    sll_real64, intent(in) :: f1  !< known function values at point 1 (relative to where we want to interpolate))
    sll_real64, intent(in) :: f2  !< known function values at point 2 (relative to where we want to interpolate)
    sll_real64, intent(in) :: f3  !< known function values at point 3 (relative to where we want to interpolate)
    sll_real64, intent(in) :: f4  !< known function values at point 4 (relative to where we want to interpolate)
    sll_real64, intent(in) :: f5  !< known function values at point 5 (relative to where we want to interpolate)
    sll_real64, intent(in) :: p   !< displacement in units of grid spacing

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


!------------------------------------------------------------------------------!
  !> vectorizable 11-pt-lagrange interpolation
  subroutine lagr_11pt_vec(fi, fp, p)
    sll_real64, intent(in)  :: fi(:) !< known function values 
    sll_real64, intent(out) :: fp(:) !< interpolated function values
    sll_real64, intent(in)  :: p     !< displacement in units of grid spacing

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


!------------------------------------------------------------------------------!
  !> @brief Lagrange interpolation, without boundary conditions
  !> @param [in] fi(:)      input array of length n
  !> @param [out] fp(:)      output array of length n
  !> @param [in]  p          offset in units of dx
  !> @param [in]  stencil    number of points in fi used for interpolation (possible values 3,5)
  subroutine sll_s_interpolate_array_disp_lagrange_fixed_no_bc(fi, fp, p, stencil)
    implicit none
    sll_real64, intent(in)   :: fi(:)   
    sll_real64, intent(out)  :: fp(:)  
    sll_real64, intent(in)   :: p
    sll_int32,  intent(in)   :: stencil
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

  end subroutine sll_s_interpolate_array_disp_lagrange_fixed_no_bc


!------------------------------------------------------------------------------!
  !> @brief Lagrange interpolation, periodic boundary conditions
  !> @param [in]  fi(:)      input array of length n
  !> @param [out] fp(:)      output array of length n
  !> @param [in]  p          offset in units of dx
  !> @param [in]  stencil    number of points in fi used for interpolation (possible values 3,5)
  subroutine sll_s_interpolate_array_disp_lagrange_fixed_periodic(fi, fp, p, stencil)
    implicit none
    sll_real64, intent(in)   :: fi(:)
    sll_real64, intent(out)  :: fp(:)
    sll_real64, intent(in)   :: p
    sll_int32,  intent(in)   :: stencil

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
  end subroutine sll_s_interpolate_array_disp_lagrange_fixed_periodic



!------------------------------------------------------------------------------!
  !> @brief Lagrange interpolation with halo cell boundaries,
  !>     where the input array already contains halo cells.
  !>     ==> This is what you would typically use for an
  !>         MPI decomposition with ghost cells.
  !>
  !> @param [in]  fi(:)      input array of length n, including the halos
  !> @param [out] fp(:)      output array of length n, only the inner part is overwritten
  !>            (ie boundaries of half stencil width are untouched)
  !> @param [in] p          offset in units of dx
  !> @param [in] stencil    number of points {3,5,7,9,11} in fi used for interpolation
  subroutine sll_s_interpolate_array_disp_lagrange_fixed_halo_cells(fi, fp, p, stencil)
    implicit none
    sll_real64, intent(in)   :: fi(:)
    sll_real64, intent(out)  :: fp(:)
    sll_real64, intent(in)   :: p
    sll_int32,  intent(in)   :: stencil

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

  end subroutine sll_s_interpolate_array_disp_lagrange_fixed_halo_cells


end module sll_m_lagrange_fast
