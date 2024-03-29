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

module sll_m_lagrange_interpolation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_working_precision.h"

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_dirichlet, &
      sll_p_periodic

   implicit none

   public :: &
      sll_s_compact_derivative_weight, &
      sll_s_compute_stencil_plus, &
      sll_f_lagrange_interpolate, &
      sll_o_weight_product_x1

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   integer, parameter :: SLL_SIZE_STENCIL_MAX = 30
   integer, parameter :: SLL_SIZE_STENCIL_MIN = -30

   interface sll_o_weight_product_x1
      module procedure weight_product1d, &
         weight_product2d_x1
   end interface

   interface weight_product_x2
      module procedure weight_product2d_x1
   end interface

contains

   ! sll_f_lagrange_interpolate returns the value of y(x), using a Lagrange
   ! interpolation based on the given array values of x(i) and y(x(i)).
   ! The interpolating polynomial is of degree 'degree'.
   function sll_f_lagrange_interpolate(x, degree, xi, yi)
      sll_real64, intent(in)               :: x
      sll_int32, intent(in)                :: degree
      sll_real64, intent(in), dimension(:) :: xi
      sll_real64, intent(in), dimension(:) :: yi
      sll_real64                           :: sll_f_lagrange_interpolate
      sll_real64, dimension(1:degree + 1)    :: p
      sll_int32                            :: i
      sll_int32                            :: deg
      sll_int32                            :: m  ! step size

      SLL_ASSERT(size(xi) >= degree + 1)
      SLL_ASSERT(size(yi) >= degree + 1)
      SLL_ASSERT(degree >= 0)

      if ((x < xi(1)) .or. (x > xi(degree + 1))) then
         !print *, 'sll_f_lagrange_interpolate() warning: x is outside of the range ', &
         !     'xi given as argument.'
         !print *, 'x = ', x, 'xmin = ', xi(1), 'xmax = ',xi(degree+1)
      end if

      ! Load the local array with the values of y(i), which are also the values
      ! of the Lagrange polynomials of order 0.
      do i = 1, degree + 1
         p(i) = yi(i)
      end do

      ! Build the Lagrange polynomials up to the desired degree. Here we use
      ! Neville's recursive relation, starting from the zeroth-order
      ! polynomials and working towards the desired degree.
      m = 1
      do deg = degree, 1, -1
         do i = 1, deg
            p(i) = (1.0_f64/(xi(i) - xi(i + m)))*((x - xi(i + m))*p(i) + (xi(i) - x)*p(i + 1))
         end do
         m = m + 1
      end do
      sll_f_lagrange_interpolate = p(1)
   end function sll_f_lagrange_interpolate

   ! compute weights w for the derivatives
   ! f'_j = sum_{k=r}^{s} w_kf_{j+k}
   ! such that the formula is exact for f polynomial of degree  <= d
   ! where d is the biggest possible
   ! we call it sometimes compact finite difference formula
   ! we suppose that r<=0 and s>=0
   subroutine sll_s_compact_derivative_weight(w, r, s)
      integer, intent(in)::r, s
      sll_real64, intent(out)::w(r:s)
      sll_int32 ::i, j
      sll_real64::tmp

      w = 0._f64

      !maple code for generation of w
      !for k from r to -1 do
      !  C[k]:=product((k-j),j=r..k-1)*product((k-j),j=k+1..s):
      !  C[k]:=1/C[k]*product((-j),j=r..k-1)*product((-j),j=k+1..-1)*product((-j),j=1..s):
      !od:
      !for k from 1 to s do
      !  C[k]:=product((k-j),j=r..k-1)*product((k-j),j=k+1..s):
      !  C[k]:=1/C[k]*product((-j),j=r..-1)*product((-j),j=1..k-1)*product((-j),j=k+1..s):
      !od:
      !C[0]:=-add(C[k],k=r..-1)-add(C[k],k=1..s):

      do i = r, -1
         tmp = 1._f64
         do j = r, i - 1
            tmp = tmp*real(i - j, f64)
         end do
         do j = i + 1, s
            tmp = tmp*real(i - j, f64)
         end do
         tmp = 1._f64/tmp
         do j = r, i - 1
            tmp = tmp*real(-j, f64)
         end do
         do j = i + 1, -1
            tmp = tmp*real(-j, f64)
         end do
         do j = 1, s
            tmp = tmp*real(-j, f64)
         end do
         w(i) = tmp
      end do

      do i = 1, s
         tmp = 1._f64
         do j = r, i - 1
            tmp = tmp*real(i - j, f64)
         end do
         do j = i + 1, s
            tmp = tmp*real(i - j, f64)
         end do
         tmp = 1._f64/tmp
         do j = r, -1
            tmp = tmp*real(-j, f64)
         end do
         do j = 1, i - 1
            tmp = tmp*real(-j, f64)
         end do
         do j = i + 1, s
            tmp = tmp*real(-j, f64)
         end do
         w(i) = tmp
      end do
      tmp = 0._f64
      do i = r, -1
         tmp = tmp + w(i)
      end do
      do i = 1, s
         tmp = tmp + w(i)
      end do
      w(0) = -tmp

   end subroutine sll_s_compact_derivative_weight

   subroutine sll_s_compute_stencil_plus(p, r, s)
      sll_int32, intent(in) :: p
      sll_int32, intent(out) :: r, s
      r = -p/2
      s = (p + 1)/2
   end subroutine sll_s_compute_stencil_plus

   subroutine compute_stencil_minus(p, r, s)
      sll_int32, intent(in) :: p
      sll_int32, intent(out) :: r, s
      r = -(p + 1)/2
      s = p/2
   end subroutine compute_stencil_minus

   ! Makes the following operation
   ! fout(j) = sum_{k=r}^s w(k)fin(j+k),j=1..N+1
   ! we suppose periodic boundary conditions fin(j+N)=fin(j)
   ! and that fin(j) ar known j=1,..,N
   subroutine weight_product1d_per(fin, fout, N, w, r, s)
      sll_real64, dimension(:), intent(in)  :: fin
      sll_real64, dimension(:), intent(out) :: fout
      sll_int32, intent(in) :: N
      sll_int32, intent(in) :: r
      sll_int32, intent(in) :: s
      sll_real64, intent(in) :: w(r:s)
      sll_int32 :: i
      sll_int32 :: i1
      sll_int32 :: j
      sll_real64 :: tmp

      do i = 1, N
         tmp = 0._f64
         do j = r, s
            i1 = i + j
            if ((i1 > N) .or. (i1 < 1)) then
               i1 = modulo(i1 - 1, N) + 1
            end if
            tmp = tmp + w(j)*fin(i1)
         end do
         fout(i) = tmp
      end do
      fout(N + 1) = fout(1)
   end subroutine weight_product1d_per

   ! Makes the following operation
   ! fout(j) = sum_{k=r}^s w(k)fin(j+k),j=1..N+1
   ! we suppose natural boundary conditions fin(j)=fin(1), j<=1
   ! and fin(j)=fin(N+1) for j>=N+1
   ! and that fin(j) ar known j=1,..,N+1
   subroutine weight_product1d_nat(fin, fout, N, w, r, s)
      sll_real64, dimension(:), intent(in)  :: fin
      sll_real64, dimension(:), intent(out) :: fout
      sll_int32, intent(in) :: N
      sll_int32, intent(in) :: r
      sll_int32, intent(in) :: s
      sll_real64, intent(in) :: w(r:s)
      sll_int32 :: i
      sll_int32 :: i1
      sll_int32 :: j
      sll_real64 :: tmp

      do i = 1, N + 1
         tmp = 0._f64
         do j = r, s
            i1 = i + j
            if (i1 >= N + 1) then
               i1 = N + 1
            end if
            if (i1 <= 1) then
               i1 = 1
            end if
            tmp = tmp + w(j)*fin(i1)
         end do
         fout(i) = tmp
      end do

   end subroutine weight_product1d_nat

   subroutine weight_product1d(fin, fout, N, w, r, s, bc_type)
      sll_real64, dimension(:), intent(in)  :: fin
      sll_real64, dimension(:), intent(out) :: fout
      sll_int32, intent(in) :: N
      sll_int32, intent(in) :: r
      sll_int32, intent(in) :: s
      sll_real64, intent(in) :: w(r:s)
      sll_int32, intent(in) :: bc_type

      select case (bc_type)
      case (sll_p_periodic)
         call weight_product1d_per(fin, fout, N, w(r:s), r, s)
      case (sll_p_dirichlet)
         call weight_product1d_nat(fin, fout, N, w(r:s), r, s)
      case default
         print *, '#bad type in weight_product1d'
         stop
      end select
   end subroutine weight_product1d

!the two following functions may be not useful

   subroutine weight_product2d_x1(fin, fout, N, w, r, s, bc_type)
      sll_real64, dimension(:, :), intent(in)  :: fin
      sll_real64, dimension(:, :), intent(out) :: fout
      sll_int32, intent(in) :: N(2)
      sll_int32, intent(in) :: r
      sll_int32, intent(in) :: s
      sll_real64, intent(in) :: w(r:s)
      sll_int32, intent(in) :: bc_type
      sll_int32 :: i
      do i = 1, N(2) + 1
         call weight_product1d(fin(:, i), fout(:, i), N(1), w(r:s), r, s, bc_type)
      end do

   end subroutine weight_product2d_x1

   subroutine weight_product2d_x2(fin, fout, N, w, r, s, bc_type)
      sll_real64, dimension(:, :), intent(in)  :: fin
      sll_real64, dimension(:, :), intent(out) :: fout
      sll_int32, intent(in) :: N(2)
      sll_int32, intent(in) :: r
      sll_int32, intent(in) :: s
      sll_real64, intent(in) :: w(r:s)
      sll_int32, intent(in) :: bc_type
      sll_int32 :: i
      sll_real64, dimension(:), allocatable :: bufin
      sll_real64, dimension(:), allocatable :: bufout
      allocate (bufin(N(2) + 1))
      allocate (bufout(N(2) + 1))

      do i = 1, N(1) + 1
         bufin(1:N(2) + 1) = fin(i, 1:N(2) + 1)
         call weight_product1d(bufin, bufout, N(2), w(r:s), r, s, bc_type)
         fout(i, 1:N(2) + 1) = bufout(1:N(2) + 1)
      end do

      deallocate (bufin)
      deallocate (bufout)

   end subroutine weight_product2d_x2

   subroutine compute_size_hermite1d(N, p, N_size, dim_size)
      sll_int32, intent(in) :: N
      sll_int32, intent(in) :: p
      sll_int32, intent(out) :: N_size
      sll_int32, intent(out) :: dim_size

      N_size = N + 1

      if (modulo(p, 2) == 0) then
         dim_size = 2
      else
         dim_size = 3
      end if

   end subroutine compute_size_hermite1d

   function initialize_interpolants_hermite1d(N, p)
      sll_real64, dimension(:, :), pointer  :: initialize_interpolants_hermite1d
      sll_int32, intent(in) :: N
      sll_int32, intent(in) :: p
      sll_int32  :: N_size
      sll_int32  :: dim_size

      call compute_size_hermite1d(N, p, N_size, dim_size)
      allocate (initialize_interpolants_hermite1d(dim_size, N_size))

   end function initialize_interpolants_hermite1d

   function initialize_interpolants_hermite2d(N, p)
      sll_real64, dimension(:, :, :, :), pointer  :: initialize_interpolants_hermite2d
      sll_int32, intent(in) :: N(2)
      sll_int32, intent(in) :: p(2)
      sll_int32 :: dim_size(2), N_size(2)
      sll_int32 :: i
      do i = 1, 2
         call compute_size_hermite1d(N(i), p(i), N_size(i), dim_size(i))
      end do

      allocate (initialize_interpolants_hermite2d(dim_size(1), &
                                                  dim_size(2), &
                                                  N_size(1), &
                                                  N_size(2)))

   end function initialize_interpolants_hermite2d

   subroutine compute_interpolants_hermite1d(fin, coef, N, p, bc_type)
      sll_real64, dimension(:), intent(in)  :: fin
      sll_real64, dimension(:, :), pointer :: coef
      sll_int32, intent(in) :: N
      sll_int32, intent(in) :: p
      sll_int32, intent(in) :: bc_type
      sll_real64 :: w(SLL_SIZE_STENCIL_MIN:SLL_SIZE_STENCIL_MIN)
      sll_int32 :: r, s

      coef(1, 1:N + 1) = fin(1:N + 1)
      call sll_s_compute_stencil_plus(p, r, s)
      call sll_s_compact_derivative_weight(w(r:s), r, s)
      call sll_o_weight_product_x1(fin, coef(2, :), N, w(r:s), r, s, bc_type)
      if (modulo(p, 2) == 1) then
         call compute_stencil_minus(p, r, s)
         call sll_s_compact_derivative_weight(w(r:s), r, s)
         call sll_o_weight_product_x1(fin, coef(3, :), N, w(r:s), r, s, bc_type)
      end if

   end subroutine compute_interpolants_hermite1d

   subroutine compute_interpolants_hermite2d(fin, coef, N, p, bc_type)
      sll_real64, dimension(:, :), intent(in)  :: fin
      sll_real64, dimension(:, :, :, :), pointer :: coef
      sll_int32, intent(in) :: N(2)
      sll_int32, intent(in) :: p(2)
      sll_int32, intent(in) :: bc_type(2)
      !sll_real64 :: w(SLL_SIZE_STENCIL_MIN:SLL_SIZE_STENCIL_MIN)
      !sll_int32 :: r,s
      sll_int32 :: i, j, k
      sll_real64, dimension(:), pointer :: bufin_x1, bufin_x2
      sll_real64, dimension(:, :), pointer :: bufout_x1, bufout_x2
      !sll_int32 :: dim(2)
      sll_int32  :: N_size(2)
      sll_int32  :: dim_size(2)

      do i = 1, 2
         call compute_size_hermite1d(N(i), p(i), N_size(i), dim_size(i))
      end do

      allocate (bufin_x1(N(1) + 1))
      allocate (bufout_x1(dim_size(1), N_size(1)))
      allocate (bufin_x2(N(2) + 1))
      allocate (bufout_x2(dim_size(2), N_size(2)))

      ! compute_interpolants in x2
      do i = 1, N(1) + 1
         bufin_x2(1:N(2) + 1) = fin(i, 1:N(2) + 1)
         call compute_interpolants_hermite1d(bufin_x2, bufout_x2, N(2), p(2), bc_type(2))
         do j = 1, dim_size(2)
            coef(1, j, i, 1:N_size(2)) = bufout_x2(j, 1:N_size(2))
         end do
      end do

      ! compute_interpolants in x1
      do i = 1, N_size(2)
         do j = 1, dim_size(2)
            bufin_x1(1:N(1) + 1) = coef(1, j, 1:N(1) + 1, i)
            call compute_interpolants_hermite1d(bufin_x1, bufout_x1, N(1), p(1), bc_type(1))
            do k = 1, dim_size(1)
               coef(k, j, 1:N_size(1), i) = bufout_x1(k, 1:N_size(1))
            end do
         end do
      end do

      deallocate (bufin_x1)
      deallocate (bufout_x1)
      deallocate (bufin_x2)
      deallocate (bufout_x2)

   end subroutine compute_interpolants_hermite2d

!  subroutine interpolate_hermite1d_c1(coef,N,i,x,fval)
!    sll_real64, dimension(:,:), pointer :: coef
!    sll_int32, intent(in) :: N
!    sll_int32, intent(in) :: i
!    sll_real64, intent(in) :: x
!    sll_real64, intent(out) :: fval
!    sll_real64 :: w(0:3)
!    sll_int32  :: i1
!
!    w(0)=(2._f64*x+1._f64)*(1._f64-x)*(1._f64-x)
!    w(1)=x*x*(3._f64-2._f64*x)
!    w(2)=x*(1._f64-x)*(1._f64-x)
!    w(3)=x*x*(x-1._f64)
!    i1=i+1
!
!    fval=w(0)*coef(1,i)+w(1)*coef(1,i1)+w(2)*coef(2,i)+w(3)*coef(2,i1)
!
!  end subroutine interpolate_hermite1d_c1

!  subroutine interpolate_hermite1d_c0(coef,N,i,x,fval)
!    sll_real64, dimension(:,:), pointer :: coef
!    sll_int32, intent(in) :: N
!    sll_int32, intent(in) :: i
!    sll_real64, intent(in) :: x
!    sll_real64, intent(out) :: fval
!    sll_real64 :: w(0:3)
!    sll_int32  :: i1
!
!    w(0)=(2._f64*x+1._f64)*(1._f64-x)*(1._f64-x);
!    w(1)=x*x*(3._f64-2._f64*x)
!    w(2)=x*(1._f64-x)*(1._f64-x)
!    w(3)=x*x*(x-1._f64)
!    i1=i+1
!
!    fval=w(0)*coef(1,i)+w(1)*coef(1,i1)+w(2)*coef(2,i)+w(3)*coef(3,i1)
!
!  end subroutine interpolate_hermite1d_c0

!  subroutine interpolate_hermite2d_c0_c1(coef,N,i,x,fval)
!
!    sll_real64, dimension(:,:,:,:), pointer :: coef
!    sll_int32, intent(in) :: N(2)
!    sll_int32, intent(in) :: i(2)
!    sll_real64, intent(in) :: x(2)
!    sll_real64, intent(out) :: fval
!
!    fval = 0
!
!
!  end subroutine interpolate_hermite2d_c0_c1

end module sll_m_lagrange_interpolation
