!> @ingroup splines
!> Implements arbitrary degree bspline interpolation on a uniform grid
!> given a B-Spline object from sll_m_arbitrary_degree_splines
module sll_m_bspline_2d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_boundary_condition_descriptors, only: &
     sll_p_periodic, &
     sll_p_hermite, &
     sll_p_greville, &
     sll_p_mirror

use sll_m_arbitrary_degree_splines, only: &
     sll_t_arbitrary_degree_spline_1d, &
     sll_s_arbitrary_degree_spline_1d_init, &
     sll_s_arbitrary_degree_spline_1d_free, &
     sll_f_find_cell, &
     sll_s_splines_at_x, &
     sll_s_spline_derivatives_at_x, &
     sll_s_splines_and_n_derivs_at_x, &
     sll_s_uniform_b_splines_at_x

use schur_complement

use sll_m_bspline_1d, only: sll_t_bspline_1d, &
  sll_s_bspline_1d_init,                      &
  sll_s_compute_bspline_1d

implicit none

public ::                                       &
     sll_t_bspline_2d,            &
     sll_s_bspline_2d_init,       &
     sll_s_bspline_2d_free,       &
     sll_s_compute_bspline_2d,                  &
     sll_f_interpolate_value_2d,                &
     sll_s_interpolate_array_values_2d,         &
     sll_f_interpolate_derivative_x1_2d,        &
     sll_f_interpolate_derivative_x2_2d,        &
     sll_s_interpolate_array_derivatives_x1_2d, &
     sll_s_interpolate_array_derivatives_x2_2d

private
  
!> @brief 
!> basic type for two-dimensional B-spline data. 
!> @details 
!> treated as an opaque type. No access to its internals is directly allowed.
type :: sll_t_bspline_2d

  type(sll_t_bspline_1d) :: bs1
  type(sll_t_bspline_1d) :: bs2
  sll_real64, pointer    :: bcoef(:,:)
  sll_real64, pointer    :: bwork(:,:)

end type sll_t_bspline_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief Initialises a 2D spline interpolation object.
!> @param[in] nx1 Number of points where the data to be interpolated are 
!>            represented.
!> @param[in] degree1 Spline degree 
!> @param[in] x1_min Minimum value of the abscissae where the data are meant 
!> to be interpolated.
!> @param[in] x1_max Maximum value of the abscissae where the data are meant 
!> to be interpolated.
!> @param[in] bc1 A boundary condition specifier. Must be one of the
!> symbols defined in the SLL_BOUNDARY_CONDITION_DESCRIPTORS module.
!> @param[in] nx2 Number of points where the data to be interpolated are 
!>            represented.
!> @param[in] degree2 Spline degree 
!> @param[in] x2_min Minimum value of the abscissae where the data are meant 
!> to be interpolated.
!> @param[in] x2_max Maximum value of the abscissae where the data are meant 
!> to be interpolated.
!> @param[in] bc2 A boundary condition specifier. Must be one of the
!> symbols defined in the sll_m_boundary_condition_descriptors module.
!> @param[in] spline_bc_type1 A boundary condition specifier (see sll_s_bspline_1d_init). 
!> @param[in] spline_bc_type2 A boundary condition specifier (see sll_s_bspline_1d_init). 
!> @param[in] bc_left1 value of function on the west boundary
!> @param[in] bc_left2 value of function on the south boundary
!> @param[in] bc_right1 value of function on the east boundary
!> @param[in] bc_right2 value of the function on the north boundary
!> @return a spline interpolation object.
subroutine sll_s_bspline_2d_init( &
  self,            &
  nx1,             &
  nx2,             &
  degree1,         &
  degree2,         &
  x1_min,          &
  x2_min,          &
  x1_max,          &
  x2_max,          &
  bc1,             &
  bc2,             &
  spline_bc_type1, &
  spline_bc_type2, &
  bc_left1,        &
  bc_left2,        &
  bc_right1,       &
  bc_right2        )

  type(sll_t_bspline_2d) :: self

  sll_int32,  intent(in)               :: nx1
  sll_int32,  intent(in)               :: degree1
  sll_real64, intent(in)               :: x1_min
  sll_real64, intent(in)               :: x1_max
  sll_int32,  intent(in)               :: bc1
  sll_int32,  intent(in)               :: nx2
  sll_int32,  intent(in)               :: degree2
  sll_real64, intent(in)               :: x2_min
  sll_real64, intent(in)               :: x2_max
  sll_int32,  intent(in)               :: bc2

  sll_int32,  optional                 :: spline_bc_type1
  sll_int32,  optional                 :: spline_bc_type2
  sll_real64, optional, pointer        :: bc_left1(:)
  sll_real64, optional, pointer        :: bc_left2(:)
  sll_real64, optional, pointer        :: bc_right1(:)
  sll_real64, optional, pointer        :: bc_right2(:)

  sll_int32                            :: n1
  sll_int32                            :: n2
  sll_int32                            :: ierr

  if (present(spline_bc_type1)) then
    call sll_s_bspline_1d_init(self%bs1,nx1,degree1, &
      x1_min,x1_max,bc1, spline_bc_type1)
  else
    call sll_s_bspline_1d_init(self%bs1,nx1,degree1, &
      x1_min,x1_max,bc1)
  end if
  if (present(spline_bc_type2)) then
    call sll_s_bspline_1d_init(self%bs2,nx2,degree2, &
      x2_min,x2_max,bc2, spline_bc_type2)
  else
    call sll_s_bspline_1d_init(self%bs2,nx2,degree2, &
      x2_min,x2_max,bc2)
  end if

  n1 = self%bs1%n
  n2 = self%bs2%n
  SLL_CLEAR_ALLOCATE(self%bwork(1:n2,1:n1), ierr)
  SLL_CLEAR_ALLOCATE(self%bcoef(1:n1,1:n2), ierr)

end subroutine sll_s_bspline_2d_init
  
subroutine sll_s_compute_bspline_2d(self, gtau, &
  val1_min, val1_max, val2_min, val2_max)

  type(sll_t_bspline_2d) :: self 
  sll_real64, intent(in)               :: gtau(:,:)
  sll_real64, intent(in), optional     :: val1_min(:,:)
  sll_real64, intent(in), optional     :: val1_max(:,:)
  sll_real64, intent(in), optional     :: val2_min(:,:)
  sll_real64, intent(in), optional     :: val2_max(:,:)

  sll_int32                            :: i
  sll_int32                            :: j

  if( present(val1_min) .and. present(val1_max)) then
    do j = 1, size(gtau,2)
      call sll_s_compute_bspline_1d( self%bs1, gtau(:,j), &
        val1_min(:,j), val1_max(:,j))
      self%bwork(j,:) = self%bs1%bcoef(:)
    end do

  else
    do j = 1, size(gtau,2)
      call sll_s_compute_bspline_1d( self%bs1, gtau(:,j))
      self%bwork(j,:) = self%bs1%bcoef(:)
    end do
  end if

  if( present(val2_min) .and. present(val2_max)) then
    do i = 1, size(self%bs1%bcoef)
      call sll_s_compute_bspline_1d( self%bs2, self%bwork(:,i), &
        val2_min(:,i), val2_max(:,i))
      self%bcoef(i,:) = self%bs2%bcoef(:)
    end do
  else
    do i = 1, size(self%bs1%bcoef)
      call sll_s_compute_bspline_1d( self%bs2, self%bwork(:,i))
      self%bcoef(i,:) = self%bs2%bcoef(:)
    end do
  end if

end subroutine sll_s_compute_bspline_2d

subroutine sll_s_interpolate_array_values_2d(self, n1, n2, x1, x2, y )

type(sll_t_bspline_2d)  :: self
sll_int32               :: n1
sll_int32               :: n2
sll_real64, intent(in)  :: x1(:,:)
sll_real64, intent(in)  :: x2(:,:)
sll_real64, intent(out) :: y(:,:)

sll_real64, allocatable :: work(:)
sll_int32               :: icell, jcell
sll_int32               :: i1, i2, k1, k2
sll_int32               :: ib, jb
sll_int32               :: j, k

SLL_ASSERT(n1 == size(x1,1) .and. n1 == size(x2,1) )
SLL_ASSERT(n2 == size(x1,2) .and. n2 == size(x2,2) )

k1 = self%bs1%deg+1
k2 = self%bs2%deg+1

allocate(work(k2))

do i2 = 1, n2
  do i1 = 1, n1
    icell = sll_f_find_cell( self%bs1%bsp, x1(i1,i2) )
    jcell = sll_f_find_cell( self%bs2%bsp, x2(i1,i2) )
    work = 0.0_f64
    do j = 1, k2
      jb = mod(jcell+j-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
      call sll_s_splines_at_x(self%bs1%bsp, icell, x1(i1,i2), self%bs1%values)
      do k=1, k1
        ib = mod(icell+k-2-self%bs1%offset+self%bs1%n,self%bs1%n) + 1
        work(j) = work(j) + self%bs1%values(k)*self%bcoef(ib,jb)
      end do
    end do
    call sll_s_splines_at_x(self%bs2%bsp, jcell, x2(i1,i2), self%bs2%values)
    y(i1,i2) = 0.0_f64
    do k=1,k2
      jb = mod(jcell+k-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
      y(i1,i2) = y(i1,i2) + self%bs2%values(k)*work(k)
    enddo
  end do
end do

deallocate(work)

end subroutine sll_s_interpolate_array_values_2d

function sll_f_interpolate_value_2d(self, xi, xj ) result (y)

type(sll_t_bspline_2d)  :: self
sll_real64, intent(in)  :: xi
sll_real64, intent(in)  :: xj
sll_real64              :: y

sll_int32               :: j
sll_int32               :: k
sll_int32               :: ib
sll_int32               :: jb
sll_int32               :: k1
sll_int32               :: k2
sll_int32               :: icell
sll_int32               :: jcell
sll_real64, allocatable :: work(:)

k1 = self%bs1%deg+1
k2 = self%bs2%deg+1
allocate(work(k2))
work = 0.0_f64

icell = sll_f_find_cell( self%bs1%bsp, xi )
jcell = sll_f_find_cell( self%bs2%bsp, xj )

do j = 1, k2
  jb = mod(jcell+j-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
  call sll_s_splines_at_x(self%bs1%bsp, icell, xi, self%bs1%values)
  do k=1, k1
    ib = mod(icell+k-2-self%bs1%offset+self%bs1%n,self%bs1%n) + 1
    work(j) = work(j) + self%bs1%values(k)*self%bcoef(ib,jb)
  end do
end do

call sll_s_splines_at_x(self%bs2%bsp, jcell, xj, self%bs2%values)
y = 0.0_f64
do k=1,k2
  jb = mod(jcell+k-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
  y = y + self%bs2%values(k)*work(k)
enddo

deallocate(work)

end function sll_f_interpolate_value_2d

function sll_f_interpolate_derivative_x1_2d(self, x1, x2 ) result(y)

type(sll_t_bspline_2d) :: self
sll_real64, intent(in) :: x1
sll_real64, intent(in) :: x2
sll_real64             :: y

sll_real64, allocatable :: work(:)
sll_int32               :: i, j
sll_int32               :: k1, k2, icell, jcell
sll_int32               :: ib, jb

k1 = self%bs1%deg+1
k2 = self%bs2%deg+1

icell =  sll_f_find_cell( self%bs1%bsp, x1 )
jcell =  sll_f_find_cell( self%bs2%bsp, x2 )

allocate(work(k2))
work = 0.0_f64
do j = 1, k2
  jb = mod(jcell+j-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
  call sll_s_spline_derivatives_at_x(self%bs1%bsp, icell, x1, self%bs1%values)
  do i = 1, k1
    ib = mod(icell+i-2-self%bs1%offset+self%bs1%n,self%bs1%n) + 1
    work(j) = work(j) + self%bs1%values(i)*self%bcoef(ib,jb)
  enddo
end do

call sll_s_splines_at_x(self%bs2%bsp, jcell, x2, self%bs2%values)
y = 0.0_f64
do j = 1, k2
  jb = mod(jcell+j-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
  y = y + self%bs2%values(j)*work(j)
enddo

deallocate(work)

end function sll_f_interpolate_derivative_x1_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function sll_f_interpolate_derivative_x2_2d(self, x1, x2 ) result(y)

type(sll_t_bspline_2d) :: self
sll_real64, intent(in) :: x1
sll_real64, intent(in) :: x2
sll_real64             :: y

sll_real64, allocatable :: work(:)
sll_int32               :: i, j
sll_int32               :: k1, k2, icell, jcell
sll_int32               :: ib, jb

k1 = self%bs1%deg+1
k2 = self%bs2%deg+1

icell =  sll_f_find_cell( self%bs1%bsp, x1 )
jcell =  sll_f_find_cell( self%bs2%bsp, x2 )

allocate(work(k2))
work = 0.0_f64
do j = 1, k2
  jb = mod(jcell+j-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
  call sll_s_splines_at_x(self%bs1%bsp, icell, x1, self%bs1%values)
  do i = 1, k1
    ib = mod(icell+i-2-self%bs1%offset+self%bs1%n,self%bs1%n) + 1
    work(j) = work(j) + self%bs1%values(i)*self%bcoef(ib,jb)
  enddo
end do

call sll_s_spline_derivatives_at_x(self%bs2%bsp, jcell, x2, self%bs2%values)
y = 0.0_f64
do j = 1, k2
  jb = mod(jcell+j-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
  y = y + self%bs2%values(j)*work(j)
enddo

deallocate(work)

end function sll_f_interpolate_derivative_x2_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sll_s_bspline_2d_free(self ) 

type(sll_t_bspline_2d) :: self


end subroutine sll_s_bspline_2d_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sll_s_interpolate_array_derivatives_x1_2d( self, n1, n2, x1, x2, y)

type(sll_t_bspline_2d) :: self
sll_int32              :: n1
sll_int32              :: n2
sll_real64             :: x1(:,:)
sll_real64             :: x2(:,:)
sll_real64             :: y(:,:)

sll_int32 :: i1, i2
sll_real64, allocatable :: work(:)
sll_int32               :: i, j
sll_int32               :: k1, k2, icell, jcell
sll_int32               :: ib, jb
sll_real64              :: xi, xj, yy

k1 = self%bs1%deg+1
k2 = self%bs2%deg+1
allocate(work(k2))

do i2 = 1, n2
  do i1 = 1, n1

    xi = x1(i1,i2)
    xj = x2(i1,i2)
    icell =  sll_f_find_cell( self%bs1%bsp, xi )
    jcell =  sll_f_find_cell( self%bs2%bsp, xj )
    
    work = 0.0_f64
    do j = 1, k2
      jb = mod(jcell+j-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
      call sll_s_spline_derivatives_at_x(self%bs1%bsp, icell, xi, self%bs1%values)
      do i = 1, k1
        ib = mod(icell+i-2-self%bs1%offset+self%bs1%n,self%bs1%n) + 1
        work(j) = work(j) + self%bs1%values(i)*self%bcoef(ib,jb)
      enddo
    end do
    call sll_s_splines_at_x(self%bs2%bsp, jcell, xj, self%bs2%values)
    yy = 0.0_f64
    do j = 1, k2
      jb = mod(jcell+j-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
      yy = yy + self%bs2%values(j)*work(j)
    enddo
    
    y(i1,i2) = yy
    
  end do
end do

deallocate(work)

end subroutine sll_s_interpolate_array_derivatives_x1_2d

subroutine sll_s_interpolate_array_derivatives_x2_2d( self, n1, n2, x1, x2, y)

type(sll_t_bspline_2d) :: self
sll_int32              :: n1
sll_int32              :: n2
sll_real64             :: x1(:,:)
sll_real64             :: x2(:,:)
sll_real64             :: y(:,:)

sll_real64, allocatable :: work(:)
sll_int32               :: i, j, i1, i2
sll_int32               :: k1, k2, icell, jcell
sll_int32               :: ib, jb
sll_real64              :: xi, xj, yy

k1 = self%bs1%deg+1
k2 = self%bs2%deg+1
allocate(work(k2))

do i2 = 1, n2
  do i1 = 1, n1
    xi = x1(i1,i2)
    xj = x2(i1,i2)
    icell =  sll_f_find_cell( self%bs1%bsp, xi )
    jcell =  sll_f_find_cell( self%bs2%bsp, xj )
    
    work = 0.0_f64
    do j = 1, k2
      jb = mod(jcell+j-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
      call sll_s_splines_at_x(self%bs1%bsp, icell, xi, self%bs1%values)
      do i = 1, k1
        ib = mod(icell+i-2-self%bs1%offset+self%bs1%n,self%bs1%n) + 1
        work(j) = work(j) + self%bs1%values(i)*self%bcoef(ib,jb)
      enddo
    end do
    
    call sll_s_spline_derivatives_at_x(self%bs2%bsp, jcell, xj, self%bs2%values)
    yy = 0.0_f64
    do j = 1, k2
      jb = mod(jcell+j-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
      yy = yy + self%bs2%values(j)*work(j)
    enddo

    y(i1,i2) = yy

  end do
end do

deallocate(work)

end subroutine sll_s_interpolate_array_derivatives_x2_2d

end module sll_m_bspline_2d
