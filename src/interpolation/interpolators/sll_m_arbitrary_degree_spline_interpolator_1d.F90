!**************************************************************
!  Copyright INRIA
!  Authors :
!     Aurore
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

!> Class interpolator and methods for arbitrary degree spline 1D interpolator
module sll_m_arbitrary_degree_spline_interpolator_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_errors.h"

use sll_m_deboor_splines_1d
use sll_m_interpolators_1d_base
use sll_m_fornberg

implicit none
private

!> Class for arbitrary degree spline 1d interpolator
type, public, extends(sll_c_interpolator_1d) :: &
  sll_arbitrary_degree_spline_interpolator_1d

  sll_int32                         :: num_pts       !< nodes number
  sll_real64                        :: eta_min       !< left boundary
  sll_real64                        :: eta_max       !< right boundary
  sll_int32                         :: bc_left       !< left boundary type
  sll_int32                         :: bc_right      !< right boundary type
  sll_int32                         :: spline_degree !< spline degree
  sll_real64, dimension(:), pointer :: eta           !< node positions
  sll_real64, dimension(:), pointer :: t             !< knots positions
  sll_real64, dimension(:), pointer :: work          !< work array to store data
  sll_int32                         :: size_t        !< knots array size
  sll_int64                         :: bc_selector   !< boundary condition type
  sll_real64, dimension(:), pointer :: coeff_splines !< coeffs array
  sll_int32                         :: size_coeffs   !< coeffs array dimension
  sll_real64                        :: slope_left    !< left boundary derivative
  sll_real64                        :: slope_right   !< right boundary derivative
  sll_real64                        :: value_left    !< left boundary value
  sll_real64                        :: value_right   !< right boundary value
  logical                           :: compute_slope_left = .TRUE. !< true
  logical                           :: compute_slope_right= .TRUE. !< true
  logical                           :: compute_value_left = .TRUE. !< true
  logical                           :: compute_value_right= .TRUE. !< true
  type(deboor_type)                 :: deboor  !< Deboor splines data object

contains

  !> Initialize the interpolator
  procedure :: initialize=>initialize_ad1d_interpolator
  !> Set spline coefficients
  procedure :: set_coefficients => set_coefficients_ad1d
  !> Compute interpolants
  procedure :: compute_interpolants => compute_interpolants_ad1d
  !> Interpolate single value
  procedure :: interpolate_value => interpolate_value_ad1d
  !> Interpolate an array 
  procedure :: interpolate_array_values => interpolate_values_ad1d
  !> Interpolate a pointer to array 
  procedure :: interpolate_pointer_values => interpolate_pointer_values_ad1d
  !> Compute derivatives
  procedure :: interpolate_derivative_eta1 => interpolate_derivative_ad1d
  !> Compute derivatives array
  procedure :: interpolate_array_derivatives => interpolate_derivatives_ad1d
  !> Compute derivatives array pointer
  procedure :: interpolate_pointer_derivatives =>interpolate_pointer_derivatives_ad1d
  !> Interpolate an array
  procedure :: interpolate_array => interpolate_array_ad1d
  !> Interpolate an array after displacement
  procedure :: interpolate_array_disp => interpolate_1d_array_disp_ad1d
  !> Get splines coefficients
  procedure :: get_coefficients => get_coefficients_ad1d
  !> Not implemented
  procedure :: reconstruct_array
  !> Destory the derived type and free memory
  procedure :: delete => delete_arbitrary_degree_1d_interpolator

end type sll_arbitrary_degree_spline_interpolator_1d

!> Deallocate
interface sll_delete
   module procedure delete_arbitrary_degree_1d_interpolator
end interface sll_delete

public sll_delete 
public new_arbitrary_degree_1d_interpolator
public set_values_at_boundary1d
public initialize_ad1d_interpolator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief delete interpolator arbitrary degree splines.
!> @details
!>
!> The parameters are
!> @param interpolator the type sll_arbitrary_degree_spline_interpolator_1d

subroutine delete_arbitrary_degree_1d_interpolator( interpolator )
  class(sll_arbitrary_degree_spline_interpolator_1d), intent(inout) :: interpolator
  sll_int32 :: ierr
  SLL_DEALLOCATE(interpolator%t,ierr)
  SLL_DEALLOCATE(interpolator%coeff_splines,ierr)
  deallocate(interpolator%eta)
  deallocate(interpolator%work)
end subroutine delete_arbitrary_degree_1d_interpolator


!> @brief Initialization of a pointer interpolator arbitrary degree splines 1d.
!> @details To have the interpolator arbitrary degree splines 1d such as a pointer
!>
!> The parameters are
!> @param[in] num_pts the number of points
!> @param[in] eta_min the minimun
!> @param[in] eta_max the maximun
!> @param[in] bc_left  the boundary condition at left
!> @param[in] bc_right the boundary condition at right
!> @param[in] spline_degree the degree of B-spline
!> @return the type interpolator arbitrary degree splines 1d

function new_arbitrary_degree_1d_interpolator( num_pts,       &
                                               eta_min,       &
                                               eta_max,       &
                                               bc_left,       &
                                               bc_right,      &
                                               spline_degree) &
result(interpolator)

class(sll_arbitrary_degree_spline_interpolator_1d),pointer :: interpolator
sll_int32,  intent(in) :: num_pts
sll_real64, intent(in) :: eta_min
sll_real64, intent(in) :: eta_max
sll_int32,  intent(in) :: bc_left
sll_int32,  intent(in) :: bc_right
sll_int32,  intent(in) :: spline_degree

sll_int32              :: ierr

SLL_ALLOCATE(interpolator,ierr)

call initialize_ad1d_interpolator( interpolator, &
                                   num_pts,      &
                                   eta_min,      &
                                   eta_max,      &
                                   bc_left,      &
                                   bc_right,     &
                                   spline_degree)

end function new_arbitrary_degree_1d_interpolator

!> @brief Initialization of interpolator arbitrary degree splines 1d.
!> @details To have the interpolator arbitrary degree splines 1d
!>
!> The parameters are
!> @param[in] num_pts the number of points
!> @param[in] eta_min the minimun
!> @param[in] eta_max the maximun
!> @param[in] bc_left  the boundary condition at left
!> @param[in] bc_right the boundary condition at right
!> @param[in] spline_degree the degree of B-spline
!> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_1d

subroutine initialize_ad1d_interpolator( interpolator, &
                                         num_pts,      &
                                         eta_min,      &
                                         eta_max,      &
                                         bc_left,      &
                                         bc_right,     &
                                         spline_degree)

class(sll_arbitrary_degree_spline_interpolator_1d), intent(inout) :: interpolator
sll_int32,       intent(in) :: num_pts
sll_real64,      intent(in) :: eta_min
sll_real64,      intent(in) :: eta_max
sll_int32,       intent(in) :: bc_left
sll_int32,       intent(in) :: bc_right
sll_int32,       intent(in) :: spline_degree

sll_int32                   :: ierr
sll_int32                   :: tmp
sll_int64                   :: bc_selector
sll_int32                   :: i, k
sll_real64                  :: delta_eta
sll_int32                   :: deriv
character(len=*), parameter :: this_sub_name = 'initialize_ad1d_interpolator'

! do some argument checking...
if(((bc_left == SLL_PERIODIC).and.(bc_right.ne. SLL_PERIODIC))) then
   print *, 'initialize_arbitrary_degree_1d_interpolator, ERROR: ', &
        'if one boundary condition is specified as periodic, then ', &
        'both must be. Error in first direction.'
end if

bc_selector = 0

if( bc_left  == SLL_DIRICHLET ) bc_selector = bc_selector + 1
if( bc_left  == SLL_NEUMANN   ) bc_selector = bc_selector + 2
if( bc_left  == SLL_HERMITE   ) bc_selector = bc_selector + 4
if( bc_right == SLL_DIRICHLET ) bc_selector = bc_selector + 8
if( bc_right == SLL_NEUMANN   ) bc_selector = bc_selector + 16
if( bc_right == SLL_HERMITE   ) bc_selector = bc_selector + 32

interpolator%spline_degree = spline_degree
interpolator%eta_min       = eta_min
interpolator%eta_max       = eta_max
interpolator%bc_left       = bc_left
interpolator%bc_right      = bc_right
interpolator%bc_selector   = bc_selector
interpolator%num_pts       = num_pts

delta_eta = (eta_max - eta_min)/(num_pts-1)
SLL_ALLOCATE(interpolator%eta(num_pts),ierr)
do i = 1, num_pts
  interpolator%eta(i) = eta_min + delta_eta*(i-1)
end do

k = spline_degree+1
deriv = 2

select case (bc_selector)

case (0) ! 1. periodic

  SLL_ALLOCATE(interpolator%coeff_splines(num_pts),ierr)
  SLL_ALLOCATE(interpolator%t(num_pts+k),ierr)
  interpolator%t(1:k) = eta_min
  interpolator%t(num_pts+1:num_pts+k) = eta_max
  
  if ( mod(k,2) == 0 ) then
    do i = k+1,num_pts
      interpolator%t(i) = interpolator%eta(i-k/2) 
    end do
  else
    do i = k+1, num_pts
      interpolator%t(i) = 0.5*(interpolator%eta(i-(k-1)/2)+    &
                               interpolator%eta(i-1-(k-1)/2))
    end do
  end if
  SLL_ALLOCATE(interpolator%work(num_pts*(2*k-1)),ierr)

case (9) ! 2. dirichlet-left, dirichlet-right

  interpolator%value_left  = 0.0_f64
  interpolator%value_right = 0.0_f64

  SLL_ALLOCATE(interpolator%coeff_splines(num_pts),ierr)
  SLL_ALLOCATE(interpolator%t(num_pts+k),ierr)
  interpolator%t(1:k) = eta_min
  interpolator%t(num_pts+1:num_pts+k) = eta_max
  
  if ( mod(k,2) == 0 ) then
    do i = k+1,num_pts
      interpolator%t(i) = interpolator%eta(i-k/2) 
    end do
  else
    do i = k+1, num_pts
      interpolator%t(i) = 0.5*(interpolator%eta(i-(k-1)/2)+    &
                               interpolator%eta(i-1-(k-1)/2))
    end do
  end if
  SLL_ALLOCATE(interpolator%work(num_pts*(2*k-1)),ierr)

case default

   tmp = num_pts*num_pts
   SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)
   interpolator%slope_right = 0.0_f64
   interpolator%slope_left = 0.0_f64
   SLL_ALLOCATE(interpolator%t(num_pts+k+deriv),ierr)
   interpolator%t = 0.0_f64
   interpolator%t(1:k) = interpolator%eta(1)
   interpolator%t(num_pts+deriv+1:num_pts+deriv+k) = interpolator%eta(num_pts)
   SLL_ALLOCATE(interpolator%work((num_pts+deriv)*(2*k-1)),ierr)
    
   if (num_pts+deriv+k == num_pts+2*(k-1)) then
     interpolator%t(k+1:num_pts+deriv) = interpolator%eta(2:num_pts-1)
   else
     SLL_ERROR( this_sub_name, 'Problem with knots settings') 
   end if

end select

interpolator%coeff_splines(:) = 0.0_f64

end subroutine initialize_ad1d_interpolator


!> Initialization of the boundary for interpolator arbitrary degree splines 1d.
!> The parameters are
!> @param[in]  value_left  contains the value in the left
!> @param[in]  value_right contains the value in the right
!> @param[in]  slope_left  contains the value in the left for derivative
!> @param[in]  slope_right contains the value in the right for derivative
!> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_1d
subroutine set_values_at_boundary1d( interpolator, &
                                     value_left,   &
                                     value_right,  &
                                     slope_left,   &
                                     slope_right)

class(sll_arbitrary_degree_spline_interpolator_1d), intent(inout) :: interpolator

sll_real64, intent(in), optional :: value_left
sll_real64, intent(in), optional :: value_right
sll_real64, intent(in), optional :: slope_left
sll_real64, intent(in), optional :: slope_right

sll_int32 :: bc_left
sll_int32 :: bc_right
sll_int64 :: bc_selector

bc_left = interpolator%bc_left
bc_right= interpolator%bc_right
bc_selector = interpolator%bc_selector

if (present(value_left)) then
  interpolator%value_left = value_left
  interpolator%compute_value_left = .FALSE.
end if

if (present(value_right)) then
  interpolator%value_right = value_right
  interpolator%compute_value_right = .FALSE.
end if

if (present(slope_left)) then
  interpolator%slope_left = slope_left
  interpolator%compute_slope_left = .FALSE.
end if

if (present(slope_right)) then
  interpolator%slope_right = slope_right
  interpolator%compute_slope_right = .FALSE.
end if

end subroutine set_values_at_boundary1d

!> @brief computing the coefficients spline with a given
!>  data_array 1D cooresponding at the values of a function
!> @details 
!>  on eta_coords of size size_eta_coords
!>  if the eta_coords and eta_coords is not given
!>  we consider that the values of the function is on the points in the mesh_1d
!>
!> The parameters are
!> @param[in]  data_array the 1d arrays corresponding at the values of a function
!> @param[in]  eta_coords the 1d arrays
!> @param[in]  size_eta_coords the size of eta_coords
!> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_1d
subroutine compute_interpolants_ad1d( interpolator,    &
                                      data_array,      &
                                      eta_coords,      &
                                      size_eta_coords)

class(sll_arbitrary_degree_spline_interpolator_1d), &
            intent(inout)           :: interpolator
sll_real64, intent(in   )           :: data_array(:)
sll_real64, intent(in   ), optional :: eta_coords(:)
sll_int32,  intent(in   ), optional :: size_eta_coords

character(len=*), parameter :: this_sub_name = 'compute_interpolants_ad1d'
sll_int32       , parameter :: sz_deriv = 2
sll_int32                   :: point_locate_eta_derivative( sz_deriv )
sll_real64                  :: data_array_derivative      ( sz_deriv )
sll_int32                   :: sz
sll_real64                  :: period
sll_int32                   :: order
sll_int32                   :: ierr

if(present(eta_coords) .or. present(size_eta_coords)) then
   SLL_ERROR( this_sub_name, 'This case is not yet implemented' )
end if

sz = interpolator%num_pts

if (interpolator%compute_slope_left) then
  interpolator%slope_left = forward_fd_5pt(data_array,interpolator%eta)
end if
if (interpolator%compute_slope_right) then
  interpolator%slope_right = backward_fd_5pt(data_array,interpolator%eta,sz)
end if

if (interpolator%compute_value_left) then
  interpolator%value_left = data_array(1)
end if
if (interpolator%compute_value_right) then
  interpolator%value_right = data_array(sz)
end if

order  = interpolator%spline_degree + 1
period = interpolator%eta_max - interpolator%eta_min


select case (interpolator%bc_selector)
case (0) ! periodic

  interpolator%size_coeffs = sz 
  interpolator%size_t = order + sz 

  call splint ( interpolator%deboor,        &
                interpolator%eta,           &
                data_array,                 &
                interpolator%t,             &
                interpolator%num_pts,       &
                order,                      &
                interpolator%work,          &
                interpolator%coeff_splines, &
                ierr)
    
case (9) ! 2. dirichlet-left, dirichlet-right

  interpolator%size_coeffs = sz
  interpolator%size_t = order + sz
  interpolator%coeff_splines  = data_array
  interpolator%deboor%ilo = 1

  call splint ( interpolator%deboor,        &
                interpolator%eta,           &
                data_array,                 &
                interpolator%t,             &
                interpolator%num_pts,       &
                order,                      &
                interpolator%work,          &
                interpolator%coeff_splines, &
                ierr)
    
  interpolator%coeff_splines(1)  = interpolator%value_left
  interpolator%coeff_splines(sz) = interpolator%value_right

case default

  ! -----------------------------------
  !!! It is only for cubic spline !!!!
  ! -----------------------------------
  interpolator%size_coeffs = sz + sz_deriv
  interpolator%size_t = order + sz + sz_deriv

  point_locate_eta_derivative(1) = 1
  if (interpolator%bc_left == SLL_NEUMANN .or. &
      interpolator%bc_left == SLL_HERMITE ) then 
    data_array_derivative(1) = interpolator%slope_left
  else
    data_array_derivative(1) = 0.0_f64
  end if

  point_locate_eta_derivative(2) = sz
  if (interpolator%bc_right == SLL_NEUMANN .or. &
      interpolator%bc_right == SLL_HERMITE ) then 
    data_array_derivative(2) = interpolator%slope_right
  else
    data_array_derivative(2) = 0.0_f64
  end if

  call splint_der(                     &
       interpolator%deboor,                       &
       interpolator%eta,                          &
       data_array,                                &
       point_locate_eta_derivative,               &
       data_array_derivative,                     &
       interpolator%t,                            &
       sz,                                        &
       sz_deriv,                                  &
       order,                                     &
       interpolator%work,                         &
       interpolator%coeff_splines(1:sz+sz_deriv), &
       ierr )

  if (interpolator%bc_left == SLL_DIRICHLET) &
     interpolator%coeff_splines(1)           = interpolator%value_left
  if (interpolator%bc_right == SLL_DIRICHLET) &
     interpolator%coeff_splines(sz+sz_deriv) = interpolator%value_right

end select

end subroutine compute_interpolants_ad1d

    
  
!> @brief Interpolation on the points eta using
!> the arbitrary degree splines interpolator 1d
!> @details computing the values with the interpolator 
!> arbitrary degree splines 1d
!> on the points eta of arbitrary degree splines 1d
!> @param[in] interpolator the type sll_arbitrary_degree_spline_interpolator_1d
!> @param[in] eta1 the point
!> @return val the values on the points eta
function interpolate_value_ad1d( interpolator, eta1) result(val)

class(sll_arbitrary_degree_spline_interpolator_1d), intent(in)  :: interpolator

sll_real64, intent(in)          :: eta1
sll_real64                      :: val
sll_int32                       :: size_coeffs
sll_real64                      :: res

size_coeffs = interpolator%size_coeffs

res = eta1

select case (interpolator%bc_selector)
case (0) ! periodic

  if( res < interpolator%eta_min ) then
     res = res+interpolator%eta_max-interpolator%eta_min
  else if( res >  interpolator%eta_max ) then
     res = res+interpolator%eta_min-interpolator%eta_max
  end if

end select

SLL_ASSERT( res >= interpolator%eta_min )
SLL_ASSERT( res <= interpolator%eta_max )

val = bvalue( interpolator%deboor,          &
              interpolator%t,               &
              interpolator%coeff_splines,   &
              size_coeffs,                  &
              interpolator%spline_degree+1, &
              res,                          &
              0)

end function interpolate_value_ad1d



!> @brief initializing the coefficients of splines.
!> @details  initializing the coefficients of splines
!>  fot the arbitrary degree splines interpolator 1d
!> The parameters are
!> @param interpolator the type sll_arbitrary_degree_spline_interpolator_1d
!> @param[in] coeffs the 1d arrays corresponding of the splines coefficients
!> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_1d
subroutine set_coefficients_ad1d( interpolator, coeffs)

class(sll_arbitrary_degree_spline_interpolator_1d), intent(inout)  :: interpolator
sll_real64, dimension(:), intent(in), optional :: coeffs

sll_int32  :: sp_deg
sll_int32  :: num_cells
sll_int32  :: tmp
sll_int32  :: i
sll_real64 :: eta_min
sll_real64 :: eta_max
sll_real64 :: delta
sll_int32  :: nb_spline_eta
sll_real64 :: eta

sp_deg    = interpolator%spline_degree
num_cells = interpolator%num_pts - 1
eta_min   = interpolator%eta_min
eta_max   = interpolator%eta_max
delta     = (eta_max - eta_min)/num_cells

tmp = (sp_deg + 1)/2

! The interpretation and further filling of the spline coefficients array
! depends on the boundary conditions.
select case (interpolator%bc_selector)
case(0) ! periodic

  interpolator%size_coeffs =  num_cells + sp_deg
  interpolator%size_t = 2*sp_deg + num_cells +1

  SLL_ASSERT(size(coeffs) == num_cells+1)

  do i = -sp_deg, num_cells + sp_deg
    interpolator%t(i+sp_deg+1)=eta_min+i*delta
  end do
  do i = 1,num_cells
    interpolator%coeff_splines(i) = coeffs( i )
  end do
  do i = 1, sp_deg
    interpolator%coeff_splines(num_cells+i ) = coeffs(i )
  end do
  do i= 1,sp_deg
    interpolator%coeff_splines(num_cells+i) = interpolator%coeff_splines(sp_deg-(i-1))
  end do

case (9) ! 2. dirichlet-left, dirichlet-right

  interpolator%size_coeffs=  num_cells + sp_deg
  interpolator%size_t = 2*sp_deg + num_cells + 1
  nb_spline_eta = num_cells + sp_deg - 2
  SLL_ASSERT( size(coeffs) == nb_spline_eta ) 

  do i = 1, sp_deg + 1
    interpolator%t(i) = eta_min
  enddo
  eta = eta_min
  do i = sp_deg + 2, num_cells + 1 + sp_deg
    eta = eta + delta
    interpolator%t(i) = eta
  enddo
  do i = num_cells + sp_deg + 2, num_cells + 1 + 2*sp_deg
    interpolator%t(i) = eta
  enddo

  do i = 1,nb_spline_eta
    interpolator%coeff_splines(i + 1 ) =  coeffs(i)
  end do

  interpolator%coeff_splines(1)               = interpolator%value_left
  interpolator%coeff_splines(nb_spline_eta+2) = interpolator%value_right

case(10) ! Neumann - Dirichlet

  interpolator%size_coeffs= num_cells + sp_deg
  interpolator%size_t     = 2*sp_deg + num_cells + 1
  nb_spline_eta = num_cells + sp_deg - 1

  SLL_ASSERT( size(coeffs) == nb_spline_eta ) 

  do i = 1, sp_deg + 1
     interpolator%t(i) = eta_min
  enddo
  eta = eta_min
  do i = sp_deg + 2, num_cells + 1 + sp_deg
     eta = eta + delta
     interpolator%t(i) = eta
  enddo
  do i = num_cells + sp_deg + 2, num_cells + 1 + 2*sp_deg
     interpolator%t(i) = eta
  enddo

  do i = 1,nb_spline_eta
     interpolator%coeff_splines(i + 1 ) =  coeffs(i)
  end do
  interpolator%coeff_splines(nb_spline_eta+2) =  interpolator%value_right

case(12) ! Hermitte- Dirichlet

  interpolator%size_coeffs= num_cells + sp_deg
  interpolator%size_t     = 2*sp_deg + num_cells + 1
  nb_spline_eta = num_cells + sp_deg - 1
  SLL_ASSERT( size(coeffs) == nb_spline_eta ) 
  do i = 1, sp_deg + 1
     interpolator%t(i) = eta_min
  enddo
  eta = eta_min
  do i = sp_deg + 2, num_cells + 1 + sp_deg
     eta = eta + delta
     interpolator%t(i) = eta
  enddo
  do i = num_cells + sp_deg + 2, num_cells + 1 + 2*sp_deg
     interpolator%t(i) = eta
  enddo

  do i = 1,nb_spline_eta
     interpolator%coeff_splines(i + 1 ) =  coeffs(i)
  end do
  interpolator%coeff_splines(nb_spline_eta+2) = interpolator%value_right

case(17) ! Dirichlet-Neumann

  interpolator%size_coeffs= num_cells + sp_deg
  interpolator%size_t     = 2*sp_deg + num_cells + 1
  nb_spline_eta = num_cells + sp_deg - 1
  SLL_ASSERT( size(coeffs) == nb_spline_eta ) 

  do i = 1, sp_deg + 1
     interpolator%t(i) = eta_min
  enddo
  eta = eta_min
  do i = sp_deg + 2, num_cells + 1 + sp_deg
     eta = eta + delta
     interpolator%t(i) = eta
  enddo
  do i = num_cells + sp_deg + 2, num_cells + 1 + 2*sp_deg
     interpolator%t(i) = eta
  enddo

  do i = 1,nb_spline_eta
     interpolator%coeff_splines(i + 1 ) =  coeffs(i)
  end do
  interpolator%coeff_splines(1) = interpolator%value_left

case default

  interpolator%size_coeffs= num_cells + sp_deg
  interpolator%size_t     = 2*sp_deg + num_cells + 1
  nb_spline_eta = num_cells + sp_deg

  SLL_ASSERT( size(coeffs) == nb_spline_eta ) 

  do i = 1, sp_deg+1
    interpolator%t(i) = eta_min
  enddo
  eta = eta_min
  do i = sp_deg+2, num_cells+1+sp_deg
    eta = eta + delta
    interpolator%t(i) = eta
  enddo
  do i = num_cells+sp_deg+2, num_cells+1+2*sp_deg
     interpolator%t(i) = eta
  enddo

  do i = 1,nb_spline_eta
     interpolator%coeff_splines(i ) = coeffs(i)
  end do

end select
end subroutine set_coefficients_ad1d

  
!> @brief First derivative interpolation on the point eta
!> @details computing the values of the first derivative
!> with the interpolator arbitrary degree splines 1d
!> on the points eta of arbitrary degree splines 1d
!>
!> The parameters are
!> @param interpolator the type sll_arbitrary_degree_spline_interpolator_1d
!> @param[in] eta1 the point
!> @return val the values on the point eta of the first derivative

function interpolate_derivative_ad1d( interpolator, eta1 ) result(val)

class(sll_arbitrary_degree_spline_interpolator_1d), intent(in)  :: interpolator

sll_real64, intent(in)           :: eta1
sll_real64                       :: val
sll_real64                       :: res

res = eta1

if (interpolator%bc_selector == 0 ) then !periodic

  if( res < interpolator%eta_min ) then
    res = res+interpolator%eta_max-interpolator%eta_min
  else if( res >  interpolator%eta_max ) then
    res = res+interpolator%eta_min-interpolator%eta_max
  end if

end if

SLL_ASSERT( res >= interpolator%eta_min )
SLL_ASSERT( res <= interpolator%eta_max )

val = bvalue( interpolator%deboor,          &
              interpolator%t,               &
              interpolator%coeff_splines,   &
              interpolator%size_coeffs,     &
              interpolator%spline_degree+1, &
              res,                          &
              1)

end function interpolate_derivative_ad1d

function interpolate_array_ad1d( this,         &
                                 num_points,   &
                                 data,         &
                                 coordinates) result(res)

class(sll_arbitrary_degree_spline_interpolator_1d), intent(in) :: this

sll_int32,  intent(in)               :: num_points
sll_real64, dimension(:), intent(in) :: coordinates
sll_real64, dimension(:), intent(in) :: data
sll_real64, dimension(num_points)    :: res
sll_int32                            :: i

SLL_ASSERT(size(data) == num_points)
SLL_ASSERT(size(coordinates) == num_points)
do i = 1, num_points

  res(i) = bvalue( this%deboor,          &
                   this%t,               &
                   this%coeff_splines,   &
                   this%size_coeffs,     &
                   this%spline_degree+1, &
                   coordinates(i),       &
                   0)
end do

end function interpolate_array_ad1d

function interpolate_1d_array_disp_ad1d( &
     this,        &
     num_points, &
     data,     &
     alpha) result(res)

class(sll_arbitrary_degree_spline_interpolator_1d), intent(in)    :: this
sll_int32, intent(in)                          :: num_points
sll_real64, dimension(:), intent(in)         :: data
sll_real64, intent(in)         :: alpha
sll_real64, dimension(num_points) :: res

print *, 'interpolate_1d_array_disp_ad1d: not implemented.'
res = -1000000._f64*alpha*data*this%spline_degree

end function interpolate_1d_array_disp_ad1d

function get_coefficients_ad1d(interpolator)

class(sll_arbitrary_degree_spline_interpolator_1d), intent(in)    :: interpolator
sll_real64, dimension(:), pointer            :: get_coefficients_ad1d

get_coefficients_ad1d => interpolator%coeff_splines

end function get_coefficients_ad1d

subroutine interpolate_values_ad1d( interpolator,        &
                                    num_pts,             &
                                    vals_to_interpolate, &
                                    output_array )

class(sll_arbitrary_degree_spline_interpolator_1d),  intent(in) :: interpolator
sll_int32,  intent(in)                 :: num_pts
sll_real64, dimension(:), intent(in)   :: vals_to_interpolate
sll_real64, dimension(:), intent(out)  :: output_array
sll_int32 :: idx

SLL_ASSERT(num_pts==size(vals_to_interpolate))
do idx=1,num_pts
  output_array(idx)=interpolate_value_ad1d( &
                            interpolator, &
                            vals_to_interpolate(idx))
enddo

end subroutine interpolate_values_ad1d

subroutine interpolate_pointer_values_ad1d( &
       interpolator, &
       num_pts, &
       vals_to_interpolate, &
       output )

class(sll_arbitrary_degree_spline_interpolator_1d),  intent(in) :: interpolator
sll_int32,  intent(in)            :: num_pts
sll_real64, dimension(:), pointer :: vals_to_interpolate
sll_real64, dimension(:), pointer :: output
sll_int32 :: idx

SLL_ASSERT(num_pts==size(vals_to_interpolate))
do idx=1,num_pts
      output(idx)=interpolate_value_ad1d( &
                            interpolator, &
                            vals_to_interpolate(idx))
enddo

end subroutine interpolate_pointer_values_ad1d

subroutine interpolate_derivatives_ad1d( &
     interpolator, &
     num_pts, &
     vals_to_interpolate, &
     output_array )

  class(sll_arbitrary_degree_spline_interpolator_1d),  intent(in) :: interpolator
  sll_int32,  intent(in)                 :: num_pts
  sll_real64, dimension(:), intent(in) :: vals_to_interpolate
  sll_real64, dimension(:), intent(out) :: output_array
  sll_int32 :: idx

  SLL_ASSERT(num_pts==size(vals_to_interpolate))
  do idx=1,num_pts
        output_array(idx)=interpolate_derivative_ad1d( &
                              interpolator, &
                              vals_to_interpolate(idx))
  enddo

end subroutine interpolate_derivatives_ad1d

subroutine interpolate_pointer_derivatives_ad1d( &
     interpolator, &
     num_pts, &
     vals_to_interpolate, &
     output )

class(sll_arbitrary_degree_spline_interpolator_1d),  intent(in) :: interpolator
sll_int32,  intent(in)              :: num_pts
sll_real64, dimension(:), pointer   :: vals_to_interpolate
sll_real64, dimension(:), pointer   :: output
sll_int32 :: idx

SLL_ASSERT(num_pts==size(vals_to_interpolate))
do idx=1,num_pts
      output(idx)=interpolate_derivative_ad1d( &
                            interpolator, &
                            vals_to_interpolate(idx))
enddo

end subroutine interpolate_pointer_derivatives_ad1d

function reconstruct_array(this, num_points, data) result(res)

  class(sll_arbitrary_degree_spline_interpolator_1d),  intent(in) :: this
  sll_int32, intent(in)                :: num_points! size of output array
  sll_real64, dimension(:), intent(in) :: data   ! data to be interpolated
  sll_real64, dimension(num_points)    :: res
  res(:) = -1000000.0_f64*data*this%spline_degree
  print*, 'reconstruct_array 1d not implemented yet'

end function reconstruct_array

! The following two functions are wrong, the stencil to compute the
! derivatives is valid only for uniform spacing between the points, the
! specific coefficients chosen here using the eta coordinates is not to
! be used. ABSOLUTE FIXME!!
function forward_fd_5pt( data,eta) 

sll_real64, dimension(:), intent(in) :: data
sll_real64, dimension(:), intent(in) :: eta
sll_real64, allocatable :: res(:,:)
sll_real64 :: forward_fd_5pt

allocate(res(0:1,size(data)))

call apply_fd(5,1,eta(1:5),data(1:5), eta(1),res(0:1,1))

forward_fd_5pt = res(1,1) 

end function forward_fd_5pt

function backward_fd_5pt( data,eta,li) result(res)

sll_real64, dimension(:), intent(in) :: data
sll_real64, dimension(:), intent(in) :: eta
sll_int32, intent(in)  :: li  ! last index of the array
sll_real64 :: res

res = (0.25_f64*data(li-4)*(eta(li-5) - eta(li-4)) -&
     (4.0_f64/3.0_f64)*  data(li-3)*(eta(li-4) - eta(li-3)) + &
      3.0_f64*           data(li-2)*(eta(li-3) - eta(li-2)) - &
      4.0_f64*           data(li-1)*(eta(li-2) - eta(li-1)) + &
     (25.0_f64/12.0_f64)*data(li)*  (eta(li-1) - eta(li)) )

end function backward_fd_5pt


end module sll_m_arbitrary_degree_spline_interpolator_1d
