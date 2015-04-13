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

!> Class of arbitrary degree version of 2d irnterpolator
module sll_module_arbitrary_degree_spline_interpolator_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h" 
use sll_module_interpolators_2d_base
use sll_utilities
use sll_module_deboor_splines_1d, only: bsplvb, bsplvd, interv
use sll_module_deboor_splines_2d, only: dvalue2d, bvalue2d

use sll_module_arbitrary_degree_spline_interpolator_1d
implicit none
private

! in what follows, the direction '1' is in the contiguous memory direction.
!> Arbitrary degree version of 2d irnterpolator
type, extends(sll_interpolator_2d_base) :: &
  sll_arbitrary_degree_spline_interpolator_2d           

  sll_int32                           :: num_pts1
  sll_int32                           :: num_pts2
  sll_real64                          :: eta1_min
  sll_real64                          :: eta1_max
  sll_real64                          :: eta2_min
  sll_real64                          :: eta2_max
  sll_int32                           :: bc_min1
  sll_int32                           :: bc_max1
  sll_int32                           :: bc_min2
  sll_int32                           :: bc_max2
  sll_int32                           :: spline_degree1
  sll_int32                           :: spline_degree2
  sll_real64, dimension(:), pointer   :: knots1
  sll_real64, dimension(:), pointer   :: knots2
  sll_real64, dimension(:), pointer   :: t1
  sll_real64, dimension(:), pointer   :: t2
  sll_int32                           :: size_t1
  sll_int32                           :: size_t2 
  sll_int64                           :: bc_selector 
  sll_real64, dimension(:,:), pointer :: coeff_splines
  sll_int32                           :: size_coeffs1
  sll_int32                           :: size_coeffs2
  logical                             :: coefficients_set = .false.
  sll_real64, dimension(:),pointer    :: slope_min1
  sll_real64, dimension(:),pointer    :: slope_max1
  sll_real64, dimension(:),pointer    :: slope_min2
  sll_real64, dimension(:),pointer    :: slope_max2
  sll_real64, dimension(:),pointer    :: value_min1
  sll_real64, dimension(:),pointer    :: value_max1
  sll_real64, dimension(:),pointer    :: value_min2
  sll_real64, dimension(:),pointer    :: value_max2

  logical :: compute_slope_min1 = .TRUE.
  logical :: compute_slope_max1 = .TRUE.
  logical :: compute_slope_max2 = .TRUE.
  logical :: compute_slope_min2 = .TRUE.
  logical :: compute_value_min1 = .TRUE.
  logical :: compute_value_max1 = .TRUE.
  logical :: compute_value_max2 = .TRUE.
  logical :: compute_value_min2 = .TRUE.

contains

  procedure :: initialize                  => initialize_ad2d_interpolator
  procedure :: set_coefficients            => set_coefficients_ad2d
  procedure :: coefficients_are_set        => coefficients_are_set_ad2d
  procedure :: compute_interpolants        => compute_interpolants_ad2d
  procedure :: interpolate_value           => interpolate_value_ad2d
  procedure :: interpolate_derivative_eta1 => interpolate_derivative1_ad2d
  procedure :: interpolate_derivative_eta2 => interpolate_derivative2_ad2d
  procedure :: interpolate_array           => interpolate_array_ad2d
  procedure :: interpolate_array_disp      => interpolate_2d_array_disp_ad2d
  procedure :: get_coefficients            => get_coefficients_ad2d
  procedure :: delete                      => delete_arbitrary_degree_2d_interpolator
  procedure :: set_values_at_boundary      => set_boundary_value2d
  procedure :: set_slopes_at_boundary      => set_slope2d

end type sll_arbitrary_degree_spline_interpolator_2d


!> Pointer to arbitrary degree version of 1d interpolator
type sll_arbitrary_degree_spline_interpolator_2d_ptr
   type(sll_arbitrary_degree_spline_interpolator_2d), pointer :: interp
end type sll_arbitrary_degree_spline_interpolator_2d_ptr


!> Deallocate the interpolator class
interface sll_delete
   module procedure delete_arbitrary_degree_2d_interpolator
end interface sll_delete

public sll_arbitrary_degree_spline_interpolator_2d           
public sll_arbitrary_degree_spline_interpolator_2d_ptr
public sll_delete
public new_arbitrary_degree_spline_interp2d
public set_slope2d
public initialize_ad2d_interpolator

contains

!> Delete interpolator arbitrary degree splines.
!> The parameters are
!> @param interpolator the type sll_arbitrary_degree_spline_interpolator_2d
!
subroutine delete_arbitrary_degree_2d_interpolator( interpolator )

class(sll_arbitrary_degree_spline_interpolator_2d), intent(inout) :: interpolator
sll_int32 :: ierr

SLL_DEALLOCATE(interpolator%knots1,ierr)
SLL_DEALLOCATE(interpolator%knots2,ierr)
SLL_DEALLOCATE(interpolator%t1,ierr)
SLL_DEALLOCATE(interpolator%t2,ierr)
SLL_DEALLOCATE(interpolator%coeff_splines,ierr)
SLL_DEALLOCATE(interpolator%value_min1,ierr)
SLL_DEALLOCATE(interpolator%value_max1,ierr)
SLL_DEALLOCATE(interpolator%value_min2,ierr)
SLL_DEALLOCATE(interpolator%value_max2,ierr)
SLL_DEALLOCATE(interpolator%slope_min1,ierr)
SLL_DEALLOCATE(interpolator%slope_max1,ierr)
SLL_DEALLOCATE(interpolator%slope_min2,ierr)
SLL_DEALLOCATE(interpolator%slope_max2,ierr)

end subroutine delete_arbitrary_degree_2d_interpolator

!> @brief Initialization of a pointer interpolator arbitrary degree splines 2d.
!> @details To have the interpolator arbitrary degree splines 2d
!> 
!> The parameters are
!> @param[in] num_pts1 the number of points in the direction eta1
!> @param[in] num_pts2 the number of points in the direction eta2
!> @param[in] eta1_min the minimun in the direction eta1
!> @param[in] eta1_max the maximun in the direction eta1
!> @param[in] eta2_min the minimun in the direction eta2
!> @param[in] eta2_max the maximun in the direction eta2
!> @param[in] bc_min1  the boundary condition at left in the direction eta1
!> @param[in] bc_max1 the boundary condition at right in the direction eta2
!> @param[in] bc_min2 the boundary condition at left in the direction eta2
!> @param[in] bc_max2 the boundary condition at right in the direction eta2
!> @param[in] spline_degree1 the degree of B-spline in the direction eta1
!> @param[in] spline_degree2 the degre of B-spline in the direction eta2
!> @return the type sll_arbitrary_degree_spline_interpolator_2d

function new_arbitrary_degree_spline_interp2d( num_pts1,       &
                                               num_pts2,       &
                                               eta1_min,       &
                                               eta1_max,       &
                                               eta2_min,       &
                                               eta2_max,       &
                                               bc_min1,        &
                                               bc_max1,        &
                                               bc_min2,        &
                                               bc_max2,        &
                                               spline_degree1, &
                                               spline_degree2) result( res )

type(sll_arbitrary_degree_spline_interpolator_2d), pointer :: res

sll_int32,  intent(in) :: num_pts1
sll_int32,  intent(in) :: num_pts2
sll_real64, intent(in) :: eta1_min
sll_real64, intent(in) :: eta1_max
sll_real64, intent(in) :: eta2_min
sll_real64, intent(in) :: eta2_max
sll_int32,  intent(in) :: bc_min1
sll_int32,  intent(in) :: bc_max1
sll_int32,  intent(in) :: bc_min2
sll_int32,  intent(in) :: bc_max2
sll_int32,  intent(in) :: spline_degree1
sll_int32,  intent(in) :: spline_degree2

sll_int32 :: ierr

SLL_ALLOCATE(res,ierr)

call initialize_ad2d_interpolator( res,            &
                                   num_pts1,       &
                                   num_pts2,       &
                                   eta1_min,       &
                                   eta1_max,       &
                                   eta2_min,       &
                                   eta2_max,       &
                                   bc_min1,        &
                                   bc_max1,        &
                                   bc_min2,        &
                                   bc_max2,        &
                                   spline_degree1, &
                                   spline_degree2)

end function new_arbitrary_degree_spline_interp2d

! -----------------------------------------------
! This subroutine allocate the type of interpolator
!    the  arbitrary_spline_interp2d
! -----------------------------------------------
!> Initialization of an interpolator arbitrary degree splines 2d.
!> The parameters are
!> @param[in] num_pts1 the number of points in the direction eta1
!> @param[in] num_pts2 the number of points in the direction eta2
!> @param[in] eta1_min the minimun in the direction eta1
!> @param[in] eta1_max the maximun in the direction eta1
!> @param[in] eta2_min the minimun in the direction eta2
!> @param[in] eta2_max the maximun in the direction eta2
!> @param[in] bc_min1  the boundary condition at left in the direction eta1
!> @param[in] bc_max1 the boundary condition at right in the direction eta2
!> @param[in] bc_min2 the boundary condition at left in the direction eta2
!> @param[in] bc_max2 the boundary condition at right in the direction eta2
!> @param[in] spline_degree1 the degree of B-spline in the direction eta1
!> @param[in] spline_degree2 the degre of B-spline in the direction eta2
!> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_2d
subroutine initialize_ad2d_interpolator( interpolator,   &
                                         num_pts1,       &
                                         num_pts2,       &
                                         eta1_min,       &
                                         eta1_max,       &
                                         eta2_min,       &
                                         eta2_max,       &
                                         bc_min1,        &
                                         bc_max1,        &
                                         bc_min2,        &
                                         bc_max2,        &
                                         spline_degree1, &
                                         spline_degree2)

class(sll_arbitrary_degree_spline_interpolator_2d):: interpolator

sll_int32,  intent(in)  :: num_pts1
sll_int32,  intent(in)  :: num_pts2
sll_real64, intent(in)  :: eta1_min
sll_real64, intent(in)  :: eta1_max
sll_real64, intent(in)  :: eta2_min
sll_real64, intent(in)  :: eta2_max
sll_int32,  intent(in)  :: bc_min1
sll_int32,  intent(in)  :: bc_max1
sll_int32,  intent(in)  :: bc_min2
sll_int32,  intent(in)  :: bc_max2
sll_int32,  intent(in)  :: spline_degree1
sll_int32,  intent(in)  :: spline_degree2

sll_int32 :: ierr
sll_int32 :: tmp1
sll_int32 :: tmp2
sll_int64 :: bc_selector


! do some argument checking...
if(((bc_min1  == SLL_PERIODIC).and.(bc_max1.ne. SLL_PERIODIC)).or.&
   ((bc_max1 == SLL_PERIODIC).and.(bc_min1 .ne. SLL_PERIODIC)))then
   print *, 'initialize_arbitrary_degree_2d_interpolator, ERROR: ', &
        'if one boundary condition is specified as periodic, then ', &
        'both must be. Error in first direction.'
end if

if(((bc_min2 == SLL_PERIODIC).and.(bc_max2.ne. SLL_PERIODIC)).or.&
   ((bc_max2 == SLL_PERIODIC).and.(bc_min2 .ne. SLL_PERIODIC)))then
   print *, 'initialize_arbitrary_degree_2d_interpolator, ERROR: ', &
        'if one boundary condition is specified as periodic, then ', &
        'both must be. Error in second direction.'
end if

bc_selector = 0

if( bc_min1 == SLL_DIRICHLET ) bc_selector = bc_selector + 1
if( bc_min1 == SLL_NEUMANN   ) bc_selector = bc_selector + 2
if( bc_min1 == SLL_HERMITE   ) bc_selector = bc_selector + 4
if( bc_max1 == SLL_DIRICHLET ) bc_selector = bc_selector + 8
if( bc_max1 == SLL_NEUMANN   ) bc_selector = bc_selector + 16
if( bc_max1 == SLL_HERMITE   ) bc_selector = bc_selector + 32
if( bc_min2 == SLL_DIRICHLET ) bc_selector = bc_selector + 64
if( bc_min2 == SLL_NEUMANN   ) bc_selector = bc_selector + 128
if( bc_min2 == SLL_HERMITE   ) bc_selector = bc_selector + 256
if( bc_max2 == SLL_DIRICHLET ) bc_selector = bc_selector + 512
if( bc_max2 == SLL_NEUMANN   ) bc_selector = bc_selector + 1024
if( bc_max2 == SLL_HERMITE   ) bc_selector = bc_selector + 2048

interpolator%spline_degree1 = spline_degree1
interpolator%spline_degree2 = spline_degree2
interpolator%eta1_min       = eta1_min
interpolator%eta1_max       = eta1_max
interpolator%eta2_min       = eta2_min
interpolator%eta2_max       = eta2_max
interpolator%bc_min1        = bc_min1
interpolator%bc_max1        = bc_max1
interpolator%bc_min2        = bc_min2
interpolator%bc_max2        = bc_max2
interpolator%bc_selector    = bc_selector
interpolator%num_pts1       = num_pts1
interpolator%num_pts2       = num_pts2

SLL_CLEAR_ALLOCATE(interpolator%value_min1(1:num_pts2  ),ierr)
SLL_CLEAR_ALLOCATE(interpolator%value_max1(1:num_pts2  ),ierr)
SLL_CLEAR_ALLOCATE(interpolator%value_min2(1:num_pts1+2),ierr)
SLL_CLEAR_ALLOCATE(interpolator%value_max2(1:num_pts1+2),ierr)
SLL_CLEAR_ALLOCATE(interpolator%slope_min1(1:num_pts2  ),ierr)
SLL_CLEAR_ALLOCATE(interpolator%slope_max1(1:num_pts2  ),ierr)
SLL_CLEAR_ALLOCATE(interpolator%slope_min2(1:num_pts1+2),ierr)
SLL_CLEAR_ALLOCATE(interpolator%slope_max2(1:num_pts1+2),ierr)

! tmp1 and tmp2 is the maximun (not absolue) for the size of coefficients
select case (bc_selector)
case (0) ! 1. periodic-periodic
   
  SLL_ALLOCATE( interpolator%knots1(2*spline_degree1+2),ierr )
  SLL_ALLOCATE( interpolator%knots2(2*spline_degree2+2),ierr )

case (9) ! 2. dirichlet-left, dirichlet-right, periodic

  SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
  SLL_ALLOCATE( interpolator%knots2(2*spline_degree2+2),ierr )

case (576) ! 3. periodic, dirichlet-bottom, dirichlet-top

  SLL_ALLOCATE( interpolator%knots1(2*spline_degree1+2),ierr )
  SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )

case default

  SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
  SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )

end select

tmp1 = num_pts1+ 4*spline_degree1
tmp2 = num_pts2+ 4*spline_degree2
SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

! knots and coeff splines allocations 
interpolator%coeff_splines(:,:) = 0.0_f64
! the minimun is to be of class C^0 everywhere on the knots
! i.e. each knot have multiplicity (spline_degree1+1) 
! so the maximun number of knots is num_pts1*(spline_degree1+1)
SLL_CLEAR_ALLOCATE( interpolator%t1(1:num_pts1*(spline_degree1+1)),ierr)
SLL_CLEAR_ALLOCATE( interpolator%t2(1:num_pts2*(spline_degree2+1)),ierr) 

end subroutine !initialize_ad2d_interpolator

subroutine set_coeff_splines_values_1d( values,        &
                                        num_pts,       &
                                        eta_min,       &
                                        eta_max,       &
                                        bc_min,        &
                                        bc_max,        &
                                        spline_degree )

sll_int32,  intent(in)    :: num_pts
sll_real64, intent(in)    :: eta_min
sll_real64, intent(in)    :: eta_max
sll_int32,  intent(in)    :: bc_min
sll_int32,  intent(in)    :: bc_max
sll_int32,  intent(in)    :: spline_degree
sll_real64, intent(inout) :: values(num_pts)

class(sll_arbitrary_degree_spline_interpolator_1d), pointer :: interp1d => null()

interp1d => new_arbitrary_degree_1d_interpolator( num_pts,      &
                                                  eta_min,      &
                                                  eta_max,      &
                                                  bc_min,       &
                                                  bc_max,       &
                                                  spline_degree )

if (bc_min == SLL_DIRICHLET .and. bc_max == SLL_DIRICHLET) then
  call set_values_at_boundary1d(interp1d, values(1), values(num_pts))
end if

call interp1d%compute_interpolants(values)

values = interp1d%coeff_splines(1:num_pts)

call sll_delete(interp1d)

end subroutine set_coeff_splines_values_1d


!> Initialization of the boundary for interpolator arbitrary degree splines 2d.
!> The parameters are
!> @param[in] value_min1 a 1d arrays contains values in the left in the direction eta1  
!> @param[in] value_max1 a 1d arrays contains values in the right in the direction eta1 
!> @param[in] value_min2 a 1d arrays contains values in the left in the direction eta2 
!> @param[in]  value_max2 a 1d arrays contains values in the right in the direction eta2
!> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_2d
subroutine set_boundary_value2d( interpolator, &
                                 value_min1,   &
                                 value_max1,   &
                                 value_min2,   &
                                 value_max2)

class(sll_arbitrary_degree_spline_interpolator_2d) :: interpolator

sll_real64, dimension(:), optional :: value_min1
sll_real64, dimension(:), optional :: value_max1
sll_real64, dimension(:), optional :: value_min2
sll_real64, dimension(:), optional :: value_max2

sll_int32  :: num_pts1
sll_int32  :: num_pts2
sll_int32  :: bc_min1
sll_int32  :: bc_max1
sll_int32  :: bc_min2
sll_int32  :: bc_max2
sll_int32  :: spline_degree1
sll_int32  :: spline_degree2
sll_real64 :: eta1_min
sll_real64 :: eta1_max
sll_real64 :: eta2_min
sll_real64 :: eta2_max

num_pts1       = interpolator%num_pts1
num_pts2       = interpolator%num_pts2
bc_min1        = interpolator%bc_min1 
bc_max1        = interpolator%bc_max1 
bc_min2        = interpolator%bc_min2  
bc_max2        = interpolator%bc_max2
spline_degree1 = interpolator%spline_degree1
spline_degree2 = interpolator%spline_degree2
eta1_min       = interpolator%eta1_min
eta1_max       = interpolator%eta1_max
eta2_min       = interpolator%eta2_min
eta2_max       = interpolator%eta2_max

if (bc_min1==SLL_DIRICHLET) then 

  if (present(value_min1)) then 
  
    call set_coeff_splines_values_1d( value_min1, &
         num_pts2,                                &
         eta2_min,                                &
         eta2_max,                                &
         bc_min2,                                 &
         bc_max2,                                 &
         spline_degree2 )
      
    interpolator%value_min1(1:num_pts2) = value_min1
    interpolator%compute_value_min1 = .FALSE.
  else
    interpolator%value_min1 = 0.0_f64
  end if

end if
  
if (bc_max1==SLL_DIRICHLET) then 

  if (present(value_max1)) then 
    call set_coeff_splines_values_1d( value_max1, &
         num_pts2,                                &
         eta2_min,                                &
         eta2_max,                                &
         bc_min2,                                 &
         bc_max2,                                 &
         spline_degree2 )
      
    interpolator%value_max1(1:num_pts2) = value_max1
    interpolator%compute_value_max1 = .FALSE.
  else
    interpolator%value_max1 = 0.0_f64
  end if

end if

if (bc_min2==SLL_DIRICHLET) then 

  if (present(value_min2)) then 
    call set_coeff_splines_values_1d( value_min2, &
                                      num_pts1,   &
                                      eta1_min,   &
                                      eta1_max,   &
                                      bc_min1,    &
                                      bc_max1,    &
                                      spline_degree1 )
      
    interpolator%value_min2(1:num_pts1) = value_min2
    interpolator%compute_value_min2 = .FALSE.
  else
    interpolator%value_min2 = 0.0_f64
  end if

end if
  
if (bc_max2==SLL_DIRICHLET) then 

  if (present(value_max2)) then 
    call set_coeff_splines_values_1d( value_max2,    &
                                      num_pts1,      &
                                      eta1_min,      &
                                      eta1_max,      &
                                      bc_min1,       &
                                      bc_max1,       &
                                      spline_degree1 )
     
    interpolator%value_max2(1:num_pts1) = value_max2
    interpolator%compute_value_max2 = .FALSE.
  else
    interpolator%value_max2 = 0.0_f64
  end if

end if
    
end subroutine set_boundary_value2d


!> @brief initializing the coefficients of splines.
!> @details  initializing the coefficients of splines
!>  in the cas of linearization of them i.e. if we have 
!>  a table in 1d corresponding of coefficients 2d.
!>  The interpretation and further filling of the spline coefficients array
!>  depends on the boundary conditions.
!> 
!> The parameters are
!> @param interpolator the type sll_arbitrary_degree_spline_interpolator_2d
!> @param[in]  coeffs_1d the 1d arrays corresponding of the splines coefficients
!> @param[in]  coeffs_2d the 2d arrays corresponding of the splines coefficients
!> @param[in]  coeff2d_size1 the number of rows of coeffs_2d
!> @param[in]  coeff2d_size2 the number of columns of coeffs_2d
!> @param[in]  knots1 the knots in the direction eta1
!> @param[in]  size_knots1 the size of knots in the direction eta1
!> @param[in]  knots2  the knots in the direction eta2
!> @param[in]  size_knots2 the size of knots in the direction eta2
!> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_2d
subroutine set_coefficients_ad2d( interpolator,  &
                                  coeffs_1d,     &
                                  coeffs_2d,     &
                                  coeff2d_size1, &
                                  coeff2d_size2, &
                                  knots1,        &
                                  size_knots1,   &
                                  knots2,        &
                                  size_knots2)

class(sll_arbitrary_degree_spline_interpolator_2d), intent(inout)  :: interpolator

sll_real64, dimension(:)  , intent(in), optional :: coeffs_1d
sll_real64, dimension(:,:), intent(in), optional :: coeffs_2d
sll_int32,                  intent(in), optional :: coeff2d_size1
sll_int32,                  intent(in), optional :: coeff2d_size2
sll_real64, dimension(:),   intent(in), optional :: knots1
sll_real64, dimension(:),   intent(in), optional :: knots2
sll_int32,                  intent(in), optional :: size_knots1
sll_int32,                  intent(in), optional :: size_knots2

sll_int32   :: sp_deg1
sll_int32   :: sp_deg2
sll_int32   :: num_cells1
sll_int32   :: num_cells2
sll_int32   :: i, j
sll_real64  :: eta1_min, eta1_max
sll_real64  :: eta2_min, eta2_max
sll_real64  :: delta1
sll_real64  :: delta2
sll_int32   :: nb_spline_eta1
sll_int32   :: nb_spline_eta2
sll_real64  :: eta1
sll_real64  :: eta2
sll_int32, parameter :: sz_derivative1=2
sll_int32, parameter :: sz_derivative2=2

sp_deg1    = interpolator%spline_degree1
sp_deg2    = interpolator%spline_degree2
num_cells1 = interpolator%num_pts1 - 1
num_cells2 = interpolator%num_pts2 - 1
eta1_min   = interpolator%eta1_min
eta2_min   = interpolator%eta2_min
eta1_max   = interpolator%eta1_max
eta2_max   = interpolator%eta2_max
delta1     = (eta1_max-eta1_min)/num_cells1
delta2     = (eta2_max-eta2_min)/num_cells2


if (present(coeffs_1d)) then 

  select case (interpolator%bc_selector)
  case(0) ! periodic-periodic
     
    interpolator%size_coeffs1 =  num_cells1 + sp_deg1 + 1
    interpolator%size_coeffs2 =  num_cells2 + sp_deg2 + 1 
    interpolator%size_t1      =  2*sp_deg1 + num_cells1 +1 +1
    interpolator%size_t2      =  2*sp_deg2 + num_cells2 +1 +1
    
    SLL_ASSERT (size(coeffs_1d,1) == num_cells1*num_cells2)
    
    do i = -sp_deg1, num_cells1 + sp_deg1 + 1
      interpolator%t1(i+sp_deg1+1) = eta1_min + i*delta1
    end do
    do i = -sp_deg2, num_cells2 + sp_deg2 + 1
      interpolator%t2(i+sp_deg2+1) = eta2_min + i*delta2
    end do
    
    do j = 1,num_cells2
    do i = 1,num_cells1
      interpolator%coeff_splines(i,j) = coeffs_1d(i+num_cells1 *(j-1))
    end do
    end do
    do j = 1,sp_deg2+1
    do i = 1,num_cells1
      interpolator%coeff_splines(i,num_cells2+j) = coeffs_1d(i+num_cells1*(j-1))
    end do
    end do
    do j = 1,num_cells2
    do i = 1,sp_deg1+1
      interpolator%coeff_splines(num_cells1+i,j) = coeffs_1d(i+num_cells1*(j-1))
    end do
    end do
    do j=1,sp_deg2+1
    do i=1,sp_deg1+1
      interpolator%coeff_splines(num_cells1+i,num_cells2+j) = interpolator%coeff_splines(i,j)
    end do
    end do

  case (9) ! 2. dirichlet-left, dirichlet-right, periodic
     
    interpolator%size_coeffs1 =  num_cells1 + sp_deg1
    interpolator%size_coeffs2 =  num_cells2 + sp_deg2 + 1
    interpolator%size_t1      =  2*sp_deg1 + num_cells1 + 1
    interpolator%size_t2      =  2*sp_deg2 + num_cells2 + 1 + 1
    nb_spline_eta1            =  num_cells1 + sp_deg1 - 2
    nb_spline_eta2            =  num_cells2
    
    SLL_ASSERT(size(coeffs_1d,1)==(num_cells1+sp_deg1-2)*num_cells2)
    do i = - sp_deg2, num_cells2 + sp_deg2 + 1
      interpolator%t2( i+ sp_deg2 + 1 ) = eta2_min + i* delta2
    end do
    do i = 1, sp_deg1 + 1
      interpolator%t1(i) = eta1_min
    enddo
    eta1 = eta1_min
    do i = sp_deg1 + 2, num_cells1 + 1 + sp_deg1
      eta1 = eta1 + delta1
      interpolator%t1(i) = eta1
    enddo
    do i = num_cells1 + sp_deg1 + 2, num_cells1 + 1 + 2*sp_deg1
      interpolator%t1(i) = eta1
    enddo
    do j = 1,nb_spline_eta2
    do i = 1 ,nb_spline_eta1
      interpolator%coeff_splines(i+1,j) = coeffs_1d(i+nb_spline_eta1*(j-1))
    end do
    end do
    do j = 1, sp_deg2 + 1
    do i = 1,nb_spline_eta1
      interpolator%coeff_splines(i+1,nb_spline_eta2+j) = coeffs_1d(i+nb_spline_eta1*(j-1))
    end do
    end do
    interpolator%coeff_splines(1,:) = 0.0_f64
    interpolator%coeff_splines(nb_spline_eta1+2,:) = 0.0_f64
    
  case(576)!3. periodic, dirichlet-bottom, dirichlet-top
     
    interpolator%size_coeffs1 =  num_cells1 + sp_deg1 + 1
    interpolator%size_coeffs2 =  num_cells2 + sp_deg2
    interpolator%size_t1      = 2*sp_deg1 + num_cells1 + 1 + 1
    interpolator%size_t2      = 2*sp_deg2 + num_cells2 + 1
    nb_spline_eta1            = num_cells1
    nb_spline_eta2            = num_cells2 + sp_deg2 - 2
   
    SLL_ASSERT(size(coeffs_1d,1)==num_cells1*(num_cells2+sp_deg2-2))
    do i = - sp_deg1, nb_spline_eta1 + sp_deg1 + 1
      interpolator%t1(i+sp_deg1+1) = eta1_min + i* delta1
    end do
    do i = 1, sp_deg2 + 1
      interpolator%t2(i) = eta2_min
    enddo
    eta2 = eta2_min
    do i = sp_deg2+2, num_cells2+1+sp_deg2
      eta2 = eta2+delta2
      interpolator%t2(i) = eta2
    enddo
    do i = num_cells2 + sp_deg2+1, num_cells2+1+2*sp_deg2
      interpolator%t2(i) = eta2_max
    enddo
    do j = 1,nb_spline_eta2
    do i = 1 , nb_spline_eta1
      interpolator%coeff_splines(i,j+1) = coeffs_1d(i+nb_spline_eta1*(j-1))
    end do
    end do
    do j = 1,nb_spline_eta2
    do i = 1, sp_deg1 + 1
      interpolator%coeff_splines(nb_spline_eta1+i,j+1) = coeffs_1d(i+nb_spline_eta1*(j-1) )
    end do
    end do
      
    interpolator%coeff_splines(:,1) = 0.0_f64
    interpolator%coeff_splines(:,nb_spline_eta2+2) = 0.0_f64
   
  case(585) ! 4. dirichlet in all sides

    interpolator%size_coeffs1=  num_cells1 + sp_deg1
    interpolator%size_coeffs2=  num_cells2 + sp_deg2
    interpolator%size_t1 = 2*sp_deg1 + num_cells1 + 1
    interpolator%size_t2 = 2*sp_deg2 + num_cells2 + 1
    nb_spline_eta1 = num_cells1 + sp_deg1 - 2
    nb_spline_eta2 = num_cells2 + sp_deg2 - 2
    
    SLL_ASSERT(size(coeffs_1d,1)==(num_cells1 + sp_deg1-2)*(num_cells2+sp_deg2-2))
    do i = 1, sp_deg1 + 1
      interpolator%t1(i) = eta1_min
    enddo
    eta1 = eta1_min
    do i = sp_deg1 + 2, num_cells1 + 1 + sp_deg1
      eta1 = eta1 + delta1
      interpolator%t1(i) = eta1
    enddo
    do i = num_cells1 + sp_deg1 + 2, num_cells1 + 1 + 2*sp_deg1
      interpolator%t1(i) = eta1
    enddo
    do i = 1, sp_deg2 + 1
      interpolator%t2(i) = eta2_min
    enddo
    eta2 = eta2_min
    do i = sp_deg2 + 2, num_cells2 + 1 + sp_deg2
      eta2 = eta2 + delta2
      interpolator%t2(i) = eta2
    enddo
    do i = num_cells2 + sp_deg2 + 2, num_cells2 + 1 + 2*sp_deg2
      interpolator%t2(i) = eta2
    enddo
    interpolator%coeff_splines(:,:) = 0.0_f64
    do j = 1,nb_spline_eta2
    do i = 1,nb_spline_eta1
      interpolator%coeff_splines(i+1,j+1) = coeffs_1d(i+nb_spline_eta1*(j-1))
    end do
    end do

  case (650) !left: Neumann, right: Dirichlet, bottom: Neumann, Top: Dirichlet

    interpolator%size_coeffs1 = num_cells1 + sp_deg1 +1
    interpolator%size_coeffs2 = num_cells2 + sp_deg2 +1
    interpolator%size_t1      = 2*sp_deg1 + num_cells1 + 1
    interpolator%size_t2      = 2*sp_deg2 + num_cells2 + 1
    nb_spline_eta1            = num_cells1 + sp_deg1 +1
    nb_spline_eta2            = num_cells2 + sp_deg2 +1
    
    SLL_ASSERT(size(coeffs_1d,1)==(num_cells1+sp_deg1+1)*(num_cells2+sp_deg2+1))
    do i = 1, sp_deg1 + 1
      interpolator%t1(i) = eta1_min
    enddo
    eta1 = eta1_min
    do i = sp_deg1 + 2, num_cells1 + 1 + sp_deg1
      eta1 = eta1 + delta1
      interpolator%t1(i) = eta1
    enddo
    do i = num_cells1 + sp_deg1 + 2, num_cells1 + 1 + 2*sp_deg1
      interpolator%t1(i) = eta1
    enddo
    do i = 1, sp_deg2 + 1
       interpolator%t2(i) = eta2_min
    enddo
    eta2 = eta2_min
    do i = sp_deg2 + 2, num_cells2 + 1 + sp_deg2
      eta2 = eta2 + delta2
      interpolator%t2(i) = eta2
    enddo
    do i = num_cells2 + sp_deg2 + 2, num_cells2 + 1 + 2*sp_deg2
      interpolator%t2(i) = eta2
    enddo
    interpolator%coeff_splines(:,:) = 0.0_f64
    do j = 1,nb_spline_eta2
    do i = 1,nb_spline_eta1
      interpolator%coeff_splines(i+1,j+1) = coeffs_1d(i+nb_spline_eta1*(j-1))
    end do
    end do
    
  case(657) !left: Dirichlet, right: Neumann, bottom: Neumann, Top: Dirichlet 

    interpolator%size_coeffs1 =  num_cells1 + sp_deg1 +1
    interpolator%size_coeffs2 =  num_cells2 + sp_deg2 +1
    interpolator%size_t1      = 2*sp_deg1 + num_cells1 + 1
    interpolator%size_t2      = 2*sp_deg2 + num_cells2 + 1
    nb_spline_eta1            = num_cells1 + sp_deg1 +1
    nb_spline_eta2            = num_cells2 + sp_deg2 +1
    
    SLL_ASSERT(size(coeffs_1d,1)==(num_cells1+sp_deg1+1)*(num_cells2+sp_deg2+1))
    do i = 1, sp_deg1 + 1
      interpolator%t1(i) = eta1_min
    enddo
    eta1 = eta1_min
    do i = sp_deg1 + 2, num_cells1 + 1 + sp_deg1
      eta1 = eta1 + delta1
      interpolator%t1(i) = eta1
    enddo
    do i = num_cells1 + sp_deg1 + 2, num_cells1 + 1 + 2*sp_deg1
      interpolator%t1(i) = eta1
    enddo
    do i = 1, sp_deg2 + 1
      interpolator%t2(i) = eta2_min
    enddo
    eta2 = eta2_min
    do i = sp_deg2 + 2, num_cells2 + 1 + sp_deg2
      eta2 = eta2 + delta2
      interpolator%t2(i) = eta2
    enddo
    do i = num_cells2 + sp_deg2 + 2, num_cells2 + 1 + 2*sp_deg2
      interpolator%t2(i) = eta2
    enddo
    interpolator%coeff_splines(:,:) = 0.0_f64
    do j = 1,nb_spline_eta2
    do i = 1,nb_spline_eta1
      interpolator%coeff_splines(i+1,j+1) = coeffs_1d(i+nb_spline_eta1*(j-1))
    end do
    end do
    
  case(780)  !left: Hermite, right: Dirichlet, bottom: Hermite, Top: Dirichlet

    interpolator%size_coeffs1 = num_cells1 + sp_deg1 +1
    interpolator%size_coeffs2 = num_cells2 + sp_deg2 +1
    interpolator%size_t1      = 2*sp_deg1 + num_cells1 + 1
    interpolator%size_t2      = 2*sp_deg2 + num_cells2 + 1
    nb_spline_eta1            = num_cells1 + sp_deg1 +1
    nb_spline_eta2            = num_cells2 + sp_deg2 +1
    
    SLL_ASSERT(size(coeffs_1d,1)==(num_cells1 + sp_deg1+1)*(num_cells2+sp_deg2+1))
    do i = 1, sp_deg1 + 1
      interpolator%t1(i) = eta1_min
    enddo
    eta1 = eta1_min
    do i = sp_deg1 + 2, num_cells1 + 1 + sp_deg1
      eta1 = eta1 + delta1
      interpolator%t1(i) = eta1
    enddo
    do i = num_cells1 + sp_deg1 + 2, num_cells1 + 1 + 2*sp_deg1
      interpolator%t1(i) = eta1
    enddo
    do i = 1, sp_deg2 + 1
      interpolator%t2(i) = eta2_min
    enddo
    eta2 = eta2_min
    do i = sp_deg2 + 2, num_cells2 + 1 + sp_deg2
      eta2 = eta2 + delta2
      interpolator%t2(i) = eta2
    enddo
    do i = num_cells2 + sp_deg2 + 2, num_cells2 + 1 + 2*sp_deg2
      interpolator%t2(i) = eta2
    enddo
    interpolator%coeff_splines(:,:) = 0.0_f64
    do j = 1,nb_spline_eta2
    do i = 1,nb_spline_eta1
      interpolator%coeff_splines(i+1,j+1) = coeffs_1d(i+nb_spline_eta1*(j-1))
    end do
    end do
    
  case(1098)  !left: Neumann, right: Dirichlet, bottom: Dirichlet, Top: Neumann

    interpolator%size_coeffs1 = num_cells1 + sp_deg1+1
    interpolator%size_coeffs2 = num_cells2 + sp_deg2+1
    interpolator%size_t1      = 2*sp_deg1  + num_cells1+1
    interpolator%size_t2      = 2*sp_deg2  + num_cells2+1
    nb_spline_eta1            = num_cells1 + sp_deg1+1
    nb_spline_eta2            = num_cells2 + sp_deg2+1
    
    SLL_ASSERT(size(coeffs_1d,1)==(num_cells1 + sp_deg1+1)*(num_cells2+sp_deg2+1))
    do i = 1, sp_deg1 + 1
      interpolator%t1(i) = eta1_min
    enddo
    eta1 = eta1_min
    do i = sp_deg1 + 2, num_cells1 + 1 + sp_deg1
      eta1 = eta1 + delta1
      interpolator%t1(i) = eta1
    enddo
    do i = num_cells1 + sp_deg1 + 2, num_cells1 + 1 + 2*sp_deg1
      interpolator%t1(i) = eta1
    enddo
    do i = 1, sp_deg2 + 1
      interpolator%t2(i) = eta2_min
    enddo
    eta2 = eta2_min
    do i = sp_deg2 + 2, num_cells2 + 1 + sp_deg2
      eta2 = eta2 + delta2
      interpolator%t2(i) = eta2
    enddo
    do i = num_cells2 + sp_deg2 + 2, num_cells2 + 1 + 2*sp_deg2
      interpolator%t2(i) = eta2
    enddo
    interpolator%coeff_splines(:,:) = 0.0_f64
    do j = 1,nb_spline_eta2
    do i = 1,nb_spline_eta1
      interpolator%coeff_splines(i+1,j+1) = coeffs_1d(i+nb_spline_eta1*(j-1))
    end do
    end do

  case(2145)  !left: Dirichlet, right: Hermite, bottom: Dirichlet, Top: Hermite  

    interpolator%size_coeffs1 = num_cells1 + sp_deg1 + 1
    interpolator%size_coeffs2 = num_cells2 + sp_deg2 + 1
    interpolator%size_t1      = 2*sp_deg1 + num_cells1 + 1
    interpolator%size_t2      = 2*sp_deg2 + num_cells2 + 1
    nb_spline_eta1            = num_cells1 + sp_deg1 + 1
    nb_spline_eta2            = num_cells2 + sp_deg2 + 1
    
    SLL_ASSERT(size(coeffs_1d,1)==(num_cells1+sp_deg1+1)*(num_cells2+sp_deg2+1))

    do i = 1, sp_deg1 + 1
      interpolator%t1(i) = eta1_min
    enddo
    eta1 = eta1_min
    do i = sp_deg1 + 2, num_cells1 + 1 + sp_deg1
      eta1 = eta1 + delta1
      interpolator%t1(i) = eta1
    enddo
    do i = num_cells1 + sp_deg1 + 2, num_cells1 + 1 + 2*sp_deg1
      interpolator%t1(i) = eta1
    enddo
    do i = 1, sp_deg2 + 1
      interpolator%t2(i) = eta2_min
    enddo
    eta2 = eta2_min
    do i = sp_deg2 + 2, num_cells2 + 1 + sp_deg2
      eta2 = eta2 + delta2
      interpolator%t2(i) = eta2
    enddo
    do i = num_cells2 + sp_deg2 + 2, num_cells2 + 1 + 2*sp_deg2
      interpolator%t2(i) = eta2
    enddo
    interpolator%coeff_splines(:,:) = 0.0_f64
    do j = 1,nb_spline_eta2
    do i = 1,nb_spline_eta1
      interpolator%coeff_splines(i+1,j+1) = coeffs_1d(i+nb_spline_eta1*(j-1))
    end do
    end do
     
 case(2124)  !left: Hermite, right: Dirichlet, bottom: Dirichlet, Top: Hermite

    interpolator%size_coeffs1 = num_cells1 + sp_deg1+1
    interpolator%size_coeffs2 = num_cells2 + sp_deg2+1
    interpolator%size_t1      = 2*sp_deg1  + num_cells1+1
    interpolator%size_t2      = 2*sp_deg2  + num_cells2+1
    nb_spline_eta1            = num_cells1 + sp_deg1+1
    nb_spline_eta2            = num_cells2 + sp_deg2+1
    
    SLL_ASSERT(size(coeffs_1d,1)==(num_cells1+sp_deg1+1)*(num_cells2+sp_deg2+1))
    do i = 1, sp_deg1 + 1
      interpolator%t1(i) = eta1_min
    enddo
    eta1 = eta1_min
    do i = sp_deg1 + 2, num_cells1 + 1 + sp_deg1
      eta1 = eta1 + delta1
      interpolator%t1(i) = eta1
    enddo
    do i = num_cells1 + sp_deg1 + 2, num_cells1 + 1 + 2*sp_deg1
      interpolator%t1(i) = eta1
    enddo
    do i = 1, sp_deg2 + 1
      interpolator%t2(i) = eta2_min
    enddo
    eta2 = eta2_min
    do i = sp_deg2 + 2, num_cells2 + 1 + sp_deg2
      eta2 = eta2 + delta2
      interpolator%t2(i) = eta2
    enddo
    do i = num_cells2 + sp_deg2 + 2, num_cells2 + 1 + 2*sp_deg2
      interpolator%t2(i) = eta2
    enddo
    interpolator%coeff_splines(:,:) = 0.0_f64
    do j = 1,nb_spline_eta2
    do i = 1,nb_spline_eta1
      interpolator%coeff_splines(i+1,j+1) = coeffs_1d(i+nb_spline_eta1 *(j-1))
    end do
    end do
  
  case(801)  !left: Dirichlet, right: Hermite, bottom: Hermite, Top: Dirichlet

    interpolator%size_coeffs1 = num_cells1 + sp_deg1+1
    interpolator%size_coeffs2 = num_cells2 + sp_deg2+1
    interpolator%size_t1      = 2*sp_deg1 + num_cells1 + 1
    interpolator%size_t2      = 2*sp_deg2 + num_cells2 + 1
    nb_spline_eta1            = num_cells1 + sp_deg1 +1
    nb_spline_eta2            = num_cells2 + sp_deg2 +1
    
    SLL_ASSERT(size(coeffs_1d,1).ne.(num_cells1+sp_deg1+1)*(num_cells2+sp_deg2+1))
    do i = 1, sp_deg1 + 1
      interpolator%t1(i) = eta1_min
    enddo
    eta1 = eta1_min
    do i = sp_deg1 + 2, num_cells1 + 1 + sp_deg1
      eta1 = eta1 + delta1
      interpolator%t1(i) = eta1
    enddo
    do i = num_cells1 + sp_deg1 + 2, num_cells1 + 1 + 2*sp_deg1
      interpolator%t1(i) = eta1
    enddo
    do i = 1, sp_deg2 + 1
      interpolator%t2(i) = eta2_min
    enddo
    eta2 = eta2_min
    do i = sp_deg2 + 2, num_cells2 + 1 + sp_deg2
      eta2 = eta2 + delta2
      interpolator%t2(i) = eta2
    enddo
    do i = num_cells2 + sp_deg2 + 2, num_cells2 + 1 + 2*sp_deg2
      interpolator%t2(i) = eta2
    enddo
    interpolator%coeff_splines(:,:) = 0.0_f64
    do j = 1,nb_spline_eta2
    do i = 1,nb_spline_eta1
      interpolator%coeff_splines(i+1,j+1) = coeffs_1d(i+nb_spline_eta1*(j-1))
    end do
    end do
    
  case(2340) ! Hermite in al sides
       
    interpolator%size_coeffs1 = num_cells1 + sp_deg1
    interpolator%size_coeffs2 = num_cells2 + sp_deg2
    interpolator%size_t1      = 2*sp_deg1 + num_cells1 + 1
    interpolator%size_t2      = 2*sp_deg2 + num_cells2 + 1
    nb_spline_eta1            = num_cells1 + sp_deg1
    nb_spline_eta2            = num_cells2 + sp_deg2
    
    SLL_ASSERT(size(coeffs_1d,1)==(num_cells1+sp_deg1)*(num_cells2+sp_deg2))
    do i = 1, sp_deg1 + 1
      interpolator%t1(i) = eta1_min
    enddo
    eta1 = eta1_min
    do i = sp_deg1 + 2, num_cells1 + 1 + sp_deg1
      eta1 = eta1 + delta1
      interpolator%t1(i) = eta1
    enddo
    do i = num_cells1 + sp_deg1 + 2, num_cells1 + 1 + 2*sp_deg1
      interpolator%t1(i) = eta1
    enddo
    do i = 1, sp_deg2 + 1
       interpolator%t2(i) = eta2_min
    enddo
    eta2 = eta2_min
    do i = sp_deg2 + 2, num_cells2 + 1 + sp_deg2
      eta2 = eta2 + delta2
      interpolator%t2(i) = eta2
    enddo
    do i = num_cells2 + sp_deg2 + 2, num_cells2 + 1 + 2*sp_deg2
      interpolator%t2(i) = eta2
    enddo
    interpolator%coeff_splines(:,:) = 0.0_f64
    do j = 1,nb_spline_eta2
    do i = 1,nb_spline_eta1
      interpolator%coeff_splines(i,j) = coeffs_1d(i+nb_spline_eta1*(j-1))
    end do
    end do
   
  case default

     stop 'arbitrary_degree_spline_2d set_spline_coefficients not recognized.'

  end select

else if (present(coeffs_2d) ) then 

  if ( present(coeff2d_size1) .and. present(coeff2d_size2)) then

    interpolator%size_coeffs1 = coeff2d_size1
    interpolator%size_coeffs2 = coeff2d_size2
    interpolator%size_t1      = sp_deg1 + coeff2d_size1 +1 
    interpolator%size_t2      = sp_deg2 + coeff2d_size2 +1
    
    if (coeff2d_size1 > num_cells1+1+4*sp_deg1) then
       stop 'size1 of coeff2d is too big'
    end if
    
    if (coeff2d_size2 > num_cells2+1+4*sp_deg2) then
       stop 'size2 of coeff2d is too big'
    end if
    
    interpolator%coeff_splines(1:coeff2d_size1,1:coeff2d_size2) = &
         coeffs_2d(1:coeff2d_size1,1:coeff2d_size2)
    
    if (present(knots1).and.present(knots2)) then 
       
      SLL_ASSERT(size_knots1==(coeff2d_size1+sp_deg1+1)) 
      SLL_ASSERT(size_knots2==(coeff2d_size2+sp_deg2+1))
       
      if ( size_knots1 > (num_cells1+1)*(sp_deg1+1)) then
        stop 'size1 of knots1 is too big'
      end if
      if ( size_knots2 >  (num_cells2+1)*(sp_deg2+1)) then
        stop 'size2 of knots2 is too big'
      end if
      
      interpolator%t1(1:interpolator%size_t1) = knots1(1:interpolator%size_t1)
      interpolator%t2(1:interpolator%size_t2) = knots2(1:interpolator%size_t2)
      
    else if ( (.not. present(knots1)).and.(.not. present(knots2))) then
       
      if ( interpolator%size_t1 > (num_cells1 + 1)*(sp_deg1+1)) then
        stop 'size1 of knots1 is too big'
      end if
      
      if ( interpolator%size_t2 >  (num_cells2 + 1)*(sp_deg2+1)) then
         stop 'size2 of knots2 is too big'
      end if

      interpolator%t1(1:sp_deg1+1)  = eta1_min
      interpolator%t1(coeff2d_size1+2:coeff2d_size1+2+sp_deg1) = eta1_max
        
      do i = 1, coeff2d_size1 -sp_deg1
        interpolator%t1(i+sp_deg1+1) = eta1_min+i*(eta1_max-eta1_min)/(coeff2d_size1-sp_deg1+1)   
      end do
        
      interpolator%t2(1:sp_deg2+1) = eta2_min
      interpolator%t2(coeff2d_size2+2:coeff2d_size2+2+sp_deg2) = eta2_max
        
      do i = 1, coeff2d_size2 -sp_deg2
        interpolator%t2(i+sp_deg2+1)=eta2_min+i*(eta2_max-eta2_min)/(coeff2d_size2-sp_deg2+1)   
      end do

    else 

      print*, 'Problem in set_coefficients in arbitrary_degree_spline_2d'
      stop 'Knots1 or Knots2 is not present'

    end if

  else 

    print*, 'Problem in set_coefficients in arbitrary_degree_spline_2d'
    print*, 'problem with the size of coeffs_2d'
    stop 'the number of coefficients must be specified'

  end if
    
else 
  stop 'Problem in set_coefficients: must be have coefficients'
end if

interpolator%coefficients_set = .true.
 
end subroutine !set_coefficients_ad2d

!> @brief computing the coefficients spline with a given 
!>  data_array 2D cooresponding at the values of a function 
!> @details computing the coefficients spline with a given 
!>  data_array 2D coorespondind at the values of a function 
!>  on eta1_coords of size size_eta1_coords in the first direction and 
!>  on eta2_coords of size size_eta2_coords in the second direction
!>  if the eta1_coords and eta2_coords is not given 
!>  we consider that the values of the function is on the points in the mesh_2d
!> 
!> The parameters are
!> @param interpolator the type sll_arbitrary_degree_spline_interpolator_2d
!> @param[in] data_array the 2d arrays corresponding at the values of a function
!> @param[in] eta1_coords the 1d arrays corresponding at the points eta1 
!> @param[in] size_eta1_coords the size of eta1_coords
!> @param[in] eta2_coords the 1d arrays corresponding at the points eta2
!> @param[in] size_eta2_coords the size of eta2_coords
!> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_2d

subroutine compute_interpolants_ad2d( interpolator,     &
                                      data_array,       &
                                      eta1_coords,      &
                                      size_eta1_coords, &
                                      eta2_coords,      &
                                      size_eta2_coords )

class(sll_arbitrary_degree_spline_interpolator_2d), intent(inout) :: interpolator

sll_real64, dimension(:,:), intent(in)           :: data_array
sll_real64, dimension(:),   intent(in), optional :: eta1_coords
sll_real64, dimension(:),   intent(in), optional :: eta2_coords
sll_int32,                  intent(in), optional :: size_eta1_coords
sll_int32,                  intent(in), optional :: size_eta2_coords

sll_real64, dimension(:),   pointer :: taux
sll_real64, dimension(:),   pointer :: tauy
sll_real64, dimension(:,:), pointer :: gtau
sll_real64, dimension(:,:), pointer :: gtau_der1
sll_real64, dimension(:,:), pointer :: gtau_der2

sll_int32, pointer :: taux_deriv(:)
sll_int32, pointer :: tauy_deriv(:)

sll_int32  :: mx,my
sll_real64 :: eta1_min, eta1_max, delta_eta1
sll_real64 :: eta2_min, eta2_max, delta_eta2
sll_int32  :: nx
sll_int32  :: ny
sll_real64 :: period1
sll_real64 :: period2
sll_int32  :: kx
sll_int32  :: ky
sll_int32  :: ierr
sll_int32  :: i

eta1_min   = interpolator%eta1_min
eta2_min   = interpolator%eta2_min
eta1_max   = interpolator%eta1_max
eta2_max   = interpolator%eta2_max
delta_eta1 = (eta1_max-eta1_min)/(interpolator%num_pts1-1)
delta_eta2 = (eta2_max-eta2_min)/(interpolator%num_pts2-1)

if(present(eta1_coords)) then 

  SLL_ASSERT(present(eta2_coords))
  SLL_ASSERT(present(size_eta2_coords))
  SLL_ASSERT(present(size_eta1_coords))

  nx = size_eta1_coords
  ny = size_eta2_coords
  
  SLL_ALLOCATE(taux(1:nx),ierr)
  SLL_ALLOCATE(tauy(1:ny),ierr)
  taux(1:nx) = eta1_coords(1:nx)
  tauy(1:ny) = eta2_coords(1:ny)

else ! size depends on BC combination, filled out at initialization.

  nx = interpolator%num_pts1
  ny = interpolator%num_pts2

  SLL_ALLOCATE(taux(1:nx),ierr)
  SLL_ALLOCATE(tauy(1:ny),ierr)
  
  do i = 1,nx
     taux(i) = interpolator%eta1_min + delta_eta1*(i-1)
  end do
  do i = 1,ny
     tauy(i) = interpolator%eta2_min + delta_eta2*(i-1)
  end do
  
end if

SLL_ASSERT(nx .le. interpolator%num_pts1+8*interpolator%spline_degree1)
SLL_ASSERT(ny .le. interpolator%num_pts2+8*interpolator%spline_degree1)
SLL_ASSERT(size(data_array,1) .ge. nx)
SLL_ASSERT(size(data_array,2) .ge. ny)
SLL_ASSERT(size(taux) .ge. nx)
SLL_ASSERT(size(tauy) .ge. ny)

kx  = interpolator%spline_degree1 + 1
ky  = interpolator%spline_degree2 + 1
period1 = interpolator%eta1_max - interpolator%eta1_min
period2 = interpolator%eta2_max - interpolator%eta2_min

! compute the knots t1 and t2
SLL_ALLOCATE(taux_deriv(2),ierr)
SLL_ALLOCATE(tauy_deriv(2),ierr)
  
select case (interpolator%bc_selector)

case (0) ! periodic-periodic

  interpolator%size_coeffs1 = nx
  interpolator%size_coeffs2 = ny
  interpolator%size_t1      = kx + nx
  interpolator%size_t2      = ky + ny

  !Apply periodic boundary conditions
  taux(nx) = taux(1)+period1
  tauy(ny) = tauy(1)+period2

  SLL_ALLOCATE(gtau(1:nx,1:ny), ierr)
  gtau(1:nx-1,1:ny-1) = data_array(1:nx-1,1:ny-1)
  gtau(nx,1:ny-1)     = data_array(1,1:ny-1 )
  gtau(1:nx-1,ny)     = data_array(1:nx-1,1)
  gtau(nx,ny)         = data_array(1,1)

  call spli2d_custom(nx, kx, taux, ny, ky, tauy, &
                     gtau,                       &
                     interpolator%coeff_splines, &
                     interpolator%t1,            &
                     interpolator%t2)
   
case (9) ! 2. dirichlet-left, dirichlet-right, periodic

  interpolator%size_coeffs1 = nx
  interpolator%size_coeffs2 = ny
  interpolator%size_t1      = kx+nx
  interpolator%size_t2      = ky+ny
       
  tauy(1:ny-1)      = tauy(1:ny-1)
  tauy(ny)          = tauy(1)+period2

  SLL_ALLOCATE(gtau(1:nx,1:ny), ierr)

  gtau(1:nx-1,1:ny) = data_array(1:nx-1,1:ny)
  gtau(nx,1:ny)     = data_array(1,1:ny)

  call spli2d_custom(nx, kx, taux, ny, ky, tauy, &
                     gtau,                       &
                     interpolator%coeff_splines, &
                     interpolator%t1,            &
                     interpolator%t2)

  interpolator%coeff_splines(1,1:ny)  = data_array(1,1:ny)
  interpolator%coeff_splines(nx,1:ny) = data_array(nx,1:ny)
  
case(576) !  3. periodic, dirichlet-bottom, dirichlet-top

  interpolator%size_coeffs1 = nx!+1
  interpolator%size_coeffs2 = ny
  interpolator%size_t1 = kx + nx !+ 1
  interpolator%size_t2 = ky + ny 

  taux(1:nx-1) = taux(1:nx-1)
  taux(nx)     = taux(1)+period1

  SLL_ALLOCATE(gtau(1:nx,1:ny), ierr)

  gtau(1:nx,1:ny-1) = data_array(1:nx,1:ny-1)
  gtau(1:nx,ny)     = data_array(1:nx,1)

  call spli2d_custom(nx, kx, taux, ny, ky, tauy, &
                     gtau,                       &
                     interpolator%coeff_splines, &
                     interpolator%t1,            &
                     interpolator%t2)

  interpolator%coeff_splines(1:nx,1)   = data_array(1:nx,1)
  interpolator%coeff_splines(1:nx,ny) = data_array(1:nx,ny)
       
case (585) ! 4. dirichlet in all sides

  interpolator%size_coeffs1 = nx
  interpolator%size_coeffs2 = ny
  interpolator%size_t1 = kx + nx 
  interpolator%size_t2 = ky + ny 
  
  SLL_ALLOCATE(gtau(1:nx,1:ny),ierr)
  gtau = data_array(1:nx,1:ny)
  call spli2d_custom( nx, kx, taux, ny, ky, tauy, &
                      gtau,                       &
                      interpolator%coeff_splines, &
                      interpolator%t1,            &
                      interpolator%t2)

  interpolator%coeff_splines(1,1:ny)  = data_array(1,1:ny)
  interpolator%coeff_splines(nx,1:ny) = data_array(nx,1:ny)
  interpolator%coeff_splines(1:nx,1)  = data_array(1:nx,1)
  interpolator%coeff_splines(1:nx,ny) = data_array(1:nx,ny)

case (650) !left: Neumann, right: Dirichlet, bottom: Neumann, Top: Dirichlet

  mx = 2
  my = 2
  interpolator%size_coeffs1 = nx + mx
  interpolator%size_coeffs2 = ny + my
  interpolator%size_t1 = kx + nx + mx
  interpolator%size_t2 = ky + ny + my
  
  SLL_ALLOCATE(gtau(1:nx,1:ny),ierr)
  SLL_CLEAR_ALLOCATE(gtau_der1(1:2,1:ny),ierr)
  SLL_CLEAR_ALLOCATE(gtau_der2(1:my,1:nx+mx),ierr)

  gtau                 = data_array(1:nx,1:ny)
  taux_deriv(1)        = 1
  taux_deriv(2)        = nx
  gtau_der1(1,1:ny)    = 0.0_f64
  gtau_der1(2,1:ny)    = interpolator%slope_max1(1:ny)
  tauy_deriv(1)        = 1
  tauy_deriv(2)        = ny
  gtau_der2(1,1:nx+mx) = 0.0_f64
  gtau_der2(2,1:nx+mx) = interpolator%slope_max2(1:nx+mx)

  call spli2d_custom_derder(nx,                         &
                            mx,                         &
                            kx,                         &
                            taux,                       &
                            taux_deriv,                 &
                            ny,                         &
                            my,                         &
                            ky,                         &
                            tauy,                       &
                            tauy_deriv,                 &
                            gtau,                       &
                            gtau_der1,                  &
                            gtau_der2,                  &
                            interpolator%coeff_splines, &
                            interpolator%t1,            &
                            interpolator%t2)

  SLL_DEALLOCATE( gtau_der1,ierr)
  SLL_DEALLOCATE( gtau_der2,ierr)

case(657) !left: Dirichlet, right: Neumann, bottom: Neumann, Top: Dirichlet 

       mx = 2
       my = 2
       interpolator%size_coeffs1 = nx + mx
       interpolator%size_coeffs2 = ny + my
       interpolator%size_t1 = kx + nx + mx
       interpolator%size_t2 = ky + ny + my
       
       !  data_array must have the same dimension than 
       !  size(  taux ) x  size(  tauy )
       !  i.e  data_array must have the dimension nx x ny
       SLL_ALLOCATE( gtau(1:nx,1:ny),ierr)
       SLL_ALLOCATE( gtau_der1(mx,1:ny),ierr)
       SLL_ALLOCATE( gtau_der2(my,1:nx+mx),ierr)
       gtau = data_array(1:nx,1:ny)
       taux_deriv(1) = 1
       taux_deriv(2) = nx
       gtau_der1(1,1:ny)     = interpolator%slope_min1(1:ny) 
       gtau_der1(2,1:ny)     = 0.0_f64
       tauy_deriv(1) = 1
       tauy_deriv(2) = ny
       gtau_der2(1,1:nx+mx)= 0.0_f64
       gtau_der2(2,1:nx+mx)=interpolator%slope_max2(1:nx+mx)
       call spli2d_custom_derder(&
            nx,&
            mx,&
            kx, &
            taux, &
            taux_deriv,&
            ny, &
            my,&
            ky, tauy, &
            tauy_deriv,&
            gtau,&
            gtau_der1,&
            gtau_der2,&
            interpolator%coeff_splines,&!(1:nx,1:ny),&
            interpolator%t1,&!(1:nx+kx), &
            interpolator%t2)!(1:ny+ky) )

       SLL_DEALLOCATE( gtau_der1,ierr)
       SLL_DEALLOCATE( gtau_der2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:ny+my)   = interpolator%value_min1(1:ny+my)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:nx,1)   = interpolator%value_min2(1:nx)
  !     interpolator%coeff_splines(1:nx,ny) = interpolator%value_max2(1:nx)


    case(780)  !left: Hermite, right: Dirichlet, bottom: Hermite, Top: Dirichlet

       mx = 2
       my = 2
       interpolator%size_coeffs1 = nx + mx
       interpolator%size_coeffs2 = ny + my
       interpolator%size_t1 = kx + nx + mx
       interpolator%size_t2 = ky + ny + my
       
       !  data_array must have the same dimension than 
       !  size(  taux ) x  size(  tauy )
       !  i.e  data_array must have the dimension nx x ny
       SLL_ALLOCATE(gtau(1:nx,1:ny),ierr)
       SLL_CLEAR_ALLOCATE(gtau_der1(1:mx,1:ny),ierr)
       SLL_CLEAR_ALLOCATE(gtau_der2(1:my,1:nx+mx),ierr)
       gtau = data_array(1:nx,1:ny)
       taux_deriv(1) = 1
       taux_deriv(2) = nx
       gtau_der1(1,1:ny) = interpolator%slope_min1(1:ny)
       gtau_der1(2,1:ny) = interpolator%slope_max1(1:ny)
       tauy_deriv(1) = 1
       tauy_deriv(2) = ny
       gtau_der2(1,1:nx+mx)= &
          interpolator%slope_min2(1:nx+mx)
       gtau_der2(2,1:nx+mx)= &
          interpolator%slope_max2(1:nx+mx)

       call spli2d_custom_derder(&
            nx,&
            mx,&
            kx, &
            taux, &
            taux_deriv,&
            ny, &
            my,&
            ky, tauy, &
            tauy_deriv,&
            gtau,&
            gtau_der1,&
            gtau_der2,&
            interpolator%coeff_splines,&!(1:nx,1:ny),&
            interpolator%t1,&!(1:nx+kx), &
            interpolator%t2)!(1:ny+ky) )


       SLL_DEALLOCATE( gtau_der1,ierr)
       SLL_DEALLOCATE( gtau_der2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:ny+my)   = interpolator%value_min1(1:ny+my)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:nx,1)   = interpolator%value_min2(1:nx)
  !     interpolator%coeff_splines(1:nx,ny) = interpolator%value_max2(1:nx)


    case(801)  !left: Dirichlet, right: Hermite, bottom: Hermite, Top: Dirichlet


       mx = 2
       my = 2
       interpolator%size_coeffs1 = nx + mx
       interpolator%size_coeffs2 = ny + my
       interpolator%size_t1 = kx + nx + mx
       interpolator%size_t2 = ky + ny + my
       
       !  data_array must have the same dimension than 
       !  size(  taux ) x  size(  tauy )
       !  i.e  data_array must have the dimension nx x ny
       SLL_ALLOCATE( gtau(1:nx,1:ny),ierr)
       SLL_ALLOCATE( gtau_der1(mx,1:ny),ierr)
       SLL_ALLOCATE( gtau_der2(my,1:nx+mx),ierr)
       gtau = data_array(1:nx,1:ny)
       taux_deriv(1) = 1
       taux_deriv(2) = nx
       gtau_der1(1,1:ny)     = interpolator%slope_min1(1:ny) 
       gtau_der1(2,1:ny)     = interpolator%slope_max1(1:ny) 
       tauy_deriv(1) = 1
       tauy_deriv(2) = ny
       gtau_der2(1,1:nx+mx)=interpolator%slope_min2(1:nx+mx)
       gtau_der2(2,1:nx+mx)=interpolator%slope_max2(1:nx+mx)
       call spli2d_custom_derder(&
            nx,&
            mx,&
            kx, &
            taux, &
            taux_deriv,&
            ny, &
            my,&
            ky, tauy, &
            tauy_deriv,&
            gtau,&
            gtau_der1,&
            gtau_der2,&
            interpolator%coeff_splines,&!(1:nx,1:ny),&
            interpolator%t1,&!(1:nx+kx), &
            interpolator%t2)!(1:ny+ky) )

       SLL_DEALLOCATE( gtau_der1,ierr)
       SLL_DEALLOCATE( gtau_der2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:ny+my)   = interpolator%value_min1(1:ny+my)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:nx,1)   = interpolator%value_min2(1:nx)
  !     interpolator%coeff_splines(1:nx,ny) = interpolator%value_max2(1:nx)

    case(804)  !left: Hermite, right: Hermite, bottom: Hermite, Top: Dirichlet
       mx = 2
       my = 2
       interpolator%size_coeffs1 = nx + mx
       interpolator%size_coeffs2 = ny + my
       interpolator%size_t1 = kx + nx + mx
       interpolator%size_t2 = ky + ny + my
       
       !  data_array must have the same dimension than 
       !  size(  taux ) x  size(  tauy )
       !  i.e  data_array must have the dimension nx x ny
       SLL_ALLOCATE( gtau(1:nx,1:ny),ierr)
       SLL_ALLOCATE( gtau_der1(mx,1:ny),ierr)
       SLL_ALLOCATE( gtau_der2(my,1:nx+mx),ierr)
       gtau = data_array(1:nx,1:ny)
       taux_deriv(1) = 1
       taux_deriv(2) = nx
       gtau_der1(1,1:ny)     = interpolator%slope_min1(1:ny) 
       gtau_der1(2,1:ny)     = interpolator%slope_max1(1:ny) 
       tauy_deriv(1) = 1
       tauy_deriv(2) = ny
       gtau_der2(1,1:nx+mx)=interpolator%slope_min2(1:nx+mx)
       gtau_der2(2,1:nx+mx)=interpolator%slope_max2(1:nx+mx)
       call spli2d_custom_derder(&
            nx,&
            mx,&
            kx, &
            taux, &
            taux_deriv,&
            ny, &
            my,&
            ky, tauy, &
            tauy_deriv,&
            gtau,&
            gtau_der1,&
            gtau_der2,&
            interpolator%coeff_splines,&!(1:nx,1:ny),&
            interpolator%t1,&!(1:nx+kx), &
            interpolator%t2)!(1:ny+ky) )

       SLL_DEALLOCATE( gtau_der1,ierr)
       SLL_DEALLOCATE( gtau_der2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:ny+my)   = interpolator%value_min1(1:ny+my)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:nx,1)   = interpolator%value_min2(1:nx)
  !     interpolator%coeff_splines(1:nx,ny) = interpolator%value_max2(1:nx)

    case(1098)  !left: Neumann, right: Dirichlet, bottom: Dirichlet, Top: Neumann
       mx = 2
       my = 2
       interpolator%size_coeffs1 = nx + mx
       interpolator%size_coeffs2 = ny + my
       interpolator%size_t1 = kx + nx + mx
       interpolator%size_t2 = ky + ny + my
       
       !  data_array must have the same dimension than 
       !  size(  taux ) x  size(  tauy )
       !  i.e  data_array must have the dimension nx x ny
       SLL_ALLOCATE( gtau(1:nx,1:ny),ierr)
       SLL_ALLOCATE( gtau_der1(mx,1:ny),ierr)
       SLL_ALLOCATE( gtau_der2(my,1:nx+mx),ierr)
       gtau = data_array(1:nx,1:ny)
       taux_deriv(1) = 1
       taux_deriv(2) = nx
       gtau_der1(1,1:ny)     = 0.0_f64
       gtau_der1(2,1:ny)     = interpolator%slope_max1(1:ny)
       tauy_deriv(1) = 1
       tauy_deriv(2) = ny
       gtau_der2(1,1:nx+mx)=interpolator%slope_min2(1:nx+mx)
       gtau_der2(2,1:nx+mx)= 0.0_f64
       call spli2d_custom_derder(&
            nx,&
            mx,&
            kx, &
            taux, &
            taux_deriv,&
            ny, &
            my,&
            ky, tauy, &
            tauy_deriv,&
            gtau,&
            gtau_der1,&
            gtau_der2,&
            interpolator%coeff_splines,&!(1:nx,1:ny),&
            interpolator%t1,&!(1:nx+kx), &
            interpolator%t2)!(1:ny+ky) )

       SLL_DEALLOCATE( gtau_der1,ierr)
       SLL_DEALLOCATE( gtau_der2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:ny+my)   = interpolator%value_min1(1:ny+my)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:nx,1)   = interpolator%value_min2(1:nx)
  !     interpolator%coeff_splines(1:nx,ny) = interpolator%value_max2(1:nx)

    case(1105)  !left: Dirichlet, right: Neumann, bottom: Dirichlet, Top: Neumann
       mx = 2
       my = 2
       interpolator%size_coeffs1 = nx + mx
       interpolator%size_coeffs2 = ny + my
       interpolator%size_t1 = kx + nx + mx
       interpolator%size_t2 = ky + ny + my
       
       !  data_array must have the same dimension than 
       !  size(  taux ) x  size(  tauy )
       !  i.e  data_array must have the dimension nx x ny
       SLL_ALLOCATE( gtau(1:nx,1:ny),ierr)
       SLL_ALLOCATE( gtau_der1(2,ny),ierr)
       SLL_ALLOCATE( gtau_der2(my,nx+mx),ierr)
       gtau = data_array(1:nx,1:ny)
       taux_deriv(1) = 1
       taux_deriv(2) = nx
       gtau_der1(1,1:ny)     = interpolator%slope_min1(1:ny)
       gtau_der1(2,1:ny)     = 0.0_f64
       tauy_deriv(1) = 1
       tauy_deriv(2) = ny
       gtau_der2(1,1:nx+mx)=interpolator%slope_min2(1:nx+mx)
       gtau_der2(2,1:nx+mx)= 0.0_f64
       call spli2d_custom_derder(&
            nx,&
            mx,&
            kx, &
            taux, &
            taux_deriv,&
            ny, &
            my,&
            ky, tauy, &
            tauy_deriv,&
            gtau,&
            gtau_der1,&
            gtau_der2,&
            interpolator%coeff_splines,&!(1:nx,1:ny),&
            interpolator%t1,&!(1:nx+kx), &
            interpolator%t2)!(1:ny+ky) )

       SLL_DEALLOCATE( gtau_der1,ierr)
       SLL_DEALLOCATE( gtau_der2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:ny+my)   = interpolator%value_min1(1:ny+my)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1:nx+mx,1)   = interpolator%value_min2(1:nx+mx)
  !     interpolator%coeff_splines(1:nx,ny) = interpolator%value_max2(1:nx)

    case(1170)  !left: Neumann, right: Neumann, bottom: Neuman, Top: Neumann

       mx = 2
       my = 2
       interpolator%size_coeffs1 = nx + mx
       interpolator%size_coeffs2 = ny + my
       interpolator%size_t1 = kx + nx + mx
       interpolator%size_t2 = ky + ny + my
       
       !  data_array must have the same dimension than 
       !  size(  taux ) x  size(  tauy )
       !  i.e  data_array must have the dimension nx x ny
       SLL_ALLOCATE( gtau(1:nx,1:ny),ierr)
       SLL_ALLOCATE( gtau_der1(2,ny),ierr)
       SLL_ALLOCATE( gtau_der2(my,nx+mx),ierr)
       gtau = data_array(1:nx,1:ny)
       taux_deriv(1) = 1
       taux_deriv(2) = nx
       gtau_der1(1,1:ny)     = 0.0_f64
       gtau_der1(2,1:ny)     = 0.0_f64
       tauy_deriv(1) = 1
       tauy_deriv(2) = ny
       gtau_der2(1,1:nx+mx)= 0.0_f64
       gtau_der2(2,1:nx+mx)= 0.0_f64
       call spli2d_custom_derder(&
            nx,&
            mx,&
            kx, &
            taux, &
            taux_deriv,&
            ny, &
            my,&
            ky, tauy, &
            tauy_deriv,&
            gtau,&
            gtau_der1,&
            gtau_der2,&
            interpolator%coeff_splines,&!(1:nx,1:ny),&
            interpolator%t1,&!(1:nx+kx), &
            interpolator%t2)!(1:ny+ky) )

       SLL_DEALLOCATE( gtau_der1,ierr)
       SLL_DEALLOCATE( gtau_der2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:ny+my)   = interpolator%value_min1(1:ny+my)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1:nx+mx,1)   = interpolator%value_min2(1:nx+mx)
  !     interpolator%coeff_splines(1:nx,ny) = interpolator%value_max2(1:nx)
    case(2338)  !left: Dirichlet, right: Hermite, bottom: Hermite, Top: Hermite
       mx = 2
       my = 2
       interpolator%size_coeffs1 = nx + mx
       interpolator%size_coeffs2 = ny + my
       interpolator%size_t1 = kx + nx + mx
       interpolator%size_t2 = ky + ny + my
       
       !  data_array must have the same dimension than 
       !  size(  taux ) x  size(  tauy )
       !  i.e  data_array must have the dimension nx x ny
       SLL_ALLOCATE( gtau(1:nx,1:ny),ierr)
       SLL_ALLOCATE( gtau_der1(2,ny),ierr)
       SLL_ALLOCATE( gtau_der2(my,nx+mx),ierr)
       gtau = data_array(1:nx,1:ny)
       taux_deriv(1) = 1
       taux_deriv(2) = nx
       gtau_der1(1,1:ny)     = interpolator%slope_min1(1:ny)
       gtau_der1(2,1:ny)     = interpolator%slope_max1(1:ny)
       tauy_deriv(1) = 1
       tauy_deriv(2) = ny
       gtau_der2(1,1:nx+mx)=interpolator%slope_min2(1:nx+mx)
       gtau_der2(2,1:nx+mx)=interpolator%slope_max2(1:nx+mx)
       call spli2d_custom_derder(&
            nx,&
            mx,&
            kx, &
            taux, &
            taux_deriv,&
            ny, &
            my,&
            ky, tauy, &
            tauy_deriv,&
            gtau,&
            gtau_der1,&
            gtau_der2,&
            interpolator%coeff_splines,&!(1:nx,1:ny),&
            interpolator%t1,&!(1:nx+kx), &
            interpolator%t2)!(1:ny+ky) )

       SLL_DEALLOCATE( gtau_der1,ierr)
       SLL_DEALLOCATE( gtau_der2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:ny+my)   = interpolator%value_min1(1:ny+my)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1:nx+mx,1)   = interpolator%value_min2(1:nx+mx)
  !     interpolator%coeff_splines(1:nx,ny) = interpolator%value_max2(1:nx)
       
    case(2145) !left: Dirichlet, right: Hermite, bottom: Dirichlet, Top: Hermite  
       mx = 2
       my = 2
       interpolator%size_coeffs1 = nx + mx
       interpolator%size_coeffs2 = ny + my
       interpolator%size_t1 = kx + nx + mx
       interpolator%size_t2 = ky + ny + my
       
       !  data_array must have the same dimension than 
       !  size(  taux ) x  size(  tauy )
       !  i.e  data_array must have the dimension nx x ny
       SLL_ALLOCATE( gtau(1:nx,1:ny),ierr)
       SLL_ALLOCATE( gtau_der1(2,ny),ierr)
       SLL_ALLOCATE( gtau_der2(my,nx+mx),ierr)
       gtau = data_array(1:nx,1:ny)
       taux_deriv(1) = 1
       taux_deriv(2) = nx
       gtau_der1(1,1:ny)     = interpolator%slope_min1(1:ny)
       gtau_der1(2,1:ny)     = interpolator%slope_max1(1:ny)
       tauy_deriv(1) = 1
       tauy_deriv(2) = ny
       gtau_der2(1,1:nx+mx)=interpolator%slope_min2(1:nx+mx)
       gtau_der2(2,1:nx+mx)=interpolator%slope_max2(1:nx+mx)
       call spli2d_custom_derder(&
            nx,&
            mx,&
            kx, &
            taux, &
            taux_deriv,&
            ny, &
            my,&
            ky, tauy, &
            tauy_deriv,&
            gtau,&
            gtau_der1,&
            gtau_der2,&
            interpolator%coeff_splines,&!(1:nx,1:ny),&
            interpolator%t1,&!(1:nx+kx), &
            interpolator%t2)!(1:ny+ky) )

       SLL_DEALLOCATE( gtau_der1,ierr)
       SLL_DEALLOCATE( gtau_der2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:ny+my)   = interpolator%value_min1(1:ny+my)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1:nx+mx,1)   = interpolator%value_min2(1:nx+mx)
  !     interpolator%coeff_splines(1:nx,ny) = interpolator%value_max2(1:nx)


    case(2124)  !left: Hermite, right: Dirichlet, bottom: Dirichlet, Top: Hermite

       mx = 2
       my = 2
       interpolator%size_coeffs1 = nx + mx
       interpolator%size_coeffs2 = ny + my
       interpolator%size_t1 = kx + nx + mx
       interpolator%size_t2 = ky + ny + my
       
       !  data_array must have the same dimension than 
       !  size(  taux ) x  size(  tauy )
       !  i.e  data_array must have the dimension nx x ny
       SLL_CLEAR_ALLOCATE( gtau(1:nx,1:ny),ierr)
       SLL_CLEAR_ALLOCATE( gtau_der1(mx,1:ny),ierr)
       SLL_CLEAR_ALLOCATE( gtau_der2(my,1:nx+mx),ierr)
       gtau = data_array(1:nx,1:ny)
       taux_deriv(1) = 1
       taux_deriv(2) = nx
       gtau_der1(1,1:ny)     = interpolator%slope_min1(1:ny)
       gtau_der1(2,1:ny)     = interpolator%slope_max1(1:ny)
       tauy_deriv(1) = 1
       tauy_deriv(2) = ny
       gtau_der2(1,1:nx+mx)=interpolator%slope_min2(1:nx+mx)
       gtau_der2(2,1:nx+mx)=interpolator%slope_max2(1:nx+mx)
       call spli2d_custom_derder(&
            nx,&
            mx,&
            kx, &
            taux, &
            taux_deriv,&
            ny, &
            my,&
            ky, tauy, &
            tauy_deriv,&
            gtau,&
            gtau_der1,&
            gtau_der2,&
            interpolator%coeff_splines,&!(1:nx,1:ny),&
            interpolator%t1,&!(1:nx+kx), &
            interpolator%t2)!(1:ny+ky) )

       SLL_DEALLOCATE( gtau_der1,ierr)
       SLL_DEALLOCATE( gtau_der2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:ny+my)   = interpolator%value_min1(1:ny+my)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:nx,1)   = interpolator%value_min2(1:nx)
  !     interpolator%coeff_splines(1:nx,ny) = interpolator%value_max2(1:nx)

    case(2148)  !left:Hermite , right: Hermite, bottom: Dirichlet, Top: Hermite
       mx = 2
       my = 2
       interpolator%size_coeffs1 = nx + mx
       interpolator%size_coeffs2 = ny + my
       interpolator%size_t1 = kx + nx + mx
       interpolator%size_t2 = ky + ny + my
       
       !  data_array must have the same dimension than 
       !  size(  taux ) x  size(  tauy )
       !  i.e  data_array must have the dimension nx x ny
       SLL_ALLOCATE( gtau(1:nx,1:ny),ierr)
       SLL_ALLOCATE( gtau_der1(mx,1:ny),ierr)
       SLL_ALLOCATE( gtau_der2(my,1:nx+mx),ierr)
       gtau = data_array(1:nx,1:ny)
       taux_deriv(1) = 1
       taux_deriv(2) = nx
       gtau_der1(1,1:ny)     = interpolator%slope_min1(1:ny)
       gtau_der1(2,1:ny)     = interpolator%slope_max1(1:ny)
       tauy_deriv(1) = 1
       tauy_deriv(2) = ny
       gtau_der2(1,1:nx+mx)=interpolator%slope_min2(1:nx+mx)
       gtau_der2(2,1:nx+mx)=interpolator%slope_max2(1:nx+mx)
       call spli2d_custom_derder(&
            nx,&
            mx,&
            kx, &
            taux, &
            taux_deriv,&
            ny, &
            my,&
            ky, tauy, &
            tauy_deriv,&
            gtau,&
            gtau_der1,&
            gtau_der2,&
            interpolator%coeff_splines,&!(1:nx,1:ny),&
            interpolator%t1,&!(1:nx+kx), &
            interpolator%t2)!(1:ny+ky) )

       SLL_DEALLOCATE( gtau_der1,ierr)
       SLL_DEALLOCATE( gtau_der2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:ny+my)   = interpolator%value_min1(1:ny+my)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:nx,1)   = interpolator%value_min2(1:nx)
  !     interpolator%coeff_splines(1:nx,ny) = interpolator%value_max2(1:nx)

    case(2316)  !left: Hermite, right: Dirichlet, bottom: Hermite, Top: Hermite
       mx = 2
       my = 2
       interpolator%size_coeffs1 = nx + mx
       interpolator%size_coeffs2 = ny + my
       interpolator%size_t1 = kx + nx + mx
       interpolator%size_t2 = ky + ny + my
       
       !  data_array must have the same dimension than 
       !  size(  taux ) x  size(  tauy )
       !  i.e  data_array must have the dimension nx x ny
       SLL_ALLOCATE( gtau(1:nx,1:ny),ierr)
       SLL_ALLOCATE( gtau_der1(mx,1:ny),ierr)
       SLL_ALLOCATE( gtau_der2(my,1:nx+mx),ierr)
       gtau = data_array(1:nx,1:ny)
       taux_deriv(1) = 1
       taux_deriv(2) = nx
       gtau_der1(1,1:ny)     = interpolator%slope_min1(1:ny)
       gtau_der1(2,1:ny)     = interpolator%slope_max1(1:ny)
       tauy_deriv(1) = 1
       tauy_deriv(2) = ny
       gtau_der2(1,1:nx+mx)=interpolator%slope_min2(1:nx+mx)
       gtau_der2(2,1:nx+mx)=interpolator%slope_max2(1:nx+mx)
       call spli2d_custom_derder(&
            nx,&
            mx,&
            kx, &
            taux, &
            taux_deriv,&
            ny, &
            my,&
            ky, tauy, &
            tauy_deriv,&
            gtau,&
            gtau_der1,&
            gtau_der2,&
            interpolator%coeff_splines,&!(1:nx,1:ny),&
            interpolator%t1,&!(1:nx+kx), &
            interpolator%t2)!(1:ny+ky) )

       SLL_DEALLOCATE( gtau_der1,ierr)
       SLL_DEALLOCATE( gtau_der2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:ny+my)   = interpolator%value_min1(1:ny+my)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:nx,1)   = interpolator%value_min2(1:nx)
  !     interpolator%coeff_splines(1:nx,ny) = interpolator%value_max2(1:nx)

    case(2340) ! Hermite in al sides
       

       mx = 2
       my = 2
       interpolator%size_coeffs1 = nx + mx
       interpolator%size_coeffs2 = ny + my
       interpolator%size_t1 = kx + nx + mx
       interpolator%size_t2 = ky + ny + my
       
       !  data_array must have the same dimension than 
       !  size(  taux ) x  size(  tauy )
       !  i.e  data_array must have the dimension nx x ny
       SLL_ALLOCATE( gtau(1:nx,1:ny),ierr)
       SLL_ALLOCATE( gtau_der1(mx,1:ny),ierr)
       SLL_ALLOCATE( gtau_der2(my,1:nx+mx),ierr)
       gtau = data_array(1:nx,1:ny)
       taux_deriv(1) = 1
       taux_deriv(2) = nx
       gtau_der1(1,1:ny)     = interpolator%slope_min1(1:ny)
       gtau_der1(2,1:ny)     = interpolator%slope_max1(1:ny)
       tauy_deriv(1) = 1
       tauy_deriv(2) = ny
       gtau_der2(1,1:nx+mx)=interpolator%slope_min2(1:nx+mx)
       gtau_der2(2,1:nx+mx)=interpolator%slope_max2(1:nx+mx)
       call spli2d_custom_derder(&
            nx,&
            mx,&
            kx, &
            taux, &
            taux_deriv,&
            ny, &
            my,&
            ky, tauy, &
            tauy_deriv,&
            gtau,&
            gtau_der1,&
            gtau_der2,&
            interpolator%coeff_splines,&!(1:nx,1:ny),&
            interpolator%t1,&!(1:nx+kx), &
            interpolator%t2)!(1:ny+ky) )

       SLL_DEALLOCATE( gtau_der1,ierr)
       SLL_DEALLOCATE( gtau_der2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:ny+my)   = interpolator%value_min1(1:ny+my)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:nx,1)   = interpolator%value_min2(1:nx)
  !     interpolator%coeff_splines(1:nx,ny) = interpolator%value_max2(1:nx)

    end select
    interpolator%coefficients_set = .true.
   
    SLL_DEALLOCATE(taux,ierr)
    SLL_DEALLOCATE(tauy,ierr)

  end subroutine !compute_interpolants_ad2d

  function coefficients_are_set_ad2d( interpolator ) result(res)
    class(sll_arbitrary_degree_spline_interpolator_2d), intent(in)  :: interpolator
    logical :: res
    res = interpolator%coefficients_set
  end function coefficients_are_set_ad2d


  !  ----------------------------------------------------------
  !  Interpolation on the points eta1 and eta2 
  !  ---------------------------------------------------------
  !> @brief Interpolation on the points eta1 and eta2 
  !> @details computing the values with the interpolator arbitrary degree splines 2d
  !>  on the points eta1 and eta2 of arbitrary degree splines 2d
  !> 
  !> The parameters are
  !> @param interpolator the type sll_arbitrary_degree_spline_interpolator_2d
  !> @param[in] eta1 the point inthe first direction
  !> @param[in] eta2 the point inthe second direction 
  !> @return val the values on the points eta1 and eta2 
  function interpolate_value_ad2d( &
    interpolator, &
    eta1, &
    eta2 ) result(val)

    use sll_timer
    class(sll_arbitrary_degree_spline_interpolator_2d), intent(in)  :: interpolator
    sll_real64, intent(in)         :: eta1
    sll_real64, intent(in)         :: eta2
    sll_real64                     :: val
    sll_int32 :: size_coeffs1
    sll_int32 :: size_coeffs2
    sll_real64 :: res1,res2
    sll_real64,dimension(:), pointer:: tmp_tx,tmp_ty
    sll_real64,dimension(:,:), pointer:: tmp_coeff

    size_coeffs1 = interpolator%size_coeffs1
    size_coeffs2 = interpolator%size_coeffs2

    res1 = eta1
    res2 = eta2

    select case (interpolator%bc_selector)
    case (0) ! periodic-periodic

       if( res1 < interpolator%eta1_min ) then
          res1 = res1 + (interpolator%eta1_max-interpolator%eta1_min)
       else if( res1 >  interpolator%eta1_max ) then
          res1 = res1 + (interpolator%eta2_min-interpolator%eta2_max)
       end if
       if( res2 < interpolator%eta2_min ) then
          res2 = res2 + (interpolator%eta2_max-interpolator%eta2_min)
       else if( res2 >  interpolator%eta2_max ) then
          res2 = res2 + (interpolator%eta2_min-interpolator%eta2_max)
       end if

    case (9) ! 2. dirichlet-left, dirichlet-right, periodic

       if( res2 < interpolator%eta2_min ) then
          res2 = res2 + (interpolator%eta2_max-interpolator%eta2_min)
       else if( res2 >  interpolator%eta2_max ) then
          res2 = res2 + (interpolator%eta2_min-interpolator%eta2_max)
       end if

       SLL_ASSERT( res1 >= interpolator%eta1_min )
       SLL_ASSERT( res1 <= interpolator%eta1_max )
  
    case(576) !  3. periodic, dirichlet-bottom, dirichlet-top

       if( res1 < interpolator%eta1_min ) then
          res1 = res1 + ( interpolator%eta1_max-interpolator%eta1_min)
       else if( res1 >  interpolator%eta1_max ) then
          res1 = res1 + (interpolator%eta2_min-interpolator%eta2_max)
       end if
       SLL_ASSERT( res2 >= interpolator%eta2_min )
       SLL_ASSERT( res2 <= interpolator%eta2_max )
 
    end select

    SLL_ASSERT( res1 >= interpolator%eta1_min )
    SLL_ASSERT( res1 <= interpolator%eta1_max )
    SLL_ASSERT( res2 >= interpolator%eta2_min )
    SLL_ASSERT( res2 <= interpolator%eta2_max )

    tmp_tx => interpolator%t1(1:interpolator%size_t1)
    tmp_ty => interpolator%t2(1:interpolator%size_t2)
    tmp_coeff =>interpolator%coeff_splines(1:size_coeffs1,1:size_coeffs2)

    call bvalue2d( &
         res1, &
         res2, &
         size_coeffs1, &
         interpolator%spline_degree1+1, &
         size_coeffs2, &
         interpolator%spline_degree2+1, &
         tmp_coeff, &
         tmp_tx, &
         tmp_ty,&
         val)

  end function interpolate_value_ad2d


  !> @brief First derivative in eta1 interpolation on the points eta1 and eta2 
  !> @details computing the values of the first derivative in eta1
  !> with the interpolator arbitrary degree splines 2d
  !> on the points eta1 and eta2 of arbitrary degree splines 2d
  !> 
  !> The parameters are
  !> @param interpolator the type sll_arbitrary_degree_spline_interpolator_2d
  !> @param[in] eta1 the point inthe first direction
  !> @param[in] eta2 the point inthe second direction 
  !> @return val the values on the points eta1 and eta2 of the first derivative in eta1
  function interpolate_derivative1_ad2d( &
    interpolator, &
    eta1, &
    eta2 ) result(val)

    class(sll_arbitrary_degree_spline_interpolator_2d), intent(in)  :: interpolator
    sll_real64, intent(in)         :: eta1
    sll_real64, intent(in)         :: eta2
    sll_real64                     :: val
    sll_int32 :: size_coeffs1
    sll_int32 :: size_coeffs2
    !sll_real64 :: dvalue2d
    sll_real64 :: res1,res2
    sll_real64, dimension(:),pointer :: knot1_tmp
    sll_real64, dimension(:),pointer :: knot2_tmp
    sll_real64, dimension(:,:),pointer :: tmp_coeff
    !sll_int32 :: ierr

    SLL_ASSERT( eta1 .ge. interpolator%eta1_min )
    SLL_ASSERT( eta1 .le. interpolator%eta1_max )
    SLL_ASSERT( eta2 .ge. interpolator%eta2_min )
    SLL_ASSERT( eta2 .le. interpolator%eta2_max )
    
    size_coeffs1 = interpolator%size_coeffs1
    size_coeffs2 = interpolator%size_coeffs2

    res1 = eta1
    res2 = eta2
    
    select case (interpolator%bc_selector)
    case (0) ! periodic-periodic

       if( res1 < interpolator%eta1_min ) then
          res1 = res1 + ( interpolator%eta1_max-interpolator%eta1_min)
       else if( res1 >  interpolator%eta1_max ) then
          res1 = res1 + (interpolator%eta2_min-interpolator%eta2_max)
       end if
       if( res2 < interpolator%eta2_min ) then
          res2 = res2 + (interpolator%eta2_max-interpolator%eta2_min)
       else if( res2 >  interpolator%eta2_max ) then
          res2 = res2 + (interpolator%eta2_min-interpolator%eta2_max)
       end if
       

    case (9) ! 2. dirichlet-left, dirichlet-right, periodic

      
       if( res2 < interpolator%eta2_min ) then
          res2 = res2 + (interpolator%eta2_max-interpolator%eta2_min)
       else if( res2 >  interpolator%eta2_max ) then
          res2 = res2 + (interpolator%eta2_min-interpolator%eta2_max)
       end if

       SLL_ASSERT( res1 >= interpolator%eta1_min )
       SLL_ASSERT( res1 <= interpolator%eta1_max )
       
    case(576) !  3. periodic, dirichlet-bottom, dirichlet-top

       if( res1 < interpolator%eta1_min ) then
          res1 = res1 + (interpolator%eta1_max-interpolator%eta1_min)
       else if( res1 >  interpolator%eta1_max ) then
          res1 = res1 + (interpolator%eta2_min-interpolator%eta2_max)
       end if

       SLL_ASSERT( res2 >= interpolator%eta2_min )
       SLL_ASSERT( res2 <= interpolator%eta2_max )
       
    end select

    SLL_ASSERT( res1 >= interpolator%eta1_min )
    SLL_ASSERT( res1 <= interpolator%eta1_max )
    SLL_ASSERT( res2 >= interpolator%eta2_min )
    SLL_ASSERT( res2 <= interpolator%eta2_max )

    knot1_tmp => interpolator%t1(1:interpolator%size_t1)
    knot2_tmp => interpolator%t2(1:interpolator%size_t2)
    tmp_coeff => interpolator%coeff_splines(1:size_coeffs1,1:size_coeffs2)

    val = dvalue2d( &
         res1, &
         res2, &
         size_coeffs1, &
         interpolator%spline_degree1+1, &
         size_coeffs2, &
         interpolator%spline_degree2+1, &
         tmp_coeff, &
         knot1_tmp, &
         knot2_tmp,&
         1,0)
    
    !SLL_DEALLOCATE(knot1_tmp,ierr)
    !SLL_DEALLOCATE(knot2_tmp,ierr)
    
  end function interpolate_derivative1_ad2d
     
  !> @brief First derivative in eta2 Interpolation on the points eta1 and eta2 
  !> using the arbitrary degree splines interpolator 2d
  !> @details computing the values of the first derivative in eta2
  !> with the interpolator arbitrary degree splines 2d
  !> on the points eta1 and eta2 of arbitrary degree splines 2d
  !> 
  !> The parameters are
  !> @param interpolator the type sll_arbitrary_degree_spline_interpolator_2d
  !> @param[in] eta1 the point inthe first direction
  !> @param[in] eta2 the point inthe second direction 
  !> @return val the values on the points eta1 and eta2 of the first derivative in eta2
  function interpolate_derivative2_ad2d( &
    interpolator, &
    eta1, &
    eta2 ) result(val)

    class(sll_arbitrary_degree_spline_interpolator_2d), intent(in)  :: interpolator
    sll_real64, intent(in)         :: eta1
    sll_real64, intent(in)         :: eta2
    sll_real64                     :: val
    sll_int32 :: size_coeffs1
    sll_int32 :: size_coeffs2
    sll_real64 :: res1,res2
    sll_real64, dimension(:),pointer :: knot1_tmp
    sll_real64, dimension(:),pointer :: knot2_tmp
    sll_real64, dimension(:,:),pointer :: tmp_coeff

    SLL_ASSERT( eta1 .ge. interpolator%eta1_min )
    SLL_ASSERT( eta1 .le. interpolator%eta1_max )
    SLL_ASSERT( eta2 .ge. interpolator%eta2_min )
    SLL_ASSERT( eta2 .le. interpolator%eta2_max )

    size_coeffs1 = interpolator%size_coeffs1
    size_coeffs2 = interpolator%size_coeffs2

    res1 = eta1
    res2 = eta2
    select case (interpolator%bc_selector)
    case (0) ! periodic-periodic
       
       if( res1 < interpolator%eta1_min ) then
          res1 = res1 + (interpolator%eta1_max-interpolator%eta1_min)
       else if( res1 >  interpolator%eta1_max ) then
          res1 = res1 + (interpolator%eta2_min-interpolator%eta2_max)
       end if
       if( res2 < interpolator%eta2_min ) then
          res2 = res2 + (interpolator%eta2_max-interpolator%eta2_min)
       else if( res2 >  interpolator%eta2_max ) then
          res2 = res2 + (interpolator%eta2_min-interpolator%eta2_max)
       end if
          
    case (9) ! 2. dirichlet-left, dirichlet-right, periodic
       
      
       if( res2 < interpolator%eta2_min ) then
          res2 = res2 + (interpolator%eta2_max-interpolator%eta2_min)
       else if( res2 >  interpolator%eta2_max ) then
          res2 = res2 + (interpolator%eta2_min-interpolator%eta2_max)
       end if
       SLL_ASSERT( res1 >= interpolator%eta1_min )
       SLL_ASSERT( res1 <= interpolator%eta1_max )
       
    case(576) !  3. periodic, dirichlet-bottom, dirichlet-top
       
       if( res1 < interpolator%eta1_min ) then
          res1 = res1+ ( interpolator%eta1_max-interpolator%eta1_min)
       else if( res1 >  interpolator%eta1_max ) then
          res1 = res1 + (interpolator%eta2_min-interpolator%eta2_max)
       end if

       SLL_ASSERT( res2 >= interpolator%eta2_min )
       SLL_ASSERT( res2 <= interpolator%eta2_max )

    end select
    SLL_ASSERT( res1 >= interpolator%eta1_min )
    SLL_ASSERT( res1 <= interpolator%eta1_max )
    SLL_ASSERT( res2 >= interpolator%eta2_min )
    SLL_ASSERT( res2 <= interpolator%eta2_max )

    knot1_tmp => interpolator%t1(1:interpolator%size_t1)
    knot2_tmp => interpolator%t2(1:interpolator%size_t2)
    tmp_coeff =>interpolator%coeff_splines(1:size_coeffs1,1:size_coeffs2)
    val = dvalue2d( &
         res1, &
         res2, &
         size_coeffs1, &
         interpolator%spline_degree1+1, &
         size_coeffs2, &
         interpolator%spline_degree2+1, &
         tmp_coeff, &
         knot1_tmp, &
         knot2_tmp,&
         0,1)
    
   ! SLL_DEALLOCATE(knot1_tmp,ierr)
   ! SLL_DEALLOCATE(knot2_tmp,ierr)
  end function interpolate_derivative2_ad2d !interpolate_derivative2_ad2d

  
  function interpolate_array_ad2d( &
  this, &
  num_points1, &
  num_points2, &
  data_in, &
  eta1, &
  eta2 ) result(res)
    
    class(sll_arbitrary_degree_spline_interpolator_2d), intent(in)  :: this
    sll_real64,  dimension(:,:), intent(in)         :: eta1
    sll_real64,  dimension(:,:), intent(in)         :: eta2
    sll_real64, dimension(:,:), intent(in)         :: data_in
    sll_int32, intent(in)         :: num_points1
    sll_int32, intent(in)         :: num_points2

    sll_real64, dimension(num_points1,num_points2) :: res
    
    print *, '#interpolate_array_ad2d: not implemented'
    res = -1000000._f64
    print *,this%num_pts1
    print *,maxval(eta1)
    print *,maxval(eta2)
    print *,maxval(data_in)
    print *,num_points1
    print *,num_points2
    stop
  end function !interpolate_array_ad2d
  
  function interpolate_2d_array_disp_ad2d( &
       this,        &
       num_points1, &
       num_points2, &
       data_in,     &
       alpha1,      &
       alpha2) result(res)
      
    class(sll_arbitrary_degree_spline_interpolator_2d), intent(in)    :: this
    sll_int32, intent(in)                          :: num_points1  
    sll_int32, intent(in)                          :: num_points2 
    sll_real64, dimension(:,:), intent(in)         :: data_in
    sll_real64, dimension(:,:), intent(in)         :: alpha1
    sll_real64, dimension(:,:), intent(in)         :: alpha2  
    sll_real64, dimension(num_points1,num_points2) :: res
    
    
    
    print *, '#interpolate_2d_array_disp_ad2d: not implemented.'
    !for preventing warning of unused objects
    print *,this%num_pts1
    print *,num_points1 
    print *,num_points2
    print *,maxval(data_in)
    print *,alpha1
    print *,alpha2     
    res = -1000000._f64
    stop
    
  end function !interpolate_2d_array_disp_ad2d
    
   
  
  function get_coefficients_ad2d(interpolator)
    class(sll_arbitrary_degree_spline_interpolator_2d), intent(in) :: interpolator
    sll_real64, dimension(:,:), pointer           :: get_coefficients_ad2d     

    get_coefficients_ad2d => interpolator%coeff_splines
  end function get_coefficients_ad2d
  
  !> Initialization of the boundary for interpolator arbitrary degree splines 2d.
  !> The parameters are
  !> @param[in] slope_min1 a 1d arrays contains values in the left in the direction eta1  
  !> @param[in] slope_max1 a 1d arrays contains values in the right in the direction eta1 
  !> @param[in] slope_min2 a 1d arrays contains values in the left in the direction eta2 
  !> @param[in] slope_max2 a 1d arrays contains values in the right in the direction eta2
  !> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_2d
  subroutine set_slope2d(&
       interpolator,&
       slope_min1,&
       slope_max1,&
       slope_min2,&
       slope_max2)

    use sll_module_arbitrary_degree_spline_interpolator_1d
    class(sll_arbitrary_degree_spline_interpolator_2d)    :: interpolator
    sll_real64, dimension(:),optional :: slope_min1
    sll_real64, dimension(:),optional :: slope_max1
    sll_real64, dimension(:),optional :: slope_min2
    sll_real64, dimension(:),optional :: slope_max2
    class(sll_arbitrary_degree_spline_interpolator_1d),pointer :: interp1d_min2=> null()
    class(sll_arbitrary_degree_spline_interpolator_1d),pointer :: interp1d_max2 => null()
    sll_int32 :: sz_slope_min2,sz_slope_max2
    sll_int64 :: bc_selector
    sll_int32 :: num_pts1
    sll_int32 :: num_pts2
    sll_int32 :: bc_min1
    sll_int32 :: bc_max1
    sll_int32 :: bc_min2
    sll_int32 :: bc_max2

    num_pts1 = interpolator%num_pts1
    num_pts2 = interpolator%num_pts2
    bc_selector = interpolator%bc_selector
    bc_min1  = interpolator%bc_min1 
    bc_max1 = interpolator%bc_max1 
    bc_min2= interpolator%bc_min2  
    bc_max2   = interpolator%bc_max2

    select case (bc_selector)
    case(0)
    case (9) ! dirichlet-left, dirichlet-right, periodic
    case (576) ! 3. periodic, dirichlet-bottom, dirichlet-top
    case (585) ! 4. dirichlet in all sides
    case (650) !left: Neumann, right: Dirichlet, bottom: Neumann, Top: Dirichlet
       !if ( present( slope_min1)) then 
       interpolator%slope_min1 = 0.0
       interpolator%compute_slope_min1= .FALSE.
       !end if
       if ( present( slope_max1)) then 
          interpolator%slope_max1 = slope_max1
          interpolator%compute_slope_max1= .FALSE.
       end if

       interpolator%slope_min2 = 0.0_f64
       interpolator%compute_slope_min2= .FALSE.
       

       if ( present( slope_max2)) then 
          interpolator%compute_slope_max2= .FALSE.

          
          sz_slope_max2 = size(slope_max2)
          if ( sz_slope_max2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_max2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_max2 => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_max2%compute_interpolants(&
               slope_max2(1:sz_slope_max2))
          
          interpolator%slope_max2(1:sz_slope_max2+2) = &
               interp1d_max2%coeff_splines(1:sz_slope_max2+2)
          call sll_delete(interp1d_max2)
       else
          print*, 'problem with slope top in case 780'
       end if
    case(657) !left: Dirichlet, right: Neumann, bottom: Neumann, Top: Dirichlet 
       if ( present( slope_min1)) then 
          interpolator%slope_min1 = slope_min1
          interpolator%compute_slope_min1= .FALSE.
       end if
       
       
       sz_slope_min2 = size(slope_min2)
       interpolator%slope_max1 = 0.0_f64
       interpolator%compute_slope_max1= .FALSE.
       
       interpolator%slope_min2(1:sz_slope_min2+2) = 0.0_f64
       interpolator%compute_slope_min2 = .FALSE.

       if ( present( slope_max2)) then 
          interpolator%compute_slope_max2= .FALSE.


          sz_slope_max2 = size(slope_max2)
          if ( sz_slope_max2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_max2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_max2 => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_max2%compute_interpolants(&
               slope_max2(1:sz_slope_max2))
          
          interpolator%slope_max2(1:sz_slope_max2+2) = &
               interp1d_max2%coeff_splines(1:sz_slope_max2+2)
          call sll_delete(interp1d_max2)
       else
          print*, 'problem with slope top in case 780'
       end if

    case(780)  !left: Hermite, right: Dirichlet, bottom: Hermite, Top: Dirichlet
       if ( present( slope_min1)) then 
          interpolator%slope_min1 = slope_min1
          interpolator%compute_slope_min1= .FALSE.
       end if
        if ( present( slope_max1)) then 
          interpolator%slope_max1 = slope_max1
          interpolator%compute_slope_max1= .FALSE.
       end if
       
       if (present(slope_min2)) then 
           sz_slope_min2 = size(slope_min2)
          if ( sz_slope_min2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_min2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_min2 =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_min2%compute_interpolants( &
               slope_min2(1:sz_slope_min2))
          
          interpolator%slope_min2(1:sz_slope_min2+2) = &
               interp1d_min2%coeff_splines(1:sz_slope_min2+2)
          call sll_delete(interp1d_min2)
          interpolator%compute_slope_min2 = .FALSE.
       else
          print*, 'problem with slope bottom in case 780'
       end if
       if ( present( slope_max2)) then 
          interpolator%compute_slope_max2= .FALSE.


          sz_slope_max2 = size(slope_max2)
          if ( sz_slope_max2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_max2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_max2 => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_max2%compute_interpolants(&
               slope_max2(1:sz_slope_max2))
          
          interpolator%slope_max2(1:sz_slope_max2+2) = &
               interp1d_max2%coeff_splines(1:sz_slope_max2+2)
          call sll_delete(interp1d_max2)
       else
          print*, 'problem with slope top in case 780'
       end if

    case(801)  !left: Dirichlet, right: Hermite, bottom: Hermite, Top: Dirichlet
        if ( present( slope_max1)) then 
           interpolator%slope_max1 = slope_max1
           interpolator%compute_slope_max1= .FALSE.
       end if

        if ( present( slope_min1)) then 
          interpolator%slope_min1 = slope_min1
          interpolator%compute_slope_min1= .FALSE.
       end if
       
       if (present(slope_min2)) then 
           sz_slope_min2 = size(slope_min2)
          if ( sz_slope_min2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_min2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_min2 =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_min2%compute_interpolants( &
               slope_min2(1:sz_slope_min2))
          
          interpolator%slope_min2(1:sz_slope_min2+2) = &
               interp1d_min2%coeff_splines(1:sz_slope_min2+2)
          call sll_delete(interp1d_min2)
          interpolator%compute_slope_min2 = .FALSE.
       else
          print*, 'problem with slope bottom in case 801'
       end if
       if ( present( slope_max2)) then 
          interpolator%compute_slope_max2= .FALSE.


          sz_slope_max2 = size(slope_max2)
          if ( sz_slope_max2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_max2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_max2 => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_max2%compute_interpolants(&
               slope_max2(1:sz_slope_max2))
          
          interpolator%slope_max2(1:sz_slope_max2+2) = &
               interp1d_max2%coeff_splines(1:sz_slope_max2+2)
          call sll_delete(interp1d_max2)
       else
          print*, 'problem with slope top in case 801'
       end if

    case(804)  !left: Hermite, right: Hermite, bottom: Hermite, Top: Dirichlet

        if ( present( slope_max1)) then 
           interpolator%slope_max1 = slope_max1
           interpolator%compute_slope_max1= .FALSE.
       end if

        if ( present( slope_min1)) then 
          interpolator%slope_min1 = slope_min1
          interpolator%compute_slope_min1= .FALSE.
       end if
       
       if (present(slope_min2)) then 
           sz_slope_min2 = size(slope_min2)
          if ( sz_slope_min2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_min2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_min2 =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_min2%compute_interpolants( &
               slope_min2(1:sz_slope_min2))
          
          interpolator%slope_min2(1:sz_slope_min2+2) = &
               interp1d_min2%coeff_splines(1:sz_slope_min2+2)
          call sll_delete(interp1d_min2)
          interpolator%compute_slope_min2 = .FALSE.
       else
          print*, 'problem with slope bottom in case 801'
       end if
       if ( present( slope_max2)) then 
          interpolator%compute_slope_max2= .FALSE.


          sz_slope_max2 = size(slope_max2)
          if ( sz_slope_max2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_max2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_max2 => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_max2%compute_interpolants(&
               slope_max2(1:sz_slope_max2))
          
          interpolator%slope_max2(1:sz_slope_max2+2) = &
               interp1d_max2%coeff_splines(1:sz_slope_max2+2)
          call sll_delete(interp1d_max2)
       else
          print*, 'problem with slope top in case 801'
       end if


    case(1098)  !left: Neumann, right: Dirichlet, bottom: Dirichlet, Top: Neumann
       interpolator%slope_min1 = 0.0_f64
       interpolator%compute_slope_min1= .FALSE.
    
       if ( present( slope_max1)) then 
          interpolator%slope_max1 = slope_max1
          interpolator%compute_slope_max1= .FALSE.
       end if
       
       sz_slope_max2= size(slope_max2)
       interpolator%slope_max2(1:sz_slope_max2+2) = 0.0_f64
       interpolator%compute_slope_max2 = .FALSE.

       if (present(slope_min2)) then 
           sz_slope_min2 = size(slope_min2)
          if ( sz_slope_min2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_min2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_min2 =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_min2%compute_interpolants( &
               slope_min2(1:sz_slope_min2))
          
          interpolator%slope_min2(1:sz_slope_min2+2) = &
               interp1d_min2%coeff_splines(1:sz_slope_min2+2)
          call sll_delete(interp1d_min2)
          interpolator%compute_slope_min2 = .FALSE.
       else
          print*, 'problem with slope bottom in case 2124'
       end if
       if ( present( slope_max2)) then 
          interpolator%compute_slope_max2= .FALSE.


          sz_slope_max2 = size(slope_max2)
          if ( sz_slope_max2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_max2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_max2 => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_max2%compute_interpolants(&
               slope_max2(1:sz_slope_max2))
          
          interpolator%slope_max2(1:sz_slope_max2+2) = &
               interp1d_max2%coeff_splines(1:sz_slope_max2+2)
          call sll_delete(interp1d_max2)
       else
          print*, 'problem with slope top in case 2124'
       end if
    case(1105)  !left: Dirichlet, right: Neumann, bottom: Dirichlet, Top: Neumann
       if ( present( slope_max1)) then 
          interpolator%slope_max1 = slope_max1
          interpolator%compute_slope_max1= .FALSE.
       end if
       interpolator%slope_min1 = 0.0_f64
       interpolator%compute_slope_min1= .FALSE.
       
       sz_slope_max2 = size(slope_max2)
       interpolator%compute_slope_max2= .FALSE.
       interpolator%slope_max2(1:sz_slope_max2+2) = 0.0_f64

       if (present(slope_min2)) then 
           sz_slope_min2 = size(slope_min2)
          if ( sz_slope_min2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_min2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_min2 =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_min2%compute_interpolants( &
               slope_min2(1:sz_slope_min2))
          
          interpolator%slope_min2(1:sz_slope_min2+2) = &
               interp1d_min2%coeff_splines(1:sz_slope_min2+2)
          call sll_delete(interp1d_min2)
          interpolator%compute_slope_min2 = .FALSE.
       else
          print*, 'problem with slope bottom in case 2145'
       end if
    case(1170)  !left: Neumann, right: Neumann, bottom: Neuman, Top: Neumann
       
       interpolator%slope_max1 = 0.0_f64
       interpolator%compute_slope_max1= .FALSE.
       
       interpolator%slope_min1 = 0.0_f64
       interpolator%compute_slope_min1= .FALSE.
       
       sz_slope_max2 = size(slope_max2)
       interpolator%compute_slope_max2= .FALSE.
       interpolator%slope_max2(1:sz_slope_max2+2) = 0.0_f64
          
       sz_slope_min2 = size(slope_min2)
       interpolator%slope_min2(1:sz_slope_min2+2) = 0.0_f64
       interpolator%compute_slope_min2 = .FALSE.

    case(2338)  !left: Dirichlet, right: Hermite, bottom: Hermite, Top: Hermite
       if ( present( slope_max1)) then 
          interpolator%slope_max1 = slope_max1
          interpolator%compute_slope_max1= .FALSE.
       end if
       if ( present( slope_min1)) then 
          interpolator%slope_min1 = slope_min1
          interpolator%compute_slope_min1= .FALSE.
       end if
       
       if ( present( slope_max2)) then 
          interpolator%compute_slope_max2= .FALSE.

          
          sz_slope_max2 = size(slope_max2)
          if ( sz_slope_max2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_max2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_max2 => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_max2%compute_interpolants(&
               slope_max2(1:sz_slope_max2))
          
          interpolator%slope_max2(1:sz_slope_max2+2) = &
               interp1d_max2%coeff_splines(1:sz_slope_max2+2)
          call sll_delete(interp1d_max2)
       else
          print*, 'problem with slope top in case 2145'
       end if

       if (present(slope_min2)) then 
           sz_slope_min2 = size(slope_min2)
          if ( sz_slope_min2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_min2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_min2 =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_min2%compute_interpolants( &
               slope_min2(1:sz_slope_min2))
          
          interpolator%slope_min2(1:sz_slope_min2+2) = &
               interp1d_min2%coeff_splines(1:sz_slope_min2+2)
          call sll_delete(interp1d_min2)
          interpolator%compute_slope_min2 = .FALSE.
       else
          print*, 'problem with slope bottom in case 2145'
       end if

    case(2145)  !left: Dirichlet, right: Hermite, bottom: Dirichlet, Top: Hermite
       if ( present( slope_max1)) then 
          interpolator%slope_max1 = slope_max1
          interpolator%compute_slope_max1= .FALSE.
       end if
        if ( present( slope_min1)) then 
          interpolator%slope_min1 = slope_min1
          interpolator%compute_slope_min1= .FALSE.
       end if
       
       if ( present( slope_max2)) then 
          interpolator%compute_slope_max2= .FALSE.

          
          sz_slope_max2 = size(slope_max2)
          if ( sz_slope_max2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_max2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_max2 => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_max2%compute_interpolants(&
               slope_max2(1:sz_slope_max2))
          
          interpolator%slope_max2(1:sz_slope_max2+2) = &
               interp1d_max2%coeff_splines(1:sz_slope_max2+2)
          call sll_delete(interp1d_max2)
       else
          print*, 'problem with slope top in case 2145'
       end if

       if (present(slope_min2)) then 
           sz_slope_min2 = size(slope_min2)
          if ( sz_slope_min2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_min2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_min2 =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_min2%compute_interpolants( &
               slope_min2(1:sz_slope_min2))
          
          interpolator%slope_min2(1:sz_slope_min2+2) = &
               interp1d_min2%coeff_splines(1:sz_slope_min2+2)
          call sll_delete(interp1d_min2)
          interpolator%compute_slope_min2 = .FALSE.
       else
          print*, 'problem with slope bottom in case 2145'
       end if

       
    case(2124)  !left: Hermite, right: Dirichlet, bottom: Dirichlet, Top: Hermite
       
       if ( present( slope_min1)) then 
          interpolator%slope_min1 = slope_min1
          interpolator%compute_slope_min1= .FALSE.
       end if
        if ( present( slope_max1)) then 
          interpolator%slope_max1 = slope_max1
          interpolator%compute_slope_max1= .FALSE.
       end if
       
       if ( present( slope_max2)) then 
          interpolator%compute_slope_max2= .FALSE.


          sz_slope_max2 = size(slope_max2)
          if ( sz_slope_max2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_max2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_max2 => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_max2%compute_interpolants(&
               slope_max2(1:sz_slope_max2))
          
          interpolator%slope_max2(1:sz_slope_max2+2) = &
               interp1d_max2%coeff_splines(1:sz_slope_max2+2)

          call sll_delete(interp1d_max2)
       else
          print*, 'problem with slope top in case 2124'
       end if

       if (present(slope_min2)) then 
           sz_slope_min2 = size(slope_min2)
          if ( sz_slope_min2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_min2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_min2 =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_min2%compute_interpolants( &
               slope_min2(1:sz_slope_min2))
          
          interpolator%slope_min2(1:sz_slope_min2+2) = &
               interp1d_min2%coeff_splines(1:sz_slope_min2+2)
          call sll_delete(interp1d_min2)
          interpolator%compute_slope_min2 = .FALSE.
       else
          print*, 'problem with slope bottom in case 2124'
       end if
       if ( present( slope_max2)) then 
          interpolator%compute_slope_max2= .FALSE.


          sz_slope_max2 = size(slope_max2)
          if ( sz_slope_max2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_max2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_max2 => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_max2%compute_interpolants(&
               slope_max2(1:sz_slope_max2))
          
          interpolator%slope_max2(1:sz_slope_max2+2) = &
               interp1d_max2%coeff_splines(1:sz_slope_max2+2)
          call sll_delete(interp1d_max2)
       else
          print*, 'problem with slope top in case 2124'
       end if
      
    case(2148)  !left:Hermite , right: Hermite, bottom: Dirichlet, Top: Hermite  
       if ( present( slope_min1)) then 
          interpolator%slope_min1 = slope_min1
          interpolator%compute_slope_min1= .FALSE.
       end if
       if ( present( slope_max1)) then 
          interpolator%slope_max1 = slope_max1
          interpolator%compute_slope_max1= .FALSE.
       end if
       
       if ( present( slope_max2)) then 
          interpolator%compute_slope_max2= .FALSE.


          sz_slope_max2 = size(slope_max2)
          if ( sz_slope_max2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_max2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_max2 => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_max2%compute_interpolants(&
               slope_max2(1:sz_slope_max2))
          
          interpolator%slope_max2(1:sz_slope_max2+2) = &
               interp1d_max2%coeff_splines(1:sz_slope_max2+2)

          call sll_delete(interp1d_max2)
       else
          print*, 'problem with slope top in case 2124'
       end if

       if (present(slope_min2)) then 
           sz_slope_min2 = size(slope_min2)
          if ( sz_slope_min2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_min2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_min2 =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_min2%compute_interpolants( &
               slope_min2(1:sz_slope_min2))
          
          interpolator%slope_min2(1:sz_slope_min2+2) = &
               interp1d_min2%coeff_splines(1:sz_slope_min2+2)
          call sll_delete(interp1d_min2)
          interpolator%compute_slope_min2 = .FALSE.
       else
          print*, 'problem with slope bottom in case 2124'
       end if
       if ( present( slope_max2)) then 
          interpolator%compute_slope_max2= .FALSE.


          sz_slope_max2 = size(slope_max2)
          if ( sz_slope_max2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_max2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_max2 => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_max2%compute_interpolants(&
               slope_max2(1:sz_slope_max2))
          
          interpolator%slope_max2(1:sz_slope_max2+2) = &
               interp1d_max2%coeff_splines(1:sz_slope_max2+2)
          call sll_delete(interp1d_max2)
       else
          print*, 'problem with slope top in case 2124'
       end if
    case(2316)  !left: Hermite, right: Dirichlet, bottom: Hermite, Top: Hermite
       if ( present( slope_max1)) then 
          interpolator%slope_max1 = slope_max1
          interpolator%compute_slope_max1= .FALSE.
       end if
       if ( present( slope_min1)) then 
          interpolator%slope_min1 = slope_min1
          interpolator%compute_slope_min1= .FALSE.
       end if
       
       if ( present( slope_max2)) then 
          interpolator%compute_slope_max2= .FALSE.


          sz_slope_max2 = size(slope_max2)
          if ( sz_slope_max2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_max2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_max2 => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_max2%compute_interpolants(&
               slope_max2(1:sz_slope_max2))
          
          interpolator%slope_max2(1:sz_slope_max2+2) = &
               interp1d_max2%coeff_splines(1:sz_slope_max2+2)
          call sll_delete(interp1d_max2)
       else
          print*, 'problem with slope top in case 2340'
       end if
       
       if (present(slope_min2)) then 
           sz_slope_min2 = size(slope_min2)
          if ( sz_slope_min2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_min2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_min2 =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_min2%compute_interpolants( &
               slope_min2(1:sz_slope_min2))
          
          interpolator%slope_min2(1:sz_slope_min2+2) = &
               interp1d_min2%coeff_splines(1:sz_slope_min2+2)
          call sll_delete(interp1d_min2)
          interpolator%compute_slope_min2 = .FALSE.
       else
          print*, 'problem with slope bottom in case 2340'
       end if
    case(2340) ! Hermite in al sides

       if ( present( slope_max1)) then 
          interpolator%slope_max1 = slope_max1
          interpolator%compute_slope_max1= .FALSE.
       end if
        if ( present( slope_min1)) then 
          interpolator%slope_min1 = slope_min1
          interpolator%compute_slope_min1= .FALSE.
       end if
       
       if ( present( slope_max2)) then 
          interpolator%compute_slope_max2= .FALSE.


          sz_slope_max2 = size(slope_max2)
          if ( sz_slope_max2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_max2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_max2 => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_max2%compute_interpolants(&
               slope_max2(1:sz_slope_max2))
          
          interpolator%slope_max2(1:sz_slope_max2+2) = &
               interp1d_max2%coeff_splines(1:sz_slope_max2+2)
          call sll_delete(interp1d_max2)
       else
          print*, 'problem with slope top in case 2340'
       end if
       
       if (present(slope_min2)) then 
           sz_slope_min2 = size(slope_min2)
          if ( sz_slope_min2 .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_min2 must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_min2 =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_min1, &
               interpolator%bc_max1, &
               interpolator%spline_degree1 )
          
          call interp1d_min2%compute_interpolants( &
               slope_min2(1:sz_slope_min2))
          
          interpolator%slope_min2(1:sz_slope_min2+2) = &
               interp1d_min2%coeff_splines(1:sz_slope_min2+2)
          call sll_delete(interp1d_min2)
          interpolator%compute_slope_min2 = .FALSE.
       else
          print*, 'problem with slope bottom in case 2340'
       end if
       
    case default
       print*,'initialize_ad2d_interpolator: BC combination not implemented.'
    end select
    
  end subroutine set_slope2d
  
subroutine spli2d_custom(nx, kx, taux, ny, ky, tauy, g, bcoef, tx, ty)

sll_int32,                           intent(in)  :: nx
sll_int32,                           intent(in)  :: kx
sll_int32,                           intent(in)  :: ny
sll_int32,                           intent(in)  :: ky
sll_real64, dimension(:),   pointer, intent(in)  :: taux
sll_real64, dimension(:),   pointer, intent(in)  :: tauy
sll_real64, dimension(:,:), pointer, intent(in)  :: g   

sll_real64, dimension(:,:), pointer, intent(out) :: bcoef
sll_real64, dimension(:),   pointer, intent(out) :: tx
sll_real64, dimension(:),   pointer, intent(out) :: ty

sll_real64, dimension(nx)                 :: work_x
sll_real64, dimension(nx*(2*kx-1))        :: qx
sll_real64, dimension((2*ky-1)*ny)        :: qy
sll_real64, dimension(ny)                 :: work_y
sll_real64, dimension(ny,nx),     target  :: bwork
sll_real64, dimension(:,:),       pointer :: pwork

sll_int32 :: i, j, flag
sll_int32 :: ierr

! *** set up knots and interpolate between knots

tx(1:kx)       = taux(1)
tx(nx+1:nx+kx) = taux(nx)
if (mod(kx,2) == 0) then
  do i = kx+1, nx
    tx(i) = taux(i-kx/2) 
  end do
else
  do i = kx+1, nx
    tx(i) = 0.5*(taux(i-(kx-1)/2)+taux(i-1-(kx-1)/2))
  end do
end if

ty = 0.0_f64
if (mod(ky,2) == 0) then
  do i = ky + 1, ny
    ty(i) = tauy(i-ky/2) 
  end do
else
  do i = ky + 1, ny
    ty(i) = 0.5*(tauy(i-(ky-1)/2)+tauy(i-1-(ky-1)/2))
   end do
end if

ty(1:ky)       = tauy(1)
ty(ny+1:ny+ky) = tauy(ny)

pwork => bwork
bcoef(1:nx,1:ny) = g

call spli2d( taux, bcoef, tx, nx, kx, ny, work_x, qx, pwork, flag)
call spli2d( tauy, pwork, ty, ny, ky, nx, work_y, qy, bcoef, flag)

end subroutine spli2d_custom

!*****************************************************************************80
!
! SPLI2D produces a interpolatory tensor product spline.
!
!  Discussion:
!
!    SPLI2D is an extended version of SPLINT.
!
!    SPLI2D produces the B-spline coefficients BCOEF(J,.) of the
!    spline of order K with knots T(1:N+K), which takes on
!    the value GTAU(I,J) at TAU(I), I=1,..., N, J=1,...,M.
!
!    The I-th equation of the linear system
!
!      A * BCOEF = B
!
!    for the B-spline coefficients of the interpolant enforces
!    interpolation at TAU(I), I=1,...,N.  Hence,  B(I) = GTAU(I),
!    for all I, and A is a band matrix with 2*K-1 bands, if it is
!    invertible.
!
!    The matrix A is generated row by row and stored, diagonal by
!    diagonal, in the rows of the array Q, with the main diagonal
!    going into row K.
!
!    The banded system is then solved by a call to BANFAC, which
!    constructs the triangular factorization for A and stores it
!    again in Q, followed by a call to BANSLV, which then obtains
!    the solution BCOEF by substitution.
!
!     The linear system to be solved is theoretically invertible if
!     and only if
!
!       T(I) < TAU(I) < TAU(I+K), for all I.
!
!     Violation of this condition is certain to lead to IFLAG = 2.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TAU(N), contains the data point abscissas.
!    TAU must be strictly increasing
!
!    Input, real ( kind = 8 ) GTAU(N,M), contains the data point ordinates.
!
!    Input, real ( kind = 8 ) T(N+K), the knot sequence.
!
!    Input, integer N, the number of data points and the
!    dimension of the spline space SPLINE(K,T)
!
!    Input, integer K, the order of the spline.
!
!    Input, integer M, the number of data sets.
!
!    Work space, real ( kind = 8 ) WORK(N).
!
!    Output, real ( kind = 8 ) Q(2*K-1)*N, the triangular
!    factorization of the coefficient matrix of the linear
!    system for the B-spline coefficients of the spline interpolant.
!    The B-spline coefficients for the interpolant of an additional
!    data set ( TAU(I), HTAU(I) ), I=1,...,N  with the same data
!    abscissae can be obtained without going through all the
!    calculations in this routine, simply by loading HTAU into
!    BCOEF and then using the statement
!      CALL BANSLV ( Q, 2*K-1, N, K-1, K-1, BCOEF )
!
!    Output, real ( kind = 8 ) BCOEF(N), the B-spline coefficients of
!    the interpolant.
!
!    Output, integer IFLAG, error indicator.
!    1, no error.
!    2, an error occurred, which may have been caused by
!       singularity of the linear system.
!
subroutine spli2d ( tau, gtau, t, n, k, m, work, q, bcoef, iflag )
    
sll_real64, dimension(:),   pointer, intent(in)  :: tau
sll_real64, dimension(:,:), pointer, intent(in)  :: gtau
sll_real64, dimension(:),   pointer, intent(in)  :: t
sll_int32                          , intent(in)  :: n
sll_int32                          , intent(in)  :: k
sll_int32                          , intent(in)  :: m
sll_real64, dimension(n)                         :: work
sll_real64, dimension((2*k-1)*n)                 :: q
sll_real64, dimension(:,:), pointer, intent(out) :: bcoef
sll_int32,                           intent(out) :: iflag

sll_int32 :: i
sll_int32 :: ilp1mx
sll_int32 :: j
sll_int32 :: jj
sll_int32 :: left
sll_real64:: taui

left = k

q(1:(2*k-1)*n) = 0.0_f64
!
!  Construct the N interpolation equations.
!

do i = 1, n
   
   taui = tau(i)
   ilp1mx = min ( i + k, n + 1 )
   !
   !  Find the index LEFT in the closed interval (I,I+K-1) such that:
   !
   !    T(LEFT) < = TAU(I) < T(LEFT+1)
   !
   !  The matrix will be singular if this is not possible.
   !
   left = max ( left, i )
   
   if ( taui < t(left) ) then
      iflag = 2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLI2D - Fatal error!'
      write ( *, '(a)' ) '  The TAU array is not strictly increasing .'
      !print*, taui, t(left),left
      stop
   end if
   
   do while ( t(left+1) <= taui )
      
      left = left + 1
      
      if ( left < ilp1mx ) then
         cycle
      end if
      
      left = left - 1
      
      if ( t(left+1) < taui ) then
         iflag = 2
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'SPLI2D - Fatal error!'
         write ( *, '(a)' ) '  The TAU array is not strictly increasing.'
         stop
      end if
      exit
   end do
   !
   !  The I-th equation enforces interpolation at TAUI, hence
   !
   !    A(I,J) = B(J,K,T)(TAUI), for all J.
   !
   !  Only the K entries with J = LEFT-K+1, ..., LEFT actually might be
   !  nonzero.  These K numbers are returned, in WORK (used for
   !  temporary storage here), by the following call:
   !
   call bsplvb ( t, k, 1, taui, left, work )
   !print*, 'achtung',taui
   ! print*, 'work', work(1:k)
   !
   !  We therefore want
   !
   !    WORK(J) = B(LEFT-K+J)(TAUI)
   !
   !  to go into
   !
   !    A(I,LEFT-K+J),
   !
   !  that is, into  Q(I-(LEFT+J)+2*K,(LEFT+J)-K) since
   !  A(I+J,J) is to go into Q(I+K,J), for all I, J, if we consider Q
   !  as a two-dimensional array, with  2*K-1 rows.  See comments in
   !  BANFAC.
   !
   !  In the present program, we treat Q as an equivalent one-dimensional
   !  array, because of fortran restrictions on dimension statements.
   !
   !  We therefore want WORK(J) to go into the entry of Q with index:
   !    I -(LEFT+J)+2*K + ((LEFT+J)-K-1)*(2*K-1)
   !    = I-LEFT+1+(LEFT -K)*(2*K-1) + (2*K-2)*J
   !
   jj = i - left + 1 + ( left - k ) * ( k + k - 1 )
   
   do j = 1, k
      jj = jj + k + k - 2
      q(jj) = work(j)
   end do
   
end do
!
!  Factor A, stored again in Q.
!
call banfac ( q, k+k-1, n, k-1, k-1, iflag )

if ( iflag == 2 ) then
   write ( *, '(a)' ) ' '
   write ( *, '(a)' ) 'SPLI2D - Fatal error!'
   write ( *, '(a)' ) '  BANFAC reports that the matrix is singular.'
   stop
end if
!  Solve
!
!    A * BCOEF = GTAU
!
!  by back substitution.

do j = 1, m
   work(1:n) = gtau(1:n,j)
   call banslv ( q, k+k-1, n, k-1, k-1, work )
   bcoef(j,1:n) = work(1:n)
end do

return
end subroutine spli2d


   subroutine spli2d_custom_derder ( &
     ai_nx,&
     ai_nx_der,&
     ai_kx,&
     apr_taux,&
     apr_taux_der,&
     ai_ny,&
     ai_ny_der,&
     ai_ky,&
     apr_tauy,&
     apr_tauy_der,&
     apr_g,&
     apr_g_der1,&
     apr_g_der2,&
     apr_Bcoef,&
     apr_tx,&
     apr_ty )
     implicit none
     ! INPUT
     sll_int32  :: ai_nx, ai_kx, ai_ny, ai_ky
     sll_int32  :: ai_nx_der,ai_ny_der
     sll_real64, dimension(:),pointer :: apr_taux !!ai_nx
     sll_real64, dimension(:),pointer :: apr_tauy !! ai_ny
     sll_int32,  dimension(:),pointer :: apr_taux_der !!ai_nx_der
     sll_int32,  dimension(:),pointer :: apr_tauy_der !!ai_ny_der
     sll_real64, dimension(:,:),pointer :: apr_g    ! ai_nx,ai_ny
     sll_real64, dimension(:,:),pointer :: apr_g_der1 ! ai_nx_der,ai_ny
     sll_real64, dimension(:,:),pointer :: apr_g_der2 !ai_ny_der,ai_nx + ai_nx_der
     ! OUTPUT
     sll_real64, dimension(:,:),pointer::apr_Bcoef!ai_nx + ai_nx_der,ai_ny+ ai_ny_der 
    sll_real64, dimension( : ),pointer:: apr_tx ! ai_nx + ai_kx + ai_nx_der
    sll_real64, dimension( : ),pointer:: apr_ty ! ai_ny + ai_ky + ai_ny_der
    ! LOCAL VARIABLES		
    sll_real64, dimension ( ai_nx + ai_nx_der , ai_ny + ai_ny_der) :: lpr_work1
    sll_real64, dimension ( ai_nx + ai_nx_der ) :: lpr_work2
    sll_real64, dimension ( (ai_nx + ai_nx_der)* (ai_ny+ai_ny_der) ) :: lpr_work3
    sll_real64, dimension ( (ai_nx+ai_nx_der) *( 2*ai_kx-1) ) :: lpr_work31
    sll_real64, dimension (( 2*ai_ky-1) * (ai_ny+ai_ny_der) ) :: lpr_work32
    sll_real64, dimension ( ai_ny +ai_ny_der) :: lpr_work4
    sll_real64, dimension (1:ai_ny,1:ai_nx+ai_nx_der),target:: lpr_work5 !  ai_ny , ai_nx
    sll_real64, dimension (:,:),pointer :: lpr_work5_ptr 
    sll_int32  :: li_iflag
    sll_int32  :: ierr
   
    
    lpr_work1(:,:) = 0.0
    lpr_work5(:,:) = 0.0
    
    ! *** set up knots
    !     interpolate between knots
    
    apr_tx = 0.0_f64
    apr_tx ( 1 : ai_kx ) = apr_taux ( 1 )
    apr_tx ( ai_nx+ ai_nx_der + 1: ai_nx + ai_nx_der + ai_kx ) = apr_taux ( ai_nx )
  
    
    if (ai_nx + ai_nx_der + ai_kx == ai_nx + 2*(ai_kx-1)) then
       apr_tx (ai_kx+1: ai_nx+ ai_nx_der) = apr_taux(2:ai_nx-1)
       
    else
       print*, 'problem with construction of knots' 
    end if
    
     !  *** construct b-coefficients of interpolant
    !
    apr_ty = 0.0_f64
    apr_ty ( 1 : ai_ky ) = apr_tauy ( 1 )
    apr_ty ( ai_ny+ ai_ny_der + 1: ai_ny + ai_ny_der + ai_ky ) = apr_tauy ( ai_ny )
    
    
    if (ai_ny + ai_ny_der + ai_ky == ai_ny + 2*(ai_ky-1)) then
       apr_ty (ai_ky+1: ai_ny+ ai_ny_der) = apr_tauy(2:ai_ny-1)
       
    else
       print*, 'problem with construction of knots' 
    end if
    
    lpr_work5_ptr => lpr_work5
    call spli2d_der ( &
         apr_taux,&
         apr_g,&
         apr_taux_der,&
         apr_g_der1,&
         apr_tx, &
         ai_nx,&
         ai_nx_der,&
         ai_kx, &
         ai_ny, &
         lpr_work2,&
         lpr_work31,&
         lpr_work5_ptr, &
         li_iflag )
    
     apr_bcoef(:,:) =0.0_8
     lpr_work4 = 0.0_8
     lpr_work3 = 0.0_8
     lpr_work32= 0.0_8
     

     call spli2d_der ( &
          apr_tauy,&
          lpr_work5_ptr,&
          apr_tauy_der,&
          apr_g_der2,&
          apr_ty,&
          ai_ny, &
          ai_ny_der,&
          ai_ky, &
          ai_nx+ai_nx_der, &
          lpr_work4, &
          lpr_work32,&
          apr_bcoef, &
          li_iflag )

   end subroutine spli2d_custom_derder

subroutine spli2d_der(&
     tau,&
     gtau,&
     tau_der,&
     gtau_der,&
     t, n,np, k, m, work, q, bcoef, iflag )

    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TAU(N), contains the data point abscissas.
    !    TAU must be strictly increasing
    !
    !    Input, real ( kind = 8 ) GTAU(N,M), contains the data point ordinates.
  
    !    Input, integer (kind= 8 )TAU_DER(Np), contains the data point abscissas.
    !    TAU must be strictly increasing
    !
    !    Input, real ( kind = 8 ) GTAU_DER(Np,M),contains the data point ordinates.
    !
    !    Input, real ( kind = 8 ) T(N+Np+K), the knot sequence.
    !
    !    Input, integer N, the number of data points and the
    !    dimension of the spline space SPLINE(K,T)
    !
    !    Input, integer K, the order of the spline.
    !
    !    Input, integer M, the number of data sets.
    !
    !    Work space, real ( kind = 8 ) WORK(N).
    !
    !    Output, real ( kind = 8 ) Q(2*K-1)*N, the triangular
    !    factorization of the coefficient matrix of the linear
    !    system for the B-spline coefficients of the spline interpolant.
    !    The B-spline coefficients for the interpolant of an additional
    !    data set ( TAU(I), HTAU(I) ), I=1,...,N  with the same data
    !    abscissae can be obtained without going through all the
    !    calculations in this routine, simply by loading HTAU into
    !    BCOEF and then using the statement
    !      CALL BANSLV ( Q, 2*K-1, N, K-1, K-1, BCOEF )
    !
    !    Output, real ( kind = 8 ) BCOEF(N), the B-spline coefficients of
    !    the interpolant.
    !
    !    Output, integer IFLAG, error indicator.
    !    1, no error.
    !    2, an error occurred, which may have been caused by
    !       singularity of the linear system.
    !
    implicit none
    
    sll_int32 :: m
    sll_int32 :: n
    sll_int32 :: np
    sll_real64,dimension(:,:),pointer:: bcoef !(m,n+np)
    sll_real64,dimension(:,:),pointer:: gtau  !(n,m)
    sll_real64,dimension(:,:),pointer:: gtau_der!(np,n)
    sll_int32 :: iflag
    sll_int32 :: j
    sll_int32 :: k
    sll_real64,dimension((2*k-1)*n):: q!((2*k-1)*n)
    sll_real64,dimension(:),pointer:: t!(n+np+k)
    sll_real64,dimension(:),pointer:: tau!(n)
    sll_int32,dimension(:),pointer:: tau_der!np
    sll_real64,dimension(n):: work!(n)
    sll_real64,dimension(np):: work_der
    sll_real64,dimension(n+np):: work_result

    
    work_result = 0.0_f64

    !print*, 'hello',gtau_der(1:np,1:m)

    do j = 1, m
       
       work(1:n) = gtau(:,j)
       work_der(1:np) = gtau_der(1:np,j)
       call splint_der(&
            tau,&
            work,&
            tau_der,&
            work_der,t,n,np,k, q, work_result, iflag )
       
       
       bcoef(j,1:n+np) = work_result(1:n+np)
       
    end do
    
    
    return
  end subroutine spli2d_der

!*************************************************************************
!
!! SPLINT_der produces the B-spline coefficients BCOEF of an 
! interpolating spline with the values of a derivative in points
!
!  Discussion:
!
!    The spline is of order K with knots T(1:N+K+M), and takes on the 
!    value GTAU(I) at TAU(I), for I = 1 to N and 
!    value of the derivative GTAU_der(I) at TAU_der(I), for I = 1 to M
!
!    The I-th equation of the linear system 
!
!      A  * BCOEF = B
!      A' * BCOEF = B'
!
!    for the B-spline coefficients of the interpolant enforces interpolation
!    at TAU(1:N) and the derivative at TAU_der(I).
!
!    Hence, B(I) = GTAU(I) and B'(I) = GTAU_der(I) , for all I,
!    and A is a band matrix with 2*K-1
!    bands, if it is invertible.
!
!    The matrix A is generated row by row and stored, diagonal by diagonal,
!    in the rows of the array Q, with the main diagonal going
!    into row K.  See comments in the program.
!
!    The banded system is then solved by a call to BANFAC, which 
!    constructs the triangular factorization for A and stores it again in
!    Q, followed by a call to BANSLV, which then obtains the solution
!    BCOEF by substitution.
!
!    BANFAC does no pivoting, since the total positivity of the matrix
!    A makes this unnecessary.
!
!    The linear system to be solved is (theoretically) invertible if
!    and only if
!      T(I) < TAU(I) < TAU(I+K), for all I.
!    Violation of this condition is certain to lead to IFLAG = 2.
!
!  Modified:
!
!    10 April 2014
!
!  Author:
!
!    Aurore Back
!
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TAU(N), the data point abscissas.The entries in
!    TAU should be strictly increasing.
!
!    Input, integer ( kind = 8 ) TAU_der(M), the node index to evaluate the derivative.
!
!    Input, real ( kind = 8 ) GTAU(N), the data ordinates.
!
!    Input, real ( kind = 8 ) GTAU_der(M), the data ordinates.
!
!    Input, real ( kind = 8 ) T(N+K+M), the knot sequence.
!
!    Input, integer ( kind = 4 ) N, the number of data points for the interpolation.
!
!    Input, integer ( kind = 4 ) M, the number of data points for the derivative.
!
!    Input, integer ( kind = 4 ) K, the order of the spline.
!
!    Output, real ( kind = 8 ) Q((2*K-1)*(N+M)), the triangular factorization
!    of the coefficient matrix of the linear system for the B-coefficients 
!    of the spline interpolant.  The B-coefficients for the interpolant 
!    of an additional data set can be obtained without going through all 
!    the calculations in this routine, simply by loading HTAU into BCOEF 
!    and then executing the call:
!      call banslv ( q, 2*k-1, n+m, k-1, k-1, bcoef )
!
!    Output, real ( kind = 8 ) BCOEF(N+M), the B-spline coefficients of 
!    the interpolant.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    1, = success.
!    2, = failure.
!
subroutine splint_der( tau,gtau,tau_der,gtau_der,t,n,m,k, q, bcoef_spline, iflag )
    
sll_int32 ::  n
sll_int32 ::  m
sll_real64,dimension(n+m):: bcoef ! (n)
sll_real64,dimension(n+m):: bcoef_spline 
sll_real64,dimension(n)::  gtau ! (n)
sll_real64,dimension(m)::  gtau_der ! (n)
sll_int32 :: i
sll_int32 :: iflag,mflag
sll_int32 :: j,l
sll_int32 :: jj
sll_int32 :: k
sll_int32 :: kpkm2
sll_int32 :: left
sll_real64, dimension((2*k-1)*(n+m)) :: q!((2*k-1)*n)
sll_real64,dimension(n+k+m) ::  t!(n+k)
sll_real64,dimension(n) ::  tau!!(n)
sll_int32,dimension(m) ::  tau_der!!(n)
sll_real64:: taui,taui_der
sll_real64, dimension(k,k):: a
sll_real64,dimension(k,2) :: bcoef_der

kpkm2 = 2 * ( k - 1 )
left = k
q(1:(2*k-1)*(n+m)) = 0.0_f64
a(1:k,1:k) = 0.0_f64
bcoef_der(1:k,1:2) = 0.0_f64

! we must suppose that m is <= than n 
if (m > n) then
   print*, 'problem m must be < = at n'
   print*, 'value m =', m, 'value n =', n
   stop
end if
l = 1 ! index for the derivative
!
!  Loop over I to construct the N interpolation equations.
!
do i = 1, n-1
   
   taui = tau(i)
   
   !
   !  Find LEFT in the closed interval (I,I+K-1) such that
   !
   !    T(LEFT) <= TAU(I) < T(LEFT+1)
   !
   !  The matrix is singular if this is not possible.
   !  With help of the Schoenberg-Whitney theorem 
   !  we can prove that if the diagonal of the 
   !  matrix B_j(x_i) is null, we have a non-inversible matrix.  

   call interv( t, n+m+k, taui, left, mflag )

   !

   !
   !  The I-th equation enforces interpolation at TAUI, hence for all J,
   !
   !    A(I,J) = B(J,K,T)(TAUI).
   !
   !Only the K entries with J = LEFT-K+1,...,LEFT actually might be nonzero.
   !
   !These K numbers are returned, in BCOEF 
   ! (used for temporary storage here),
   !  by the following.
   !
   
   call bsplvb ( t, k, 1, taui, left, bcoef )
   
      
   !
   !  We therefore want BCOEF(J) = B(LEFT-K+J)(TAUI) to go into
   !  A(I,LEFT-K+J), that is, into Q(I-(LEFT+J)+2*K,(LEFT+J)-K) since
   !  A(I+J,J) is to go into Q(I+K,J), for all I, J, if we consider Q
   !  as a two-dimensional array, with  2*K-1 rows.  See comments in
   !  BANFAC.
   !
   !  In the present program, we treat Q as an equivalent
   !  one-dimensional array, because of fortran restrictions on
   !  dimension statements.
   !
   !  We therefore want  BCOEF(J) to go into the entry of Q with index:
   !
   !    I -(LEFT+J)+2*K + ((LEFT+J)-K-1)*(2*K-1)
   !    =  begin_ligne +  (begin_col -1) * number_coef_different_0
   !   = I-LEFT+1+(LEFT -K)*(2*K-1) + (2*K-2)*J
   !
   jj = i - left + 1 + ( left - k ) * ( k + k - 1 ) + l - 1
   
   do j = 1, k
      jj = jj + kpkm2  ! kpkm2 = 2*(k-1)
      q(jj) = bcoef(j)
   end do

   bcoef_spline(i+ l-1) = gtau(i)
   if ( tau_der(l) == i ) then   
      taui_der = taui
      
      call bsplvd( t, k, taui_der, left, a, bcoef_der, 2)

      l = l + 1
      jj = i - left + 1 + ( left - k ) * ( k + k - 1 ) + l - 1
   
      do j = 1, k
         jj = jj + kpkm2  ! kpkm2 = 2*(k-1)
         q(jj) = bcoef_der(j,2)
      end do
   bcoef_spline(i+ l-1) = gtau_der(l-1)
   end if
   

end do


taui = tau(n)
call interv( t, n+m+k, taui, left, mflag )
if ( tau_der(l)== n ) then   
      taui_der = taui
      
      call bsplvd( t, k, taui_der, left, a, bcoef_der, 2)

      
      jj = n - left + 1 + ( left - k ) * ( k + k - 1 ) + l - 1
   
      do j = 1, k
         jj = jj + kpkm2  ! kpkm2 = 2*(k-1)
         q(jj) = bcoef_der(j,2)
      end do
      bcoef_spline(n+ l-1) = gtau_der(l)
      l = l + 1
      
   end if
   

 call bsplvb ( t, k, 1, taui, left, bcoef )
 jj = n - left + 1 + ( left - k ) * ( k + k - 1 ) + l - 1
   
 do j = 1, k
    jj = jj + kpkm2  ! kpkm2 = 2*(k-1)
    q(jj) = bcoef(j)
 end do
 bcoef_spline(n+l-1) = gtau(n)

!
!  Obtain factorization of A, stored again in Q.
!
call banfac ( q, k+k-1, n+m, k-1, k-1, iflag )

if ( iflag == 2 ) then
   write ( *, '(a)' ) ' '
   write ( *, '(a)' ) 'SPLINT - Fatal Error!'
   write ( *, '(a)' ) '  The linear system is not invertible!'
   return
end if
!
!  Solve 
!
!    A * BCOEF = GTAU
!
!  by back substitution.
!
call banslv ( q, k+k-1, n+m, k-1, k-1, bcoef_spline )

end subroutine splint_der
  
end module sll_module_arbitrary_degree_spline_interpolator_2d
