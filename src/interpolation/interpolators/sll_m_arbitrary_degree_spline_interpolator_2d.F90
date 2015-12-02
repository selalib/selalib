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
module sll_m_arbitrary_degree_spline_interpolator_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h" 
use sll_m_interpolators_2d_base
use sll_m_utilities
use sll_m_deboor_splines_1d, only: &
  deboor_type, bsplvb, bsplvd, interv, bvalue, splint_der
use sll_m_arbitrary_degree_spline_interpolator_1d

implicit none
private

! in what follows, the direction '1' is in the contiguous memory direction.
!> Arbitrary degree version of 2d irnterpolator
type, extends(sll_c_interpolator_2d) :: &
  sll_arbitrary_degree_spline_interpolator_2d           

  sll_int32                           :: num_pts1
  sll_int32                           :: num_pts2
  sll_real64                          :: eta1_min
  sll_real64                          :: eta1_max
  sll_real64                          :: eta2_min
  sll_real64                          :: eta2_max
  sll_int32                           :: bc1_min
  sll_int32                           :: bc1_max
  sll_int32                           :: bc2_min
  sll_int32                           :: bc2_max
  sll_int32                           :: spline_degree1
  sll_int32                           :: spline_degree2
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

  sll_real64, dimension(:,:), pointer :: gtau
  sll_real64, dimension(:,:), pointer :: gtau_der1
  sll_real64, dimension(:,:), pointer :: gtau_der2

  type(deboor_type)                   :: deboor(2)  !< Deboor splines data object

contains

  procedure :: initialize                  => initialize_ad2d_interpolator
  procedure :: set_coefficients            => set_coefficients_ad2d
  procedure :: coefficients_are_set        => coefficients_are_set_ad2d
  procedure :: compute_interpolants        => compute_interpolants_ad2d
  procedure :: interpolate_from_interpolant_value           => interpolate_value_ad2d
  procedure :: interpolate_from_interpolant_derivative_eta1 => interpolate_derivative1_ad2d
  procedure :: interpolate_from_interpolant_derivative_eta2 => interpolate_derivative2_ad2d
  procedure :: interpolate_array           => interpolate_array_ad2d
  procedure :: interpolate_array_disp      => interpolate_2d_array_disp_ad2d
  procedure :: get_coefficients            => get_coefficients_ad2d
  procedure :: delete                      => delete_arbitrary_degree_2d_interpolator
  procedure :: set_values_at_boundary      => set_boundary_value2d
  procedure :: set_slopes_at_boundary      => set_slope2d

end type sll_arbitrary_degree_spline_interpolator_2d


!> Pointer to arbitrary degree version of 2d interpolator
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
public set_coeff_splines_values_1d

contains

!> Delete interpolator arbitrary degree splines.
!> The parameters are
!> @param interpolator the type sll_arbitrary_degree_spline_interpolator_2d
!
subroutine delete_arbitrary_degree_2d_interpolator( interpolator )

class(sll_arbitrary_degree_spline_interpolator_2d), intent(inout) :: interpolator
sll_int32 :: ierr

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
!> @param[in] bc1_min  the boundary condition at left in the direction eta1
!> @param[in] bc1_max  the boundary condition at right in the direction eta2
!> @param[in] bc2_min  the boundary condition at left in the direction eta2
!> @param[in] bc2_max  the boundary condition at right in the direction eta2
!> @param[in] spline_degree1 the degree of B-spline in the direction eta1
!> @param[in] spline_degree2 the degre of B-spline in the direction eta2
!> @return the type sll_arbitrary_degree_spline_interpolator_2d

function new_arbitrary_degree_spline_interp2d( num_pts1,       &
                                               num_pts2,       &
                                               eta1_min,       &
                                               eta1_max,       &
                                               eta2_min,       &
                                               eta2_max,       &
                                               bc1_min,        &
                                               bc1_max,        &
                                               bc2_min,        &
                                               bc2_max,        &
                                               spline_degree1, &
                                               spline_degree2) result( res )

type(sll_arbitrary_degree_spline_interpolator_2d), pointer :: res

sll_int32,  intent(in) :: num_pts1
sll_int32,  intent(in) :: num_pts2
sll_real64, intent(in) :: eta1_min
sll_real64, intent(in) :: eta1_max
sll_real64, intent(in) :: eta2_min
sll_real64, intent(in) :: eta2_max
sll_int32,  intent(in) :: bc1_min
sll_int32,  intent(in) :: bc1_max
sll_int32,  intent(in) :: bc2_min
sll_int32,  intent(in) :: bc2_max
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
                                   bc1_min,        &
                                   bc1_max,        &
                                   bc2_min,        &
                                   bc2_max,        &
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
!> @param[in] bc1_min  the boundary condition at left in the direction eta1
!> @param[in] bc1_max the boundary condition at right in the direction eta2
!> @param[in] bc2_min the boundary condition at left in the direction eta2
!> @param[in] bc2_max the boundary condition at right in the direction eta2
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
                                         bc1_min,        &
                                         bc1_max,        &
                                         bc2_min,        &
                                         bc2_max,        &
                                         spline_degree1, &
                                         spline_degree2)

class(sll_arbitrary_degree_spline_interpolator_2d):: interpolator

sll_int32,  intent(in)  :: num_pts1
sll_int32,  intent(in)  :: num_pts2
sll_real64, intent(in)  :: eta1_min
sll_real64, intent(in)  :: eta1_max
sll_real64, intent(in)  :: eta2_min
sll_real64, intent(in)  :: eta2_max
sll_int32,  intent(in)  :: bc1_min
sll_int32,  intent(in)  :: bc1_max
sll_int32,  intent(in)  :: bc2_min
sll_int32,  intent(in)  :: bc2_max
sll_int32,  intent(in)  :: spline_degree1
sll_int32,  intent(in)  :: spline_degree2

sll_int32 :: ierr
sll_int32 :: tmp1
sll_int32 :: tmp2
sll_int64 :: bc_selector


! do some argument checking...
if( bc1_min == SLL_PERIODIC .and. bc1_max .ne. SLL_PERIODIC .or.&
    bc1_max == SLL_PERIODIC .and. bc1_min .ne. SLL_PERIODIC )then
   print *, 'initialize_arbitrary_degree_2d_interpolator, ERROR: ', &
        'if one boundary condition is specified as periodic, then ', &
        'both must be. Error in first direction.'
end if

if(((bc2_min == SLL_PERIODIC).and.(bc2_max .ne. SLL_PERIODIC)).or.&
   ((bc2_max == SLL_PERIODIC).and.(bc2_min .ne. SLL_PERIODIC)))then
   print *, 'initialize_arbitrary_degree_2d_interpolator, ERROR: ', &
        'if one boundary condition is specified as periodic, then ', &
        'both must be. Error in second direction.'
end if

bc_selector = 0

if( bc1_min == SLL_DIRICHLET ) bc_selector = bc_selector + 1
if( bc1_min == SLL_NEUMANN   ) bc_selector = bc_selector + 2
if( bc1_min == SLL_HERMITE   ) bc_selector = bc_selector + 4
if( bc1_max == SLL_DIRICHLET ) bc_selector = bc_selector + 8
if( bc1_max == SLL_NEUMANN   ) bc_selector = bc_selector + 16
if( bc1_max == SLL_HERMITE   ) bc_selector = bc_selector + 32
if( bc2_min == SLL_DIRICHLET ) bc_selector = bc_selector + 64
if( bc2_min == SLL_NEUMANN   ) bc_selector = bc_selector + 128
if( bc2_min == SLL_HERMITE   ) bc_selector = bc_selector + 256
if( bc2_max == SLL_DIRICHLET ) bc_selector = bc_selector + 512
if( bc2_max == SLL_NEUMANN   ) bc_selector = bc_selector + 1024
if( bc2_max == SLL_HERMITE   ) bc_selector = bc_selector + 2048

interpolator%spline_degree1 = spline_degree1
interpolator%spline_degree2 = spline_degree2
interpolator%eta1_min       = eta1_min
interpolator%eta1_max       = eta1_max
interpolator%eta2_min       = eta2_min
interpolator%eta2_max       = eta2_max
interpolator%bc1_min        = bc1_min
interpolator%bc1_max        = bc1_max
interpolator%bc2_min        = bc2_min
interpolator%bc2_max        = bc2_max
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

tmp1 = num_pts1+ 4*spline_degree1
tmp2 = num_pts2+ 4*spline_degree2
SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

interpolator%coeff_splines(:,:) = 0.0_f64
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
deallocate(interp1d)
nullify(interp1d)

end subroutine set_coeff_splines_values_1d


!> Initialization of the boundary for interpolator arbitrary degree splines 2d.
!> The parameters are
!> @param[in]  value_min1 a 1d array contains values in the left  in the direction 1  
!> @param[in]  value_max1 a 1d array contains values in the right in the direction 1 
!> @param[in]  value_min2 a 1d array contains values in the left  in the direction 2 
!> @param[in]  value_max2 a 1d array contains values in the right in the direction 2
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
sll_int32  :: bc1_min
sll_int32  :: bc1_max
sll_int32  :: bc2_min
sll_int32  :: bc2_max
sll_int32  :: spline_degree1
sll_int32  :: spline_degree2
sll_real64 :: eta1_min
sll_real64 :: eta1_max
sll_real64 :: eta2_min
sll_real64 :: eta2_max

num_pts1       = interpolator%num_pts1
num_pts2       = interpolator%num_pts2
bc1_min        = interpolator%bc1_min 
bc1_max        = interpolator%bc1_max 
bc2_min        = interpolator%bc2_min  
bc2_max        = interpolator%bc2_max
spline_degree1 = interpolator%spline_degree1
spline_degree2 = interpolator%spline_degree2
eta1_min       = interpolator%eta1_min
eta1_max       = interpolator%eta1_max
eta2_min       = interpolator%eta2_min
eta2_max       = interpolator%eta2_max

if (bc1_min==SLL_DIRICHLET) then 

  if (present(value_min1)) then 
  
    call set_coeff_splines_values_1d( value_min1, &
         num_pts2,                                &
         eta2_min,                                &
         eta2_max,                                &
         bc2_min,                                 &
         bc2_max,                                 &
         spline_degree2 )
      
    interpolator%value_min1(1:num_pts2) = value_min1
    interpolator%compute_value_min1 = .FALSE.
  else
    interpolator%value_min1 = 0.0_f64
  end if

end if
  
if (bc1_max==SLL_DIRICHLET) then 

  if (present(value_max1)) then 
    call set_coeff_splines_values_1d( value_max1, &
         num_pts2,                                &
         eta2_min,                                &
         eta2_max,                                &
         bc2_min,                                 &
         bc2_max,                                 &
         spline_degree2 )
      
    interpolator%value_max1(1:num_pts2) = value_max1
    interpolator%compute_value_max1 = .FALSE.
  else
    interpolator%value_max1 = 0.0_f64
  end if

end if

if (bc2_min==SLL_DIRICHLET) then 

  if (present(value_min2)) then 
    call set_coeff_splines_values_1d( value_min2, &
                                      num_pts1,   &
                                      eta1_min,   &
                                      eta1_max,   &
                                      bc1_min,    &
                                      bc1_max,    &
                                      spline_degree1 )
      
    interpolator%value_min2(1:num_pts1) = value_min2
    interpolator%compute_value_min2 = .FALSE.
  else
    interpolator%value_min2 = 0.0_f64
  end if

end if
  
if (bc2_max==SLL_DIRICHLET) then 

  if (present(value_max2)) then 
    call set_coeff_splines_values_1d( value_max2,    &
                                      num_pts1,      &
                                      eta1_min,      &
                                      eta1_max,      &
                                      bc1_min,       &
                                      bc1_max,       &
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


if (present(coeffs_1d)) then  !This case is used to set the solution from
                              !the general coordinate elliptic solver

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
    
    SLL_ASSERT(size(coeffs_1d,1)==(num_cells1+sp_deg1-2)*(num_cells2+sp_deg2-2))
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
    interpolator%size_t1      = 2*sp_deg1  + num_cells1 + 1
    interpolator%size_t2      = 2*sp_deg2  + num_cells2 + 1
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
    interpolator%size_t1      = 2*sp_deg1  + num_cells1+1
    interpolator%size_t2      = 2*sp_deg2  + num_cells2+1
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
    
  case(2340) ! Hermite in all sides
       
    interpolator%size_coeffs1 = num_cells1 + sp_deg1
    interpolator%size_coeffs2 = num_cells2 + sp_deg2
    interpolator%size_t1      = 2*sp_deg1  + num_cells1+1
    interpolator%size_t2      = 2*sp_deg2  + num_cells2+1
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
 
end subroutine set_coefficients_ad2d

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

sll_int32 :: taux_deriv(2)
sll_int32 :: tauy_deriv(2)

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

  
select case (interpolator%bc_selector)

case (0) ! periodic-periodic

  interpolator%size_coeffs1 = nx
  interpolator%size_coeffs2 = ny
  interpolator%size_t1      = kx + nx
  interpolator%size_t2      = ky + ny

  !Apply periodic boundary conditions
  taux(nx) = taux(1)+period1
  tauy(ny) = tauy(1)+period2

  SLL_ALLOCATE(interpolator%gtau(1:nx,1:ny), ierr)
  interpolator%gtau(1:nx-1,1:ny-1) = data_array(1:nx-1,1:ny-1)
  interpolator%gtau(nx,1:ny-1)     = data_array(1,1:ny-1 )
  interpolator%gtau(1:nx-1,ny)     = data_array(1:nx-1,1)
  interpolator%gtau(nx,ny)         = data_array(1,1)

  call spli2d_custom(interpolator%deboor,        &
                     nx, kx, taux, ny, ky, tauy, &
                     interpolator%gtau,          &
                     interpolator%coeff_splines, &
                     interpolator%t1,            &
                     interpolator%t2)
  deallocate(interpolator%gtau)
   
case (9) ! 2. dirichlet-left, dirichlet-right, periodic

  interpolator%size_coeffs1 = nx
  interpolator%size_coeffs2 = ny
  interpolator%size_t1      = kx+nx
  interpolator%size_t2      = ky+ny
       
  tauy(1:ny-1)      = tauy(1:ny-1)
  tauy(ny)          = tauy(1)+period2

  SLL_ALLOCATE(interpolator%gtau(1:nx,1:ny), ierr)

  interpolator%gtau(1:nx-1,1:ny) = data_array(1:nx-1,1:ny)
  interpolator%gtau(nx,1:ny)     = data_array(1,1:ny)

  call spli2d_custom(interpolator%deboor,        &
                     nx, kx, taux, ny, ky, tauy, &
                     interpolator%gtau,          &
                     interpolator%coeff_splines, &
                     interpolator%t1,            &
                     interpolator%t2)

  interpolator%coeff_splines(1,1:ny)  = data_array(1,1:ny)
  interpolator%coeff_splines(nx,1:ny) = data_array(nx,1:ny)
  deallocate(interpolator%gtau)
  
case(576) !  3. periodic, dirichlet-bottom, dirichlet-top

  interpolator%size_coeffs1 = nx!+1
  interpolator%size_coeffs2 = ny
  interpolator%size_t1 = kx + nx !+ 1
  interpolator%size_t2 = ky + ny 

  taux(1:nx-1) = taux(1:nx-1)
  taux(nx)     = taux(1)+period1

  SLL_ALLOCATE(interpolator%gtau(1:nx,1:ny), ierr)

  interpolator%gtau(1:nx,1:ny-1) = data_array(1:nx,1:ny-1)
  interpolator%gtau(1:nx,ny)     = data_array(1:nx,1)

  call spli2d_custom(interpolator%deboor,        &
                     nx, kx, taux, ny, ky, tauy, &
                     interpolator%gtau,          &
                     interpolator%coeff_splines, &
                     interpolator%t1,            &
                     interpolator%t2)

  interpolator%coeff_splines(1:nx,1)   = data_array(1:nx,1)
  interpolator%coeff_splines(1:nx,ny) = data_array(1:nx,ny)
  deallocate(interpolator%gtau)
       
case (585) ! 4. dirichlet in all sides

  interpolator%size_coeffs1 = nx
  interpolator%size_coeffs2 = ny
  interpolator%size_t1 = kx + nx 
  interpolator%size_t2 = ky + ny 
  
  SLL_ALLOCATE(interpolator%gtau(1:nx,1:ny),ierr)
  interpolator%gtau = data_array(1:nx,1:ny)
  call spli2d_custom( interpolator%deboor,        &
                      nx, kx, taux, ny, ky, tauy, &
                      interpolator%gtau,          &
                      interpolator%coeff_splines, &
                      interpolator%t1,            &
                      interpolator%t2)

  interpolator%coeff_splines(1,1:ny)  = data_array(1,1:ny)
  interpolator%coeff_splines(nx,1:ny) = data_array(nx,1:ny)
  interpolator%coeff_splines(1:nx,1)  = data_array(1:nx,1)
  interpolator%coeff_splines(1:nx,ny) = data_array(1:nx,ny)
  deallocate(interpolator%gtau)

case (650) !left: Neumann, right: Dirichlet, bottom: Neumann, Top: Dirichlet

  mx = 2
  my = 2
  interpolator%size_coeffs1 = nx + mx
  interpolator%size_coeffs2 = ny + my
  interpolator%size_t1 = kx + nx + mx
  interpolator%size_t2 = ky + ny + my
  
  SLL_ALLOCATE(interpolator%gtau(1:nx,1:ny),ierr)
  SLL_CLEAR_ALLOCATE(interpolator%gtau_der1(1:2,1:ny),ierr)
  SLL_CLEAR_ALLOCATE(interpolator%gtau_der2(1:my,1:nx+mx),ierr)

  interpolator%gtau                 = data_array(1:nx,1:ny)
  taux_deriv(1)        = 1
  taux_deriv(2)        = nx
  interpolator%gtau_der1(1,1:ny)    = 0.0_f64
  interpolator%gtau_der1(2,1:ny)    = interpolator%slope_max1(1:ny)
  tauy_deriv(1)        = 1
  tauy_deriv(2)        = ny
  interpolator%gtau_der2(1,1:nx+mx) = 0.0_f64
  interpolator%gtau_der2(2,1:nx+mx) = interpolator%slope_max2(1:nx+mx)

  call spli2d_custom_derder(interpolator%deboor,        &
                            nx,                         &
                            mx,                         &
                            kx,                         &
                            taux,                       &
                            taux_deriv,                 &
                            ny,                         &
                            my,                         &
                            ky,                         &
                            tauy,                       &
                            tauy_deriv,                 &
                            interpolator%gtau,          &
                            interpolator%gtau_der1,     &
                            interpolator%gtau_der2,     &
                            interpolator%coeff_splines, &
                            interpolator%t1,            &
                            interpolator%t2)

  SLL_DEALLOCATE( interpolator%gtau_der1,ierr)
  SLL_DEALLOCATE( interpolator%gtau_der2,ierr)
  deallocate(interpolator%gtau)

case(657) !left: Dirichlet, right: Neumann, bottom: Neumann, Top: Dirichlet 

  mx = 2
  my = 2
  interpolator%size_coeffs1 = nx + mx
  interpolator%size_coeffs2 = ny + my
  interpolator%size_t1 = kx + nx + mx
  interpolator%size_t2 = ky + ny + my
  
  SLL_ALLOCATE( interpolator%gtau(1:nx,1:ny),ierr)
  SLL_ALLOCATE( interpolator%gtau_der1(mx,1:ny),ierr)
  SLL_ALLOCATE( interpolator%gtau_der2(my,1:nx+mx),ierr)

  interpolator%gtau                 = data_array(1:nx,1:ny)
  taux_deriv(1)        = 1
  taux_deriv(2)        = nx
  interpolator%gtau_der1(1,1:ny)    = interpolator%slope_min1(1:ny) 
  interpolator%gtau_der1(2,1:ny)    = 0.0_f64
  tauy_deriv(1)        = 1
  tauy_deriv(2)        = ny
  interpolator%gtau_der2(1,1:nx+mx) = 0.0_f64
  interpolator%gtau_der2(2,1:nx+mx) = interpolator%slope_max2(1:nx+mx)

  call spli2d_custom_derder(interpolator%deboor,        &
                            nx,                         &
                            mx,                         &
                            kx,                         &
                            taux,                       &
                            taux_deriv,                 &
                            ny,                         &
                            my,                         &
                            ky, tauy,                   &
                            tauy_deriv,                 &
                            interpolator%gtau,          &
                            interpolator%gtau_der1,     &
                            interpolator%gtau_der2,     &
                            interpolator%coeff_splines, &
                            interpolator%t1,            &
                            interpolator%t2)

  SLL_DEALLOCATE( interpolator%gtau_der1,ierr)
  SLL_DEALLOCATE( interpolator%gtau_der2,ierr)
  deallocate(interpolator%gtau)

case(780)  !left: Hermite, right: Dirichlet, bottom: Hermite, Top: Dirichlet

  mx = 2
  my = 2
  interpolator%size_coeffs1 = nx + mx
  interpolator%size_coeffs2 = ny + my
  interpolator%size_t1 = kx + nx + mx
  interpolator%size_t2 = ky + ny + my
  
  SLL_ALLOCATE(interpolator%gtau(1:nx,1:ny),ierr)
  SLL_CLEAR_ALLOCATE(interpolator%gtau_der1(1:mx,1:ny),ierr)
  SLL_CLEAR_ALLOCATE(interpolator%gtau_der2(1:my,1:nx+mx),ierr)

  interpolator%gtau                  = data_array(1:nx,1:ny)
  taux_deriv(1)         = 1
  taux_deriv(2)         = nx
  interpolator%gtau_der1(1,1:ny)     = interpolator%slope_min1(1:ny)
  interpolator%gtau_der1(2,1:ny)     = interpolator%slope_max1(1:ny)
  tauy_deriv(1)         = 1
  tauy_deriv(2)         = ny 
  interpolator%gtau_der2(1,1:nx+mx)  = interpolator%slope_min2(1:nx+mx)
  interpolator%gtau_der2(2,1:nx+mx)  = interpolator%slope_max2(1:nx+mx)

  call spli2d_custom_derder(interpolator%deboor,        &
                            nx,                         &
                            mx,                         &
                            kx,                         &
                            taux,                       &
                            taux_deriv,                 &
                            ny,                         &
                            my,                         & 
                            ky, tauy,                   &
                            tauy_deriv,                 &
                            interpolator%gtau,          &
                            interpolator%gtau_der1,     &
                            interpolator%gtau_der2,     &
                            interpolator%coeff_splines, &
                            interpolator%t1,            &
                            interpolator%t2)


  SLL_DEALLOCATE( interpolator%gtau_der1,ierr)
  SLL_DEALLOCATE( interpolator%gtau_der2,ierr)
  deallocate(interpolator%gtau)

case(801)  !left: Dirichlet, right: Hermite, bottom: Hermite, Top: Dirichlet


  mx = 2
  my = 2
  interpolator%size_coeffs1 = nx + mx
  interpolator%size_coeffs2 = ny + my
  interpolator%size_t1 = kx + nx + mx
  interpolator%size_t2 = ky + ny + my
  
  SLL_ALLOCATE( interpolator%gtau(1:nx,1:ny),ierr)
  SLL_ALLOCATE( interpolator%gtau_der1(mx,1:ny),ierr)
  SLL_ALLOCATE( interpolator%gtau_der2(my,1:nx+mx),ierr)

  interpolator%gtau                 = data_array(1:nx,1:ny)
  taux_deriv(1)                     = 1
  taux_deriv(2)                     = nx
  interpolator%gtau_der1(1,1:ny)    = interpolator%slope_min1(1:ny) 
  interpolator%gtau_der1(2,1:ny)    = interpolator%slope_max1(1:ny) 
  tauy_deriv(1)                     = 1
  tauy_deriv(2)                     = ny
  interpolator%gtau_der2(1,1:nx+mx) = interpolator%slope_min2(1:nx+mx)
  interpolator%gtau_der2(2,1:nx+mx) = interpolator%slope_max2(1:nx+mx)

  call spli2d_custom_derder(interpolator%deboor,        &
                            nx,                         &
                            mx,                         &
                            kx,                         &
                            taux,                       &
                            taux_deriv,                 &
                            ny,                         &
                            my,                         &
                            ky, tauy,                   &
                            tauy_deriv,                 &
                            interpolator%gtau,          &
                            interpolator%gtau_der1,     &
                            interpolator%gtau_der2,     &
                            interpolator%coeff_splines, &
                            interpolator%t1,            &
                            interpolator%t2)

  SLL_DEALLOCATE(interpolator%gtau_der1,ierr)
  SLL_DEALLOCATE(interpolator%gtau_der2,ierr)
  deallocate(interpolator%gtau)

case(804)  !left: Hermite, right: Hermite, bottom: Hermite, Top: Dirichlet

  mx = 2
  my = 2
  interpolator%size_coeffs1 = nx + mx
  interpolator%size_coeffs2 = ny + my
  interpolator%size_t1 = kx + nx + mx
  interpolator%size_t2 = ky + ny + my
  
  SLL_ALLOCATE(interpolator%gtau(1:nx,1:ny),ierr)
  SLL_ALLOCATE(interpolator%gtau_der1(mx,1:ny),ierr)
  SLL_ALLOCATE(interpolator%gtau_der2(my,1:nx+mx),ierr)

  interpolator%gtau                 = data_array(1:nx,1:ny)
  taux_deriv(1)        = 1
  taux_deriv(2)        = nx
  interpolator%gtau_der1(1,1:ny)    = interpolator%slope_min1(1:ny) 
  interpolator%gtau_der1(2,1:ny)    = interpolator%slope_max1(1:ny) 
  tauy_deriv(1)        = 1
  tauy_deriv(2)        = ny
  interpolator%gtau_der2(1,1:nx+mx) = interpolator%slope_min2(1:nx+mx)
  interpolator%gtau_der2(2,1:nx+mx) = interpolator%slope_max2(1:nx+mx)

  call spli2d_custom_derder( interpolator%deboor,       &
                             nx,                        &
                             mx,                        &
                             kx,                        &
                             taux,                      &
                             taux_deriv,                &
                             ny,                        &
                             my,                        &
                             ky, tauy,                  &
                             tauy_deriv,                &
                             interpolator%gtau,         &
                             interpolator%gtau_der1,    &
                             interpolator%gtau_der2,    &
                             interpolator%coeff_splines,&
                             interpolator%t1,           &
                             interpolator%t2)

  SLL_DEALLOCATE( interpolator%gtau_der1,ierr)
  SLL_DEALLOCATE( interpolator%gtau_der2,ierr)
  deallocate(interpolator%gtau)

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
  SLL_ALLOCATE(interpolator%gtau(1:nx,1:ny),ierr)
  SLL_ALLOCATE(interpolator%gtau_der1(mx,1:ny),ierr)
  SLL_ALLOCATE(interpolator%gtau_der2(my,1:nx+mx),ierr)

  interpolator%gtau                 = data_array(1:nx,1:ny)
  taux_deriv(1)        = 1
  taux_deriv(2)        = nx
  interpolator%gtau_der1(1,1:ny)    = 0.0_f64
  interpolator%gtau_der1(2,1:ny)    = interpolator%slope_max1(1:ny)
  tauy_deriv(1)        = 1
  tauy_deriv(2)        = ny
  interpolator%gtau_der2(1,1:nx+mx) = interpolator%slope_min2(1:nx+mx)
  interpolator%gtau_der2(2,1:nx+mx) = 0.0_f64

  call spli2d_custom_derder(interpolator%deboor,       &
                            nx,                        &
                            mx,                        &
                            kx,                        &
                            taux,                      &
                            taux_deriv,                &
                            ny,                        &
                            my,                        &
                            ky, tauy,                  &
                            tauy_deriv,                &
                            interpolator%gtau,         &
                            interpolator%gtau_der1,    &
                            interpolator%gtau_der2,    &
                            interpolator%coeff_splines,&
                            interpolator%t1,           &
                            interpolator%t2)

  SLL_DEALLOCATE( interpolator%gtau_der1,ierr)
  SLL_DEALLOCATE( interpolator%gtau_der2,ierr)
  deallocate(interpolator%gtau)

case(1105)  !left: Dirichlet, right: Neumann, bottom: Dirichlet, Top: Neumann

   mx = 2
   my = 2
   interpolator%size_coeffs1 = nx + mx
   interpolator%size_coeffs2 = ny + my
   interpolator%size_t1 = kx + nx + mx
   interpolator%size_t2 = ky + ny + my
   
   SLL_ALLOCATE(interpolator%gtau(1:nx,1:ny),ierr)
   SLL_ALLOCATE(interpolator%gtau_der1(2,ny),ierr)
   SLL_ALLOCATE(interpolator%gtau_der2(my,nx+mx),ierr)

   interpolator%gtau                 = data_array(1:nx,1:ny)
   taux_deriv(1)        = 1
   taux_deriv(2)        = nx
   interpolator%gtau_der1(1,1:ny)    = interpolator%slope_min1(1:ny)
   interpolator%gtau_der1(2,1:ny)    = 0.0_f64
   tauy_deriv(1)        = 1
   tauy_deriv(2)        = ny
   interpolator%gtau_der2(1,1:nx+mx) = interpolator%slope_min2(1:nx+mx)
   interpolator%gtau_der2(2,1:nx+mx) = 0.0_f64

  call spli2d_custom_derder(interpolator%deboor,        &
                            nx,                         &
                            mx,                         &
                            kx,                         &
                            taux,                       &
                            taux_deriv,                 &
                            ny,                         &
                            my,                         &
                            ky, tauy,                   &
                            tauy_deriv,                 &
                            interpolator%gtau,          &
                            interpolator%gtau_der1,     &
                            interpolator%gtau_der2,     &
                            interpolator%coeff_splines, &
                            interpolator%t1,            &
                            interpolator%t2)

   SLL_DEALLOCATE(interpolator%gtau_der1,ierr)
   SLL_DEALLOCATE(interpolator%gtau_der2,ierr)
  deallocate(interpolator%gtau)

case(1170)  !left: Neumann, right: Neumann, bottom: Neuman, Top: Neumann

  mx = 2
  my = 2
  interpolator%size_coeffs1 = nx + mx
  interpolator%size_coeffs2 = ny + my
  interpolator%size_t1 = kx + nx + mx
  interpolator%size_t2 = ky + ny + my
  
  SLL_ALLOCATE(interpolator%gtau(1:nx,1:ny),ierr)
  SLL_ALLOCATE(interpolator%gtau_der1(2,ny),ierr)
  SLL_ALLOCATE(interpolator%gtau_der2(my,nx+mx),ierr)

  interpolator%gtau                 = data_array(1:nx,1:ny)
  taux_deriv(1)                     = 1
  taux_deriv(2)                     = nx
  interpolator%gtau_der1(1,1:ny)    = 0.0_f64
  interpolator%gtau_der1(2,1:ny)    = 0.0_f64
  tauy_deriv(1)                     = 1
  tauy_deriv(2)                     = ny
  interpolator%gtau_der2(1,1:nx+mx) = 0.0_f64
  interpolator%gtau_der2(2,1:nx+mx) = 0.0_f64

  call spli2d_custom_derder(interpolator%deboor,        &
                            nx,                         &
                            mx,                         &
                            kx,                         &
                            taux,                       &
                            taux_deriv,                 &
                            ny,                         &
                            my,                         &
                            ky, tauy,                   &
                            tauy_deriv,                 &
                            interpolator%gtau,          &
                            interpolator%gtau_der1,     &
                            interpolator%gtau_der2,     &
                            interpolator%coeff_splines, &
                            interpolator%t1,            &
                            interpolator%t2)

  SLL_DEALLOCATE( interpolator%gtau_der1,ierr)
  SLL_DEALLOCATE( interpolator%gtau_der2,ierr)
  deallocate(interpolator%gtau)

case(2338)  !left: Dirichlet, right: Hermite, bottom: Hermite, Top: Hermite

  mx = 2
  my = 2
  interpolator%size_coeffs1 = nx + mx
  interpolator%size_coeffs2 = ny + my
  interpolator%size_t1 = kx + nx + mx
  interpolator%size_t2 = ky + ny + my
  
  SLL_ALLOCATE(interpolator%gtau(1:nx,1:ny),ierr)
  SLL_ALLOCATE(interpolator%gtau_der1(2,ny),ierr)
  SLL_ALLOCATE(interpolator%gtau_der2(my,nx+mx),ierr)

  interpolator%gtau = data_array(1:nx,1:ny)
  taux_deriv(1) = 1
  taux_deriv(2) = nx
  interpolator%gtau_der1(1,1:ny)     = interpolator%slope_min1(1:ny)
  interpolator%gtau_der1(2,1:ny)     = interpolator%slope_max1(1:ny)
  tauy_deriv(1) = 1
  tauy_deriv(2) = ny
  interpolator%gtau_der2(1,1:nx+mx)=interpolator%slope_min2(1:nx+mx)
  interpolator%gtau_der2(2,1:nx+mx)=interpolator%slope_max2(1:nx+mx)

  call spli2d_custom_derder(interpolator%deboor,       &
                            nx,                        &
                            mx,                        &
                            kx,                        &
                            taux,                      &
                            taux_deriv,                &
                            ny,                        &
                            my,                        &
                            ky, tauy,                  &
                            tauy_deriv,                &
                            interpolator%gtau,         &
                            interpolator%gtau_der1,    &
                            interpolator%gtau_der2,    &
                            interpolator%coeff_splines,&
                            interpolator%t1,           &
                            interpolator%t2)

  SLL_DEALLOCATE( interpolator%gtau_der1,ierr)
  SLL_DEALLOCATE( interpolator%gtau_der2,ierr)
  deallocate(interpolator%gtau)
       
case(2145) !left: Dirichlet, right: Hermite, bottom: Dirichlet, Top: Hermite  

  mx = 2
  my = 2
  interpolator%size_coeffs1 = nx + mx
  interpolator%size_coeffs2 = ny + my
  interpolator%size_t1 = kx + nx + mx
  interpolator%size_t2 = ky + ny + my
  
  SLL_ALLOCATE( interpolator%gtau(1:nx,1:ny),ierr)
  SLL_ALLOCATE( interpolator%gtau_der1(2,ny),ierr)
  SLL_ALLOCATE( interpolator%gtau_der2(my,nx+mx),ierr)

  interpolator%gtau                  = data_array(1:nx,1:ny)
  taux_deriv(1)         = 1
  taux_deriv(2)         = nx
  interpolator%gtau_der1(1,1:ny)     = interpolator%slope_min1(1:ny)
  interpolator%gtau_der1(2,1:ny)     = interpolator%slope_max1(1:ny)
  tauy_deriv(1)         = 1
  tauy_deriv(2)         = ny
  interpolator%gtau_der2(1,1:nx+mx)  = interpolator%slope_min2(1:nx+mx)
  interpolator%gtau_der2(2,1:nx+mx)  = interpolator%slope_max2(1:nx+mx)

  call spli2d_custom_derder(interpolator%deboor,        &
                            nx,                         &
                            mx,                         &
                            kx,                         &
                            taux,                       &
                            taux_deriv,                 &
                            ny,                         &
                            my,                         &
                            ky, tauy,                   &
                            tauy_deriv,                 &
                            interpolator%gtau,          &
                            interpolator%gtau_der1,     &
                            interpolator%gtau_der2,     &
                            interpolator%coeff_splines, &
                            interpolator%t1,            &
                            interpolator%t2)

  SLL_DEALLOCATE( interpolator%gtau_der1,ierr)
  SLL_DEALLOCATE( interpolator%gtau_der2,ierr)
  deallocate(interpolator%gtau)

case(2124)  !left: Hermite, right: Dirichlet, bottom: Dirichlet, Top: Hermite

  mx = 2
  my = 2
  interpolator%size_coeffs1 = nx + mx
  interpolator%size_coeffs2 = ny + my
  interpolator%size_t1 = kx + nx + mx
  interpolator%size_t2 = ky + ny + my
  
  SLL_CLEAR_ALLOCATE(interpolator%gtau(1:nx,1:ny),ierr)
  SLL_CLEAR_ALLOCATE(interpolator%gtau_der1(mx,1:ny),ierr)
  SLL_CLEAR_ALLOCATE(interpolator%gtau_der2(my,1:nx+mx),ierr)

  interpolator%gtau                  = data_array(1:nx,1:ny)
  taux_deriv(1)         = 1
  taux_deriv(2)         = nx
  interpolator%gtau_der1(1,1:ny)     = interpolator%slope_min1(1:ny)
  interpolator%gtau_der1(2,1:ny)     = interpolator%slope_max1(1:ny)
  tauy_deriv(1)         = 1
  tauy_deriv(2)         = ny
  interpolator%gtau_der2(1,1:nx+mx)  = interpolator%slope_min2(1:nx+mx)
  interpolator%gtau_der2(2,1:nx+mx)  = interpolator%slope_max2(1:nx+mx)

  call spli2d_custom_derder(interpolator%deboor,       &
                            nx,                        &
                            mx,                        &
                            kx,                        &
                            taux,                      &
                            taux_deriv,                &
                            ny,                        &
                            my,                        &
                            ky, tauy,                  &
                            tauy_deriv,                &
                            interpolator%gtau,         &
                            interpolator%gtau_der1,    &
                            interpolator%gtau_der2,    &
                            interpolator%coeff_splines,&
                            interpolator%t1,           &
                            interpolator%t2)

  SLL_DEALLOCATE( interpolator%gtau_der1,ierr)
  SLL_DEALLOCATE( interpolator%gtau_der2,ierr)
  deallocate(interpolator%gtau)

case(2148)  !left:Hermite , right: Hermite, bottom: Dirichlet, Top: Hermite

  mx = 2
  my = 2
  interpolator%size_coeffs1 = nx + mx
  interpolator%size_coeffs2 = ny + my
  interpolator%size_t1 = kx + nx + mx
  interpolator%size_t2 = ky + ny + my
  
  SLL_ALLOCATE(interpolator%gtau(1:nx,1:ny),ierr)
  SLL_ALLOCATE(interpolator%gtau_der1(mx,1:ny),ierr)
  SLL_ALLOCATE(interpolator%gtau_der2(my,1:nx+mx),ierr)

  interpolator%gtau                  = data_array(1:nx,1:ny)
  taux_deriv(1)         = 1
  taux_deriv(2)         = nx
  interpolator%gtau_der1(1,1:ny)     = interpolator%slope_min1(1:ny)
  interpolator%gtau_der1(2,1:ny)     = interpolator%slope_max1(1:ny)
  tauy_deriv(1)         = 1
  tauy_deriv(2)         = ny
  interpolator%gtau_der2(1,1:nx+mx)  = interpolator%slope_min2(1:nx+mx)
  interpolator%gtau_der2(2,1:nx+mx)  = interpolator%slope_max2(1:nx+mx)

  call spli2d_custom_derder(interpolator%deboor,        &
                            nx,                         &
                            mx,                         &
                            kx,                         &
                            taux,                       &
                            taux_deriv,                 &
                            ny,                         &
                            my,                         &
                            ky, tauy,                   &
                            tauy_deriv,                 &
                            interpolator%gtau,          &
                            interpolator%gtau_der1,     &
                            interpolator%gtau_der2,     &
                            interpolator%coeff_splines, &
                            interpolator%t1,            &
                            interpolator%t2)

  SLL_DEALLOCATE( interpolator%gtau_der1,ierr)
  SLL_DEALLOCATE( interpolator%gtau_der2,ierr)
  deallocate(interpolator%gtau)

case(2316)  !left: Hermite, right: Dirichlet, bottom: Hermite, Top: Hermite

  mx = 2
  my = 2
  interpolator%size_coeffs1 = nx + mx
  interpolator%size_coeffs2 = ny + my
  interpolator%size_t1 = kx + nx + mx
  interpolator%size_t2 = ky + ny + my
  
  SLL_ALLOCATE(interpolator%gtau(1:nx,1:ny),ierr)
  SLL_ALLOCATE(interpolator%gtau_der1(mx,1:ny),ierr)
  SLL_ALLOCATE(interpolator%gtau_der2(my,1:nx+mx),ierr)

  interpolator%gtau                  = data_array(1:nx,1:ny)
  taux_deriv(1)         = 1
  taux_deriv(2)         = nx
  interpolator%gtau_der1(1,1:ny)     = interpolator%slope_min1(1:ny)
  interpolator%gtau_der1(2,1:ny)     = interpolator%slope_max1(1:ny)
  tauy_deriv(1)         = 1
  tauy_deriv(2)         = ny
  interpolator%gtau_der2(1,1:nx+mx)  = interpolator%slope_min2(1:nx+mx)
  interpolator%gtau_der2(2,1:nx+mx)  = interpolator%slope_max2(1:nx+mx)

  call spli2d_custom_derder(interpolator%deboor,        &
                            nx,                         &
                            mx,                         &
                            kx,                         &
                            taux,                       &
                            taux_deriv,                 &
                            ny,                         &
                            my,                         &
                            ky, tauy,                   &
                            tauy_deriv,                 &
                            interpolator%gtau,          &
                            interpolator%gtau_der1,     &
                            interpolator%gtau_der2,     &
                            interpolator%coeff_splines, &
                            interpolator%t1,            &
                            interpolator%t2)

  SLL_DEALLOCATE( interpolator%gtau_der1,ierr)
  SLL_DEALLOCATE( interpolator%gtau_der2,ierr)
  deallocate(interpolator%gtau)

case(2340) ! Hermite in al sides
       
  mx = 2
  my = 2
  interpolator%size_coeffs1 = nx + mx
  interpolator%size_coeffs2 = ny + my
  interpolator%size_t1 = kx + nx + mx
  interpolator%size_t2 = ky + ny + my
  
  SLL_ALLOCATE(interpolator%gtau(1:nx,1:ny),ierr)
  SLL_ALLOCATE(interpolator%gtau_der1(mx,1:ny),ierr)
  SLL_ALLOCATE(interpolator%gtau_der2(my,1:nx+mx),ierr)

  interpolator%gtau                  = data_array(1:nx,1:ny)
  taux_deriv(1)         = 1
  taux_deriv(2)         = nx
  interpolator%gtau_der1(1,1:ny)     = interpolator%slope_min1(1:ny)
  interpolator%gtau_der1(2,1:ny)     = interpolator%slope_max1(1:ny)
  tauy_deriv(1)         = 1
  tauy_deriv(2)         = ny
  interpolator%gtau_der2(1,1:nx+mx)  = interpolator%slope_min2(1:nx+mx)
  interpolator%gtau_der2(2,1:nx+mx)  = interpolator%slope_max2(1:nx+mx)

  call spli2d_custom_derder(interpolator%deboor,        &
                            nx,                         &
                            mx,                         &
                            kx,                         &
                            taux,                       &
                            taux_deriv,                 &
                            ny,                         &
                            my,                         &
                            ky, tauy,                   &
                            tauy_deriv,                 &
                            interpolator%gtau,          &
                            interpolator%gtau_der1,     &
                            interpolator%gtau_der2,     &
                            interpolator%coeff_splines, &
                            interpolator%t1,            &
                            interpolator%t2)

  SLL_DEALLOCATE( interpolator%gtau_der1,ierr)
  SLL_DEALLOCATE( interpolator%gtau_der2,ierr)
  deallocate(interpolator%gtau)

end select

interpolator%coefficients_set = .true.

SLL_DEALLOCATE(taux,ierr)
SLL_DEALLOCATE(tauy,ierr)

end subroutine compute_interpolants_ad2d

function coefficients_are_set_ad2d( interpolator ) result(res)

  class(sll_arbitrary_degree_spline_interpolator_2d), intent(in)  :: interpolator
  logical :: res

  res = interpolator%coefficients_set

end function coefficients_are_set_ad2d

!> @brief Interpolation on the points eta1 and eta2 
!> @details computing the values with the interpolator arbitrary degree splines 2d
!>  on the points eta1 and eta2 of arbitrary degree splines 2d
!> 
!> The parameters are
!> @param interpolator the type sll_arbitrary_degree_spline_interpolator_2d
!> @param[in] eta1 points in the first direction
!> @param[in] eta2 points in the second direction 
!> @return val the values on the points eta1 and eta2 
function interpolate_value_ad2d( interpolator, eta1, eta2 ) result(val)

class(sll_arbitrary_degree_spline_interpolator_2d), intent(in)  :: interpolator

sll_real64, intent(in) :: eta1
sll_real64, intent(in) :: eta2
sll_real64             :: val
sll_int32              :: nx
sll_int32              :: ny
sll_int32              :: kx
sll_int32              :: ky

sll_real64                        :: x
sll_real64                        :: y
sll_real64                        :: length1
sll_real64                        :: length2
sll_real64, dimension(:), pointer :: ty
sll_real64, dimension(:), pointer :: tab
sll_real64, dimension(:), pointer :: coef

sll_int32 :: j
sll_int32 :: lefty
sll_int32 :: ierr

sll_real64, dimension(:),   pointer :: t1
sll_real64, dimension(:),   pointer :: t2
sll_real64, dimension(:,:), pointer :: coeff

nx = interpolator%size_coeffs1
ny = interpolator%size_coeffs2
kx = interpolator%spline_degree1+1
ky = interpolator%spline_degree2+1

SLL_ALLOCATE(coef(ky),ierr)
SLL_ALLOCATE(tab(nx),ierr)
SLL_ALLOCATE(ty(2*ky),ierr)

length1 = interpolator%eta1_max-interpolator%eta1_min
length2 = interpolator%eta2_max-interpolator%eta2_min

x   = eta1
y   = eta2
val = 0.0_f64

!if (interpolator%bc1_min==SLL_PERIODIC .and. eta1<interpolator%eta1_min) x = eta1 + length1 
!if (interpolator%bc1_max==SLL_PERIODIC .and. eta1>interpolator%eta1_max) x = eta1 - length1
!if (interpolator%bc2_min==SLL_PERIODIC .and. eta2<interpolator%eta2_min) y = eta2 + length2
!if (interpolator%bc2_max==SLL_PERIODIC .and. eta2>interpolator%eta2_max) y = eta2 - length2

SLL_ASSERT( x >= interpolator%eta1_min )
SLL_ASSERT( x <= interpolator%eta1_max )
SLL_ASSERT( y >= interpolator%eta2_min )
SLL_ASSERT( y <= interpolator%eta2_max )

t1    => interpolator%t1(1:interpolator%size_t1)
t2    => interpolator%t2(1:interpolator%size_t2)
coeff => interpolator%coeff_splines(1:nx,1:ny)

call interv(interpolator%deboor(2), t2(1:ny+ky), ny+ky, y, lefty, ierr)

if (ierr .ne. 0) return 

do j = 1, ky
  tab     = interpolator%coeff_splines(1:nx,lefty-ky+j)
  coef(j) = bvalue(interpolator%deboor(1), t1(1:nx+kx), tab, nx, kx, x, 0)
end do

ty = interpolator%t2(lefty-ky+1:lefty+ky)

val = bvalue(interpolator%deboor(2), ty, coef, ky, ky, y, 0)

deallocate(tab)
deallocate(coef)
deallocate(ty)

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
function interpolate_derivative1_ad2d( interpolator, eta1, eta2 ) result(val)

class(sll_arbitrary_degree_spline_interpolator_2d), intent(in)  :: interpolator

sll_real64, intent(in)         :: eta1
sll_real64, intent(in)         :: eta2
sll_real64                     :: val

sll_int32                      :: nx
sll_int32                      :: ny
sll_int32                      :: kx
sll_int32                      :: ky
sll_real64                     :: x
sll_real64                     :: y
sll_real64                     :: length1
sll_real64                     :: length2

sll_int32, parameter  :: deriv1 = 1
sll_int32, parameter  :: deriv2 = 0

sll_real64, dimension(:), pointer :: coef 
sll_real64, dimension(:), pointer :: tab
sll_real64, dimension(:), pointer :: ty

sll_real64, dimension(:),   pointer :: t1
sll_real64, dimension(:),   pointer :: t2
sll_real64, dimension(:,:), pointer :: coeff

sll_int32 :: j
sll_int32 :: lefty
sll_int32 :: ierr

val = 0.0_f64

length1 = interpolator%eta1_max-interpolator%eta1_min
length2 = interpolator%eta2_max-interpolator%eta2_min

nx = interpolator%size_coeffs1
ny = interpolator%size_coeffs2
kx = interpolator%spline_degree1+1
ky = interpolator%spline_degree2+1

x = eta1
y = eta2

if (interpolator%bc1_min==SLL_PERIODIC .and. eta1<interpolator%eta1_min) x = eta1+length1 
if (interpolator%bc1_max==SLL_PERIODIC .and. eta1>interpolator%eta1_max) x = eta1-length1
if (interpolator%bc2_min==SLL_PERIODIC .and. eta2<interpolator%eta2_min) y = eta2+length2
if (interpolator%bc2_max==SLL_PERIODIC .and. eta2>interpolator%eta2_max) y = eta2-length2

SLL_ASSERT( x >= interpolator%eta1_min )
SLL_ASSERT( x <= interpolator%eta1_max )
SLL_ASSERT( y >= interpolator%eta2_min )
SLL_ASSERT( y <= interpolator%eta2_max )

t1 => interpolator%t1(1:interpolator%size_t1)
t2 => interpolator%t2(1:interpolator%size_t2)
coeff => interpolator%coeff_splines(1:nx,1:nx)

call interv(interpolator%deboor(2),t2,ny+ky,y,lefty,ierr)
    
if (ierr .ne. 0)  return 

SLL_ALLOCATE(coef(1:ky), ierr) 
SLL_ALLOCATE(tab(1:nx),  ierr)
SLL_ALLOCATE(ty(1:2*ky), ierr)
    
do j = 1, ky
       
  tab     = interpolator%coeff_splines(1:nx,lefty-ky+j)
  coef(j) = bvalue(interpolator%deboor(1), t1, tab, nx, kx, x, deriv1 )
       
end do

ty = t2(lefty-ky+1:lefty+ky)

val = bvalue(interpolator%deboor(2), ty, coef, ky, ky, y, deriv2 )

deallocate(ty)
deallocate(tab)
deallocate(coef)

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
function interpolate_derivative2_ad2d( interpolator, eta1, eta2 ) result(val)

class(sll_arbitrary_degree_spline_interpolator_2d), intent(in)  :: interpolator

sll_real64, intent(in)         :: eta1
sll_real64, intent(in)         :: eta2
sll_real64                     :: val
sll_int32                      :: nx
sll_int32                      :: ny
sll_int32                      :: kx
sll_int32                      :: ky
sll_real64                     :: x
sll_real64                     :: y
sll_real64                     :: length1
sll_real64                     :: length2

sll_real64, dimension(:), pointer :: coef
sll_real64, dimension(:), pointer :: tab
sll_real64, dimension(:), pointer :: ty
sll_int32, parameter  :: deriv1 = 0
sll_int32, parameter  :: deriv2 = 1
sll_int32 :: j
sll_int32 :: lefty
sll_int32 :: ierr = 0
sll_real64, dimension(:),   pointer :: t1
sll_real64, dimension(:),   pointer :: t2
sll_real64, dimension(:,:), pointer :: coeff

val = 0.0_f64 

nx = interpolator%size_coeffs1
ny = interpolator%size_coeffs2
kx = interpolator%spline_degree1+1
ky = interpolator%spline_degree2+1
length1 = interpolator%eta1_max-interpolator%eta1_min
length2 = interpolator%eta2_max-interpolator%eta2_min

x = eta1
y = eta2
if (interpolator%bc1_min==SLL_PERIODIC .and. eta1<interpolator%eta1_min) x = eta1+length1 
if (interpolator%bc1_max==SLL_PERIODIC .and. eta1>interpolator%eta1_max) x = eta1-length1
if (interpolator%bc2_min==SLL_PERIODIC .and. eta2<interpolator%eta2_min) y = eta2+length2
if (interpolator%bc2_max==SLL_PERIODIC .and. eta2>interpolator%eta2_max) y = eta2-length2

SLL_ASSERT( x >= interpolator%eta1_min )
SLL_ASSERT( x <= interpolator%eta1_max )
SLL_ASSERT( y >= interpolator%eta2_min )
SLL_ASSERT( y <= interpolator%eta2_max )

t1 => interpolator%t1(1:interpolator%size_t1)
t2 => interpolator%t2(1:interpolator%size_t2)
coeff => interpolator%coeff_splines(1:nx,1:ny)

call interv(interpolator%deboor(2), t2(1:ny+ky), ny+ky, y, lefty, ierr )
    
if ( ierr .ne. 0 ) return

SLL_ALLOCATE(coef(1:ky), ierr) 
SLL_ALLOCATE(tab(1:nx),  ierr)
SLL_ALLOCATE(ty(1:2*ky), ierr)
    
do j = 1, ky
       
  tab = interpolator%coeff_splines(1:nx, lefty-ky+j)
  coef(j) = bvalue(interpolator%deboor(1), t1, tab, nx, kx, x, deriv1 )
       
end do

ty =  t2(lefty-ky+1:lefty+ky)
val = bvalue(interpolator%deboor(2), ty, coef, ky, ky, y, deriv2 )

deallocate(ty)
deallocate(tab)
deallocate(coef)

end function interpolate_derivative2_ad2d

  
subroutine interpolate_array_ad2d( this,            &
                                 num_points1,     &
                                 num_points2,     &
                                 data_in,         &
                                 eta1,            &
                                 eta2,            &
                                 data_out)
  
class(sll_arbitrary_degree_spline_interpolator_2d), intent(in)  :: this

sll_real64, dimension(:,:), intent(in) :: eta1
sll_real64, dimension(:,:), intent(in) :: eta2
sll_real64, dimension(:,:), intent(in) :: data_in
sll_int32,                  intent(in) :: num_points1
sll_int32,                  intent(in) :: num_points2

sll_real64,                 intent(out):: data_out(num_points1, num_points2)

print *, '#interpolate_array_ad2d: not implemented'
data_out = -1000000._f64
print *,this%num_pts1
print *,maxval(eta1)
print *,maxval(eta2)
print *,maxval(data_in)
print *,num_points1
print *,num_points2
stop
end subroutine interpolate_array_ad2d  !interpolate_array_ad2d
  
subroutine interpolate_2d_array_disp_ad2d( this,        &
                                         num_points1, &
                                         num_points2, &
                                         data_in,     &
                                         alpha1,      &
                                         alpha2,      &
                                         data_out)
    
class(sll_arbitrary_degree_spline_interpolator_2d), intent(in)    :: this

sll_int32,                  intent(in)         :: num_points1  
sll_int32,                  intent(in)         :: num_points2 
sll_real64, dimension(:,:), intent(in)         :: data_in
sll_real64, dimension(:,:), intent(in)         :: alpha1
sll_real64, dimension(:,:), intent(in)         :: alpha2  
sll_real64,                 intent(out)        :: data_out(num_points1,num_points2)

print *, '#interpolate_2d_array_disp_ad2d: not implemented.'
!for preventing warning of unused objects
print *,this%num_pts1
print *,num_points1 
print *,num_points2
print *,maxval(data_in)
print *,alpha1
print *,alpha2     
data_out = -1000000._f64
stop
  
end subroutine interpolate_2d_array_disp_ad2d  !interpolate_2d_array_disp_ad2d
    
function get_coefficients_ad2d(interpolator)
class(sll_arbitrary_degree_spline_interpolator_2d), intent(in) :: interpolator
sll_real64, dimension(:,:), pointer                            :: get_coefficients_ad2d     

get_coefficients_ad2d => interpolator%coeff_splines

end function get_coefficients_ad2d
  
!> Initialization of the boundary for interpolator arbitrary degree splines 2d.
!> The parameters are
!> @param[in] slope_min1 a 1d arrays contains values in the left in the direction eta1  
!> @param[in] slope_max1 a 1d arrays contains values in the right in the direction eta1 
!> @param[in] slope_min2 a 1d arrays contains values in the left in the direction eta2 
!> @param[in] slope_max2 a 1d arrays contains values in the right in the direction eta2
!> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_2d
subroutine set_slope2d( interpolator, &
                        slope_min1,   &
                        slope_max1,   &
                        slope_min2,   &
                        slope_max2)

class(sll_arbitrary_degree_spline_interpolator_2d)    :: interpolator
sll_real64, dimension(:),optional :: slope_min1
sll_real64, dimension(:),optional :: slope_max1
sll_real64, dimension(:),optional :: slope_min2
sll_real64, dimension(:),optional :: slope_max2
class(sll_arbitrary_degree_spline_interpolator_1d),pointer :: interp1d_min2=> null()
class(sll_arbitrary_degree_spline_interpolator_1d),pointer :: interp1d_max2=> null()
sll_int32 :: sz_slope_min2,sz_slope_max2
sll_int64 :: bc_selector
sll_int32 :: num_pts1
sll_int32 :: num_pts2
sll_int32 :: bc1_min
sll_int32 :: bc1_max
sll_int32 :: bc2_min
sll_int32 :: bc2_max

num_pts1    = interpolator%num_pts1
num_pts2    = interpolator%num_pts2
bc_selector = interpolator%bc_selector
bc1_min     = interpolator%bc1_min 
bc1_max     = interpolator%bc1_max 
bc2_min     = interpolator%bc2_min  
bc2_max     = interpolator%bc2_max

select case (bc_selector)
case(0)
case (9) ! dirichlet-left, dirichlet-right, periodic
case (576) ! 3. periodic, dirichlet-bottom, dirichlet-top
case (585) ! 4. dirichlet in all sides
case (650) !left: Neumann, right: Dirichlet, bottom: Neumann, Top: Dirichlet

  interpolator%slope_min1 = 0.0_f64
  interpolator%compute_slope_min1= .FALSE.
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
          interpolator%bc1_min, &
          interpolator%bc1_max, &
          interpolator%spline_degree1 )
     
     call interp1d_max2%compute_interpolants(slope_max2(1:sz_slope_max2))
     
     interpolator%slope_max2(1:sz_slope_max2+2) = &
          interp1d_max2%coeff_splines(1:sz_slope_max2+2)
     call sll_delete(interp1d_max2)
     deallocate(interp1d_max2)
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
          interpolator%bc1_min, &
          interpolator%bc1_max, &
          interpolator%spline_degree1 )
     
     call interp1d_max2%compute_interpolants(&
          slope_max2(1:sz_slope_max2))
     
     interpolator%slope_max2(1:sz_slope_max2+2) = &
          interp1d_max2%coeff_splines(1:sz_slope_max2+2)
     call sll_delete(interp1d_max2)
     deallocate(interp1d_max2)
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
          interpolator%bc1_min, &
          interpolator%bc1_max, &
          interpolator%spline_degree1 )
     
     call interp1d_min2%compute_interpolants( &
          slope_min2(1:sz_slope_min2))
     
     interpolator%slope_min2(1:sz_slope_min2+2) = &
          interp1d_min2%coeff_splines(1:sz_slope_min2+2)
     call sll_delete(interp1d_min2)
     deallocate(interp1d_min2)
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
          interpolator%bc1_min, &
          interpolator%bc1_max, &
          interpolator%spline_degree1 )
     
     call interp1d_max2%compute_interpolants(&
          slope_max2(1:sz_slope_max2))
     
     interpolator%slope_max2(1:sz_slope_max2+2) = &
          interp1d_max2%coeff_splines(1:sz_slope_max2+2)
     call sll_delete(interp1d_max2)
     deallocate(interp1d_max2)
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
           interpolator%spline_degree1 )
      
      call interp1d_max2%compute_interpolants(&
           slope_max2(1:sz_slope_max2))
      
      interpolator%slope_max2(1:sz_slope_max2+2) = &
           interp1d_max2%coeff_splines(1:sz_slope_max2+2)
      call sll_delete(interp1d_max2)
      deallocate(interp1d_max2)
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
           interpolator%spline_degree1 )
      
      call interp1d_min2%compute_interpolants( &
           slope_min2(1:sz_slope_min2))
      
      interpolator%slope_min2(1:sz_slope_min2+2) = &
           interp1d_min2%coeff_splines(1:sz_slope_min2+2)
      call sll_delete(interp1d_min2)
      deallocate(interp1d_min2)
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
           interpolator%spline_degree1 )
      
      call interp1d_max2%compute_interpolants(&
           slope_max2(1:sz_slope_max2))
      
      interpolator%slope_max2(1:sz_slope_max2+2) = &
           interp1d_max2%coeff_splines(1:sz_slope_max2+2)
      call sll_delete(interp1d_max2)
      deallocate(interp1d_max2)
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
           interpolator%spline_degree1 )
      
      call interp1d_min2%compute_interpolants( &
           slope_min2(1:sz_slope_min2))
      
      interpolator%slope_min2(1:sz_slope_min2+2) = &
           interp1d_min2%coeff_splines(1:sz_slope_min2+2)
      call sll_delete(interp1d_min2)
      deallocate(interp1d_min2)
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
           interpolator%spline_degree1 )
      
      call interp1d_max2%compute_interpolants(&
           slope_max2(1:sz_slope_max2))
      
      interpolator%slope_max2(1:sz_slope_max2+2) = &
           interp1d_max2%coeff_splines(1:sz_slope_max2+2)
      call sll_delete(interp1d_max2)
      deallocate(interp1d_max2)
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
           interpolator%spline_degree1 )
      
      call interp1d_min2%compute_interpolants( &
           slope_min2(1:sz_slope_min2))
      
      interpolator%slope_min2(1:sz_slope_min2+2) = &
           interp1d_min2%coeff_splines(1:sz_slope_min2+2)
      call sll_delete(interp1d_min2)
      deallocate(interp1d_min2)
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
           interpolator%spline_degree1 )
      
      call interp1d_max2%compute_interpolants(&
           slope_max2(1:sz_slope_max2))
      
      interpolator%slope_max2(1:sz_slope_max2+2) = &
           interp1d_max2%coeff_splines(1:sz_slope_max2+2)
      call sll_delete(interp1d_max2)
      deallocate(interp1d_max2)
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
           interpolator%spline_degree1 )
      
      call interp1d_min2%compute_interpolants( &
           slope_min2(1:sz_slope_min2))
      
      interpolator%slope_min2(1:sz_slope_min2+2) = &
           interp1d_min2%coeff_splines(1:sz_slope_min2+2)
      call sll_delete(interp1d_min2)
      deallocate(interp1d_min2)
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
           interpolator%spline_degree1 )
      
      call interp1d_max2%compute_interpolants(&
           slope_max2(1:sz_slope_max2))
      
      interpolator%slope_max2(1:sz_slope_max2+2) = &
           interp1d_max2%coeff_splines(1:sz_slope_max2+2)
      call sll_delete(interp1d_max2)
      deallocate(interp1d_max2)
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
           interpolator%spline_degree1 )
      
      call interp1d_min2%compute_interpolants( &
           slope_min2(1:sz_slope_min2))
      
      interpolator%slope_min2(1:sz_slope_min2+2) = &
           interp1d_min2%coeff_splines(1:sz_slope_min2+2)
      call sll_delete(interp1d_min2)
      deallocate(interp1d_min2)
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
           interpolator%spline_degree1 )
      
      call interp1d_max2%compute_interpolants(&
           slope_max2(1:sz_slope_max2))
      
      interpolator%slope_max2(1:sz_slope_max2+2) = &
           interp1d_max2%coeff_splines(1:sz_slope_max2+2)

      call sll_delete(interp1d_max2)
      deallocate(interp1d_max2)
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
           interpolator%spline_degree1 )
      
      call interp1d_min2%compute_interpolants( &
           slope_min2(1:sz_slope_min2))
      
      interpolator%slope_min2(1:sz_slope_min2+2) = &
           interp1d_min2%coeff_splines(1:sz_slope_min2+2)
      call sll_delete(interp1d_min2)
      deallocate(interp1d_min2)
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
           interpolator%spline_degree1 )
      
      call interp1d_max2%compute_interpolants(&
           slope_max2(1:sz_slope_max2))
      
      interpolator%slope_max2(1:sz_slope_max2+2) = &
           interp1d_max2%coeff_splines(1:sz_slope_max2+2)
      call sll_delete(interp1d_max2)
      deallocate(interp1d_max2)
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
           interpolator%spline_degree1 )
      
      call interp1d_max2%compute_interpolants(&
           slope_max2(1:sz_slope_max2))
      
      interpolator%slope_max2(1:sz_slope_max2+2) = &
           interp1d_max2%coeff_splines(1:sz_slope_max2+2)

      call sll_delete(interp1d_max2)
      deallocate(interp1d_max2)
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
           interpolator%spline_degree1 )
      
      call interp1d_min2%compute_interpolants( &
           slope_min2(1:sz_slope_min2))
      
      interpolator%slope_min2(1:sz_slope_min2+2) = &
           interp1d_min2%coeff_splines(1:sz_slope_min2+2)
      call sll_delete(interp1d_min2)
      deallocate(interp1d_min2)
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
           interpolator%spline_degree1 )
      
      call interp1d_max2%compute_interpolants(&
           slope_max2(1:sz_slope_max2))
      
      interpolator%slope_max2(1:sz_slope_max2+2) = &
           interp1d_max2%coeff_splines(1:sz_slope_max2+2)
      call sll_delete(interp1d_max2)
      deallocate(interp1d_max2)
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
           interpolator%spline_degree1 )
      
      call interp1d_max2%compute_interpolants(&
           slope_max2(1:sz_slope_max2))
      
      interpolator%slope_max2(1:sz_slope_max2+2) = &
           interp1d_max2%coeff_splines(1:sz_slope_max2+2)
      call sll_delete(interp1d_max2)
      deallocate(interp1d_max2)
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
           interpolator%spline_degree1 )
      
      call interp1d_min2%compute_interpolants( &
           slope_min2(1:sz_slope_min2))
      
      interpolator%slope_min2(1:sz_slope_min2+2) = &
           interp1d_min2%coeff_splines(1:sz_slope_min2+2)
      call sll_delete(interp1d_min2)
      deallocate(interp1d_min2)
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
           interpolator%spline_degree1 )
      
      call interp1d_max2%compute_interpolants(&
           slope_max2(1:sz_slope_max2))
      
      interpolator%slope_max2(1:sz_slope_max2+2) = &
           interp1d_max2%coeff_splines(1:sz_slope_max2+2)
      call sll_delete(interp1d_max2)
      deallocate(interp1d_max2)
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
           interpolator%bc1_min, &
           interpolator%bc1_max, &
           interpolator%spline_degree1 )
      
      call interp1d_min2%compute_interpolants( &
           slope_min2(1:sz_slope_min2))
      
      interpolator%slope_min2(1:sz_slope_min2+2) = &
           interp1d_min2%coeff_splines(1:sz_slope_min2+2)
      call sll_delete(interp1d_min2)
      deallocate(interp1d_min2)
      interpolator%compute_slope_min2 = .FALSE.
   else
      print*, 'problem with slope bottom in case 2340'
   end if
   
case default
   print*,'initialize_ad2d_interpolator: BC combination not implemented.'
end select
  
end subroutine set_slope2d
  
subroutine spli2d_custom(db, nx, kx, taux, ny, ky, tauy, g, bcoef, tx, ty)

type(deboor_type)                                :: db(2)
sll_int32,                           intent(in)  :: nx
sll_int32,                           intent(in)  :: kx
sll_int32,                           intent(in)  :: ny
sll_int32,                           intent(in)  :: ky
sll_real64, dimension(:),            intent(in)  :: taux
sll_real64, dimension(:),            intent(in)  :: tauy
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

sll_int32 :: i, flag

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

call spli2d( db(1), taux, bcoef, tx, nx, kx, ny, work_x, qx, pwork, flag)
call spli2d( db(2), tauy, pwork, ty, ny, ky, nx, work_y, qy, bcoef, flag)

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
subroutine spli2d ( db, tau, gtau, t, n, k, m, work, q, bcoef, iflag )
    
type(deboor_type)                                :: db
sll_real64, dimension(:),            intent(in)  :: tau
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
   call bsplvb ( db, t, k, 1, taui, left, work )
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


subroutine spli2d_custom_derder ( db,         &
                                  nx,         &
                                  mx,         &
                                  kx,         &
                                  taux,       &
                                  taux_der,   &
                                  ny,         &
                                  my,         &
                                  ky,         &
                                  tauy,       &
                                  tauy_der,   &
                                  gtau,       &
                                  gtau_der1,  &
                                  gtau_der2,  &
                                  bcoef,      &
                                  tx,         &
                                  ty          )

type(deboor_type)                               :: db(2)
sll_int32,                           intent(in) :: nx
sll_int32,                           intent(in) :: kx
sll_int32,                           intent(in) :: ny
sll_int32,                           intent(in) :: ky
sll_int32,                           intent(in) :: mx
sll_int32,                           intent(in) :: my
sll_real64, dimension(:),   pointer, intent(in) :: taux 
sll_real64, dimension(:),   pointer, intent(in) :: tauy 
sll_int32,  dimension(:),            intent(in) :: taux_der 
sll_int32,  dimension(:),            intent(in) :: tauy_der
sll_real64, dimension(:,:), pointer, intent(in) :: gtau    
sll_real64, dimension(:,:), pointer, intent(in) :: gtau_der1 
sll_real64, dimension(:,:), pointer, intent(in) :: gtau_der2

sll_real64, dimension(:,:), pointer, intent(out) :: bcoef
sll_real64, dimension( : ), pointer, intent(out) :: tx 
sll_real64, dimension( : ), pointer, intent(out) :: ty 


sll_real64, dimension(nx+mx)                :: wx
sll_real64, dimension(ny+my)                :: wy
sll_real64, dimension((nx+mx)*(2*kx-1))     :: qx
sll_real64, dimension((ny+my)*(2*ky-1))     :: qy
sll_real64, dimension(1:ny,1:nx+mx)         :: tmp 

sll_int32 :: i, j
sll_int32 :: ierr

tx(1:kx)             = taux(1)
tx(nx+mx+1:nx+mx+kx) = taux(nx)

if (nx+mx+kx == nx+2*(kx-1)) then
   tx(kx+1:nx+mx) = taux(2:nx-1)
else
   stop 'problem with construction of knots' 
end if

ty(1:ky)             = tauy(1)
ty(ny+my+1:ny+my+ky) = tauy(ny)

if (ny+my+ky == ny+2*(ky-1)) then
   ty(ky+1:ny+my) = tauy(2:ny-1)
else
   stop 'problem with construction of knots' 
end if

do j = 1, ny
   
   call splint_der( db(1),          &
                    taux,           &
                    gtau(1:nx,j),   &
                    taux_der,       &
                    gtau_der1(:,j), &
                    tx,             &
                    nx,             &
                    mx,             &
                    kx,             &
                    qx,             &
                    wx,             &
                    ierr )

   tmp(j,1:nx+mx) = wx

end do

do i = 1, nx+mx
   
    call splint_der( db(2),          &
                     tauy,           &
                     tmp(1:ny,i),    &
                     tauy_der,       &
                     gtau_der2(:,i), &
                     ty,             &
                     ny,             &
                     my,             &
                     ky,             &
                     qy,             &
                     wy,             &
                     ierr )

   bcoef(i,1:ny+my) = wy

end do


end subroutine spli2d_custom_derder

end module sll_m_arbitrary_degree_spline_interpolator_2d
