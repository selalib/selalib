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

#define sll_interpolator class(sll_arbitrary_degree_spline_interpolator_2d)

!> @ingroup interpolators
!> @brief
!> Class of arbitrary degree version of 2d interpolator
!> @details
module sll_module_arbitrary_degree_spline_interpolator_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h" 
#include "sll_utilities.h"

use sll_module_interpolators_2d_base
use sll_module_arbitrary_degree_spline_interpolator_1d
use sll_module_deboor_splines_2d

implicit none
private

integer, parameter :: N = 1
integer, parameter :: E = 2
integer, parameter :: S = 3
integer, parameter :: W = 4

type :: edge
   sll_int32           :: bc  
   sll_real64, pointer :: values(:)
   sll_real64, pointer :: slopes(:)
end type edge


! in what follows, the direction '1' is in the contiguous memory direction.
!> Arbitrary degree version of 2d irnterpolator
type, public, extends(sll_interpolator_2d_base) :: &
     sll_arbitrary_degree_spline_interpolator_2d           

   private

   sll_int32,  public          :: num_pts1          !< PLEASE ADD DOCUMENTATION
   sll_int32,  public          :: num_pts2          !< PLEASE ADD DOCUMENTATION
   sll_real64, public, pointer :: coeff_splines(:,:)!< PLEASE ADD DOCUMENTATION
   sll_int32,  public          :: size_coeffs1      !< PLEASE ADD DOCUMENTATION
   sll_int32,  public          :: size_coeffs2      !< PLEASE ADD DOCUMENTATION

   sll_real64 :: eta1_min           !< PLEASE ADD DOCUMENTATION
   sll_real64 :: eta1_max           !< PLEASE ADD DOCUMENTATION
   sll_real64 :: eta2_min           !< PLEASE ADD DOCUMENTATION
   sll_real64 :: eta2_max           !< PLEASE ADD DOCUMENTATION
   sll_int32  :: spline_degree1     !< PLEASE ADD DOCUMENTATION
   sll_int32  :: spline_degree2     !< PLEASE ADD DOCUMENTATION
   sll_real64, pointer :: knots1(:) !< PLEASE ADD DOCUMENTATION
   sll_real64, pointer :: knots2(:) !< PLEASE ADD DOCUMENTATION
   sll_real64, pointer :: t1(:)     !< PLEASE ADD DOCUMENTATION
   sll_real64, pointer :: t2(:)     !< PLEASE ADD DOCUMENTATION
   sll_int32  :: size_t1            !< PLEASE ADD DOCUMENTATION
   sll_int32  :: size_t2            !< PLEASE ADD DOCUMENTATION
   sll_int64  :: bc_selector        !< this is set in initializ:wation
   logical    :: coefficients_set = .false.   !< PLEASE ADD DOCUMENTATION

   ! table contains the coeff spline of the function in boundary 
   ! in the case of dirichlet boundary condition non homogene 

   type(edge), dimension(4) :: side

   logical, dimension(4) :: compute_slope = .TRUE. !< PLEASE ADD DOCUMENTATION
   logical, dimension(4) :: compute_value = .TRUE. !< PLEASE ADD DOCUMENTATION

contains

   !> PLEASE ADD DOCUMENTATION
   procedure :: initialize => initialize_ad2d_interpolator
   !> PLEASE ADD DOCUMENTATION
   procedure, pass(interpolator) :: set_coefficients => set_coefficients_ad2d
   !> PLEASE ADD DOCUMENTATION
   procedure, pass(interpolator) :: coefficients_are_set => coefficients_are_set_ad2d
   !> PLEASE ADD DOCUMENTATION
   procedure :: compute_interpolants        => compute_interpolants_ad2d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_value           => interpolate_value_ad2d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_derivative_eta1 => interpolate_derivative1_ad2d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_derivative_eta2 => interpolate_derivative2_ad2d
   !> PLEASE ADD DOCUMENTATION
   procedure, pass:: interpolate_array      => interpolate_array_ad2d
   !> PLEASE ADD DOCUMENTATION
   procedure, pass:: interpolate_array_disp => interpolate_2d_array_disp_ad2d
   !> PLEASE ADD DOCUMENTATION
   procedure, pass:: get_coefficients       => get_coefficients_ad2d
   !> PLEASE ADD DOCUMENTATION
   procedure, pass:: set_values_at_boundary => set_boundary_value2d
   !> PLEASE ADD DOCUMENTATION
   procedure, pass:: set_slopes_at_boundary => set_slope2d
   !> PLEASE ADD DOCUMENTATION
   procedure, pass:: delete => delete_arbitrary_degree_2d_interpolator

end type sll_arbitrary_degree_spline_interpolator_2d


!> Pointer to arbitrary degree version of 1d interpolator
type, public :: sll_arbitrary_degree_spline_interpolator_2d_ptr
   type(sll_arbitrary_degree_spline_interpolator_2d), pointer :: interp
end type sll_arbitrary_degree_spline_interpolator_2d_ptr


!> Deallocate the interpolator class
interface sll_delete
   module procedure delete_arbitrary_degree_2d_interpolator
end interface sll_delete

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
  sll_interpolator, intent(inout) :: interpolator
  sll_int32 :: ierr
  SLL_DEALLOCATE(interpolator%knots1,ierr)
  SLL_DEALLOCATE(interpolator%knots2,ierr)
  SLL_DEALLOCATE(interpolator%t1,ierr)
  SLL_DEALLOCATE(interpolator%t2,ierr)
  SLL_DEALLOCATE(interpolator%coeff_splines,ierr)
  SLL_DEALLOCATE(interpolator%value,ierr)
  SLL_DEALLOCATE(interpolator%slope,ierr)
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
!> @param[in] bc_w the boundary condition at left in the direction eta1
!> @param[in] bc_e the boundary condition at right in the direction eta2
!> @param[in] bc_s the boundary condition at left in the direction eta2
!> @param[in] bc_n the boundary condition at right in the direction eta2
!> @param[in] spline_degree1 the degree of B-spline in the direction eta1
!> @param[in] spline_degree2 the degre of B-spline in the direction eta2
!> @return the type sll_arbitrary_degree_spline_interpolator_2d

function new_arbitrary_degree_spline_interp2d( &
  num_pts1,                                    &
  num_pts2,                                    &
  eta1_min,                                    & 
  eta1_max,                                    &
  eta2_min,                                    &
  eta2_max,                                    &
  bc_w,                                        &
  bc_e,                                        &
  bc_s,                                        &
  bc_n,                                        &
  spline_degree1,                              &
  spline_degree2) result( res )

  type(sll_arbitrary_degree_spline_interpolator_2d), pointer :: res

  sll_int32,  intent(in) :: num_pts1
  sll_int32,  intent(in) :: num_pts2
  sll_real64, intent(in) :: eta1_min
  sll_real64, intent(in) :: eta1_max
  sll_real64, intent(in) :: eta2_min
  sll_real64, intent(in) :: eta2_max
  sll_int32,  intent(in) :: bc_w
  sll_int32,  intent(in) :: bc_e
  sll_int32,  intent(in) :: bc_s
  sll_int32,  intent(in) :: bc_n
  sll_int32,  intent(in) :: spline_degree1
  sll_int32,  intent(in) :: spline_degree2

  sll_int32 :: ierr
    
  SLL_ALLOCATE(res,ierr)

  call initialize_ad2d_interpolator( &
       res,                          &
       num_pts1,                     &
       num_pts2,                     &
       eta1_min,                     &
       eta1_max,                     &
       eta2_min,                     &
       eta2_max,                     &
       bc_w,                         &
       bc_e,                         &
       bc_s,                         &
       bc_n,                         &
       spline_degree1,               &
       spline_degree2)


end function new_arbitrary_degree_spline_interp2d

! -----------------------------------------------
! This subroutine allocate the type of interpolator
!    the  arbitrary_spline_interp2d
! -----------------------------------------------
!> Initialization of an interpolator arbitrary degree splines 2d.
!> The parameters are
!> @param[in]  num_pts1 the number of points in the direction eta1
!> @param[in]  num_pts2 the number of points in the direction eta2
!> @param[in]  eta1_min the minimun in the direction eta1
!> @param[in]  eta1_max the maximun in the direction eta1
!> @param[in]  eta2_min the minimun in the direction eta2
!> @param[in]  eta2_max the maximun in the direction eta2
!> @param[in]  bc_w  the boundary condition at left in the direction eta1
!> @param[in]  bc_e the boundary condition at right in the direction eta2
!> @param[in]  bc_s the boundary condition at left in the direction eta2
!> @param[in]  bc_n the boundary condition at right in the direction eta2
!> @param[in]  spline_degree1 the degree of B-spline in the direction eta1
!> @param[in]  spline_degree2 the degre of B-spline in the direction eta2
!> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_2d
subroutine initialize_ad2d_interpolator( &
  interpolator,                          &
  num_pts1,                              &
  num_pts2,                              &
  eta1_min,                              &
  eta1_max,                              &
  eta2_min,                              &
  eta2_max,                              &
  bc_w,                                  &
  bc_e,                                  &
  bc_s,                                  &
  bc_n,                                  &
  spline_degree1,                        &
  spline_degree2)

  sll_interpolator:: interpolator

  sll_int32,  intent(in) :: num_pts1
  sll_int32,  intent(in) :: num_pts2
  sll_real64, intent(in) :: eta1_min
  sll_real64, intent(in) :: eta1_max
  sll_real64, intent(in) :: eta2_min
  sll_real64, intent(in) :: eta2_max
  sll_int32,  intent(in) :: bc_w
  sll_int32,  intent(in) :: bc_e
  sll_int32,  intent(in) :: bc_s
  sll_int32,  intent(in) :: bc_n
  sll_int32,  intent(in) :: spline_degree1
  sll_int32,  intent(in) :: spline_degree2

  sll_int32 :: ierr
  sll_int32 :: tmp1
  sll_int32 :: tmp2
  sll_int64 :: bc_selector
   
  if(((bc_w == SLL_PERIODIC).and.(bc_e .ne. SLL_PERIODIC)).or.&
     ((bc_e == SLL_PERIODIC).and.(bc_w .ne. SLL_PERIODIC)))then
     SLL_ERROR("if one boundary condition is specified as periodic, then &
               & both must be. Error in first direction.")
  end if

  if(((bc_s == SLL_PERIODIC).and.(bc_n .ne. SLL_PERIODIC)).or.&
     ((bc_n == SLL_PERIODIC).and.(bc_s .ne. SLL_PERIODIC)))then
     SLL_ERROR("if one boundary condition is specified as periodic, then &
               & both must be. Error in second direction.")
  end if

  bc_selector = 0

  if( bc_w == SLL_DIRICHLET ) bc_selector = bc_selector + 1
  if( bc_w == SLL_NEUMANN   ) bc_selector = bc_selector + 2
  if( bc_w == SLL_HERMITE   ) bc_selector = bc_selector + 4
  if( bc_e == SLL_DIRICHLET ) bc_selector = bc_selector + 8
  if( bc_e == SLL_NEUMANN   ) bc_selector = bc_selector + 16
  if( bc_e == SLL_HERMITE   ) bc_selector = bc_selector + 32
  if( bc_s == SLL_DIRICHLET ) bc_selector = bc_selector + 64
  if( bc_s == SLL_NEUMANN   ) bc_selector = bc_selector + 128
  if( bc_s == SLL_HERMITE   ) bc_selector = bc_selector + 256
  if( bc_n == SLL_DIRICHLET ) bc_selector = bc_selector + 512
  if( bc_n == SLL_NEUMANN   ) bc_selector = bc_selector + 1024
  if( bc_n == SLL_HERMITE   ) bc_selector = bc_selector + 2048

  ! Initialization in the type of interpolator
  interpolator%spline_degree1 = spline_degree1
  interpolator%spline_degree2 = spline_degree2
  interpolator%eta1_min       = eta1_min
  interpolator%eta1_max       = eta1_max
  interpolator%eta2_min       = eta2_min
  interpolator%eta2_max       = eta2_max
  interpolator%side(N)%bc     = bc_n
  interpolator%side(E)%bc     = bc_e
  interpolator%side(S)%bc     = bc_s
  interpolator%side(W)%bc     = bc_w
  interpolator%bc_selector    = bc_selector
  interpolator%num_pts1       = num_pts1
  interpolator%num_pts2       = num_pts2

  SLL_CLEAR_ALLOCATE(interpolator%side(N)%values(1:num_pts1),ierr)
  SLL_CLEAR_ALLOCATE(interpolator%side(N)%slopes(1:num_pts1),ierr)
  SLL_CLEAR_ALLOCATE(interpolator%side(E)%values(1:num_pts2),ierr)
  SLL_CLEAR_ALLOCATE(interpolator%side(E)%slopes(1:num_pts2),ierr)
  SLL_CLEAR_ALLOCATE(interpolator%side(S)%values(1:num_pts1),ierr)
  SLL_CLEAR_ALLOCATE(interpolator%side(S)%slopes(1:num_pts1),ierr)
  SLL_CLEAR_ALLOCATE(interpolator%side(W)%values(1:num_pts2),ierr)
  SLL_CLEAR_ALLOCATE(interpolator%side(W)%slopes(1:num_pts2),ierr)

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

  case (585) ! 4. dirichlet in all sides

    SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
    SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )

  case(650) !left: Neumann, right: Dirichlet, bottom: Neumann, Top: Dirichlet 

    SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
    SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )
       
  case(657) !left: Dirichlet, right: Neumann, bottom: Neumann, Top: Dirichlet 

    SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
    SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )
       
  case(780)  !left: Hermite, right: Dirichlet, bottom: Hermite, Top: Dirichlet

    SLL_ALLOCATE( interpolator%knots1(1:num_pts1+2*spline_degree1),ierr )
    SLL_ALLOCATE( interpolator%knots2(1:num_pts2+2*spline_degree2),ierr )

  case(801)  !left: Dirichlet, right: Hermite, bottom: Hermite, Top: Dirichlet

    SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
    SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )

  case(804)  !left: Hermite, right: Hermite, bottom: Hermite, Top: Dirichlet
       
    SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
    SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )
       
  case(1098)  !left: Neumann, right: Dirichlet, bottom: Dirichlet, Top: Neumann
       
    SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
    SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )

  case(1105)  !left: Dirichlet, right: Neumann, bottom: Dirichlet, Top: Neumann
       
    SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
    SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )

  case(1170)  !left: Neumann, right: Neumann, bottom: Neuman, Top: Neumann
       
    SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
    SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )

  case(2124)  !left: Hermite, right: Dirichlet, bottom: Dirichlet, Top: Hermite

    SLL_ALLOCATE( interpolator%knots1(1:num_pts1+2*spline_degree1),ierr )
    SLL_ALLOCATE( interpolator%knots2(1:num_pts2+2*spline_degree2),ierr )

  case(2145)  !left: Dirichlet, right: Hermite, bottom: Dirichlet, Top: Hermite  

    SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
    SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )

  case(2148)  !left:Hermite , right: Hermite, bottom: Dirichlet, Top: Hermite  

    SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
    SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )

  case(2316)  !left: Hermite, right: Dirichlet, bottom: Hermite, Top: Hermite
       
    SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
    SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )
    
  case(2338)  !left: Dirichlet, right: Hermite, bottom: Hermite, Top: Hermite
       
    SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
    SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )
       
  case(2340) ! Hermite in all sides

    SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
    SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )

  case default
    SLL_WARNING(" BC combination not implemented ")
  end select

  ! tmp1 and tmp2 is the maximun (not absolue) for the size of coefficients
  tmp1 = num_pts1+4*spline_degree1
  tmp2 = num_pts2+4*spline_degree2
  ! knots and coeff splines allocations 
  SLL_CLEAR_ALLOCATE( interpolator%coeff_splines(1:tmp1,1:tmp2),ierr)
  ! the minimun is to be of class C^0 everywhere on the knots
  ! i.e. each knot have multiplicity (spline_degree1+1) 
  ! so the maximun number of knots is num_pts1*(spline_degree1+1)
  SLL_CLEAR_ALLOCATE( interpolator%t1(1:num_pts1*(spline_degree1+1)),ierr)
  SLL_CLEAR_ALLOCATE( interpolator%t2(1:num_pts2*(spline_degree2+1)),ierr) 

end subroutine !initialize_ad2d_interpolator


!> @brief
!> Initialization of the boundary for interpolator arbitrary degree splines 2d.
!> @details
!> @param[in] slope_w a 1d array contains values in the left in the direction eta1  
!> @param[in] slope_e a 1d array contains values in the right in the direction eta1 
!> @param[in] slope_s a 1d array contains values in the left in the direction eta2 
!> @param[in] slope_n a 1d array contains values in the right in the direction eta2
!> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_2d
subroutine set_slope2d( interpolator,  &
                        slope_w,       &
                        slope_e,       &
                        slope_s,       &
                        slope_n)

  sll_interpolator :: interpolator

  sll_real64, dimension(:),optional :: slope_w
  sll_real64, dimension(:),optional :: slope_e
  sll_real64, dimension(:),optional :: slope_s
  sll_real64, dimension(:),optional :: slope_n

  class(sll_arbitrary_degree_spline_interpolator_1d),pointer :: interp_eta1 => null()
  class(sll_arbitrary_degree_spline_interpolator_1d),pointer :: interp_eta2 => null()

  sll_int32 :: sz_slope_s
  sll_int32 :: sz_slope_n
  sll_int32 :: sz_slope_w
  sll_int32 :: sz_slope_e
  sll_int32 :: num_pts1
  sll_int32 :: num_pts2
  sll_int32 :: side

  num_pts1 = interpolator%num_pts1
  num_pts2 = interpolator%num_pts2

  if (present(slope_n) .or. present(slope_s)) then 
    interp_eta1 => new_arbitrary_degree_1d_interpolator( &
                   interpolator%num_pts1,                &
                   interpolator%eta1_min,                &
                   interpolator%eta1_max,                &
                   interpolator%bc(W),                   &
                   interpolator%bc(E),                   &
                   interpolator%spline_degree1 )
  end if

  if (present(slope_e) .or. present(slope_w)) then 
    interp_eta2 => new_arbitrary_degree_1d_interpolator( &
                   interpolator%num_pts2,                &
                   interpolator%eta2_min,                &
                   interpolator%eta2_max,                &
                   interpolator%bc(S),                   &
                   interpolator%bc(N),                   &
                   interpolator%spline_degree2 )
  end if

  interpolator%slope = 0.0
  interpolator%compute_slope = .FALSE.

  if (present(slope_e)) then 
    sz_slope_e = size(slope_e)
    SLL_ASSERT( sz_slope_e == interpolator%num_pts2 ) 
    interpolator%slope(:,E) = slope_e
    interpolator%compute_slope(E)= .FALSE.
  end if

  if (present(slope_w)) then 
    sz_slope_w = size(slope_w)
    SLL_ASSERT( sz_slope_w == interpolator%num_pts2 ) 
    interpolator%slope(:,W) = slope_w
    interpolator%compute_slope(W)= .FALSE.
  end if

  if (present(slope_n)) then 
    sz_slope_n = size(slope_n)
    SLL_ASSERT( sz_slope_n == interpolator%num_pts1 ) 
    interpolator%slope(:,N) = slope_n
    interpolator%compute_slope(N)= .FALSE.
  end if

  if (present(slope_s)) then 
    sz_slope_s = size(slope_s)
    SLL_ASSERT( sz_slope_s == interpolator%num_pts1 ) 
    interpolator%slope(:,S) = slope_s
    interpolator%compute_slope(S)= .FALSE.
  end if

  do side = 1, 4

    if (interpolator%side(i)%bc == SLL_HERMITE) then
      if      (present(slope_n) .and. side == N) then 
        call interp_eta1%compute_interpolants(slope_n)
        interpolator%slope(1:sz_slope_n+2,side) = interp_eta1%coeff_splines(1:sz_slope_n+2)
      else if (present(slope_s) .and. side == S) then 
        call interp_eta1%compute_interpolants(slope_s)
        interpolator%slope(1:sz_slope_s+2,side) = interp_eta1%coeff_splines(1:sz_slope_s+2)
      else if (present(slope_w) .and. side == W) then 
        call interp_eta2%compute_interpolants(slope_w)
        interpolator%slope(1:sz_slope_w,side) = interp_eta2%coeff_splines(1:sz_slope_w)
      else if (present(slope_e) .and. side == E) then 
        call interp_eta2%compute_interpolants(slope_e)
        interpolator%slope(1:sz_slope_e,side) = interp_eta2%coeff_splines(1:sz_slope_e)
      end if
    end if

  end do

  if (present(slope_n) .or. present(slope_s))     &
    call sll_delete(interp_eta1)

  if (present(slope_w) .or. present(slope_e))     &
    call sll_delete(interp_eta2)

    
end subroutine set_slope2d
  
  
!> @brief
!> Initialization of the boundary for interpolator arbitrary degree splines 2d.
!> @details
!> @param[in]  value_w a 1d array contains values in the left direction eta1  
!> @param[in]  value_e a 1d array contains values in the right direction eta1 
!> @param[in]  value_s a 1d array contains values in the left direction eta2 
!> @param[in]  value_n a 1d array contains values in the right direction eta2
!> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_2d

subroutine set_boundary_value2d( this,  &
                                 value_w,    &
                                 value_e,   &
                                 value_s,  &
                                 value_n)

sll_interpolator :: this
sll_real64, optional :: value_w(:)
sll_real64, optional :: value_e(:)
sll_real64, optional :: value_s(:)
sll_real64, optional :: value_n(:)

class(sll_arbitrary_degree_spline_interpolator_1d),pointer :: interp_eta1 => null()
class(sll_arbitrary_degree_spline_interpolator_1d),pointer :: interp_eta2 => null()

sll_int32 :: sz_value_w,sz_value_e,sz_value_s,sz_value_n
sll_int32 :: num_pts1
sll_int32 :: num_pts2
sll_int32 :: side

num_pts1    = this%num_pts1
num_pts2    = this%num_pts2

if (present(value_n) .or. present(value_s))     &
  interp_eta1 => new_arbitrary_degree_1d_interpolator( &
                   this%num_pts1,                      &
                   this%eta1_min,                      &
                   this%eta1_max,                      &
                   this%bc(W),                         &
                   this%bc(E),                         &
                   this%spline_degree1 )

if (present(value_w) .or. present(value_e))     &
  interp_eta2 => new_arbitrary_degree_1d_interpolator( &
                   this%num_pts2,                      &
                   this%eta2_min,                      &
                   this%eta2_max,                      &
                   this%bc(S),                         &
                   this%bc(N),                         &
                   this%spline_degree2 )

this%value = 0.0_f64

do side = 1, 4
       
  if (this%(side) == SLL_DIRICHLET) then

    if (side == W .and. present(value_w)) then 
      sz_value_w = size(value_w)
      SLL_ASSERT(sz_value_w == this%num_pts2)
      call interp_eta2%compute_interpolants(value_w)
      this%side(W)%values = interp_eta2%coeff_splines(1:sz_value_w)
      this%compute_value(side) = .FALSE.
    end if

    if (side == E .and. present(value_e)) then 
      sz_value_e = size(value_e)
      SLL_ASSERT(sz_value_e == this%num_pts2)
      call interp_eta2%compute_interpolants(value_e)
      this%side(E)%values = interp_eta2%coeff_splines(1:sz_value_e)
      this%compute_value(side) = .FALSE.
    end if

    if (side == S .and. present(value_s)) then 
      sz_value_s = size(value_s)
      SLL_ASSERT(sz_value_s == this%num_pts1)
      call interp_eta1%compute_interpolants(value_s)
      this%side(S)%values = interp_eta1%coeff_splines(1:sz_value_s)
      this%compute_value(side) = .FALSE.
    end if
       
    if (side == N .and. present(value_n)) then 
      sz_value_n = size(value_n)
      SLL_ASSERT(sz_value_n == this%num_pts1)
      call interp_eta1%compute_interpolants(value_n)
      this%side(N)%values = interp_eta1%coeff_splines(1:sz_value_n)
      this%compute_value(side) = .FALSE.
    end if

  end if

end do

if (present(value_n) .or. present(value_s)) call sll_delete(interp_eta1)

if (present(value_w) .or. present(value_e)) call sll_delete(interp_eta2)

end subroutine set_boundary_value2d


!> @brief initializing the coefficients of splines.
!> @details  initializing the coefficients of splines
!>  in the cas of linearization of them i.e. if we have 
!>  a table in 1d corresponding of coefficients 2d
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
subroutine set_coefficients_ad2d( &
  interpolator,                   &
  coeffs_1d,                      &
  coeffs_2d,                      &
  coeff2d_size1,                  &
  coeff2d_size2,                  &
  knots1,                         &
  size_knots1,                    &
  knots2,                         &
  size_knots2)

sll_interpolator,           intent(inout)        :: interpolator
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
sll_int32   :: sz_derivative1
sll_int32   :: sz_derivative2

sp_deg1    = interpolator%spline_degree1
sp_deg2    = interpolator%spline_degree2
num_cells1 = interpolator%num_pts1 - 1
num_cells2 = interpolator%num_pts2 - 1
eta1_min   = interpolator%eta1_min
eta2_min   = interpolator%eta2_min
eta1_max   = interpolator%eta1_max
eta2_max   = interpolator%eta2_max
delta1     = (eta1_max - eta1_min)/num_cells1
delta2     = (eta2_max - eta2_min)/num_cells2


if (present(coeffs_1d) ) then 


  ! The interpretation and further filling of the spline coefficients array
  ! depends on the boundary conditions.
  select case (interpolator%bc_selector)

  case(0) ! periodic-periodic
         
    SLL_ASSERT( size(coeffs_1d,1) == num_cells1*num_cells2)

    interpolator%size_coeffs1 =  num_cells1 + sp_deg1 + 1
    interpolator%size_coeffs2 =  num_cells2 + sp_deg2 + 1 
    interpolator%size_t1      =  2*sp_deg1 + num_cells1 +1 +1
    interpolator%size_t2      =  2*sp_deg2 + num_cells2 +1 +1
         
    ! ------------------------------------------------------------
    ! allocation and definition of knots
    ! ------------------------------------------------------------
         
    do i = -sp_deg1, num_cells1 + sp_deg1 + 1
      interpolator%t1( i + sp_deg1 + 1 ) = eta1_min + i*delta1
    end do
         
    do i = -sp_deg2, num_cells2 + sp_deg2 + 1
      interpolator%t2( i + sp_deg2 + 1 ) = eta2_min + i*delta2
    end do
         
    ! ------------------------------------------------------------
    !   reorganization of spline coefficients 1D in coefficients 2D 
    ! ------------------------------------------------------------
         
    do i = 1,num_cells1
      do j = 1,num_cells2
        interpolator%coeff_splines(i,j) = coeffs_1d( i + num_cells1 *(j-1) )
      end do
    end do
         
    do j = 1, sp_deg2 + 1
      do i = 1,num_cells1
        interpolator%coeff_splines(i ,num_cells2+j) = coeffs_1d(i+num_cells1*(j-1))
      end do
    end do

    do i = 1, sp_deg1 + 1
      do j = 1,num_cells2
        interpolator%coeff_splines(num_cells1+i,j) = coeffs_1d(i+num_cells1*(j-1))
      end do
    end do

    do i= 1,sp_deg1 + 1
      do j=1,sp_deg2 + 1
        interpolator%coeff_splines(num_cells1 +  i ,num_cells2 + j) = &
          interpolator%coeff_splines(i,j)
      end do
    end do

  case (9) ! 2. dirichlet-left, dirichlet-right, periodic
         
    interpolator%size_coeffs1 =  num_cells1 + sp_deg1
    interpolator%size_coeffs2 =  num_cells2 + sp_deg2 + 1
    interpolator%size_t1      =  2*sp_deg1 + num_cells1 + 1
    interpolator%size_t2      =  2*sp_deg2 + num_cells2 + 1 + 1
    nb_spline_eta1            =  num_cells1 + sp_deg1 - 2
    nb_spline_eta2            =  num_cells2
         
    SLL_ASSERT(size(coeffs_1d,1) == (num_cells1 + sp_deg1 - 2)*num_cells2)

    ! ------------------------------------------------------------
    ! allocation and definition of knots
    ! ------------------------------------------------------------
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
         
    ! ------------------------------------------------------------
    ! reorganization of spline coefficients 1D in coefficients 2D 
    ! ------------------------------------------------------------
    do i = 1 ,nb_spline_eta1
      do j = 1,nb_spline_eta2
        interpolator%coeff_splines(i+1,j) = coeffs_1d(i+nb_spline_eta1*(j-1))
      end do
    end do
         
    do j = 1, sp_deg2 + 1
      do i = 1,nb_spline_eta1
               
        interpolator%coeff_splines(i + 1 ,nb_spline_eta2 + j ) = &
                    coeffs_1d(i+nb_spline_eta1*(j-1))
      end do
    end do
         
    interpolator%coeff_splines(1,:) = 0.0_8
    interpolator%coeff_splines(nb_spline_eta1+2,:) = 0.0_8
         
  case(576)!3. periodic, dirichlet-bottom, dirichlet-top
       
    interpolator%size_coeffs1 =  num_cells1 + sp_deg1 + 1
    interpolator%size_coeffs2 =  num_cells2 + sp_deg2
    interpolator%size_t1      = 2*sp_deg1 + num_cells1 + 1 + 1
    interpolator%size_t2      = 2*sp_deg2 + num_cells2 + 1
    nb_spline_eta1            = num_cells1
    nb_spline_eta2            = num_cells2 + sp_deg2 - 2
       
    SLL_ASSERT(size(coeffs_1d,1) == num_cells1*( num_cells2 + sp_deg2 - 2))

    ! ------------------------------------------------------------
    ! allocation and definition of knots
    ! ------------------------------------------------------------
    do i = - sp_deg1, nb_spline_eta1 + sp_deg1 + 1
      interpolator%t1( i+ sp_deg1 + 1 ) = eta1_min + i* delta1
    end do
       
    do i = 1, sp_deg2 + 1
      interpolator%t2(i) = eta2_min
    enddo
    eta2 = eta2_min
    do i = sp_deg2 + 2, num_cells2 + 1 + sp_deg2
      eta2 = eta2 + delta2
      interpolator%t2(i) = eta2
    enddo
    do i = num_cells2 + sp_deg2 + 1, num_cells2 + 1 + 2*sp_deg2
      interpolator%t2(i) = eta2_max
    enddo
       
    ! ------------------------------------------------------------
    ! reorganization of spline coefficients 1D in coefficients 2D 
    ! -----------------------------------------------------------
    do i = 1 , nb_spline_eta1
      do j = 1,nb_spline_eta2
        interpolator%coeff_splines(i ,j+1) = &
                  coeffs_1d(i+nb_spline_eta1 *(j-1) )
      end do
    end do
       
    do i = 1, sp_deg1 + 1
      do j = 1,nb_spline_eta2
        interpolator%coeff_splines(nb_spline_eta1 + i ,j+1) = &
                  coeffs_1d(i+nb_spline_eta1 *(j-1) )
             
      end do
    end do
         
    interpolator%coeff_splines(:,1) = 0.0_8
    interpolator%coeff_splines(:,nb_spline_eta2+2) = 0.0_8
    ! ------------------------------------------------------------
       
  case(585) ! 4. dirichlet in all sides

    interpolator%size_coeffs1=  num_cells1 + sp_deg1
    interpolator%size_coeffs2=  num_cells2 + sp_deg2
    interpolator%size_t1 = 2*sp_deg1 + num_cells1 + 1
    interpolator%size_t2 = 2*sp_deg2 + num_cells2 + 1
    nb_spline_eta1 = num_cells1 + sp_deg1 - 2
    nb_spline_eta2 = num_cells2 + sp_deg2 - 2
         
    SLL_ASSERT(size(coeffs_1d,1)==(num_cells1 + sp_deg1-2)*(num_cells2+sp_deg2-2))

    ! ------------------------------------------------------------
    ! allocation and definition of knots
    ! ------------------------------------------------------------
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
       
    ! ------------------------------------------------------------
    ! reorganization of spline coefficients 1D in coefficients 2D 
    ! ------------------------------------------------------------
    ! achtung ! normaly interpolator%slope(:,W) and interpolator%value_e(:)
    ! achtung ! normaly interpolator%value_s(:) and interpolator%value_n(:)

    interpolator%coeff_splines(:,:) = 0.0_8
    ! allocation coefficient spline
    do i = 1,nb_spline_eta1
      do j = 1,nb_spline_eta2
             
        interpolator%coeff_splines(i+1,j+1) = coeffs_1d(i+nb_spline_eta1*(j-1))
      end do
    end do
    
  case (650) !left: Neumann, right: Dirichlet, bottom: Neumann, Top: Dirichlet

    interpolator%size_coeffs1=  num_cells1 + sp_deg1 +1
    interpolator%size_coeffs2=  num_cells2 + sp_deg2 +1
    interpolator%size_t1 = 2*sp_deg1 + num_cells1 + 1
    interpolator%size_t2 = 2*sp_deg2 + num_cells2 + 1
    nb_spline_eta1 = num_cells1 + sp_deg1 +1
    nb_spline_eta2 = num_cells2 + sp_deg2 +1
         
    SLL_ASSERT(size(coeffs_1d,1) == (num_cells1+sp_deg1+1)*(num_cells2+sp_deg2+1))

    ! ------------------------------------------------------------
    ! allocation and definition of knots
    ! ------------------------------------------------------------
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
        
    ! ------------------------------------------------------------
    ! reorganization of spline coefficients 1D in coefficients 2D 
    ! ------------------------------------------------------------
    ! achtung ! normaly interpolator%slope(:,W) and interpolator%value_e(:)
    ! achtung ! normaly interpolator%value_s(:) and interpolator%value_n(:)

    interpolator%coeff_splines(:,:) = 0.0_8
    ! allocation coefficient spline
    do i = 1,nb_spline_eta1
      do j = 1,nb_spline_eta2
        interpolator%coeff_splines(i+1,j+1) = coeffs_1d(i+nb_spline_eta1*(j-1))
      end do
    end do

  case(657) !left: Dirichlet, right: Neumann, bottom: Neumann, Top: Dirichlet 

    interpolator%size_coeffs1=  num_cells1 + sp_deg1 +1
    interpolator%size_coeffs2=  num_cells2 + sp_deg2 +1
    interpolator%size_t1 = 2*sp_deg1 + num_cells1 + 1
    interpolator%size_t2 = 2*sp_deg2 + num_cells2 + 1
    nb_spline_eta1 = num_cells1 + sp_deg1 +1
    nb_spline_eta2 = num_cells2 + sp_deg2 +1
        
    SLL_ASSERT(size(coeffs_1d,1) == (num_cells1+sp_deg1+1)*(num_cells2+sp_deg2+1))

    ! ------------------------------------------------------------
    ! allocation and definition of knots
    ! ------------------------------------------------------------
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
       
    ! ------------------------------------------------------------
    ! reorganization of spline coefficients 1D in coefficients 2D 
    ! ------------------------------------------------------------
    ! achtung ! normaly interpolator%slope(:,W) and interpolator%value_e(:)
    ! achtung ! normaly interpolator%value_s(:) and interpolator%value_n(:)

    interpolator%coeff_splines(:,:) = 0.0_8
    ! allocation coefficient spline
    do i = 1,nb_spline_eta1
       do j = 1,nb_spline_eta2
          
          interpolator%coeff_splines(i+1,j+1) = &
               coeffs_1d( i + nb_spline_eta1 *(j-1))
       end do
    end do
 
  case(780)  !left: Hermite, right: Dirichlet, bottom: Hermite, Top: Dirichlet

    interpolator%size_coeffs1=  num_cells1 + sp_deg1 +1
    interpolator%size_coeffs2=  num_cells2 + sp_deg2 +1
    interpolator%size_t1 = 2*sp_deg1 + num_cells1 + 1
    interpolator%size_t2 = 2*sp_deg2 + num_cells2 + 1
    nb_spline_eta1 = num_cells1 + sp_deg1 +1
    nb_spline_eta2 = num_cells2 + sp_deg2 +1
         
    SLL_ASSERT(size(coeffs_1d,1) == (num_cells1+sp_deg1+1)*(num_cells2+sp_deg2+1))

    ! ------------------------------------------------------------
    ! allocation and definition of knots
    ! ------------------------------------------------------------
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
        
    ! ------------------------------------------------------------
    ! reorganization of spline coefficients 1D in coefficients 2D 
    ! ------------------------------------------------------------
    ! achtung ! normaly interpolator%slope(:,W) and interpolator%value_e(:)
    ! achtung ! normaly interpolator%value_s(:) and interpolator%value_n(:)

    interpolator%coeff_splines(:,:) = 0.0_8
    ! allocation coefficient spline
    do i = 1,nb_spline_eta1
      do j = 1,nb_spline_eta2
        interpolator%coeff_splines(i+1,j+1) = coeffs_1d( i + nb_spline_eta1 *(j-1))
      end do
    end do
  
  case(1098)  !left: Neumann, right: Dirichlet, bottom: Dirichlet, Top: Neumann

    interpolator%size_coeffs1=  num_cells1 + sp_deg1 +1
    interpolator%size_coeffs2=  num_cells2 + sp_deg2 +1
    interpolator%size_t1 = 2*sp_deg1 + num_cells1 + 1
    interpolator%size_t2 = 2*sp_deg2 + num_cells2 + 1
    nb_spline_eta1 = num_cells1 + sp_deg1 +1
    nb_spline_eta2 = num_cells2 + sp_deg2 +1
         
    SLL_ASSERT(size(coeffs_1d,1) == (num_cells1+sp_deg1+1)*(num_cells2+sp_deg2+1))

    ! ------------------------------------------------------------
    ! allocation and definition of knots
    ! ------------------------------------------------------------
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
    
    ! ------------------------------------------------------------
    ! reorganization of spline coefficients 1D in coefficients 2D 
    ! ------------------------------------------------------------
    ! achtung ! normaly interpolator%slope(:,W) and interpolator%value_e(:)
    ! achtung ! normaly interpolator%value_s(:) and interpolator%value_n(:)

    interpolator%coeff_splines(:,:) = 0.0_8
    ! allocation coefficient spline
    do i = 1,nb_spline_eta1
       do j = 1,nb_spline_eta2
          interpolator%coeff_splines(i+1,j+1) = coeffs_1d( i + nb_spline_eta1 *(j-1))
       end do
    end do

  case(2145)  !left: Dirichlet, right: Hermite, bottom: Dirichlet, Top: Hermite  

    sz_derivative1 = 2
    sz_derivative2 = 2
    interpolator%size_coeffs1=  num_cells1 + sp_deg1 + 1
    interpolator%size_coeffs2=  num_cells2 + sp_deg2 + 1
    interpolator%size_t1 = 2*sp_deg1 + num_cells1 + 1
    interpolator%size_t2 = 2*sp_deg2 + num_cells2 + 1
    nb_spline_eta1 = num_cells1 + sp_deg1 + 1
    nb_spline_eta2 = num_cells2 + sp_deg2 + 1
    
    SLL_ASSERT(size(coeffs_1d,1) == (num_cells1+sp_deg1+1)*(num_cells2+sp_deg2+1))

    ! ------------------------------------------------------------
    ! allocation and definition of knots
    ! ------------------------------------------------------------
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
         
    ! ------------------------------------------------------------
    ! reorganization of spline coefficients 1D in coefficients 2D 
    ! ------------------------------------------------------------
    ! achtung ! normaly interpolator%slope(:,W) and interpolator%value_e(:)
    ! achtung ! normaly interpolator%value_s(:) and interpolator%value_n(:)

    interpolator%coeff_splines(:,:) = 0.0_8
    ! allocation coefficient spline
    do i = 1,nb_spline_eta1
       do j = 1,nb_spline_eta2
          
          interpolator%coeff_splines(i+1,j+1) = &
               coeffs_1d( i + nb_spline_eta1 *(j-1))
       end do
    end do
    ! ------------------------------------------------------------
    
  case(2124)  !left: Hermite, right: Dirichlet, bottom: Dirichlet, Top: Hermite

    interpolator%size_coeffs1=  num_cells1 + sp_deg1 +1
    interpolator%size_coeffs2=  num_cells2 + sp_deg2 +1
    interpolator%size_t1 = 2*sp_deg1 + num_cells1 + 1
    interpolator%size_t2 = 2*sp_deg2 + num_cells2 + 1
    nb_spline_eta1 = num_cells1 + sp_deg1 +1
    nb_spline_eta2 = num_cells2 + sp_deg2 +1
    
    SLL_ASSERT(size(coeffs_1d,1).ne.(num_cells1+sp_deg1+1)*(num_cells2+sp_deg2+1))

    ! ------------------------------------------------------------
    ! allocation and definition of knots
    ! ------------------------------------------------------------
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
    
    ! ------------------------------------------------------------
    ! reorganization of spline coefficients 1D in coefficients 2D 
    ! ------------------------------------------------------------
    ! achtung ! normaly interpolator%slope(:,W) and interpolator%value_e(:)
    ! achtung ! normaly interpolator%value_s(:) and interpolator%value_n(:)

    interpolator%coeff_splines(:,:) = 0.0_8
    ! allocation coefficient spline
    do i = 1,nb_spline_eta1
      do j = 1,nb_spline_eta2
        interpolator%coeff_splines(i+1,j+1) = coeffs_1d(i+nb_spline_eta1*(j-1))
      end do
    end do
    ! ------------------------------------------------------------
 
  case(801)  !left: Dirichlet, right: Hermite, bottom: Hermite, Top: Dirichlet

    interpolator%size_coeffs1=  num_cells1 + sp_deg1+1
    interpolator%size_coeffs2=  num_cells2 + sp_deg2+1
    interpolator%size_t1 = 2*sp_deg1 + num_cells1 + 1
    interpolator%size_t2 = 2*sp_deg2 + num_cells2 + 1
    nb_spline_eta1 = num_cells1 + sp_deg1 +1
    nb_spline_eta2 = num_cells2 + sp_deg2 +1
       
    SLL_ASSERT(size(coeffs_1d,1) == (num_cells1+sp_deg1+1)*(num_cells2+sp_deg2+1))

    ! ------------------------------------------------------------
    ! allocation and definition of knots
    ! ------------------------------------------------------------
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
    
    ! ------------------------------------------------------------
    ! reorganization of spline coefficients 1D in coefficients 2D 
    ! ------------------------------------------------------------
    ! achtung ! normaly interpolator%slope(:,W) and interpolator%value_e(:)
    ! achtung ! normaly interpolator%value_s(:) and interpolator%value_n(:)

    interpolator%coeff_splines(:,:) = 0.0_8
    ! allocation coefficient spline
    do i = 1,nb_spline_eta1
      do j = 1,nb_spline_eta2
        interpolator%coeff_splines(i+1,j+1) = coeffs_1d( i + nb_spline_eta1 *(j-1))
      end do
    end do
    
  case(2340) ! Hermite in al sides
      
    interpolator%size_coeffs1=  num_cells1 + sp_deg1
    interpolator%size_coeffs2=  num_cells2 + sp_deg2
    interpolator%size_t1 = 2*sp_deg1 + num_cells1 + 1
    interpolator%size_t2 = 2*sp_deg2 + num_cells2 + 1
    nb_spline_eta1 = num_cells1 + sp_deg1
    nb_spline_eta2 = num_cells2 + sp_deg2
    
    SLL_ASSERT(size(coeffs_1d,1) == (num_cells1 + sp_deg1)*(num_cells2+sp_deg2))
    ! ------------------------------------------------------------
    ! allocation and definition of knots
    ! ------------------------------------------------------------
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
       
    ! ------------------------------------------------------------
    ! reorganization of spline coefficients 1D in coefficients 2D 
    ! ------------------------------------------------------------
    ! achtung ! normaly interpolator%slope(:,W) and interpolator%value_e(:)
    ! achtung ! normaly interpolator%value_s(:) and interpolator%value_n(:)

    interpolator%coeff_splines(:,:) = 0.0_8
    ! allocation coefficient spline
    do i = 1,nb_spline_eta1
      do j = 1,nb_spline_eta2
        interpolator%coeff_splines(i,j) = coeffs_1d(i+nb_spline_eta1*(j-1))
      end do
    end do
    
       
  case default

    SLL_ERROR("set_spline_coefficients not recognized.")

  end select

else if (present(coeffs_2d) ) then 

  if ( present(coeff2d_size1) .and. present(coeff2d_size2)) then

    interpolator%size_coeffs1 = coeff2d_size1
    interpolator%size_coeffs2 = coeff2d_size2
    interpolator%size_t1      = sp_deg1 + coeff2d_size1 +1 
    interpolator%size_t2      = sp_deg2 + coeff2d_size2 +1
       
    SLL_ASSERT( coeff2d_size1 <= num_cells1+1+4*sp_deg1) 
    SLL_ASSERT( coeff2d_size2 <= num_cells2+1+4*sp_deg2) 
         
    interpolator%coeff_splines(1:coeff2d_size1,1:coeff2d_size2) = &
      coeffs_2d(1:coeff2d_size1,1:coeff2d_size2)

         
    if ( present(knots1) .and. present(knots2) ) then 
            
      SLL_ASSERT( size_knots1 == coeff2d_size1+sp_deg1+1)
      SLL_ASSERT( size_knots2 == coeff2d_size2+sp_deg2+1)
      SLL_ASSERT( size_knots1 <= (num_cells1+1)*(sp_deg1+1)) 
      SLL_ASSERT( size_knots2 <= (num_cells2+1)*(sp_deg2+1))
            
      interpolator%t1(1:interpolator%size_t1) = knots1(1:interpolator%size_t1)
      interpolator%t2(1:interpolator%size_t2) = knots2(1:interpolator%size_t2)
            
    else if ( (.not. present(knots1)).and.(.not. present(knots2))) then
            
      SLL_ASSERT(interpolator%size_t1 <= (num_cells1+1)*(sp_deg1+1)) 
      SLL_ASSERT(interpolator%size_t2 <= (num_cells2+1)*(sp_deg2+1))

      interpolator%t1 ( 1 : sp_deg1 + 1 )  = eta1_min
      interpolator%t1 ( coeff2d_size1 + 2: coeff2d_size1 + 2 + sp_deg1) = eta1_max
            
      do i = 1, coeff2d_size1 -sp_deg1
        interpolator%t1 ( i + sp_deg1 + 1 ) = eta1_min + &
          i * (eta1_max - eta1_min) / (coeff2d_size1-sp_deg1 + 1)   
      end do
            
      interpolator%t2 ( 1 : sp_deg2 + 1 )  = eta2_min
      interpolator%t2 ( coeff2d_size2 + 2: coeff2d_size2 + 2 + sp_deg2) = eta2_max
            
      do i = 1, coeff2d_size2 -sp_deg2
        interpolator%t2 ( i + sp_deg2 + 1 ) = eta2_min + &
          i * (eta2_max - eta2_min) / (coeff2d_size2-sp_deg2 + 1)   
      end do
            
    else 

      SLL_ERROR(" Knots1 or Knots2 is not present ")
          
    end if
         
  else 

    SLL_ERROR("the number of coefficients must be specified")
  end if
      
else 

  SLL_ERROR('Problem in set_coefficients: must be have coefficients')

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

subroutine compute_interpolants_ad2d( &
  interpolator,                       &
  data_array,                         &
  eta1_coords,                        &
  size_eta1_coords,                   &
  eta2_coords,                        &
  size_eta2_coords )

  sll_interpolator, intent(inout)  :: interpolator

  sll_real64, dimension(:,:), intent(in)         :: data_array
  sll_real64, dimension(:), intent(in),optional  :: eta1_coords
  sll_real64, dimension(:), intent(in),optional  :: eta2_coords
  sll_int32, intent(in),optional                 :: size_eta1_coords
  sll_int32, intent(in),optional                 :: size_eta2_coords

  sll_real64, dimension(:),pointer               :: point_location_eta1
  sll_real64, dimension(:),pointer               :: point_location_eta2
  sll_real64, dimension(:),pointer               :: point_location_eta1_tmp
  sll_real64, dimension(:),pointer               :: point_location_eta2_tmp
  sll_real64, dimension(:,:),pointer             :: data_array_tmp
  sll_real64, dimension(:,:),pointer             :: data_array_deriv_eta1
  sll_real64, dimension(:,:),pointer             :: data_array_deriv_eta2

  sll_int32, pointer :: point_location_eta1_deriv(:)
  sll_int32, pointer :: point_location_eta2_deriv(:)

  sll_int32  :: sz_derivative_eta1,sz_derivative_eta2
  sll_real64 :: delta_eta1
  sll_real64 :: delta_eta2
  sll_int32  :: sz1
  sll_int32  :: sz2
  sll_real64 :: period1
  sll_real64 :: period2
  sll_int32  :: order1
  sll_int32  :: order2
  sll_int32  :: ierr
  sll_int32  :: i
  logical    :: user_coords

  user_coords = .false.

  if(present(eta1_coords) .and. (.not. present(size_eta1_coords))) then
    SLL_ERROR("eta1_coords size must be specifieda")
  end if
    
  if(present(eta2_coords) .and. (.not. present(size_eta2_coords))) then
    SLL_ERROR("eta2_coords size must be specifieda")
  end if
    
  if ((present(eta1_coords) .and. (.not. present(eta2_coords))) .or.&
      (present(eta2_coords) .and. (.not. present(eta1_coords))) ) then
    SLL_ERROR('eta1_coords or eta2_coords must be specified')
  end if
    
  if( present(eta1_coords) .and. present(eta2_coords) ) then
    user_coords = .true.
  end if
    
  if (user_coords) then

    sz1 = size_eta1_coords
    sz2 = size_eta2_coords
       
    SLL_ALLOCATE(point_location_eta1(1:sz1),ierr)
    SLL_ALLOCATE(point_location_eta2(1:sz2),ierr)
    point_location_eta1 = eta1_coords
    point_location_eta2 = eta2_coords

  else ! size depends on BC combination, filled out at initialization.

    sz1 = interpolator%num_pts1
    sz2 = interpolator%num_pts2

    delta_eta1 = (interpolator%eta1_max - interpolator%eta1_min)&
            /(interpolator%num_pts1 -1)
    delta_eta2 = (interpolator%eta2_max - interpolator%eta2_min)&
            /(interpolator%num_pts2 -1)
    SLL_ALLOCATE(point_location_eta1(1:sz1),ierr)
    SLL_ALLOCATE(point_location_eta2(1:sz2),ierr)
      
    do i = 1,sz1
      point_location_eta1(i) = interpolator%eta1_min + delta_eta1*(i-1)
    end do
    do i = 1,sz2
      point_location_eta2(i) = interpolator%eta2_min + delta_eta2*(i-1)
    end do

  end if

  SLL_ALLOCATE(point_location_eta1_tmp(1:sz1-1),ierr)
  SLL_ALLOCATE(point_location_eta2_tmp(1:sz2-1),ierr)
  point_location_eta1_tmp = point_location_eta1(1:sz1-1)
  point_location_eta2_tmp = point_location_eta2(1:sz2-1)
    
  SLL_ASSERT(sz1 <= interpolator%num_pts1 + 8*interpolator%spline_degree1)
  SLL_ASSERT(sz2 <= interpolator%num_pts2 + 8*interpolator%spline_degree1)
  SLL_ASSERT(size(data_array,1) >= sz1)
  SLL_ASSERT(size(data_array,2) >= sz2)
  SLL_ASSERT(size(point_location_eta1) >= sz1)
  SLL_ASSERT(size(point_location_eta2) >= sz2)
    
  order1  = interpolator%spline_degree1 + 1
  order2  = interpolator%spline_degree2 + 1
  period1 = interpolator%eta1_max - interpolator%eta1_min
  period2 = interpolator%eta2_max - interpolator%eta2_min
    
  ! we compute the coefficients spline associate to the values 
  ! data_array and we compute also the knots t1 and t2 using to 
  ! construct the spline to have a good interpolation
    
  SLL_ALLOCATE(point_location_eta1_deriv(2),ierr)
  SLL_ALLOCATE(point_location_eta2_deriv(2),ierr)

  select case (interpolator%bc_selector)
  case (0) ! periodic-periodic
    interpolator%size_coeffs1 = sz1!+1
    interpolator%size_coeffs2 = sz2!+1
    interpolator%size_t1 = order1 + sz1 !+ 1
    interpolator%size_t2 = order2 + sz2 !+ 1 

    !  data_array must have the same dimension than 
    !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
    !  i.e  data_array must have the dimension sz1 x sz2
    SLL_ALLOCATE( data_array_tmp(1:sz1-1,1:sz2-1),ierr)
 
    data_array_tmp = data_array(1:sz1-1,1:sz2-1)
    if ( .not. associated(point_location_eta1_tmp)) &
      SLL_ALLOCATE(point_location_eta1_tmp(sz1-1),ierr)
    if ( .not. associated(point_location_eta2_tmp)) &
      SLL_ALLOCATE(point_location_eta2_tmp(sz2-1),ierr)
    call spli2d_perper( &
      period1, sz1, order1, point_location_eta1_tmp,&!(1:sz1-1), & !+1
      period2, sz2, order2, point_location_eta2_tmp,&!(1:sz2-1), & !+1
      data_array_tmp, interpolator%coeff_splines,&!(1:sz1,1:sz2),&
      interpolator%t1,&!(1:order1 + sz1 ), &!+ 1), &
      interpolator%t2)!(1:order2 + sz2 ))!+ 1) )
       
  case (9) ! 2. dirichlet-left, dirichlet-right, periodic

    interpolator%size_coeffs1 = sz1
    interpolator%size_coeffs2 = sz2!+1
    interpolator%size_t1 = order1 + sz1
    interpolator%size_t2 = order2 + sz2 !+ 1
       
       !print*, size(data_array(1:sz1,1:sz2-1),1)
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2

       SLL_ALLOCATE( data_array_tmp(1:sz1,1:sz2-1),ierr)
       data_array_tmp = data_array(1:sz1,1:sz2-1)
       call spli2d_dirper( sz1, order1, point_location_eta1,&!(1:sz1), &
            period2, sz2, order2, point_location_eta2_tmp,&!(1:sz2-1), & !+1
            data_array_tmp, interpolator%coeff_splines,&!(1:sz1,1:sz2),&!+1
            interpolator%t1,&!(1:sz1+order1), &
            interpolator%t2)!(1:sz2+order2) ) !+1


      ! print*, 'oulala'
       ! boundary condition non homogene  a revoir !!!!! 
       !print*,'zarrrr', interpolator%value_w(1:sz2)
       interpolator%coeff_splines(1,1:sz2)   = data_array(1,1:sz2)!interpolator%value_w(1:sz2)
       interpolator%coeff_splines(sz1,1:sz2) = data_array(sz1,1:sz2)!interpolator%value_e(1:sz2)
  
    case(576) !  3. periodic, dirichlet-bottom, dirichlet-top
       interpolator%size_coeffs1 = sz1!+1
       interpolator%size_coeffs2 = sz2
       interpolator%size_t1 = order1 + sz1 !+ 1
       interpolator%size_t2 = order2 + sz2 
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       SLL_ALLOCATE( data_array_tmp(1:sz1-1,1:sz2),ierr)
       data_array_tmp = data_array(1:sz1-1,1:sz2)
       call spli2d_perdir( period1, sz1, order1, point_location_eta1_tmp,&!(1:sz1-1), & !+ 1
            sz2, order2, point_location_eta2, &
            data_array_tmp, interpolator%coeff_splines,&!(1:sz1,1:sz2),& !+ 1
            interpolator%t1,&!(1:sz1+order1), & ! + 1
            interpolator%t2)!)(1:sz2+order2) )

       ! boundary condition non homogene
       interpolator%coeff_splines(1:sz1,1)   = data_array(1:sz1,1)
       interpolator%coeff_splines(1:sz1,sz2) = data_array(1:sz1,sz2)
       
    case (585) ! 4. dirichlet in all sides
       !print*, 'her'
       interpolator%size_coeffs1 = sz1
       interpolator%size_coeffs2 = sz2
       interpolator%size_t1 = order1 + sz1 
       interpolator%size_t2 = order2 + sz2 
       
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       SLL_ALLOCATE( data_array_tmp(1:sz1,1:sz2),ierr)
       data_array_tmp = data_array(1:sz1,1:sz2)
       call spli2d_custom( sz1, order1, point_location_eta1, &
            sz2, order2, point_location_eta2, &
            data_array_tmp, interpolator%coeff_splines,&!(1:sz1,1:sz2),&
            interpolator%t1,&!(1:sz1+order1), &
            interpolator%t2)!(1:sz2+order2) )

       ! boundary condition non homogene
       interpolator%coeff_splines(1,1:sz2)   = data_array(1,1:sz2)
       interpolator%coeff_splines(sz1,1:sz2) = data_array(sz1,1:sz2)
       ! boundary condition non homogene
       interpolator%coeff_splines(1:sz1,1)   = data_array(1:sz1,1)
       interpolator%coeff_splines(1:sz1,sz2) = data_array(1:sz1,sz2)

    case (650) !left: Neumann, right: Dirichlet, bottom: Neumann, Top: Dirichlet
       sz_derivative_eta1 = 2
       sz_derivative_eta2 = 2
       interpolator%size_coeffs1 = sz1 + sz_derivative_eta1
       interpolator%size_coeffs2 = sz2 + sz_derivative_eta2
       interpolator%size_t1 = order1 + sz1 + sz_derivative_eta1
       interpolator%size_t2 = order2 + sz2 + sz_derivative_eta2
       
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       SLL_ALLOCATE( data_array_tmp(1:sz1,1:sz2),ierr)
       SLL_CLEAR_ALLOCATE( data_array_deriv_eta1(1:2,1:sz2),ierr)
       SLL_CLEAR_ALLOCATE( data_array_deriv_eta2(1:sz_derivative_eta2,1:sz1+sz_derivative_eta1),ierr)
       data_array_tmp = data_array(1:sz1,1:sz2)
       point_location_eta1_deriv(1) = 1
       point_location_eta1_deriv(2) = sz1
       data_array_deriv_eta1(1,1:sz2)     = 0.0_f64
       data_array_deriv_eta1(2,1:sz2)     = interpolator%slope(1:sz2,E)
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=0.0_f64
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)=interpolator%slope(1:sz1+sz_derivative_eta1,N)
       call spli2d_custom_derder(&
            sz1,&
            sz_derivative_eta1,&
            order1, &
            point_location_eta1, &
            point_location_eta1_deriv,&
            sz2, &
            sz_derivative_eta2,&
            order2, point_location_eta2, &
            point_location_eta2_deriv,&
            data_array_tmp,&
            data_array_deriv_eta1,&
            data_array_deriv_eta2,&
            interpolator%coeff_splines,&!(1:sz1,1:sz2),&
            interpolator%t1,&!(1:sz1+order1), &
            interpolator%t2)!(1:sz2+order2) )

       SLL_DEALLOCATE( data_array_deriv_eta1,ierr)
       SLL_DEALLOCATE( data_array_deriv_eta2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_w(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1:sz1+sz_derivative_eta1,1)   = interpolator%value_s(1:sz1+sz_derivative_eta1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_n(1:sz1)

    case(657) !left: Dirichlet, right: Neumann, bottom: Neumann, Top: Dirichlet 
       sz_derivative_eta1 = 2
       sz_derivative_eta2 = 2
       interpolator%size_coeffs1 = sz1 + sz_derivative_eta1
       interpolator%size_coeffs2 = sz2 + sz_derivative_eta2
       interpolator%size_t1 = order1 + sz1 + sz_derivative_eta1
       interpolator%size_t2 = order2 + sz2 + sz_derivative_eta2
       
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       SLL_ALLOCATE( data_array_tmp(1:sz1,1:sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta1(sz_derivative_eta1,1:sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta2(sz_derivative_eta2,1:sz1+sz_derivative_eta1),ierr)
       data_array_tmp = data_array(1:sz1,1:sz2)
       point_location_eta1_deriv(1) = 1
       point_location_eta1_deriv(2) = sz1
       data_array_deriv_eta1(1,1:sz2)     = interpolator%slope(1:sz2,W) 
       data_array_deriv_eta1(2,1:sz2)     = 0.0_f64
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)= 0.0_f64
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)=interpolator%slope(1:sz1+sz_derivative_eta1,N)
       call spli2d_custom_derder(&
            sz1,&
            sz_derivative_eta1,&
            order1, &
            point_location_eta1, &
            point_location_eta1_deriv,&
            sz2, &
            sz_derivative_eta2,&
            order2, point_location_eta2, &
            point_location_eta2_deriv,&
            data_array_tmp,&
            data_array_deriv_eta1,&
            data_array_deriv_eta2,&
            interpolator%coeff_splines,&!(1:sz1,1:sz2),&
            interpolator%t1,&!(1:sz1+order1), &
            interpolator%t2)!(1:sz2+order2) )

       SLL_DEALLOCATE( data_array_deriv_eta1,ierr)
       SLL_DEALLOCATE( data_array_deriv_eta2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_w(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:sz1,1)   = interpolator%value_s(1:sz1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_n(1:sz1)


    case(780)  !left: Hermite, right: Dirichlet, bottom: Hermite, Top: Dirichlet

       sz_derivative_eta1 = 2
       sz_derivative_eta2 = 2
       interpolator%size_coeffs1 = sz1 + sz_derivative_eta1
       interpolator%size_coeffs2 = sz2 + sz_derivative_eta2
       interpolator%size_t1 = order1 + sz1 + sz_derivative_eta1
       interpolator%size_t2 = order2 + sz2 + sz_derivative_eta2
       
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       SLL_ALLOCATE(data_array_tmp(1:sz1,1:sz2),ierr)
       SLL_CLEAR_ALLOCATE(data_array_deriv_eta1(1:sz_derivative_eta1,1:sz2),ierr)
       SLL_CLEAR_ALLOCATE(data_array_deriv_eta2(1:sz_derivative_eta2,1:sz1+sz_derivative_eta1),ierr)
       data_array_tmp = data_array(1:sz1,1:sz2)
       point_location_eta1_deriv(1) = 1
       point_location_eta1_deriv(2) = sz1
       data_array_deriv_eta1(1,1:sz2) = interpolator%slope(1:sz2,W)
       data_array_deriv_eta1(2,1:sz2) = interpolator%slope(1:sz2,E)
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)= &
          interpolator%slope(1:sz1+sz_derivative_eta1,S)
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)= &
          interpolator%slope(1:sz1+sz_derivative_eta1,N)

       call spli2d_custom_derder(&
            sz1,&
            sz_derivative_eta1,&
            order1, &
            point_location_eta1, &
            point_location_eta1_deriv,&
            sz2, &
            sz_derivative_eta2,&
            order2, point_location_eta2, &
            point_location_eta2_deriv,&
            data_array_tmp,&
            data_array_deriv_eta1,&
            data_array_deriv_eta2,&
            interpolator%coeff_splines,&!(1:sz1,1:sz2),&
            interpolator%t1,&!(1:sz1+order1), &
            interpolator%t2)!(1:sz2+order2) )


       SLL_DEALLOCATE( data_array_deriv_eta1,ierr)
       SLL_DEALLOCATE( data_array_deriv_eta2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_w(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:sz1,1)   = interpolator%value_s(1:sz1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_n(1:sz1)


    case(801)  !left: Dirichlet, right: Hermite, bottom: Hermite, Top: Dirichlet


       sz_derivative_eta1 = 2
       sz_derivative_eta2 = 2
       interpolator%size_coeffs1 = sz1 + sz_derivative_eta1
       interpolator%size_coeffs2 = sz2 + sz_derivative_eta2
       interpolator%size_t1 = order1 + sz1 + sz_derivative_eta1
       interpolator%size_t2 = order2 + sz2 + sz_derivative_eta2
       
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       SLL_ALLOCATE( data_array_tmp(1:sz1,1:sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta1(sz_derivative_eta1,1:sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta2(sz_derivative_eta2,1:sz1+sz_derivative_eta1),ierr)
       data_array_tmp = data_array(1:sz1,1:sz2)
       point_location_eta1_deriv(1) = 1
       point_location_eta1_deriv(2) = sz1
       data_array_deriv_eta1(1,1:sz2)     = interpolator%slope(1:sz2,W) 
       data_array_deriv_eta1(2,1:sz2)     = interpolator%slope(1:sz2,E) 
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=interpolator%slope(1:sz1+sz_derivative_eta1,S)
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)=interpolator%slope(1:sz1+sz_derivative_eta1,N)
       call spli2d_custom_derder(&
            sz1,&
            sz_derivative_eta1,&
            order1, &
            point_location_eta1, &
            point_location_eta1_deriv,&
            sz2, &
            sz_derivative_eta2,&
            order2, point_location_eta2, &
            point_location_eta2_deriv,&
            data_array_tmp,&
            data_array_deriv_eta1,&
            data_array_deriv_eta2,&
            interpolator%coeff_splines,&!(1:sz1,1:sz2),&
            interpolator%t1,&!(1:sz1+order1), &
            interpolator%t2)!(1:sz2+order2) )

       SLL_DEALLOCATE( data_array_deriv_eta1,ierr)
       SLL_DEALLOCATE( data_array_deriv_eta2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_w(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:sz1,1)   = interpolator%value_s(1:sz1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_n(1:sz1)

    case(804)  !left: Hermite, right: Hermite, bottom: Hermite, Top: Dirichlet
       sz_derivative_eta1 = 2
       sz_derivative_eta2 = 2
       interpolator%size_coeffs1 = sz1 + sz_derivative_eta1
       interpolator%size_coeffs2 = sz2 + sz_derivative_eta2
       interpolator%size_t1 = order1 + sz1 + sz_derivative_eta1
       interpolator%size_t2 = order2 + sz2 + sz_derivative_eta2
       
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       SLL_ALLOCATE( data_array_tmp(1:sz1,1:sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta1(sz_derivative_eta1,1:sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta2(sz_derivative_eta2,1:sz1+sz_derivative_eta1),ierr)
       data_array_tmp = data_array(1:sz1,1:sz2)
       point_location_eta1_deriv(1) = 1
       point_location_eta1_deriv(2) = sz1
       data_array_deriv_eta1(1,1:sz2)     = interpolator%slope(1:sz2,W) 
       data_array_deriv_eta1(2,1:sz2)     = interpolator%slope(1:sz2,E) 
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=interpolator%slope(1:sz1+sz_derivative_eta1,S)
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)=interpolator%slope(1:sz1+sz_derivative_eta1,N)
       call spli2d_custom_derder(&
            sz1,&
            sz_derivative_eta1,&
            order1, &
            point_location_eta1, &
            point_location_eta1_deriv,&
            sz2, &
            sz_derivative_eta2,&
            order2, point_location_eta2, &
            point_location_eta2_deriv,&
            data_array_tmp,&
            data_array_deriv_eta1,&
            data_array_deriv_eta2,&
            interpolator%coeff_splines,&!(1:sz1,1:sz2),&
            interpolator%t1,&!(1:sz1+order1), &
            interpolator%t2)!(1:sz2+order2) )

       SLL_DEALLOCATE( data_array_deriv_eta1,ierr)
       SLL_DEALLOCATE( data_array_deriv_eta2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_w(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:sz1,1)   = interpolator%value_s(1:sz1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_n(1:sz1)

    case(1098)  !left: Neumann, right: Dirichlet, bottom: Dirichlet, Top: Neumann
       sz_derivative_eta1 = 2
       sz_derivative_eta2 = 2
       interpolator%size_coeffs1 = sz1 + sz_derivative_eta1
       interpolator%size_coeffs2 = sz2 + sz_derivative_eta2
       interpolator%size_t1 = order1 + sz1 + sz_derivative_eta1
       interpolator%size_t2 = order2 + sz2 + sz_derivative_eta2
       
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       SLL_ALLOCATE( data_array_tmp(1:sz1,1:sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta1(sz_derivative_eta1,1:sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta2(sz_derivative_eta2,1:sz1+sz_derivative_eta1),ierr)
       data_array_tmp = data_array(1:sz1,1:sz2)
       point_location_eta1_deriv(1) = 1
       point_location_eta1_deriv(2) = sz1
       data_array_deriv_eta1(1,1:sz2)     = 0.0_f64
       data_array_deriv_eta1(2,1:sz2)     = interpolator%slope(1:sz2,E)
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=interpolator%slope(1:sz1+sz_derivative_eta1,S)
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)= 0.0_f64
       call spli2d_custom_derder(&
            sz1,&
            sz_derivative_eta1,&
            order1, &
            point_location_eta1, &
            point_location_eta1_deriv,&
            sz2, &
            sz_derivative_eta2,&
            order2, point_location_eta2, &
            point_location_eta2_deriv,&
            data_array_tmp,&
            data_array_deriv_eta1,&
            data_array_deriv_eta2,&
            interpolator%coeff_splines,&!(1:sz1,1:sz2),&
            interpolator%t1,&!(1:sz1+order1), &
            interpolator%t2)!(1:sz2+order2) )

       SLL_DEALLOCATE( data_array_deriv_eta1,ierr)
       SLL_DEALLOCATE( data_array_deriv_eta2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_w(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:sz1,1)   = interpolator%value_s(1:sz1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_n(1:sz1)

    case(1105)  !left: Dirichlet, right: Neumann, bottom: Dirichlet, Top: Neumann
       sz_derivative_eta1 = 2
       sz_derivative_eta2 = 2
       interpolator%size_coeffs1 = sz1 + sz_derivative_eta1
       interpolator%size_coeffs2 = sz2 + sz_derivative_eta2
       interpolator%size_t1 = order1 + sz1 + sz_derivative_eta1
       interpolator%size_t2 = order2 + sz2 + sz_derivative_eta2
       
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       SLL_ALLOCATE( data_array_tmp(1:sz1,1:sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta1(2,sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta2(sz_derivative_eta2,sz1+sz_derivative_eta1),ierr)
       data_array_tmp = data_array(1:sz1,1:sz2)
       point_location_eta1_deriv(1) = 1
       point_location_eta1_deriv(2) = sz1
       data_array_deriv_eta1(1,1:sz2)     = interpolator%slope(1:sz2,W)
       data_array_deriv_eta1(2,1:sz2)     = 0.0_f64
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=interpolator%slope(1:sz1+sz_derivative_eta1,S)
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)= 0.0_f64
       call spli2d_custom_derder(&
            sz1,&
            sz_derivative_eta1,&
            order1, &
            point_location_eta1, &
            point_location_eta1_deriv,&
            sz2, &
            sz_derivative_eta2,&
            order2, point_location_eta2, &
            point_location_eta2_deriv,&
            data_array_tmp,&
            data_array_deriv_eta1,&
            data_array_deriv_eta2,&
            interpolator%coeff_splines,&!(1:sz1,1:sz2),&
            interpolator%t1,&!(1:sz1+order1), &
            interpolator%t2)!(1:sz2+order2) )

       SLL_DEALLOCATE( data_array_deriv_eta1,ierr)
       SLL_DEALLOCATE( data_array_deriv_eta2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_w(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1:sz1+sz_derivative_eta1,1)   = interpolator%value_s(1:sz1+sz_derivative_eta1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_n(1:sz1)

    case(1170)  !left: Neumann, right: Neumann, bottom: Neuman, Top: Neumann

       sz_derivative_eta1 = 2
       sz_derivative_eta2 = 2
       interpolator%size_coeffs1 = sz1 + sz_derivative_eta1
       interpolator%size_coeffs2 = sz2 + sz_derivative_eta2
       interpolator%size_t1 = order1 + sz1 + sz_derivative_eta1
       interpolator%size_t2 = order2 + sz2 + sz_derivative_eta2
       
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       SLL_ALLOCATE( data_array_tmp(1:sz1,1:sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta1(2,sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta2(sz_derivative_eta2,sz1+sz_derivative_eta1),ierr)
       data_array_tmp = data_array(1:sz1,1:sz2)
       point_location_eta1_deriv(1) = 1
       point_location_eta1_deriv(2) = sz1
       data_array_deriv_eta1(1,1:sz2)     = 0.0_f64
       data_array_deriv_eta1(2,1:sz2)     = 0.0_f64
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)= 0.0_f64
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)= 0.0_f64
       call spli2d_custom_derder(&
            sz1,&
            sz_derivative_eta1,&
            order1, &
            point_location_eta1, &
            point_location_eta1_deriv,&
            sz2, &
            sz_derivative_eta2,&
            order2, point_location_eta2, &
            point_location_eta2_deriv,&
            data_array_tmp,&
            data_array_deriv_eta1,&
            data_array_deriv_eta2,&
            interpolator%coeff_splines,&!(1:sz1,1:sz2),&
            interpolator%t1,&!(1:sz1+order1), &
            interpolator%t2)!(1:sz2+order2) )

       SLL_DEALLOCATE( data_array_deriv_eta1,ierr)
       SLL_DEALLOCATE( data_array_deriv_eta2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_w(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1:sz1+sz_derivative_eta1,1)   = interpolator%value_s(1:sz1+sz_derivative_eta1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_n(1:sz1)
    case(2338)  !left: Dirichlet, right: Hermite, bottom: Hermite, Top: Hermite
       sz_derivative_eta1 = 2
       sz_derivative_eta2 = 2
       interpolator%size_coeffs1 = sz1 + sz_derivative_eta1
       interpolator%size_coeffs2 = sz2 + sz_derivative_eta2
       interpolator%size_t1 = order1 + sz1 + sz_derivative_eta1
       interpolator%size_t2 = order2 + sz2 + sz_derivative_eta2
       
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       SLL_ALLOCATE( data_array_tmp(1:sz1,1:sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta1(2,sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta2(sz_derivative_eta2,sz1+sz_derivative_eta1),ierr)
       data_array_tmp = data_array(1:sz1,1:sz2)
       point_location_eta1_deriv(1) = 1
       point_location_eta1_deriv(2) = sz1
       data_array_deriv_eta1(1,1:sz2)     = interpolator%slope(1:sz2,W)
       data_array_deriv_eta1(2,1:sz2)     = interpolator%slope(1:sz2,E)
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=interpolator%slope(1:sz1+sz_derivative_eta1,S)
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)=interpolator%slope(1:sz1+sz_derivative_eta1,N)
       call spli2d_custom_derder(&
            sz1,&
            sz_derivative_eta1,&
            order1, &
            point_location_eta1, &
            point_location_eta1_deriv,&
            sz2, &
            sz_derivative_eta2,&
            order2, point_location_eta2, &
            point_location_eta2_deriv,&
            data_array_tmp,&
            data_array_deriv_eta1,&
            data_array_deriv_eta2,&
            interpolator%coeff_splines,&!(1:sz1,1:sz2),&
            interpolator%t1,&!(1:sz1+order1), &
            interpolator%t2)!(1:sz2+order2) )

       SLL_DEALLOCATE( data_array_deriv_eta1,ierr)
       SLL_DEALLOCATE( data_array_deriv_eta2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_w(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1:sz1+sz_derivative_eta1,1)   = interpolator%value_s(1:sz1+sz_derivative_eta1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_n(1:sz1)
       
    case(2145) !left: Dirichlet, right: Hermite, bottom: Dirichlet, Top: Hermite  
       sz_derivative_eta1 = 2
       sz_derivative_eta2 = 2
       interpolator%size_coeffs1 = sz1 + sz_derivative_eta1
       interpolator%size_coeffs2 = sz2 + sz_derivative_eta2
       interpolator%size_t1 = order1 + sz1 + sz_derivative_eta1
       interpolator%size_t2 = order2 + sz2 + sz_derivative_eta2
       
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       SLL_ALLOCATE( data_array_tmp(1:sz1,1:sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta1(2,sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta2(sz_derivative_eta2,sz1+sz_derivative_eta1),ierr)
       data_array_tmp = data_array(1:sz1,1:sz2)
       point_location_eta1_deriv(1) = 1
       point_location_eta1_deriv(2) = sz1
       data_array_deriv_eta1(1,1:sz2)     = interpolator%slope(1:sz2,W)
       data_array_deriv_eta1(2,1:sz2)     = interpolator%slope(1:sz2,E)
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=interpolator%slope(1:sz1+sz_derivative_eta1,S)
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)=interpolator%slope(1:sz1+sz_derivative_eta1,N)
       call spli2d_custom_derder(&
            sz1,&
            sz_derivative_eta1,&
            order1, &
            point_location_eta1, &
            point_location_eta1_deriv,&
            sz2, &
            sz_derivative_eta2,&
            order2, point_location_eta2, &
            point_location_eta2_deriv,&
            data_array_tmp,&
            data_array_deriv_eta1,&
            data_array_deriv_eta2,&
            interpolator%coeff_splines,&!(1:sz1,1:sz2),&
            interpolator%t1,&!(1:sz1+order1), &
            interpolator%t2)!(1:sz2+order2) )

       SLL_DEALLOCATE( data_array_deriv_eta1,ierr)
       SLL_DEALLOCATE( data_array_deriv_eta2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_w(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1:sz1+sz_derivative_eta1,1)   = interpolator%value_s(1:sz1+sz_derivative_eta1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_n(1:sz1)


    case(2124)  !left: Hermite, right: Dirichlet, bottom: Dirichlet, Top: Hermite

       sz_derivative_eta1 = 2
       sz_derivative_eta2 = 2
       interpolator%size_coeffs1 = sz1 + sz_derivative_eta1
       interpolator%size_coeffs2 = sz2 + sz_derivative_eta2
       interpolator%size_t1 = order1 + sz1 + sz_derivative_eta1
       interpolator%size_t2 = order2 + sz2 + sz_derivative_eta2
       
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       SLL_CLEAR_ALLOCATE( data_array_tmp(1:sz1,1:sz2),ierr)
       SLL_CLEAR_ALLOCATE( data_array_deriv_eta1(sz_derivative_eta1,1:sz2),ierr)
       SLL_CLEAR_ALLOCATE( data_array_deriv_eta2(sz_derivative_eta2,1:sz1+sz_derivative_eta1),ierr)
       data_array_tmp = data_array(1:sz1,1:sz2)
       point_location_eta1_deriv(1) = 1
       point_location_eta1_deriv(2) = sz1
       data_array_deriv_eta1(1,1:sz2)     = interpolator%slope(1:sz2,W)
       data_array_deriv_eta1(2,1:sz2)     = interpolator%slope(1:sz2,E)
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=interpolator%slope(1:sz1+sz_derivative_eta1,S)
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)=interpolator%slope(1:sz1+sz_derivative_eta1,N)
       call spli2d_custom_derder(&
            sz1,&
            sz_derivative_eta1,&
            order1, &
            point_location_eta1, &
            point_location_eta1_deriv,&
            sz2, &
            sz_derivative_eta2,&
            order2, point_location_eta2, &
            point_location_eta2_deriv,&
            data_array_tmp,&
            data_array_deriv_eta1,&
            data_array_deriv_eta2,&
            interpolator%coeff_splines,&!(1:sz1,1:sz2),&
            interpolator%t1,&!(1:sz1+order1), &
            interpolator%t2)!(1:sz2+order2) )

       SLL_DEALLOCATE( data_array_deriv_eta1,ierr)
       SLL_DEALLOCATE( data_array_deriv_eta2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_w(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:sz1,1)   = interpolator%value_s(1:sz1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_n(1:sz1)

    case(2148)  !left:Hermite , right: Hermite, bottom: Dirichlet, Top: Hermite
       sz_derivative_eta1 = 2
       sz_derivative_eta2 = 2
       interpolator%size_coeffs1 = sz1 + sz_derivative_eta1
       interpolator%size_coeffs2 = sz2 + sz_derivative_eta2
       interpolator%size_t1 = order1 + sz1 + sz_derivative_eta1
       interpolator%size_t2 = order2 + sz2 + sz_derivative_eta2
       
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       SLL_ALLOCATE( data_array_tmp(1:sz1,1:sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta1(sz_derivative_eta1,1:sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta2(sz_derivative_eta2,1:sz1+sz_derivative_eta1),ierr)
       data_array_tmp = data_array(1:sz1,1:sz2)
       point_location_eta1_deriv(1) = 1
       point_location_eta1_deriv(2) = sz1
       data_array_deriv_eta1(1,1:sz2)     = interpolator%slope(1:sz2,W)
       data_array_deriv_eta1(2,1:sz2)     = interpolator%slope(1:sz2,E)
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=interpolator%slope(1:sz1+sz_derivative_eta1,S)
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)=interpolator%slope(1:sz1+sz_derivative_eta1,N)
       call spli2d_custom_derder(&
            sz1,&
            sz_derivative_eta1,&
            order1, &
            point_location_eta1, &
            point_location_eta1_deriv,&
            sz2, &
            sz_derivative_eta2,&
            order2, point_location_eta2, &
            point_location_eta2_deriv,&
            data_array_tmp,&
            data_array_deriv_eta1,&
            data_array_deriv_eta2,&
            interpolator%coeff_splines,&!(1:sz1,1:sz2),&
            interpolator%t1,&!(1:sz1+order1), &
            interpolator%t2)!(1:sz2+order2) )

       SLL_DEALLOCATE( data_array_deriv_eta1,ierr)
       SLL_DEALLOCATE( data_array_deriv_eta2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_w(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:sz1,1)   = interpolator%value_s(1:sz1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_n(1:sz1)

    case(2316)  !left: Hermite, right: Dirichlet, bottom: Hermite, Top: Hermite
       sz_derivative_eta1 = 2
       sz_derivative_eta2 = 2
       interpolator%size_coeffs1 = sz1 + sz_derivative_eta1
       interpolator%size_coeffs2 = sz2 + sz_derivative_eta2
       interpolator%size_t1 = order1 + sz1 + sz_derivative_eta1
       interpolator%size_t2 = order2 + sz2 + sz_derivative_eta2
       
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       SLL_ALLOCATE( data_array_tmp(1:sz1,1:sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta1(sz_derivative_eta1,1:sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta2(sz_derivative_eta2,1:sz1+sz_derivative_eta1),ierr)
       data_array_tmp = data_array(1:sz1,1:sz2)
       point_location_eta1_deriv(1) = 1
       point_location_eta1_deriv(2) = sz1
       data_array_deriv_eta1(1,1:sz2)     = interpolator%slope(1:sz2,W)
       data_array_deriv_eta1(2,1:sz2)     = interpolator%slope(1:sz2,E)
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=interpolator%slope(1:sz1+sz_derivative_eta1,S)
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)=interpolator%slope(1:sz1+sz_derivative_eta1,N)
       call spli2d_custom_derder(&
            sz1,&
            sz_derivative_eta1,&
            order1, &
            point_location_eta1, &
            point_location_eta1_deriv,&
            sz2, &
            sz_derivative_eta2,&
            order2, point_location_eta2, &
            point_location_eta2_deriv,&
            data_array_tmp,&
            data_array_deriv_eta1,&
            data_array_deriv_eta2,&
            interpolator%coeff_splines,&!(1:sz1,1:sz2),&
            interpolator%t1,&!(1:sz1+order1), &
            interpolator%t2)!(1:sz2+order2) )

       SLL_DEALLOCATE( data_array_deriv_eta1,ierr)
       SLL_DEALLOCATE( data_array_deriv_eta2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_w(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:sz1,1)   = interpolator%value_s(1:sz1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_n(1:sz1)

    case(2340) ! Hermite in al sides
       

       sz_derivative_eta1 = 2
       sz_derivative_eta2 = 2
       interpolator%size_coeffs1 = sz1 + sz_derivative_eta1
       interpolator%size_coeffs2 = sz2 + sz_derivative_eta2
       interpolator%size_t1 = order1 + sz1 + sz_derivative_eta1
       interpolator%size_t2 = order2 + sz2 + sz_derivative_eta2
       
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       SLL_ALLOCATE( data_array_tmp(1:sz1,1:sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta1(sz_derivative_eta1,1:sz2),ierr)
       SLL_ALLOCATE( data_array_deriv_eta2(sz_derivative_eta2,1:sz1+sz_derivative_eta1),ierr)
       data_array_tmp = data_array(1:sz1,1:sz2)
       point_location_eta1_deriv(1) = 1
       point_location_eta1_deriv(2) = sz1
       data_array_deriv_eta1(1,1:sz2)     = interpolator%slope(1:sz2,W)
       data_array_deriv_eta1(2,1:sz2)     = interpolator%slope(1:sz2,E)
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=interpolator%slope(1:sz1+sz_derivative_eta1,S)
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)=interpolator%slope(1:sz1+sz_derivative_eta1,N)
       call spli2d_custom_derder(&
            sz1,&
            sz_derivative_eta1,&
            order1, &
            point_location_eta1, &
            point_location_eta1_deriv,&
            sz2, &
            sz_derivative_eta2,&
            order2, point_location_eta2, &
            point_location_eta2_deriv,&
            data_array_tmp,&
            data_array_deriv_eta1,&
            data_array_deriv_eta2,&
            interpolator%coeff_splines,&!(1:sz1,1:sz2),&
            interpolator%t1,&!(1:sz1+order1), &
            interpolator%t2)!(1:sz2+order2) )

       SLL_DEALLOCATE( data_array_deriv_eta1,ierr)
       SLL_DEALLOCATE( data_array_deriv_eta2,ierr)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_w(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:sz1,1)   = interpolator%value_s(1:sz1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_n(1:sz1)

    end select
    interpolator%coefficients_set = .true.
   
    SLL_DEALLOCATE(point_location_eta2,ierr)
    SLL_DEALLOCATE(point_location_eta1,ierr)
    SLL_DEALLOCATE(point_location_eta1_tmp,ierr)
    SLL_DEALLOCATE(point_location_eta2_tmp,ierr)
    SLL_DEALLOCATE(data_array_tmp,ierr)

  end subroutine !compute_interpolants_ad2d

  function coefficients_are_set_ad2d( interpolator ) result(res)
    sll_interpolator, intent(in)  :: interpolator
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
    sll_interpolator, intent(in)  :: interpolator
    sll_real64, intent(in)         :: eta1
    sll_real64, intent(in)         :: eta2
    sll_real64                     :: val
    sll_int32 :: size_coeffs1
    sll_int32 :: size_coeffs2
    sll_int32 :: partint
    !sll_real64 :: bvalue2d
    !integer  :: li_i, li_j, li_mflag, li_wy
    sll_real64 :: res1,res2
    !sll_int32  :: ierr
    sll_real64,dimension(:), pointer:: tmp_tx,tmp_ty
    sll_real64,dimension(:,:), pointer:: tmp_coeff
    !type(sll_time_mark)  :: t0 
    !double precision :: time

    size_coeffs1 = interpolator%size_coeffs1
    size_coeffs2 = interpolator%size_coeffs2

    res1 = eta1
    res2 = eta2

    select case (interpolator%bc_selector)
    case (0) ! periodic-periodic

       if( res1 < interpolator%eta1_min ) then
          partint = &
               IDINT(abs(res1)/(interpolator%eta1_max-interpolator%eta1_min))
          res1 = res1 + &
               (partint+1) *( interpolator%eta1_max-interpolator%eta1_min)
       else if( res1 >  interpolator%eta1_max ) then
          partint = &
               IDINT(abs(res1)/(interpolator%eta2_max-interpolator%eta2_min))
          res1 = res1 + partint*(interpolator%eta2_min-interpolator%eta2_max)
       end if
       if( res2 < interpolator%eta2_min ) then
          partint = &
               IDINT(abs(res2)/(interpolator%eta2_max-interpolator%eta2_min))
          res2 = res2 + &
               (partint+1) *(interpolator%eta2_max-interpolator%eta2_min)
       else if( res2 >  interpolator%eta2_max ) then
          partint = &
               IDINT(abs(res2)/(interpolator%eta2_max-interpolator%eta2_min))
          res2 = res2 + partint*(interpolator%eta2_min-interpolator%eta2_max)
       end if
       

    case (9) ! 2. dirichlet-left, dirichlet-right, periodic


       if( res2 < interpolator%eta2_min ) then
          partint = &
               IDINT(abs(res2)/(interpolator%eta2_max-interpolator%eta2_min))
          res2 = res2 + &
               (partint+1) *(interpolator%eta2_max-interpolator%eta2_min)
       else if( res2 >  interpolator%eta2_max ) then
          partint = &
               IDINT(abs(res2)/(interpolator%eta2_max-interpolator%eta2_min))
          res2 = res2 + partint*(interpolator%eta2_min-interpolator%eta2_max)
       end if

       SLL_ASSERT( res1 >= interpolator%eta1_min )
       SLL_ASSERT( res1 <= interpolator%eta1_max )
       if ( res1 > interpolator%eta1_max) then 
          print*, 'problem  x > eta1_max'
          !res1 =  interpolator%eta1_max
          stop
       end if
       if ( res1 < interpolator%eta1_min) then 
          print*, 'problem  x < eta1_min'
          !res1 =  interpolator%eta1_min
          !print*, res1,interpolator%eta1_min 
          stop
       end if
      
  
    case(576) !  3. periodic, dirichlet-bottom, dirichlet-top



       if( res1 < interpolator%eta1_min ) then
          partint = &
               IDINT(abs(res1)/(interpolator%eta1_max-interpolator%eta1_min))
          res1 = res1 + &
               (partint+1) *( interpolator%eta1_max-interpolator%eta1_min)
       else if( res1 >  interpolator%eta1_max ) then
          partint = &
               IDINT(abs(res1)/(interpolator%eta2_max-interpolator%eta2_min))
          res1 = res1 + partint*(interpolator%eta2_min-interpolator%eta2_max)
       end if
       SLL_ASSERT( res2 >= interpolator%eta2_min )
       SLL_ASSERT( res2 <= interpolator%eta2_max )
       if ( res2 > interpolator%eta2_max) then 
          print*, 'problem  y > eta2_max'
          stop
       end if
       if ( res2 < interpolator%eta2_min) then 
          print*, 'problem  y < eta2_min'
  
          stop
       end if
       
 
    end select

    SLL_ASSERT( res1 >= interpolator%eta1_min )
    SLL_ASSERT( res1 <= interpolator%eta1_max )
    SLL_ASSERT( res2 >= interpolator%eta2_min )
    SLL_ASSERT( res2 <= interpolator%eta2_max )
    if ( res1 > interpolator%eta1_max) then 
       print*, 'ERROR, interpolate_value_ad2d(): problem  x > eta1_max'
       stop
    end if
    if ( res1 < interpolator%eta1_min) then 
       print*, 'ERROR, interpolate_value_ad2d(): problem  x < eta1_min'
    end if
    if ( res2 > interpolator%eta2_max) then 
       print*, 'ERROR, interpolate_value_ad2d(): problem  y > eta2_max'
       print*, res2, interpolator%eta2_min,interpolator%eta2_max
    end if
    if ( res2 < interpolator%eta2_min) then 
       print*, 'ERROR, interpolate_value_ad2d(): problem  y < eta2_min'
    end if

    tmp_tx => interpolator%t1(1:interpolator%size_t1)
    !print*, 't1',tmp_tx
    tmp_ty => interpolator%t2(1:interpolator%size_t2)
    !print*, 't2',tmp_ty
    tmp_coeff =>interpolator%coeff_splines(1:size_coeffs1,1:size_coeffs2)
    !print*, 'coef',tmp_coeff
    !call interv( tmp_ty, interpolator%size_t2, res2, li_wy, li_mflag )
  !  call set_time_mark(t0)
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
    ! time = time_elapsed_since(t0)
   !  print *, 'time elapsed since t0 : ',time
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

    sll_interpolator, intent(in)  :: interpolator
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
    sll_int32 :: partint

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
          partint = IDINT(abs(res1)/(interpolator%eta1_max-interpolator%eta1_min))
          res1 = res1+ (partint+1) *( interpolator%eta1_max-interpolator%eta1_min)
       else if( res1 >  interpolator%eta1_max ) then
          partint = IDINT(abs(res1)/(interpolator%eta2_max-interpolator%eta2_min))
          res1 = res1 + partint*(interpolator%eta2_min-interpolator%eta2_max)
       end if
       if( res2 < interpolator%eta2_min ) then
          partint = IDINT(abs(res2)/(interpolator%eta2_max-interpolator%eta2_min))
          res2 = res2 + (partint+1) *(interpolator%eta2_max-interpolator%eta2_min)
       else if( res2 >  interpolator%eta2_max ) then
          partint = IDINT(abs(res2)/(interpolator%eta2_max-interpolator%eta2_min))
          res2 = res2 + partint*(interpolator%eta2_min-interpolator%eta2_max)
       end if
       

    case (9) ! 2. dirichlet-left, dirichlet-right, periodic

      
       if( res2 < interpolator%eta2_min ) then
          partint = IDINT(abs(res2)/(interpolator%eta2_max-interpolator%eta2_min))
          res2 = res2 + (partint+1) *(interpolator%eta2_max-interpolator%eta2_min)
       else if( res2 >  interpolator%eta2_max ) then
          partint = IDINT(abs(res2)/(interpolator%eta2_max-interpolator%eta2_min))
          res2 = res2 + partint*(interpolator%eta2_min-interpolator%eta2_max)
       end if

       SLL_ASSERT( res1 >= interpolator%eta1_min )
       SLL_ASSERT( res1 <= interpolator%eta1_max )
       if ( res1 > interpolator%eta1_max) then 
          print*, 'problem  x > eta1_max'
          stop
       end if
       if ( res1 < interpolator%eta1_min) then 
          print*, 'problem  x < eta1_min'
          stop
       end if
       
    case(576) !  3. periodic, dirichlet-bottom, dirichlet-top

       if( res1 < interpolator%eta1_min ) then
          partint = IDINT(abs(res1)/(interpolator%eta1_max-interpolator%eta1_min))
          res1 = res1+ (partint+1) *( interpolator%eta1_max-interpolator%eta1_min)
       else if( res1 >  interpolator%eta1_max ) then
          partint = IDINT(abs(res1)/(interpolator%eta2_max-interpolator%eta2_min))
          res1 = res1 + partint*(interpolator%eta2_min-interpolator%eta2_max)
       end if


       SLL_ASSERT( res2 >= interpolator%eta2_min )
       SLL_ASSERT( res2 <= interpolator%eta2_max )
       if ( res2 > interpolator%eta2_max) then 
          print*, 'problem  y > eta2_max'
          stop
       end if
       if ( res2 < interpolator%eta2_min) then 
          print*, 'problem  y < eta2_min'
          stop
       end if
       
    end select

    SLL_ASSERT( res1 >= interpolator%eta1_min )
    SLL_ASSERT( res1 <= interpolator%eta1_max )
    SLL_ASSERT( res2 >= interpolator%eta2_min )
    SLL_ASSERT( res2 <= interpolator%eta2_max )
    if ( res1 > interpolator%eta1_max) then 
       print*, 'problem  x > eta1_max'
       stop
    end if
    if ( res1 < interpolator%eta1_min) then 
       print*, 'problem  x < eta1_min'
       stop
    end if
    if ( res2 > interpolator%eta2_max) then 
       print*, 'problem  y > eta2_max'
       stop
    end if
    if ( res2 < interpolator%eta2_min) then 
       print*, 'problem  y < eta2_min'
       stop
    end if
    !SLL_ALLOCATE(knot1_tmp(interpolator%size_t1),ierr)
    !SLL_ALLOCATE(knot2_tmp(interpolator%size_t2),ierr)
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

    sll_interpolator, intent(in)  :: interpolator
    sll_real64, intent(in)         :: eta1
    sll_real64, intent(in)         :: eta2
    sll_real64                     :: val
    sll_int32 :: size_coeffs1
    sll_int32 :: size_coeffs2
    sll_int32 :: partint
    !sll_real64 :: dvalue2d
    sll_real64 :: res1,res2
    !sll_int32 :: ierr
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
          partint = IDINT(abs(res1)/(interpolator%eta1_max-interpolator%eta1_min))
          res1 = res1+ (partint+1) *( interpolator%eta1_max-interpolator%eta1_min)
       else if( res1 >  interpolator%eta1_max ) then
          partint = IDINT(abs(res1)/(interpolator%eta2_max-interpolator%eta2_min))
          res1 = res1 + partint*(interpolator%eta2_min-interpolator%eta2_max)
       end if
       if( res2 < interpolator%eta2_min ) then
          partint = IDINT(abs(res2)/(interpolator%eta2_max-interpolator%eta2_min))
          res2 = res2 + (partint+1) *(interpolator%eta2_max-interpolator%eta2_min)
       else if( res2 >  interpolator%eta2_max ) then
          partint = IDINT(abs(res2)/(interpolator%eta2_max-interpolator%eta2_min))
          res2 = res2 + partint*(interpolator%eta2_min-interpolator%eta2_max)
       end if
          
    case (9) ! 2. dirichlet-left, dirichlet-right, periodic
       
      
       if( res2 < interpolator%eta2_min ) then
          partint = IDINT(abs(res2)/(interpolator%eta2_max-interpolator%eta2_min))
          res2 = res2 + (partint+1) *(interpolator%eta2_max-interpolator%eta2_min)
       else if( res2 >  interpolator%eta2_max ) then
          partint = IDINT(abs(res2)/(interpolator%eta2_max-interpolator%eta2_min))
          res2 = res2 + partint*(interpolator%eta2_min-interpolator%eta2_max)
       end if
       SLL_ASSERT( res1 >= interpolator%eta1_min )
       SLL_ASSERT( res1 <= interpolator%eta1_max )
       if ( res1 > interpolator%eta1_max) then 
          print*,'problem  x > eta1_max'
          stop
       end if
       if ( res1 < interpolator%eta1_min) then 
          print*, 'problem  x < eta1_min'
          stop
       end if
       
    case(576) !  3. periodic, dirichlet-bottom, dirichlet-top
       
       if( res1 < interpolator%eta1_min ) then
          partint = IDINT(abs(res1)/(interpolator%eta1_max-interpolator%eta1_min))
          res1 = res1+ (partint+1) *( interpolator%eta1_max-interpolator%eta1_min)
       else if( res1 >  interpolator%eta1_max ) then
          partint = IDINT(abs(res1)/(interpolator%eta2_max-interpolator%eta2_min))
          res1 = res1 + partint*(interpolator%eta2_min-interpolator%eta2_max)
       end if

       
       SLL_ASSERT( res2 >= interpolator%eta2_min )
       SLL_ASSERT( res2 <= interpolator%eta2_max )
       if ( res2 > interpolator%eta2_max) then 
          print*, 'problem  y > eta2_max'
          stop
       end if
       if ( res2 < interpolator%eta2_min) then 
          print*, 'problem  y < eta2_min'
          stop
       end if
       

    end select
    SLL_ASSERT( res1 >= interpolator%eta1_min )
    SLL_ASSERT( res1 <= interpolator%eta1_max )
    SLL_ASSERT( res2 >= interpolator%eta2_min )
    SLL_ASSERT( res2 <= interpolator%eta2_max )
    if ( res1 > interpolator%eta1_max) then 
       print*, 'problem  x > eta1_max'
       stop
    end if
    if ( res1 < interpolator%eta1_min) then 
       print*, 'problem  x < eta1_min'
       stop
    end if
    if ( res2 > interpolator%eta2_max) then 
       print*, 'problem  y > eta2_max'
       stop
    end if
    if ( res2 < interpolator%eta2_min) then 
       print*, 'problem  y < eta2_min'
       stop
    end if
   ! SLL_ALLOCATE(knot1_tmp(interpolator%size_t1),ierr)
   ! SLL_ALLOCATE(knot2_tmp(interpolator%size_t2),ierr)
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

  
function interpolate_array_ad2d( this,              &
                                 num_points1,       &
                                 num_points2,       &
                                 data_in,           &
                                 eta1,              &
                                 eta2 ) result(res)
    
  sll_interpolator, intent(in)  :: this
  sll_real64,  dimension(:,:), intent(in)         :: eta1
  sll_real64,  dimension(:,:), intent(in)         :: eta2
  sll_real64, dimension(:,:), intent(in)         :: data_in
  sll_int32, intent(in)         :: num_points1
  sll_int32, intent(in)         :: num_points2

  sll_real64, dimension(num_points1,num_points2) :: res
    
  
  res = -1000000._f64
  print *,this%num_pts1
  print *,maxval(eta1)
  print *,maxval(eta2)
  print *,maxval(data_in)
  print *,num_points1
  print *,num_points2
  SLL_ERROR( " interpolate_array_ad2d: not implemented ")

end function 
  
function interpolate_2d_array_disp_ad2d( this,        &
                                         num_points1, &
                                         num_points2, &
                                         data_in,     &
                                         alpha1,      &
                                         alpha2) result(res)
      
  sll_interpolator, intent(in)                   :: this
  sll_int32, intent(in)                          :: num_points1  
  sll_int32, intent(in)                          :: num_points2 
  sll_real64, dimension(:,:), intent(in)         :: data_in
  sll_real64, dimension(:,:), intent(in)         :: alpha1
  sll_real64, dimension(:,:), intent(in)         :: alpha2  
  sll_real64, dimension(num_points1,num_points2) :: res
    
  print *,this%num_pts1
  print *,num_points1 
  print *,num_points2
  print *,maxval(data_in)
  print *,alpha1
  print *,alpha2     
  res = -1000000._f64
  SLL_ERROR("interpolate_2d_array_disp_ad2d: not implemented.")
    
end function !interpolate_2d_array_disp_ad2d
    
function get_coefficients_ad2d(interpolator)
  sll_interpolator, intent(in)        :: interpolator
  sll_real64, dimension(:,:), pointer :: get_coefficients_ad2d     

  get_coefficients_ad2d => interpolator%coeff_splines

end function get_coefficients_ad2d
  
end module sll_module_arbitrary_degree_spline_interpolator_2d
