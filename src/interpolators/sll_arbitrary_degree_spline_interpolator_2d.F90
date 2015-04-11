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

implicit none
private

! in what follows, the direction '1' is in the contiguous memory direction.
!> Arbitrary degree version of 2d irnterpolator
type, extends(sll_interpolator_2d_base) :: sll_arbitrary_degree_spline_interpolator_2d           

   sll_int32                         :: num_pts1
   sll_int32                         :: num_pts2
   sll_real64                        :: eta1_min
   sll_real64                        :: eta1_max
   sll_real64                        :: eta2_min
   sll_real64                        :: eta2_max
   sll_int32                         :: bc_left
   sll_int32                         :: bc_right
   sll_int32                         :: bc_bottom
   sll_int32                         :: bc_top
   sll_int32                         :: spline_degree1
   sll_int32                         :: spline_degree2
   sll_real64, dimension(:), pointer :: knots1
   sll_real64, dimension(:), pointer :: knots2

   ! some knot-like arrays needed by the spli2d_per routine
   sll_real64, dimension(:), pointer :: t1
   sll_real64, dimension(:), pointer :: t2
   sll_int32  :: size_t1
   sll_int32  :: size_t2 
   sll_int64  :: bc_selector ! this is set in initialization
   sll_real64, dimension(:,:), pointer :: coeff_splines
   sll_int32                  :: size_coeffs1
   sll_int32                  :: size_coeffs2
   logical    :: coefficients_set = .false.
   ! table contains the coeff spline of the function in boundary 
   ! in the case of dirichlet boundary condition non homogene 
   sll_real64, dimension(:),pointer :: slope_left
   sll_real64, dimension(:),pointer :: slope_right
   sll_real64, dimension(:),pointer :: slope_bottom
   sll_real64, dimension(:),pointer :: slope_top
   sll_real64, dimension(:),pointer :: value_left
   sll_real64, dimension(:),pointer :: value_right
   sll_real64, dimension(:),pointer :: value_bottom
   sll_real64, dimension(:),pointer :: value_top
   logical    :: compute_slope_left = .TRUE.
   logical    :: compute_slope_right= .TRUE.
   logical    :: compute_slope_top = .TRUE.
   logical    :: compute_slope_bottom= .TRUE.
   logical    :: compute_value_left = .TRUE.
   logical    :: compute_value_right= .TRUE.
   logical    :: compute_value_top = .TRUE.
   logical    :: compute_value_bottom= .TRUE.
contains
   procedure, pass(interpolator) :: initialize=>initialize_ad2d_interpolator
   procedure, pass(interpolator) :: set_coefficients => set_coefficients_ad2d
   procedure, pass(interpolator) :: coefficients_are_set => &
         coefficients_are_set_ad2d
   ! better: pre-compute-interpolation-information or something...
   procedure :: compute_interpolants => compute_interpolants_ad2d
   ! procedure,  pass(interpolator) :: compute_spline_coefficients => &
   !     compute_spline_coefficients_ad2d
   !procedure, pass:: compute_spline_coefficients =>compute_spline_coefficients_ad2d
   procedure :: interpolate_value => interpolate_value_ad2d
   procedure :: interpolate_derivative_eta1 => interpolate_derivative1_ad2d
   procedure :: interpolate_derivative_eta2 => interpolate_derivative2_ad2d
   procedure, pass:: interpolate_array => interpolate_array_ad2d
   procedure, pass:: interpolate_array_disp => interpolate_2d_array_disp_ad2d
   procedure, pass:: get_coefficients => get_coefficients_ad2d
   procedure, pass:: delete => delete_arbitrary_degree_2d_interpolator
   procedure, pass:: set_values_at_boundary => set_boundary_value2d
   procedure, pass:: set_slopes_at_boundary => set_slope2d
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
    SLL_DEALLOCATE(interpolator%value_left,ierr)
    SLL_DEALLOCATE(interpolator%value_right,ierr)
    SLL_DEALLOCATE(interpolator%value_bottom,ierr)
    SLL_DEALLOCATE(interpolator%value_top,ierr)
    SLL_DEALLOCATE(interpolator%slope_left,ierr)
    SLL_DEALLOCATE(interpolator%slope_right,ierr)
    SLL_DEALLOCATE(interpolator%slope_bottom,ierr)
    SLL_DEALLOCATE(interpolator%slope_top,ierr)
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
  !> @param[in] bc_left  the boundary condition at left in the direction eta1
  !> @param[in] bc_right the boundary condition at right in the direction eta2
  !> @param[in] bc_bottom the boundary condition at left in the direction eta2
  !> @param[in] bc_top the boundary condition at right in the direction eta2
  !> @param[in] spline_degree1 the degree of B-spline in the direction eta1
  !> @param[in] spline_degree2 the degre of B-spline in the direction eta2
  !> @return the type sll_arbitrary_degree_spline_interpolator_2d

  function new_arbitrary_degree_spline_interp2d( &
    num_pts1, &
    num_pts2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    bc_left, &
    bc_right, &
    bc_bottom, &
    bc_top, &
    spline_degree1, &
    spline_degree2) result( res )

    type(sll_arbitrary_degree_spline_interpolator_2d), pointer :: res
    sll_int32, intent(in) :: num_pts1
    sll_int32, intent(in) :: num_pts2
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_real64, intent(in) :: eta2_min
    sll_real64, intent(in) :: eta2_max
    sll_int32, intent(in) :: bc_left
    sll_int32, intent(in) :: bc_right
    sll_int32, intent(in) :: bc_bottom
    sll_int32, intent(in) :: bc_top
    sll_int32, intent(in) :: spline_degree1
    sll_int32, intent(in) :: spline_degree2
    sll_int32 :: ierr
    
    SLL_ALLOCATE(res,ierr)

    call initialize_ad2d_interpolator( &
         res,&
         num_pts1, &
         num_pts2, &
         eta1_min, &
         eta1_max, &
         eta2_min, &
         eta2_max, &
         bc_left, &
         bc_right, &
         bc_bottom, &
         bc_top, &
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
  !> @param[in] bc_left  the boundary condition at left in the direction eta1
  !> @param[in] bc_right the boundary condition at right in the direction eta2
  !> @param[in] bc_bottom the boundary condition at left in the direction eta2
  !> @param[in] bc_top the boundary condition at right in the direction eta2
  !> @param[in] spline_degree1 the degree of B-spline in the direction eta1
  !> @param[in] spline_degree2 the degre of B-spline in the direction eta2
  !> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_2d
       subroutine initialize_ad2d_interpolator( &
    interpolator, &
    num_pts1, &
    num_pts2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    bc_left, &
    bc_right, &
    bc_bottom, &
    bc_top, &
    spline_degree1, &
    spline_degree2)

    class(sll_arbitrary_degree_spline_interpolator_2d):: interpolator
    sll_int32, intent(in) :: num_pts1
    sll_int32, intent(in) :: num_pts2
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_real64, intent(in) :: eta2_min
    sll_real64, intent(in) :: eta2_max
    sll_int32, intent(in) :: bc_left
    sll_int32, intent(in) :: bc_right
    sll_int32, intent(in) :: bc_bottom
    sll_int32, intent(in) :: bc_top
    sll_int32, intent(in) :: spline_degree1
    sll_int32, intent(in) :: spline_degree2
    sll_int32 :: ierr
    sll_int32 :: tmp1
    sll_int32 :: tmp2
    sll_int64 :: bc_selector
   

    ! do some argument checking...
    if(((bc_left  == SLL_PERIODIC).and.(bc_right.ne. SLL_PERIODIC)).or.&
       ((bc_right == SLL_PERIODIC).and.(bc_left .ne. SLL_PERIODIC)))then
       print *, 'initialize_arbitrary_degree_2d_interpolator, ERROR: ', &
            'if one boundary condition is specified as periodic, then ', &
            'both must be. Error in first direction.'
    end if

    if(((bc_bottom == SLL_PERIODIC).and.(bc_top.ne. SLL_PERIODIC)).or.&
       ((bc_top == SLL_PERIODIC).and.(bc_bottom .ne. SLL_PERIODIC)))then
       print *, 'initialize_arbitrary_degree_2d_interpolator, ERROR: ', &
            'if one boundary condition is specified as periodic, then ', &
            'both must be. Error in second direction.'
    end if

    bc_selector = 0

    if( bc_left == SLL_DIRICHLET ) then
       bc_selector = bc_selector + 1
    end if

    if( bc_left == SLL_NEUMANN ) then
       bc_selector = bc_selector + 2
    end if

    if( bc_left == SLL_HERMITE ) then
       bc_selector = bc_selector + 4
    end if

    if( bc_right == SLL_DIRICHLET ) then
       bc_selector = bc_selector + 8
    end if

    if( bc_right == SLL_NEUMANN ) then
       bc_selector = bc_selector + 16
    end if

   if( bc_right == SLL_HERMITE ) then
       bc_selector = bc_selector + 32
    end if

    if( bc_bottom == SLL_DIRICHLET ) then
       bc_selector = bc_selector + 64
    end if

    if( bc_bottom == SLL_NEUMANN ) then
       bc_selector = bc_selector + 128
    end if

    if( bc_bottom == SLL_HERMITE ) then
       bc_selector = bc_selector + 256
    end if

    if( bc_top == SLL_DIRICHLET ) then
       bc_selector = bc_selector + 512
    end if

    if( bc_top == SLL_NEUMANN ) then
       bc_selector = bc_selector + 1024
    end if

   if( bc_top == SLL_HERMITE ) then
       bc_selector = bc_selector + 2048
    end if

    ! Initialization in the type of interpolator
    interpolator%spline_degree1 = spline_degree1
    interpolator%spline_degree2 = spline_degree2
    interpolator%eta1_min = eta1_min
    interpolator%eta1_max = eta1_max
    interpolator%eta2_min = eta2_min
    interpolator%eta2_max = eta2_max
    interpolator%bc_left  = bc_left
    interpolator%bc_right = bc_right
    interpolator%bc_bottom= bc_bottom
    interpolator%bc_top   = bc_top
    interpolator%bc_selector = bc_selector
    interpolator%num_pts1 = num_pts1
    interpolator%num_pts2 = num_pts2
   

    SLL_CLEAR_ALLOCATE(interpolator%value_left  (1:num_pts2),ierr)
    SLL_CLEAR_ALLOCATE(interpolator%value_right (1:num_pts2),ierr)
    SLL_CLEAR_ALLOCATE(interpolator%value_bottom(1:num_pts1+2),ierr)
    SLL_CLEAR_ALLOCATE(interpolator%value_top   (1:num_pts1+2),ierr)

    SLL_CLEAR_ALLOCATE(interpolator%slope_left  (1:num_pts2),ierr)
    SLL_CLEAR_ALLOCATE(interpolator%slope_right (1:num_pts2),ierr)
    SLL_CLEAR_ALLOCATE(interpolator%slope_bottom(1:num_pts1+2),ierr)
    SLL_CLEAR_ALLOCATE(interpolator%slope_top   (1:num_pts1+2),ierr)
  
    ! tmp1 and tmp2 is the maximun (not absolue) for the size of coefficients
    select case (bc_selector)
    case (0) ! 1. periodic-periodic
       
       ! Allocate the knots in each direction 
       SLL_ALLOCATE( interpolator%knots1(2*spline_degree1+2),ierr )
       SLL_ALLOCATE( interpolator%knots2(2*spline_degree2+2),ierr )

       ! Allocate the coefficients spline
       tmp1 = num_pts1+ 4*spline_degree1! *num_pts1 !+ 2*spline_degree1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2 !+ 2*spline_degree2
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

    case (9) ! 2. dirichlet-left, dirichlet-right, periodic
       ! Allocate the knots in each direction 
       SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
       SLL_ALLOCATE( interpolator%knots2(2*spline_degree2+2),ierr )

       ! Allocate the coefficients spline
       tmp1 = num_pts1+ 4*spline_degree1!*num_pts1! + spline_degree1 !- 1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2! + 2*spline_degree2
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

    case (576) ! 3. periodic, dirichlet-bottom, dirichlet-top

        ! Allocate the knots in each direction 
       SLL_ALLOCATE( interpolator%knots1(2*spline_degree1+2),ierr )
       SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )

       ! Allocate the coefficients spline
       tmp1 = num_pts1+ 4*spline_degree1!*num_pts1! + 2*spline_degree1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2 + spline_degree2 !- 1
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

    case (585) ! 4. dirichlet in all sides
        ! Allocate the knots in each direction
       SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
       SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )

       ! Allocate the coefficients spline
       tmp1 = num_pts1+ 4*spline_degree1!*num_pts1! + spline_degree1 !- 1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2! + spline_degree2 !- 1
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)


    case(650) !left: Neumann, right: Dirichlet, bottom: Neumann, Top: Dirichlet 
       ! Allocate the knots in each direction
       SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
       SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )
       
       ! Allocate the coefficients spline
       tmp1 = num_pts1+ 4*spline_degree1!*num_pts1! + spline_degree1 !- 1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2! + spline_degree2 !- 1
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

    case(657) !left: Dirichlet, right: Neumann, bottom: Neumann, Top: Dirichlet 
       ! Allocate the knots in each direction
       SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
       SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )
       
       ! Allocate the coefficients spline
       tmp1 = num_pts1+ 4*spline_degree1!*num_pts1! + spline_degree1 !- 1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2! + spline_degree2 !- 1
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

    case(780)  !left: Hermite, right: Dirichlet, bottom: Hermite, Top: Dirichlet

       SLL_CLEAR_ALLOCATE( interpolator%knots1(1:num_pts1+2*spline_degree1),ierr )
       SLL_CLEAR_ALLOCATE( interpolator%knots2(1:num_pts2+2*spline_degree2),ierr )

       ! Allocate the coefficients spline
       tmp1 = num_pts1+ 4*spline_degree1!*num_pts1! + spline_degree1 !- 1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2! + spline_degree2 !- 1
       SLL_CLEAR_ALLOCATE( interpolator%coeff_splines(1:tmp1,1:tmp2),ierr)

    case(801)  !left: Dirichlet, right: Hermite, bottom: Hermite, Top: Dirichlet

       SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
       SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )

       ! Allocate the coefficients spline
       tmp1 = num_pts1+ 4*spline_degree1!*num_pts1! + spline_degree1 !- 1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2! + spline_degree2 !- 1
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

    case(804)  !left: Hermite, right: Hermite, bottom: Hermite, Top: Dirichlet
       
       SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
       SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )
       
       ! Allocate the coefficients spline
       tmp1 = num_pts1+ 4*spline_degree1!*num_pts1! + spline_degree1 !- 1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2! + spline_degree2 !- 1
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

    case(1098)  !left: Neumann, right: Dirichlet, bottom: Dirichlet, Top: Neumann
       
       SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
       SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )

       ! Allocate the coefficients spline
       tmp1 = num_pts1+ 4*spline_degree1!*num_pts1! + spline_degree1 !- 1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2! + spline_degree2 !- 1
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

    case(1105)  !left: Dirichlet, right: Neumann, bottom: Dirichlet, Top: Neumann
       
       SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
       SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )

       ! Allocate the coefficients spline
       tmp1 = num_pts1+ 4*spline_degree1!*num_pts1! + spline_degree1 !- 1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2! + spline_degree2 !- 1
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

    case(1170)  !left: Neumann, right: Neumann, bottom: Neuman, Top: Neumann
       
       SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
       SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )

       ! Allocate the coefficients spline
       tmp1 = num_pts1+ 4*spline_degree1!*num_pts1! + spline_degree1 !- 1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2! + spline_degree2 !- 1
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

    case(2124)  !left: Hermite, right: Dirichlet, bottom: Dirichlet, Top: Hermite

       SLL_CLEAR_ALLOCATE( interpolator%knots1(1:num_pts1+2*spline_degree1),ierr )
       SLL_CLEAR_ALLOCATE( interpolator%knots2(1:num_pts2+2*spline_degree2),ierr )

       ! Allocate the coefficients spline
       tmp1 = num_pts1+ 4*spline_degree1!*num_pts1! + spline_degree1 !- 1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2! + spline_degree2 !- 1
       SLL_CLEAR_ALLOCATE( interpolator%coeff_splines(1:tmp1,1:tmp2),ierr)

    case(2145)  !left: Dirichlet, right: Hermite, bottom: Dirichlet, Top: Hermite  

       SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
       SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )

       ! Allocate the coefficients spline
       tmp1 = num_pts1+ 4*spline_degree1!*num_pts1! + spline_degree1 !- 1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2! + spline_degree2 !- 1
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

    case(2148)  !left:Hermite , right: Hermite, bottom: Dirichlet, Top: Hermite  

       SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
       SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )

       ! Allocate the coefficients spline
       tmp1 = num_pts1+ 4*spline_degree1!*num_pts1! + spline_degree1 !- 1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2! + spline_degree2 !- 1
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

    case(2316)  !left: Hermite, right: Dirichlet, bottom: Hermite, Top: Hermite
       
       SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
       SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )
       
       ! Allocate the coefficients spline
       tmp1 = num_pts1+ 4*spline_degree1!*num_pts1! + spline_degree1 !- 1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2! + spline_degree2 !- 1
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

    case(2338)  !left: Dirichlet, right: Hermite, bottom: Hermite, Top: Hermite
       
       SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
       SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )
       
       ! Allocate the coefficients spline
       tmp1 = num_pts1+ 4*spline_degree1!*num_pts1! + spline_degree1 !- 1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2! + spline_degree2 !- 1
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

    case(2340) ! Hermite in all sides

       SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
       SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )

       ! Allocate the coefficients spline
       tmp1 = num_pts1+ 4*spline_degree1!*num_pts1! + spline_degree1 !- 1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2! + spline_degree2 !- 1
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)
       
    case default
       print*,'initialize_ad2d_interpolator: BC combination not implemented.'
    end select

    ! knots and coeff splines allocations 
    interpolator%coeff_splines(:,:) = 0.0_f64
    ! the minimun is to be of class C^0 everywhere on the knots
    ! i.e. each knot have multiplicity (spline_degree1+1) 
    ! so the maximun number of knots is num_pts1*(spline_degree1+1)
    SLL_CLEAR_ALLOCATE( interpolator%t1(1:num_pts1*(spline_degree1+1)),ierr)
    SLL_CLEAR_ALLOCATE( interpolator%t2(1:num_pts2*(spline_degree2+1)),ierr) 

  end subroutine !initialize_ad2d_interpolator

  !> Initialization of the boundary for interpolator arbitrary degree splines 2d.
  !> The parameters are
  !> @param[in] slope_left a 1d arrays contains values in the left in the direction eta1  
  !> @param[in] slope_right a 1d arrays contains values in the right in the direction eta1 
  !> @param[in] slope_bottom a 1d arrays contains values in the left in the direction eta2 
  !> @param[in] slope_top a 1d arrays contains values in the right in the direction eta2
  !> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_2d
  subroutine set_slope2d(&
       interpolator,&
       slope_left,&
       slope_right,&
       slope_bottom,&
       slope_top)

    use sll_module_arbitrary_degree_spline_interpolator_1d
    class(sll_arbitrary_degree_spline_interpolator_2d)    :: interpolator
    sll_real64, dimension(:),optional :: slope_left
    sll_real64, dimension(:),optional :: slope_right
    sll_real64, dimension(:),optional :: slope_bottom
    sll_real64, dimension(:),optional :: slope_top
    class(sll_arbitrary_degree_spline_interpolator_1d),pointer :: interp1d_bottom=> null()
    class(sll_arbitrary_degree_spline_interpolator_1d),pointer :: interp1d_top => null()
    sll_int32 :: sz_slope_bottom,sz_slope_top
    sll_int64 :: bc_selector
    sll_int32 :: num_pts1
    sll_int32 :: num_pts2
    sll_int32 :: bc_left
    sll_int32 :: bc_right
    sll_int32 :: bc_bottom
    sll_int32 :: bc_top

    num_pts1 = interpolator%num_pts1
    num_pts2 = interpolator%num_pts2
    bc_selector = interpolator%bc_selector
    bc_left  = interpolator%bc_left 
    bc_right = interpolator%bc_right 
    bc_bottom= interpolator%bc_bottom  
    bc_top   = interpolator%bc_top

    select case (bc_selector)
    case(0)
    case (9) ! dirichlet-left, dirichlet-right, periodic
    case (576) ! 3. periodic, dirichlet-bottom, dirichlet-top
    case (585) ! 4. dirichlet in all sides
    case (650) !left: Neumann, right: Dirichlet, bottom: Neumann, Top: Dirichlet
       !if ( present( slope_left)) then 
       interpolator%slope_left = 0.0
       interpolator%compute_slope_left= .FALSE.
       !end if
       if ( present( slope_right)) then 
          interpolator%slope_right = slope_right
          interpolator%compute_slope_right= .FALSE.
       end if

       interpolator%slope_bottom = 0.0_f64
       interpolator%compute_slope_bottom= .FALSE.
       

       if ( present( slope_top)) then 
          interpolator%compute_slope_top= .FALSE.

          
          sz_slope_top = size(slope_top)
          if ( sz_slope_top .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_top must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_top => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_top%compute_interpolants(&
               slope_top(1:sz_slope_top))
          
          interpolator%slope_top(1:sz_slope_top+2) = &
               interp1d_top%coeff_splines(1:sz_slope_top+2)
          call sll_delete(interp1d_top)
       else
          print*, 'problem with slope top in case 780'
       end if
    case(657) !left: Dirichlet, right: Neumann, bottom: Neumann, Top: Dirichlet 
       if ( present( slope_left)) then 
          interpolator%slope_left = slope_left
          interpolator%compute_slope_left= .FALSE.
       end if
       
       
       sz_slope_bottom = size(slope_bottom)
       interpolator%slope_right = 0.0_f64
       interpolator%compute_slope_right= .FALSE.
       
       interpolator%slope_bottom(1:sz_slope_bottom+2) = 0.0_f64
       interpolator%compute_slope_bottom = .FALSE.

       if ( present( slope_top)) then 
          interpolator%compute_slope_top= .FALSE.


          sz_slope_top = size(slope_top)
          if ( sz_slope_top .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_top must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_top => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_top%compute_interpolants(&
               slope_top(1:sz_slope_top))
          
          interpolator%slope_top(1:sz_slope_top+2) = &
               interp1d_top%coeff_splines(1:sz_slope_top+2)
          call sll_delete(interp1d_top)
       else
          print*, 'problem with slope top in case 780'
       end if

    case(780)  !left: Hermite, right: Dirichlet, bottom: Hermite, Top: Dirichlet
       if ( present( slope_left)) then 
          interpolator%slope_left = slope_left
          interpolator%compute_slope_left= .FALSE.
       end if
        if ( present( slope_right)) then 
          interpolator%slope_right = slope_right
          interpolator%compute_slope_right= .FALSE.
       end if
       
       if (present(slope_bottom)) then 
           sz_slope_bottom = size(slope_bottom)
          if ( sz_slope_bottom .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_bottom must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_bottom =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_bottom%compute_interpolants( &
               slope_bottom(1:sz_slope_bottom))
          
          interpolator%slope_bottom(1:sz_slope_bottom+2) = &
               interp1d_bottom%coeff_splines(1:sz_slope_bottom+2)
          call sll_delete(interp1d_bottom)
          interpolator%compute_slope_bottom = .FALSE.
       else
          print*, 'problem with slope bottom in case 780'
       end if
       if ( present( slope_top)) then 
          interpolator%compute_slope_top= .FALSE.


          sz_slope_top = size(slope_top)
          if ( sz_slope_top .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_top must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_top => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_top%compute_interpolants(&
               slope_top(1:sz_slope_top))
          
          interpolator%slope_top(1:sz_slope_top+2) = &
               interp1d_top%coeff_splines(1:sz_slope_top+2)
          call sll_delete(interp1d_top)
       else
          print*, 'problem with slope top in case 780'
       end if

    case(801)  !left: Dirichlet, right: Hermite, bottom: Hermite, Top: Dirichlet
        if ( present( slope_right)) then 
           interpolator%slope_right = slope_right
           interpolator%compute_slope_right= .FALSE.
       end if

        if ( present( slope_left)) then 
          interpolator%slope_left = slope_left
          interpolator%compute_slope_left= .FALSE.
       end if
       
       if (present(slope_bottom)) then 
           sz_slope_bottom = size(slope_bottom)
          if ( sz_slope_bottom .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_bottom must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_bottom =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_bottom%compute_interpolants( &
               slope_bottom(1:sz_slope_bottom))
          
          interpolator%slope_bottom(1:sz_slope_bottom+2) = &
               interp1d_bottom%coeff_splines(1:sz_slope_bottom+2)
          call sll_delete(interp1d_bottom)
          interpolator%compute_slope_bottom = .FALSE.
       else
          print*, 'problem with slope bottom in case 801'
       end if
       if ( present( slope_top)) then 
          interpolator%compute_slope_top= .FALSE.


          sz_slope_top = size(slope_top)
          if ( sz_slope_top .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_top must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_top => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_top%compute_interpolants(&
               slope_top(1:sz_slope_top))
          
          interpolator%slope_top(1:sz_slope_top+2) = &
               interp1d_top%coeff_splines(1:sz_slope_top+2)
          call sll_delete(interp1d_top)
       else
          print*, 'problem with slope top in case 801'
       end if

    case(804)  !left: Hermite, right: Hermite, bottom: Hermite, Top: Dirichlet

        if ( present( slope_right)) then 
           interpolator%slope_right = slope_right
           interpolator%compute_slope_right= .FALSE.
       end if

        if ( present( slope_left)) then 
          interpolator%slope_left = slope_left
          interpolator%compute_slope_left= .FALSE.
       end if
       
       if (present(slope_bottom)) then 
           sz_slope_bottom = size(slope_bottom)
          if ( sz_slope_bottom .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_bottom must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_bottom =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_bottom%compute_interpolants( &
               slope_bottom(1:sz_slope_bottom))
          
          interpolator%slope_bottom(1:sz_slope_bottom+2) = &
               interp1d_bottom%coeff_splines(1:sz_slope_bottom+2)
          call sll_delete(interp1d_bottom)
          interpolator%compute_slope_bottom = .FALSE.
       else
          print*, 'problem with slope bottom in case 801'
       end if
       if ( present( slope_top)) then 
          interpolator%compute_slope_top= .FALSE.


          sz_slope_top = size(slope_top)
          if ( sz_slope_top .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_top must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_top => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_top%compute_interpolants(&
               slope_top(1:sz_slope_top))
          
          interpolator%slope_top(1:sz_slope_top+2) = &
               interp1d_top%coeff_splines(1:sz_slope_top+2)
          call sll_delete(interp1d_top)
       else
          print*, 'problem with slope top in case 801'
       end if


    case(1098)  !left: Neumann, right: Dirichlet, bottom: Dirichlet, Top: Neumann
       interpolator%slope_left = 0.0_f64
       interpolator%compute_slope_left= .FALSE.
    
       if ( present( slope_right)) then 
          interpolator%slope_right = slope_right
          interpolator%compute_slope_right= .FALSE.
       end if
       
       sz_slope_top= size(slope_top)
       interpolator%slope_top(1:sz_slope_top+2) = 0.0_f64
       interpolator%compute_slope_top = .FALSE.

       if (present(slope_bottom)) then 
           sz_slope_bottom = size(slope_bottom)
          if ( sz_slope_bottom .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_bottom must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_bottom =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_bottom%compute_interpolants( &
               slope_bottom(1:sz_slope_bottom))
          
          interpolator%slope_bottom(1:sz_slope_bottom+2) = &
               interp1d_bottom%coeff_splines(1:sz_slope_bottom+2)
          call sll_delete(interp1d_bottom)
          interpolator%compute_slope_bottom = .FALSE.
       else
          print*, 'problem with slope bottom in case 2124'
       end if
       if ( present( slope_top)) then 
          interpolator%compute_slope_top= .FALSE.


          sz_slope_top = size(slope_top)
          if ( sz_slope_top .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_top must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_top => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_top%compute_interpolants(&
               slope_top(1:sz_slope_top))
          
          interpolator%slope_top(1:sz_slope_top+2) = &
               interp1d_top%coeff_splines(1:sz_slope_top+2)
          call sll_delete(interp1d_top)
       else
          print*, 'problem with slope top in case 2124'
       end if
    case(1105)  !left: Dirichlet, right: Neumann, bottom: Dirichlet, Top: Neumann
       if ( present( slope_right)) then 
          interpolator%slope_right = slope_right
          interpolator%compute_slope_right= .FALSE.
       end if
       interpolator%slope_left = 0.0_f64
       interpolator%compute_slope_left= .FALSE.
       
       sz_slope_top = size(slope_top)
       interpolator%compute_slope_top= .FALSE.
       interpolator%slope_top(1:sz_slope_top+2) = 0.0_f64

       if (present(slope_bottom)) then 
           sz_slope_bottom = size(slope_bottom)
          if ( sz_slope_bottom .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_bottom must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_bottom =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_bottom%compute_interpolants( &
               slope_bottom(1:sz_slope_bottom))
          
          interpolator%slope_bottom(1:sz_slope_bottom+2) = &
               interp1d_bottom%coeff_splines(1:sz_slope_bottom+2)
          call sll_delete(interp1d_bottom)
          interpolator%compute_slope_bottom = .FALSE.
       else
          print*, 'problem with slope bottom in case 2145'
       end if
    case(1170)  !left: Neumann, right: Neumann, bottom: Neuman, Top: Neumann
       
       interpolator%slope_right = 0.0_f64
       interpolator%compute_slope_right= .FALSE.
       
       interpolator%slope_left = 0.0_f64
       interpolator%compute_slope_left= .FALSE.
       
       sz_slope_top = size(slope_top)
       interpolator%compute_slope_top= .FALSE.
       interpolator%slope_top(1:sz_slope_top+2) = 0.0_f64
          
       sz_slope_bottom = size(slope_bottom)
       interpolator%slope_bottom(1:sz_slope_bottom+2) = 0.0_f64
       interpolator%compute_slope_bottom = .FALSE.

    case(2338)  !left: Dirichlet, right: Hermite, bottom: Hermite, Top: Hermite
       if ( present( slope_right)) then 
          interpolator%slope_right = slope_right
          interpolator%compute_slope_right= .FALSE.
       end if
       if ( present( slope_left)) then 
          interpolator%slope_left = slope_left
          interpolator%compute_slope_left= .FALSE.
       end if
       
       if ( present( slope_top)) then 
          interpolator%compute_slope_top= .FALSE.

          
          sz_slope_top = size(slope_top)
          if ( sz_slope_top .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_top must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_top => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_top%compute_interpolants(&
               slope_top(1:sz_slope_top))
          
          interpolator%slope_top(1:sz_slope_top+2) = &
               interp1d_top%coeff_splines(1:sz_slope_top+2)
          call sll_delete(interp1d_top)
       else
          print*, 'problem with slope top in case 2145'
       end if

       if (present(slope_bottom)) then 
           sz_slope_bottom = size(slope_bottom)
          if ( sz_slope_bottom .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_bottom must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_bottom =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_bottom%compute_interpolants( &
               slope_bottom(1:sz_slope_bottom))
          
          interpolator%slope_bottom(1:sz_slope_bottom+2) = &
               interp1d_bottom%coeff_splines(1:sz_slope_bottom+2)
          call sll_delete(interp1d_bottom)
          interpolator%compute_slope_bottom = .FALSE.
       else
          print*, 'problem with slope bottom in case 2145'
       end if

    case(2145)  !left: Dirichlet, right: Hermite, bottom: Dirichlet, Top: Hermite
       if ( present( slope_right)) then 
          interpolator%slope_right = slope_right
          interpolator%compute_slope_right= .FALSE.
       end if
        if ( present( slope_left)) then 
          interpolator%slope_left = slope_left
          interpolator%compute_slope_left= .FALSE.
       end if
       
       if ( present( slope_top)) then 
          interpolator%compute_slope_top= .FALSE.

          
          sz_slope_top = size(slope_top)
          if ( sz_slope_top .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_top must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_top => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_top%compute_interpolants(&
               slope_top(1:sz_slope_top))
          
          interpolator%slope_top(1:sz_slope_top+2) = &
               interp1d_top%coeff_splines(1:sz_slope_top+2)
          call sll_delete(interp1d_top)
       else
          print*, 'problem with slope top in case 2145'
       end if

       if (present(slope_bottom)) then 
           sz_slope_bottom = size(slope_bottom)
          if ( sz_slope_bottom .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_bottom must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_bottom =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_bottom%compute_interpolants( &
               slope_bottom(1:sz_slope_bottom))
          
          interpolator%slope_bottom(1:sz_slope_bottom+2) = &
               interp1d_bottom%coeff_splines(1:sz_slope_bottom+2)
          call sll_delete(interp1d_bottom)
          interpolator%compute_slope_bottom = .FALSE.
       else
          print*, 'problem with slope bottom in case 2145'
       end if

       
    case(2124)  !left: Hermite, right: Dirichlet, bottom: Dirichlet, Top: Hermite
       
       if ( present( slope_left)) then 
          interpolator%slope_left = slope_left
          interpolator%compute_slope_left= .FALSE.
       end if
        if ( present( slope_right)) then 
          interpolator%slope_right = slope_right
          interpolator%compute_slope_right= .FALSE.
       end if
       
       if ( present( slope_top)) then 
          interpolator%compute_slope_top= .FALSE.


          sz_slope_top = size(slope_top)
          if ( sz_slope_top .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_top must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_top => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_top%compute_interpolants(&
               slope_top(1:sz_slope_top))
          
          interpolator%slope_top(1:sz_slope_top+2) = &
               interp1d_top%coeff_splines(1:sz_slope_top+2)

          call sll_delete(interp1d_top)
       else
          print*, 'problem with slope top in case 2124'
       end if

       if (present(slope_bottom)) then 
           sz_slope_bottom = size(slope_bottom)
          if ( sz_slope_bottom .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_bottom must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_bottom =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_bottom%compute_interpolants( &
               slope_bottom(1:sz_slope_bottom))
          
          interpolator%slope_bottom(1:sz_slope_bottom+2) = &
               interp1d_bottom%coeff_splines(1:sz_slope_bottom+2)
          call sll_delete(interp1d_bottom)
          interpolator%compute_slope_bottom = .FALSE.
       else
          print*, 'problem with slope bottom in case 2124'
       end if
       if ( present( slope_top)) then 
          interpolator%compute_slope_top= .FALSE.


          sz_slope_top = size(slope_top)
          if ( sz_slope_top .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_top must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_top => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_top%compute_interpolants(&
               slope_top(1:sz_slope_top))
          
          interpolator%slope_top(1:sz_slope_top+2) = &
               interp1d_top%coeff_splines(1:sz_slope_top+2)
          call sll_delete(interp1d_top)
       else
          print*, 'problem with slope top in case 2124'
       end if
      
    case(2148)  !left:Hermite , right: Hermite, bottom: Dirichlet, Top: Hermite  
       if ( present( slope_left)) then 
          interpolator%slope_left = slope_left
          interpolator%compute_slope_left= .FALSE.
       end if
       if ( present( slope_right)) then 
          interpolator%slope_right = slope_right
          interpolator%compute_slope_right= .FALSE.
       end if
       
       if ( present( slope_top)) then 
          interpolator%compute_slope_top= .FALSE.


          sz_slope_top = size(slope_top)
          if ( sz_slope_top .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_top must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_top => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_top%compute_interpolants(&
               slope_top(1:sz_slope_top))
          
          interpolator%slope_top(1:sz_slope_top+2) = &
               interp1d_top%coeff_splines(1:sz_slope_top+2)

          call sll_delete(interp1d_top)
       else
          print*, 'problem with slope top in case 2124'
       end if

       if (present(slope_bottom)) then 
           sz_slope_bottom = size(slope_bottom)
          if ( sz_slope_bottom .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_bottom must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_bottom =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_bottom%compute_interpolants( &
               slope_bottom(1:sz_slope_bottom))
          
          interpolator%slope_bottom(1:sz_slope_bottom+2) = &
               interp1d_bottom%coeff_splines(1:sz_slope_bottom+2)
          call sll_delete(interp1d_bottom)
          interpolator%compute_slope_bottom = .FALSE.
       else
          print*, 'problem with slope bottom in case 2124'
       end if
       if ( present( slope_top)) then 
          interpolator%compute_slope_top= .FALSE.


          sz_slope_top = size(slope_top)
          if ( sz_slope_top .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_top must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_top => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_top%compute_interpolants(&
               slope_top(1:sz_slope_top))
          
          interpolator%slope_top(1:sz_slope_top+2) = &
               interp1d_top%coeff_splines(1:sz_slope_top+2)
          call sll_delete(interp1d_top)
       else
          print*, 'problem with slope top in case 2124'
       end if
    case(2316)  !left: Hermite, right: Dirichlet, bottom: Hermite, Top: Hermite
       if ( present( slope_right)) then 
          interpolator%slope_right = slope_right
          interpolator%compute_slope_right= .FALSE.
       end if
       if ( present( slope_left)) then 
          interpolator%slope_left = slope_left
          interpolator%compute_slope_left= .FALSE.
       end if
       
       if ( present( slope_top)) then 
          interpolator%compute_slope_top= .FALSE.


          sz_slope_top = size(slope_top)
          if ( sz_slope_top .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_top must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_top => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_top%compute_interpolants(&
               slope_top(1:sz_slope_top))
          
          interpolator%slope_top(1:sz_slope_top+2) = &
               interp1d_top%coeff_splines(1:sz_slope_top+2)
          call sll_delete(interp1d_top)
       else
          print*, 'problem with slope top in case 2340'
       end if
       
       if (present(slope_bottom)) then 
           sz_slope_bottom = size(slope_bottom)
          if ( sz_slope_bottom .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_bottom must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_bottom =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_bottom%compute_interpolants( &
               slope_bottom(1:sz_slope_bottom))
          
          interpolator%slope_bottom(1:sz_slope_bottom+2) = &
               interp1d_bottom%coeff_splines(1:sz_slope_bottom+2)
          call sll_delete(interp1d_bottom)
          interpolator%compute_slope_bottom = .FALSE.
       else
          print*, 'problem with slope bottom in case 2340'
       end if
    case(2340) ! Hermite in al sides

       if ( present( slope_right)) then 
          interpolator%slope_right = slope_right
          interpolator%compute_slope_right= .FALSE.
       end if
        if ( present( slope_left)) then 
          interpolator%slope_left = slope_left
          interpolator%compute_slope_left= .FALSE.
       end if
       
       if ( present( slope_top)) then 
          interpolator%compute_slope_top= .FALSE.


          sz_slope_top = size(slope_top)
          if ( sz_slope_top .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_top must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_top => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_top%compute_interpolants(&
               slope_top(1:sz_slope_top))
          
          interpolator%slope_top(1:sz_slope_top+2) = &
               interp1d_top%coeff_splines(1:sz_slope_top+2)
          call sll_delete(interp1d_top)
       else
          print*, 'problem with slope top in case 2340'
       end if
       
       if (present(slope_bottom)) then 
           sz_slope_bottom = size(slope_bottom)
          if ( sz_slope_bottom .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_bottom must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_bottom =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_bottom%compute_interpolants( &
               slope_bottom(1:sz_slope_bottom))
          
          interpolator%slope_bottom(1:sz_slope_bottom+2) = &
               interp1d_bottom%coeff_splines(1:sz_slope_bottom+2)
          call sll_delete(interp1d_bottom)
          interpolator%compute_slope_bottom = .FALSE.
       else
          print*, 'problem with slope bottom in case 2340'
       end if
       
    case default
       print*,'initialize_ad2d_interpolator: BC combination not implemented.'
    end select
    
  end subroutine set_slope2d
  
  
  !> Initialization of the boundary for interpolator arbitrary degree splines 2d.
  !> The parameters are
  !> @param[in] value_left a 1d arrays contains values in the left in the direction eta1  
  !> @param[in] value_right a 1d arrays contains values in the right in the direction eta1 
  !> @param[in] value_bottom a 1d arrays contains values in the left in the direction eta2 
  !> @param[in]  value_top a 1d arrays contains values in the right in the direction eta2
  !> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_2d

  subroutine set_boundary_value2d(&
       interpolator,&
       value_left,&
       value_right,&
       value_bottom,&
       value_top)

    use sll_module_arbitrary_degree_spline_interpolator_1d
    class(sll_arbitrary_degree_spline_interpolator_2d)    :: interpolator
    sll_real64, dimension(:),optional :: value_left
    sll_real64, dimension(:),optional :: value_right
    sll_real64, dimension(:),optional :: value_bottom
    sll_real64, dimension(:),optional :: value_top
    class(sll_arbitrary_degree_spline_interpolator_1d),pointer :: interp1d_left => null()
    class(sll_arbitrary_degree_spline_interpolator_1d),pointer :: interp1d_right => null()
    class(sll_arbitrary_degree_spline_interpolator_1d),pointer :: interp1d_bottom=> null()
    class(sll_arbitrary_degree_spline_interpolator_1d),pointer :: interp1d_top => null()
    sll_int32 :: sz_value_left,sz_value_right,sz_value_bottom,sz_value_top
    sll_int64 :: bc_selector
    sll_int32 :: num_pts1
    sll_int32 :: num_pts2
    sll_int32 :: bc_left
    sll_int32 :: bc_right
    sll_int32 :: bc_bottom
    sll_int32 :: bc_top

    num_pts1 = interpolator%num_pts1
    num_pts2 = interpolator%num_pts2
    bc_selector = interpolator%bc_selector
    bc_left  = interpolator%bc_left 
    bc_right = interpolator%bc_right 
    bc_bottom= interpolator%bc_bottom  
    bc_top   = interpolator%bc_top

    select case (bc_selector)
    case(0)
       
    case (9) ! dirichlet-left, dirichlet-right, periodic
       if (present(value_left)) then 
          sz_value_left = size(value_left)
          if ( sz_value_left .ne. interpolator%num_pts2 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, 'value_left must have the size of numbers of pts in direction 2 '
             stop
          end if
          interp1d_left => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts2, &
               interpolator%eta2_min, &
               interpolator%eta2_max, &
               interpolator%bc_bottom, &
               interpolator%bc_top, &
               interpolator%spline_degree2 )
          
          call interp1d_left%compute_interpolants( &
               value_left(1:sz_value_left))
          
          interpolator%value_left(1:sz_value_left) = &
               interp1d_left%coeff_splines(1:sz_value_left)
          call sll_delete(interp1d_left)

          interpolator%compute_value_left = .FALSE.
       else
          interpolator%value_left(:) = 0.0_f64
       end if

       if (present(value_right)) then 
          sz_value_right = size(value_right)
          if ( sz_value_right .ne. interpolator%num_pts2 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' value_right must have the size of numbers of pts in direction 2 '
             stop
          end if
          
          interp1d_right => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts2, &
               interpolator%eta2_min, &
               interpolator%eta2_max, &
               interpolator%bc_bottom, &
               interpolator%bc_top, &
               interpolator%spline_degree2 )
          
          call interp1d_right%compute_interpolants( &
               value_right(1:sz_value_right))
          
          interpolator%value_right(1:sz_value_right) = &
               interp1d_right%coeff_splines(1:sz_value_right)
          call sll_delete(interp1d_right)
          interpolator%compute_value_right = .FALSE.
       else
          interpolator%value_right(:) = 0.0_f64
       end if
    case (576) ! 3. periodic, dirichlet-bottom, dirichlet-top
       
       if (present(value_bottom)) then 
          sz_value_bottom = size(value_bottom)
          if ( sz_value_bottom .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' value_bottom must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_bottom =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_bottom%compute_interpolants( &
               value_bottom(1:sz_value_bottom))
          
          interpolator%value_bottom(1:sz_value_bottom) = &
               interp1d_bottom%coeff_splines(1:sz_value_bottom)
          call sll_delete(interp1d_bottom)
          interpolator%compute_value_bottom = .FALSE.
       else
          interpolator%value_bottom(:) = 0.0_f64
       end if
       
       if (present(value_top)) then 
          sz_value_top = size(value_top)
          if ( sz_value_top .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' value_top must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_top => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_top%compute_interpolants(&
               value_top(1:sz_value_top))
          
          interpolator%value_top(1:sz_value_top) = &
               interp1d_top%coeff_splines(1:sz_value_top)
          call sll_delete(interp1d_top)
          interpolator%compute_value_top = .FALSE.
       else
          interpolator%value_top(:) = 0.0_f64
       end if
    case (585) ! 4. dirichlet in all sides
       
       if (present(value_left)) then 
          sz_value_left = size(value_left)
          if ( sz_value_left .ne. interpolator%num_pts2 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' value_left must have the size of numbers of pts in direction 2 '
             stop
          end if

          interp1d_left =>  new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts2, &
               interpolator%eta2_min, &
               interpolator%eta2_max, &
               interpolator%bc_bottom, &
               interpolator%bc_top, &
               interpolator%spline_degree2)

          call set_values_at_boundary1d(&
               interp1d_left,&
               value_left(1),&
               value_left(sz_value_left))
          
          call interp1d_left%compute_interpolants( &
               value_left(1:sz_value_left))
          
          interpolator%value_left(1:sz_value_left) = &
               interp1d_left%coeff_splines(1:sz_value_left)
          call sll_delete(interp1d_left)
          interpolator%compute_value_left = .FALSE.
       else
          interpolator%value_left(:) = 0.0_f64
       end if
       
       if (present(value_right)) then 
          sz_value_right = size(value_right)
          if ( sz_value_right .ne. interpolator%num_pts2 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' value_right must have the size of numbers of pts in direction 2 '
             stop
          end if
          
          interp1d_right => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts2, &
               interpolator%eta2_min, &
               interpolator%eta2_max, &
               interpolator%bc_bottom, &
               interpolator%bc_top, &
               interpolator%spline_degree2)
          
          call set_values_at_boundary1d(&
               interp1d_right,&
               value_right(1),&
               value_right(sz_value_right))
          
          call interp1d_right%compute_interpolants( &
               value_right(1:sz_value_right))
          
          interpolator%value_right(1:sz_value_right) = &
               interp1d_right%coeff_splines(1:sz_value_right)
          call sll_delete(interp1d_right)
          interpolator%compute_value_right = .FALSE.
       else
          interpolator%value_right(:) = 0.0_f64
       end if
       
       
       if (present(value_bottom)) then 
          sz_value_bottom = size(value_bottom)
          if ( sz_value_bottom .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' value_bottom must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_bottom=> new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1)

          call set_values_at_boundary1d(&
               interp1d_bottom,&
               value_bottom(1),&
               value_bottom(sz_value_bottom))
          
          call interp1d_bottom%compute_interpolants( &
               value_bottom(1:sz_value_bottom))
          
          interpolator%value_bottom(1:sz_value_bottom) = &
               interp1d_bottom%coeff_splines(1:sz_value_bottom)
          call sll_delete(interp1d_bottom)
          interpolator%compute_value_bottom = .FALSE.
       else
          interpolator%value_bottom(:) = 0.0_f64
       end if
       
       if (present(value_top)) then 
          sz_value_top = size(value_top)
          if ( sz_value_top .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' value_top must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          interp1d_top => new_arbitrary_degree_1d_interpolator(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1)
          
          call set_values_at_boundary1d(&
               interp1d_top,&
               value_top(1),&
               value_top(sz_value_top))
          
          call interp1d_top%compute_interpolants( &
               value_top(1:sz_value_top))
          
          interpolator%value_top(1:sz_value_top) = &
               interp1d_top%coeff_splines(1:sz_value_top)
          call sll_delete(interp1d_top)
          interpolator%compute_value_top = .FALSE.
       else
          interpolator%value_top(:) = 0.0_f64
       end if
    case(650) !left: Neumann, right: Dirichlet, bottom: Neumann, Top: Dirichlet
    case(657) !left: Dirichlet, right: Neumann, bottom: Neumann, Top: Dirichlet 
    case(780)  !left: Hermite, right: Dirichlet, bottom: Hermite, Top: Dirichlet
    case(801)  !left: Dirichlet, right: Hermite, bottom: Hermite, Top: Dirichlet
    case(804)  !left: Hermite, right: Hermite, bottom: Hermite, Top: Dirichlet
    case(1098)  !left: Neumann, right: Dirichlet, bottom: Dirichlet, Top: Neumann
    case(1105)  !left: Dirichlet, right: Neumann, bottom: Dirichlet, Top: Neumann
    case(1170)  !left: Neumann, right: Neumann, bottom: Neuman, Top: Neumann
    case(2124)  !left: Hermite, right: Dirichlet, bottom: Dirichlet, Top: Hermite
    case(2145)  !left: Dirichlet, right: Hermite, bottom: Dirichlet, Top: Hermite
    case(2148)  !left:Hermite , right: Hermite, bottom: Dirichlet, Top: Hermite
    case(2316)  !left: Hermite, right: Dirichlet, bottom: Hermite, Top: Hermite
    case(2338)  !left: Dirichlet, right: Hermite, bottom: Hermite, Top: Hermite
    case(2340) ! Hermite in all sides
    case default
       print*,'initialize_ad2d_interpolator: BC combination not implemented.'
    end select

  end subroutine set_boundary_value2d


  ! -------------------------------------------------------------
  !  subroutine initializing the coefficients of splines
  !  in the cas of linearization of them i.e. if we have 
  !  a table in 1d corresponding of coefficients 2d
  ! 
  ! -------------------------------------------------------------
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
   interpolator, &
   coeffs_1d, &
   coeffs_2d,&
   coeff2d_size1,&
   coeff2d_size2,&
   knots1,&
   size_knots1,&
   knots2,&
   size_knots2)

   class(sll_arbitrary_degree_spline_interpolator_2d), intent(inout)  :: interpolator
   sll_real64, dimension(:)  , intent(in), optional :: coeffs_1d
   sll_real64, dimension(:,:), intent(in), optional :: coeffs_2d
   ! size coeffs 2D 
   sll_int32, intent(in), optional :: coeff2d_size1
   sll_int32, intent(in), optional :: coeff2d_size2
   sll_real64, dimension(:), intent(in), optional   :: knots1
   sll_real64, dimension(:), intent(in), optional   :: knots2
   sll_int32, intent(in), optional :: size_knots1
   sll_int32, intent(in), optional :: size_knots2

   ! Local variables
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
   sll_int32   :: sz_derivative1,sz_derivative2

   
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
         
         interpolator%size_coeffs1 =  num_cells1 + sp_deg1 + 1
         interpolator%size_coeffs2 =  num_cells2 + sp_deg2 + 1 
         interpolator%size_t1      =  2*sp_deg1 + num_cells1 +1 +1
         interpolator%size_t2      =  2*sp_deg2 + num_cells2 +1 +1
         
         if ( size( coeffs_1d,1) .ne. num_cells1*num_cells2) then
            print*, 'Problem in set_coefficients in arbitrary_degree_spline_2d'
            print*, ' Problem with the size coeffs_1d must have the size equal to '
            print*, ' num_cells1*num_cells2=', num_cells1*num_cells2
            stop
         end if
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
         
         
         ! ------------------------------------------------------------
         !   reorganization of spline coefficients 1D in coefficients 2D 
         ! ------------------------------------------------------------
         
         
         do i = 1,num_cells1
            do j = 1,num_cells2
               interpolator%coeff_splines(i,j) = &
                    coeffs_1d( i + num_cells1 *(j-1) )
            end do
         end do
         
         do j = 1, sp_deg2 + 1
            do i = 1,num_cells1
               
               interpolator%coeff_splines(i ,num_cells2 + j ) = &
                    coeffs_1d(i+num_cells1*(j-1))
            end do
         end do
         do i = 1, sp_deg1 + 1
            do j = 1,num_cells2
               
               interpolator%coeff_splines(num_cells1 + i ,j) = &
                    coeffs_1d(i+num_cells1 *(j-1) )
            end do
         end do
         do i= 1,sp_deg1 + 1
            do j=1,sp_deg2 + 1
               
               interpolator%coeff_splines(num_cells1 +  i ,num_cells2 + j) = &
                    interpolator%coeff_splines(i,j)
            end do
         end do

        
      ! ------------------------------------------------------------
      case (9) ! 2. dirichlet-left, dirichlet-right, periodic
         
         
         
         interpolator%size_coeffs1 =  num_cells1 + sp_deg1
         interpolator%size_coeffs2 =  num_cells2 + sp_deg2 + 1
         interpolator%size_t1      =  2*sp_deg1 + num_cells1 + 1
         interpolator%size_t2      =  2*sp_deg2 + num_cells2 + 1 + 1
         nb_spline_eta1            =  num_cells1 + sp_deg1 - 2
         nb_spline_eta2            =  num_cells2
         
         if ( size( coeffs_1d,1) .ne. (num_cells1 + sp_deg1 - 2)*num_cells2) then
            print*, 'Problem in set_coefficients in arbitrary_degree_spline_2d'
            print*, ' Problem with the size coeffs_1d must have the size equal to '
            print*, ' (num_cells1 + sp_deg1 - 2)*num_cells2=', &
                 (num_cells1 + sp_deg1 - 2)*num_cells2
            stop
         end if
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
               interpolator%coeff_splines(i+1,j) = &
                    coeffs_1d(i+nb_spline_eta1*(j-1))
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
         ! ------------------------------------------------------------
      case(576)!3. periodic, dirichlet-bottom, dirichlet-top
       
         
         interpolator%size_coeffs1 =  num_cells1 + sp_deg1 + 1
         interpolator%size_coeffs2 =  num_cells2 + sp_deg2
         interpolator%size_t1      = 2*sp_deg1 + num_cells1 + 1 + 1
         interpolator%size_t2      = 2*sp_deg2 + num_cells2 + 1
         nb_spline_eta1            = num_cells1
         nb_spline_eta2            = num_cells2 + sp_deg2 - 2
       

       if ( size( coeffs_1d,1) .ne. num_cells1*( num_cells2 + sp_deg2 - 2)) then
          print*, 'Problem in set_coefficients in arbitrary_degree_spline_2d'
          print*, ' Problem with the size coeffs_1d must have the size equal to '
          print*, ' num_cells1*( num_cells2 + sp_deg2 - 2)=',&
               num_cells1*( num_cells2 + sp_deg2 - 2)
          stop
       end if
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
         
         if(size(coeffs_1d,1).ne.(num_cells1 + sp_deg1-2)*(num_cells2+sp_deg2-2))then
            print*, 'Problem in set_coefficients in arbitrary_degree_spline_2d'
            print*, ' Problem with the size coeffs_1d must have the size equal to '
            print*, ' (num_cells1 + sp_deg1 - 2)*( num_cells2 + sp_deg2 - 2)=',&
                 (num_cells1 + sp_deg1 - 2)*( num_cells2 + sp_deg2 - 2)
            stop
         end if
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
         ! achtung ! normaly interpolator%slope_left(:) and interpolator%value_right(:)
         ! achtung ! normaly interpolator%value_bottom(:) and interpolator%value_top(:)

         interpolator%coeff_splines(:,:) = 0.0_8
         ! allocation coefficient spline
         do i = 1,nb_spline_eta1
            do j = 1,nb_spline_eta2
               
               interpolator%coeff_splines(i+1,j+1) = &
                    coeffs_1d( i + nb_spline_eta1 *(j-1))
            end do
         end do
         ! ------------------------------------------------------------

      case (650) !left: Neumann, right: Dirichlet, bottom: Neumann, Top: Dirichlet
         interpolator%size_coeffs1=  num_cells1 + sp_deg1 +1
         interpolator%size_coeffs2=  num_cells2 + sp_deg2 +1
         interpolator%size_t1 = 2*sp_deg1 + num_cells1 + 1
         interpolator%size_t2 = 2*sp_deg2 + num_cells2 + 1
         nb_spline_eta1 = num_cells1 + sp_deg1 +1
         nb_spline_eta2 = num_cells2 + sp_deg2 +1
         
         if(size(coeffs_1d,1).ne.(num_cells1 + sp_deg1+1)*(num_cells2+sp_deg2+1))then
            print*, 'Problem in set_coefficients in arbitrary_degree_spline_2d'
            print*, ' Problem with the size coeffs_1d must have the size equal to '
            print*, ' (num_cells1 + sp_deg1 +1)*( num_cells2 + sp_deg2 +1)=',&
                 (num_cells1 + sp_deg1 +1)*( num_cells2 + sp_deg2+1)
            stop
         end if
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
         ! achtung ! normaly interpolator%slope_left(:) and interpolator%value_right(:)
         ! achtung ! normaly interpolator%value_bottom(:) and interpolator%value_top(:)

         interpolator%coeff_splines(:,:) = 0.0_8
         ! allocation coefficient spline
         do i = 1,nb_spline_eta1
            do j = 1,nb_spline_eta2
               
               interpolator%coeff_splines(i+1,j+1) = &
                    coeffs_1d( i + nb_spline_eta1 *(j-1))
            end do
         end do
         ! ------------------------------------------------------------
      case(657) !left: Dirichlet, right: Neumann, bottom: Neumann, Top: Dirichlet 
         interpolator%size_coeffs1=  num_cells1 + sp_deg1 +1
         interpolator%size_coeffs2=  num_cells2 + sp_deg2 +1
         interpolator%size_t1 = 2*sp_deg1 + num_cells1 + 1
         interpolator%size_t2 = 2*sp_deg2 + num_cells2 + 1
         nb_spline_eta1 = num_cells1 + sp_deg1 +1
         nb_spline_eta2 = num_cells2 + sp_deg2 +1
         
         if(size(coeffs_1d,1).ne.(num_cells1 + sp_deg1+1)*(num_cells2+sp_deg2+1))then
            print*, 'Problem in set_coefficients in arbitrary_degree_spline_2d'
            print*, ' Problem with the size coeffs_1d must have the size equal to '
            print*, ' (num_cells1 + sp_deg1 +1)*( num_cells2 + sp_deg2 +1)=',&
                 (num_cells1 + sp_deg1 +1)*( num_cells2 + sp_deg2+1)
            stop
         end if
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
         ! achtung ! normaly interpolator%slope_left(:) and interpolator%value_right(:)
         ! achtung ! normaly interpolator%value_bottom(:) and interpolator%value_top(:)

         interpolator%coeff_splines(:,:) = 0.0_8
         ! allocation coefficient spline
         do i = 1,nb_spline_eta1
            do j = 1,nb_spline_eta2
               
               interpolator%coeff_splines(i+1,j+1) = &
                    coeffs_1d( i + nb_spline_eta1 *(j-1))
            end do
         end do
         ! ------------------------------------------------------------
      case(780)  !left: Hermite, right: Dirichlet, bottom: Hermite, Top: Dirichlet
         interpolator%size_coeffs1=  num_cells1 + sp_deg1 +1
         interpolator%size_coeffs2=  num_cells2 + sp_deg2 +1
         interpolator%size_t1 = 2*sp_deg1 + num_cells1 + 1
         interpolator%size_t2 = 2*sp_deg2 + num_cells2 + 1
         nb_spline_eta1 = num_cells1 + sp_deg1 +1
         nb_spline_eta2 = num_cells2 + sp_deg2 +1
         
         if(size(coeffs_1d,1).ne.(num_cells1 + sp_deg1+1)*(num_cells2+sp_deg2+1))then
            print*, 'Problem in set_coefficients in arbitrary_degree_spline_2d'
            print*, ' Problem with the size coeffs_1d must have the size equal to '
            print*, ' (num_cells1 + sp_deg1 +1)*( num_cells2 + sp_deg2 +1)=',&
                 (num_cells1 + sp_deg1 +1)*( num_cells2 + sp_deg2+1)
            stop
         end if
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
         ! achtung ! normaly interpolator%slope_left(:) and interpolator%value_right(:)
         ! achtung ! normaly interpolator%value_bottom(:) and interpolator%value_top(:)

         interpolator%coeff_splines(:,:) = 0.0_8
         ! allocation coefficient spline
         do i = 1,nb_spline_eta1
            do j = 1,nb_spline_eta2
               
               interpolator%coeff_splines(i+1,j+1) = &
                    coeffs_1d( i + nb_spline_eta1 *(j-1))
            end do
         end do
         ! ------------------------------------------------------------
      case(1098)  !left: Neumann, right: Dirichlet, bottom: Dirichlet, Top: Neumann
         interpolator%size_coeffs1=  num_cells1 + sp_deg1 +1
         interpolator%size_coeffs2=  num_cells2 + sp_deg2 +1
         interpolator%size_t1 = 2*sp_deg1 + num_cells1 + 1
         interpolator%size_t2 = 2*sp_deg2 + num_cells2 + 1
         nb_spline_eta1 = num_cells1 + sp_deg1 +1
         nb_spline_eta2 = num_cells2 + sp_deg2 +1
         
         if(size(coeffs_1d,1).ne.(num_cells1 + sp_deg1+1)*(num_cells2+sp_deg2+1))then
            print*, 'Problem in set_coefficients in arbitrary_degree_spline_2d'
            print*, ' Problem with the size coeffs_1d must have the size equal to '
            print*, ' (num_cells1 + sp_deg1 +1)*( num_cells2 + sp_deg2 +1)=',&
                 (num_cells1 + sp_deg1 +1)*( num_cells2 + sp_deg2+1)
            stop
         end if
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
         ! achtung ! normaly interpolator%slope_left(:) and interpolator%value_right(:)
         ! achtung ! normaly interpolator%value_bottom(:) and interpolator%value_top(:)

         interpolator%coeff_splines(:,:) = 0.0_8
         ! allocation coefficient spline
         do i = 1,nb_spline_eta1
            do j = 1,nb_spline_eta2
               
               interpolator%coeff_splines(i+1,j+1) = &
                    coeffs_1d( i + nb_spline_eta1 *(j-1))
            end do
         end do
         ! ------------------------------------------------------------

      case(2145)  !left: Dirichlet, right: Hermite, bottom: Dirichlet, Top: Hermite  

         sz_derivative1 = 2
         sz_derivative2 = 2
         interpolator%size_coeffs1=  num_cells1 + sp_deg1 + 1
         interpolator%size_coeffs2=  num_cells2 + sp_deg2 + 1
         interpolator%size_t1 = 2*sp_deg1 + num_cells1 + 1
         interpolator%size_t2 = 2*sp_deg2 + num_cells2 + 1
         nb_spline_eta1 = num_cells1 + sp_deg1 + 1
         nb_spline_eta2 = num_cells2 + sp_deg2 + 1
         
         if(size(coeffs_1d,1).ne.(num_cells1+sp_deg1+1)*(num_cells2+sp_deg2+1))then
            print*, 'Problem in set_coefficients in arbitrary_degree_spline_2d'
            print*, ' Problem with the size coeffs_1d must have the size equal to '
            print*, ' (num_cells1 + sp_deg1 +1)*( num_cells2 + sp_deg2 +1)=',&
                 (num_cells1 + sp_deg1 +1)*( num_cells2 + sp_deg2 +1)
            stop
         end if
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
         ! achtung ! normaly interpolator%slope_left(:) and interpolator%value_right(:)
         ! achtung ! normaly interpolator%value_bottom(:) and interpolator%value_top(:)

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
         
         if(size(coeffs_1d,1).ne.(num_cells1 + sp_deg1+1)*(num_cells2+sp_deg2+1))then
            print*, 'Problem in set_coefficients in arbitrary_degree_spline_2d'
            print*, ' Problem with the size coeffs_1d must have the size equal to '
            print*, ' (num_cells1 + sp_deg1 +1)*( num_cells2 + sp_deg2 +1)=',&
                 (num_cells1 + sp_deg1 +1)*( num_cells2 + sp_deg2 +1)
            stop
         end if
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
         ! achtung ! normaly interpolator%slope_left(:) and interpolator%value_right(:)
         ! achtung ! normaly interpolator%value_bottom(:) and interpolator%value_top(:)

         interpolator%coeff_splines(:,:) = 0.0_8
         ! allocation coefficient spline
         do i = 1,nb_spline_eta1
            do j = 1,nb_spline_eta2
               
               interpolator%coeff_splines(i+1,j+1) = &
                    coeffs_1d( i + nb_spline_eta1 *(j-1))
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
         
         if(size(coeffs_1d,1).ne.(num_cells1 + sp_deg1+1)*(num_cells2+sp_deg2+1))then
            print*, 'Problem in set_coefficients in arbitrary_degree_spline_2d'
            print*, ' Problem with the size coeffs_1d must have the size equal to '
            print*, ' (num_cells1 + sp_deg1 +1)*( num_cells2 + sp_deg2 +1)=',&
                 (num_cells1 + sp_deg1 +1)*( num_cells2 + sp_deg2 +1)
            stop
         end if
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
         ! achtung ! normaly interpolator%slope_left(:) and interpolator%value_right(:)
         ! achtung ! normaly interpolator%value_bottom(:) and interpolator%value_top(:)

         interpolator%coeff_splines(:,:) = 0.0_8
         ! allocation coefficient spline
         do i = 1,nb_spline_eta1
            do j = 1,nb_spline_eta2
               
               interpolator%coeff_splines(i+1,j+1) = &
                    coeffs_1d( i + nb_spline_eta1 *(j-1))
            end do
         end do
         ! ------------------------------------------------------------
        case(2340) ! Hermite in al sides
           
         interpolator%size_coeffs1=  num_cells1 + sp_deg1
         interpolator%size_coeffs2=  num_cells2 + sp_deg2
         interpolator%size_t1 = 2*sp_deg1 + num_cells1 + 1
         interpolator%size_t2 = 2*sp_deg2 + num_cells2 + 1
         nb_spline_eta1 = num_cells1 + sp_deg1
         nb_spline_eta2 = num_cells2 + sp_deg2
         
         if(size(coeffs_1d,1).ne.(num_cells1 + sp_deg1)*(num_cells2+sp_deg2))then
            print*, 'Problem in set_coefficients in arbitrary_degree_spline_2d'
            print*, ' Problem with the size coeffs_1d must have the size equal to '
            print*, ' (num_cells1 + sp_deg1)*( num_cells2 + sp_deg2)=',&
                 (num_cells1 + sp_deg1)*( num_cells2 + sp_deg2)

            print*, 'but it is equal to =', size(coeffs_1d,1)
            stop
         end if
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
         ! achtung ! normaly interpolator%slope_left(:) and interpolator%value_right(:)
         ! achtung ! normaly interpolator%value_bottom(:) and interpolator%value_top(:)

         interpolator%coeff_splines(:,:) = 0.0_8
         ! allocation coefficient spline
         do i = 1,nb_spline_eta1
            do j = 1,nb_spline_eta2
               
               interpolator%coeff_splines(i,j) = &
                    coeffs_1d( i + nb_spline_eta1 *(j-1))
            end do
         end do
         ! ------------------------------------------------------------
       
      case default
         print *, 'arbitrary_degree_spline_2d() error: set_spline_coefficients ',&
              'not recognized.'
         stop
      end select
   else if (present(coeffs_2d) ) then 

      if ( present(coeff2d_size1) .and. present(coeff2d_size2)) then

         interpolator%size_coeffs1 = coeff2d_size1
         interpolator%size_coeffs2 = coeff2d_size2
         interpolator%size_t1      = sp_deg1 + coeff2d_size1 +1 
         interpolator%size_t2      = sp_deg2 + coeff2d_size2 +1
         
         if ( coeff2d_size1 > num_cells1 + 1 + 4*sp_deg1) then
            print*, 'size1 of coeff2d is too big'
            stop
         end if
         
         if ( coeff2d_size2 > num_cells2 + 1 + 4*sp_deg2) then
            print*, 'size2 of coeff2d is too big'
            stop
         end if
         
         interpolator%coeff_splines(1:coeff2d_size1,1:coeff2d_size2) = &
              coeffs_2d(1:coeff2d_size1,1:coeff2d_size2)

         
         if ( present(knots1) .and. present(knots2) ) then 
            
            if ( ( size_knots1 .ne. (coeff2d_size1 + sp_deg1 + 1)  ) .OR.&
                 ( size_knots2 .ne. (coeff2d_size2 + sp_deg2 + 1)  ))  then
               print*, 'Problem in set_coefficients in arbitrary_degree_spline_2d'
               print*, 'problem with the size of knots'
               print*, 'size(knots1) must be equal to',coeff2d_size1 + sp_deg1 + 1
               print*, 'size(knots2) must be equal to',coeff2d_size2 + sp_deg2 + 1
               stop
            end if
             
            if ( size_knots1 > (num_cells1 + 1)*(sp_deg1+1)) then
               print*, 'size1 of knots1 is too big'
               stop
            end if
            
            if ( size_knots2 >  (num_cells2 + 1)*(sp_deg2+1)) then
               print*, 'size2 of knots2 is too big'
               stop
            end if
            
            
            
            interpolator%t1(1:interpolator%size_t1 ) = &
                 knots1(1:interpolator%size_t1 )
            interpolator%t2(1:interpolator%size_t2 ) =&
                 knots2(1:interpolator%size_t2 )
            
         else if ( (.not. present(knots1)).and.(.not. present(knots2))) then
            
            
            if ( interpolator%size_t1 > (num_cells1 + 1)*(sp_deg1+1)) then
               print*, 'size1 of knots1 is too big'
               stop
            end if
            
            if ( interpolator%size_t2 >  (num_cells2 + 1)*(sp_deg2+1)) then
               print*, 'size2 of knots2 is too big'
               stop
            end if

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
            print*, 'Problem in set_coefficients in arbitrary_degree_spline_2d'
            print*, 'Knots1 or Knots2 is not present'
            stop
            
         end if
         
      else 
         print*, 'Problem in set_coefficients in arbitrary_degree_spline_2d'
         print*, 'problem with the size of coeffs_2d'
         print*, 'the number of coefficients must be specified'
         stop
         
      end if
      
   else 
      print*, 'Problem in set_coefficients: must be have coefficients'
      stop
   end if
   interpolator%coefficients_set = .true.
    
 end subroutine !set_coefficients_ad2d


 ! ----------------------------------------------------------------
 ! subroutine computing the coefficients spline with a given 
 !  data_array 2D coorespondind at the values of a function 
 !  on eta1_coords of size size_eta1_coords in the first direction and 
 !  on eta2_coords of size size_eta2_coords in the second direction
 !  if the eta1_coords and eta2_coords is not given 
 !  we consider that the values of the function is on the points in the mesh_2d
 !   ----------------------------------------------------------------

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
    interpolator, &
    data_array, &
    eta1_coords, &
    size_eta1_coords, &
    eta2_coords, &
    size_eta2_coords )

    class(sll_arbitrary_degree_spline_interpolator_2d), intent(inout)  :: interpolator

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

    sll_int32 :: sz_derivative_eta1,sz_derivative_eta2
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

    
    !print*, data_array
    if(present(eta1_coords) .and. (.not. present(size_eta1_coords))) then
       print *, 'compute_interpolants_ad2d(), ERROR: if eta1_coords is ', &
            'passed, its size must be specified as well through ', &
            'size_eta1_coords.'
       stop
    end if
    
    if(present(eta2_coords) .and. (.not. present(size_eta2_coords))) then
       print *, 'compute_interpolants_ad2d(), ERROR: if eta2_coords is ', &
            'passed, its size must be specified as well through ', &
            'size_eta2_coords.'
       stop
    end if
    
    if ( (present(eta1_coords) .and. (.not. present(eta2_coords))) .or.&
       (present(eta2_coords) .and. (.not. present(eta1_coords))) ) then
       print *, 'compute_interpolants_ad2d(), ERROR: if either, ', &
            'eta1_coords or eta2_coords is specified, the other must be also.'
       stop
    end if
    
    if( present(eta1_coords) .and. present(eta2_coords) ) then
       user_coords = .true.
    else
       user_coords = .false.
    end if
    
    if (user_coords .eqv. .true.) then
       sz1 = size_eta1_coords
       sz2 = size_eta2_coords
       
       SLL_ALLOCATE(point_location_eta1(1:sz1),ierr)
       SLL_ALLOCATE(point_location_eta2(1:sz2),ierr)
       point_location_eta1(1:sz1) = eta1_coords(1:sz1)
       point_location_eta2(1:sz2) = eta2_coords(1:sz2)

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
    
    
    ! the size of data_array  must be <= interpolator%num_pts1 + 4*interpolator%spline_degree1
    ! because we have not need more !! 
    SLL_ASSERT(sz1 .le. interpolator%num_pts1 + 8*interpolator%spline_degree1)
    SLL_ASSERT(sz2 .le. interpolator%num_pts2 + 8*interpolator%spline_degree1)
    SLL_ASSERT(size(data_array,1) .ge. sz1)
    SLL_ASSERT(size(data_array,2) .ge. sz2)
    SLL_ASSERT(size(point_location_eta1)  .ge. sz1)
    SLL_ASSERT(size(point_location_eta2)  .ge. sz2)
    
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
       !print*,'zarrrr', interpolator%value_left(1:sz2)
       interpolator%coeff_splines(1,1:sz2)   = data_array(1,1:sz2)!interpolator%value_left(1:sz2)
       interpolator%coeff_splines(sz1,1:sz2) = data_array(sz1,1:sz2)!interpolator%value_right(1:sz2)
  
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
       data_array_deriv_eta1(2,1:sz2)     = interpolator%slope_right(1:sz2)
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=0.0_f64
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)=interpolator%slope_top(1:sz1+sz_derivative_eta1)
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
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_left(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1:sz1+sz_derivative_eta1,1)   = interpolator%value_bottom(1:sz1+sz_derivative_eta1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_top(1:sz1)

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
       data_array_deriv_eta1(1,1:sz2)     = interpolator%slope_left(1:sz2) 
       data_array_deriv_eta1(2,1:sz2)     = 0.0_f64
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)= 0.0_f64
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)=interpolator%slope_top(1:sz1+sz_derivative_eta1)
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
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_left(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:sz1,1)   = interpolator%value_bottom(1:sz1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_top(1:sz1)


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
       data_array_deriv_eta1(1,1:sz2) = interpolator%slope_left(1:sz2)
       data_array_deriv_eta1(2,1:sz2) = interpolator%slope_right(1:sz2)
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)= &
          interpolator%slope_bottom(1:sz1+sz_derivative_eta1)
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)= &
          interpolator%slope_top(1:sz1+sz_derivative_eta1)

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
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_left(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:sz1,1)   = interpolator%value_bottom(1:sz1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_top(1:sz1)


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
       data_array_deriv_eta1(1,1:sz2)     = interpolator%slope_left(1:sz2) 
       data_array_deriv_eta1(2,1:sz2)     = interpolator%slope_right(1:sz2) 
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=interpolator%slope_bottom(1:sz1+sz_derivative_eta1)
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)=interpolator%slope_top(1:sz1+sz_derivative_eta1)
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
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_left(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:sz1,1)   = interpolator%value_bottom(1:sz1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_top(1:sz1)

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
       data_array_deriv_eta1(1,1:sz2)     = interpolator%slope_left(1:sz2) 
       data_array_deriv_eta1(2,1:sz2)     = interpolator%slope_right(1:sz2) 
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=interpolator%slope_bottom(1:sz1+sz_derivative_eta1)
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)=interpolator%slope_top(1:sz1+sz_derivative_eta1)
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
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_left(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:sz1,1)   = interpolator%value_bottom(1:sz1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_top(1:sz1)

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
       data_array_deriv_eta1(2,1:sz2)     = interpolator%slope_right(1:sz2)
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=interpolator%slope_bottom(1:sz1+sz_derivative_eta1)
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
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_left(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:sz1,1)   = interpolator%value_bottom(1:sz1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_top(1:sz1)

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
       data_array_deriv_eta1(1,1:sz2)     = interpolator%slope_left(1:sz2)
       data_array_deriv_eta1(2,1:sz2)     = 0.0_f64
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=interpolator%slope_bottom(1:sz1+sz_derivative_eta1)
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
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_left(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1:sz1+sz_derivative_eta1,1)   = interpolator%value_bottom(1:sz1+sz_derivative_eta1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_top(1:sz1)

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
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_left(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1:sz1+sz_derivative_eta1,1)   = interpolator%value_bottom(1:sz1+sz_derivative_eta1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_top(1:sz1)
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
       data_array_deriv_eta1(1,1:sz2)     = interpolator%slope_left(1:sz2)
       data_array_deriv_eta1(2,1:sz2)     = interpolator%slope_right(1:sz2)
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=interpolator%slope_bottom(1:sz1+sz_derivative_eta1)
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)=interpolator%slope_top(1:sz1+sz_derivative_eta1)
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
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_left(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1:sz1+sz_derivative_eta1,1)   = interpolator%value_bottom(1:sz1+sz_derivative_eta1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_top(1:sz1)
       
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
       data_array_deriv_eta1(1,1:sz2)     = interpolator%slope_left(1:sz2)
       data_array_deriv_eta1(2,1:sz2)     = interpolator%slope_right(1:sz2)
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=interpolator%slope_bottom(1:sz1+sz_derivative_eta1)
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)=interpolator%slope_top(1:sz1+sz_derivative_eta1)
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
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_left(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
       !interpolator%coeff_splines(1:sz1+sz_derivative_eta1,1)   = interpolator%value_bottom(1:sz1+sz_derivative_eta1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_top(1:sz1)


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
       data_array_deriv_eta1(1,1:sz2)     = interpolator%slope_left(1:sz2)
       data_array_deriv_eta1(2,1:sz2)     = interpolator%slope_right(1:sz2)
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=interpolator%slope_bottom(1:sz1+sz_derivative_eta1)
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)=interpolator%slope_top(1:sz1+sz_derivative_eta1)
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
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_left(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:sz1,1)   = interpolator%value_bottom(1:sz1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_top(1:sz1)

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
       data_array_deriv_eta1(1,1:sz2)     = interpolator%slope_left(1:sz2)
       data_array_deriv_eta1(2,1:sz2)     = interpolator%slope_right(1:sz2)
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=interpolator%slope_bottom(1:sz1+sz_derivative_eta1)
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)=interpolator%slope_top(1:sz1+sz_derivative_eta1)
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
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_left(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:sz1,1)   = interpolator%value_bottom(1:sz1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_top(1:sz1)

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
       data_array_deriv_eta1(1,1:sz2)     = interpolator%slope_left(1:sz2)
       data_array_deriv_eta1(2,1:sz2)     = interpolator%slope_right(1:sz2)
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=interpolator%slope_bottom(1:sz1+sz_derivative_eta1)
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)=interpolator%slope_top(1:sz1+sz_derivative_eta1)
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
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_left(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:sz1,1)   = interpolator%value_bottom(1:sz1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_top(1:sz1)

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
       data_array_deriv_eta1(1,1:sz2)     = interpolator%slope_left(1:sz2)
       data_array_deriv_eta1(2,1:sz2)     = interpolator%slope_right(1:sz2)
       point_location_eta2_deriv(1) = 1
       point_location_eta2_deriv(2) = sz2
       data_array_deriv_eta2(1,1:sz1+sz_derivative_eta1)=interpolator%slope_bottom(1:sz1+sz_derivative_eta1)
       data_array_deriv_eta2(2,1:sz1+sz_derivative_eta1)=interpolator%slope_top(1:sz1+sz_derivative_eta1)
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
       !interpolator%coeff_splines(1,1:sz2+sz_derivative_eta2)   = interpolator%value_left(1:sz2+sz_derivative_eta2)
       ! boundary condition non homogene
 !      interpolator%coeff_splines(1:sz1,1)   = interpolator%value_bottom(1:sz1)
  !     interpolator%coeff_splines(1:sz1,sz2) = interpolator%value_top(1:sz1)

    end select
    interpolator%coefficients_set = .true.
   
    SLL_DEALLOCATE(point_location_eta2,ierr)
    SLL_DEALLOCATE(point_location_eta1,ierr)
    SLL_DEALLOCATE(point_location_eta1_tmp,ierr)
    SLL_DEALLOCATE(point_location_eta2_tmp,ierr)
    SLL_DEALLOCATE(data_array_tmp,ierr)

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
    !sll_real64 :: bvalue2d
    !integer  :: li_i, li_j, li_mflag, li_lefty
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


  ! Just for the record: what is ar_x, ar_y, ai_nx, etc., etc.? What is this
  ! function supposed to do? This is a deeply frustrating file.
  subroutine bvalue2d(&
       ar_x,&
       ar_y,&
       ai_nx,&
       ai_kx,&
       ai_ny,&
       ai_ky,&
       apr_Bcoef,&
       apr_tx,&
       apr_ty,&
       val )
    implicit none
    ! INPUT
    sll_real64,intent(in) :: ar_x, ar_y
    sll_int32,intent(in)  :: ai_nx, ai_kx, ai_ny, ai_ky
    sll_real64, dimension(:), pointer :: apr_tx !  ai_nx + ai_kx 
    sll_real64, dimension(:), pointer :: apr_ty !  ai_ny + ai_ky	
    sll_real64, dimension(:,:),pointer :: apr_Bcoef!( ai_nx,ai_ny)
    !OUTPUT
    sll_real64,intent(out) ::val
    ! LOCAL VARIABLES
    sll_int32  :: li_j, li_mflag, li_lefty
    sll_real64, dimension(1:ai_ky),target :: lpr_coef ! ai_ky
    sll_real64, dimension(:),pointer :: lpr_coef_ptr ! ai_ky
    sll_real64, dimension(1:ai_nx),target :: tmp_tab
    sll_real64, dimension(:),pointer :: tmp_tab_ptr
    sll_real64, dimension(1:2*ai_ky),target :: tmp_ty
    sll_real64, dimension(:),pointer :: tmp_ty_ptr
    sll_int32 :: ierr
   

    call interv ( apr_ty,ai_ny + ai_ky, ar_y, li_lefty, li_mflag )

    if ( li_mflag .NE. 0 ) then
       val = 0.0_8
       return 
    end if
    
    do li_j = 1, ai_ky
       
       
       tmp_tab = apr_bcoef ( 1:ai_nx , li_lefty - ai_ky + li_j )
       tmp_tab_ptr => tmp_tab
       lpr_coef ( li_j ) = bvalue(&
            apr_tx,&
            tmp_tab_ptr,&
            ai_nx,&
            ai_kx,&
            ar_x,&
            0 )
       
       
    end do
   
    lpr_coef_ptr => lpr_coef
    tmp_ty =  apr_ty ( li_lefty - ai_ky + 1 : li_lefty + ai_ky)
    tmp_ty_ptr => tmp_ty
    val = bvalue(&
         tmp_ty_ptr,&
         lpr_coef_ptr,&
         ai_ky,&
         ai_ky,&
         ar_y,&
         0 )
    end subroutine bvalue2d


  function dvalue2d(&
       ar_x,&
       ar_y,&
       ai_nx,&
       ai_kx,&
       ai_ny,&
       ai_ky,&
       apr_Bcoef,&
       apr_tx,&
       apr_ty,deriv1,deriv2 ) result(res)
    implicit none
    ! INPUT
    sll_real64 :: ar_x, ar_y
    sll_real64 :: res
    sll_int32  :: ai_nx, ai_kx, ai_ny, ai_ky
    sll_int32  :: deriv1,deriv2
    sll_real64, dimension ( : ),pointer :: apr_tx ! ai_nx + ai_kx
    sll_real64, dimension ( : ),pointer :: apr_ty ! ai_ny + ai_ky
    sll_real64, dimension ( : , : ),pointer :: apr_Bcoef !(ai_nx,ai_ny)
    ! LOCAL VARIABLES
    sll_int32  :: li_j, li_mflag, li_lefty
    sll_real64, dimension (1:ai_ky),target:: lpr_coef ! ai_ky
    sll_real64, dimension (:),pointer:: lpr_coef_ptr
    sll_real64, dimension (1:ai_nx),target :: tmp_coef
    sll_real64,dimension(:),pointer::  tmp_coef_ptr
    sll_real64, dimension (1:2*ai_ky),target :: tmp_ty
    sll_real64, dimension (:),pointer :: tmp_ty_ptr
    sll_int32:: ierr
    
    call interv ( apr_ty, ai_ny + ai_ky, ar_y, li_lefty, li_mflag )
    
    if ( li_mflag .NE. 0 ) then
       res = 0.0_8
       return 
    end if
    
    do li_j = 1, ai_ky
       
       tmp_coef = apr_bcoef ( 1:ai_nx , li_lefty - ai_ky + li_j )
       tmp_coef_ptr => tmp_coef
       lpr_coef ( li_j ) = bvalue(&
            apr_tx,&
            tmp_coef_ptr,&
            ai_nx,&
            ai_kx,&
            ar_x,&
            deriv1 )
       
    end do
    lpr_coef_ptr => lpr_coef
    tmp_ty =  apr_ty ( li_lefty - ai_ky + 1 : li_lefty + ai_ky)
    tmp_ty_ptr => tmp_ty
    res = bvalue(&
         tmp_ty_ptr,&
         lpr_coef_ptr,&
         ai_ky,&
         ai_ky,&
         ar_y,&
         deriv2 )

  end function dvalue2d
  


  subroutine spli2d ( tau, gtau, t, n, k, m, work, q, bcoef, iflag )

    !*****************************************************************************80
    !
    !! SPLI2D produces a interpolatory tensor product spline.
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
    implicit none
    
    sll_int32 :: m
    sll_int32 :: n
    
    sll_real64,dimension(:,:),pointer:: bcoef !(m,n)
    sll_real64,dimension(:,:),pointer:: gtau  !(n,m)
    sll_int32 :: i
    sll_int32 :: iflag
    sll_int32 :: ilp1mx
    sll_int32 :: j
    sll_int32 :: jj
    sll_int32 :: k
    sll_int32 :: left
    sll_real64,dimension((2*k-1)*n):: q!((2*k-1)*n)
    sll_real64,dimension(:),pointer:: t!(n+k)
    sll_real64,dimension(:),pointer:: tau!(n)
    sll_real64:: taui
    sll_real64,dimension(n):: work!(n)
    
    left = k
    
    !print*, t
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
             !print*, taui, t(left+1),left
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

  subroutine spli2d_custom ( &
     ai_nx,&
     ai_kx,&
     apr_taux,&
     ai_ny,&
     ai_ky,&
     apr_tauy,&
     apr_g,&
     apr_Bcoef,&
     apr_tx,&
     apr_ty )
    implicit none
    ! INPUT
    sll_int32  :: ai_nx, ai_kx, ai_ny, ai_ky
    sll_real64, dimension (:),pointer :: apr_taux !!ai_nx
    sll_real64, dimension (:),pointer :: apr_tauy !! ai_ny	
    sll_real64, dimension (:,:),pointer :: apr_g    ! ai_nx,ai_ny	
    ! OUTPUT
    sll_real64, dimension (:,:),pointer :: apr_Bcoef !ai_nx , ai_ny 
    sll_real64, dimension ( : ),pointer :: apr_tx ! ai_nx + ai_kx
    sll_real64, dimension ( : ),pointer :: apr_ty ! ai_ny + ai_ky 
    ! LOCAL VARIABLES		
    sll_real64, dimension ( ai_nx , ai_ny ) :: lpr_work1
    sll_real64, dimension ( ai_nx         ) :: lpr_work2
    sll_real64, dimension ( ai_nx * ai_ny ) :: lpr_work3
    sll_real64, dimension ( ai_nx *( 2*ai_kx-1) ) :: lpr_work31
    sll_real64, dimension (( 2*ai_ky-1) * ai_ny ) :: lpr_work32
    sll_real64, dimension ( ai_ny         ) :: lpr_work4
    sll_real64, dimension (1:ai_ny,1:ai_nx),target :: lpr_work5 !  ai_ny , ai_nx 
    sll_real64, dimension (:,:),pointer :: lpr_work5_ptr !  ai_ny , ai_nx 
    sll_real64, dimension (1:ai_ny),target:: apr_ty_bis
    sll_real64, dimension (:),pointer:: apr_ty_bis_ptr
    sll_int32  :: li_i, li_j, li_iflag
    sll_int32 :: ierr
    
    
    lpr_work1(:,:) = 0.0
    
    ! *** set up knots
    !     interpolate between knots
    
    apr_tx ( 1 : ai_kx ) = apr_taux ( 1 )
    apr_tx ( ai_nx + 1 : ai_nx + ai_kx ) = apr_taux ( ai_nx )
    
    if ( mod(ai_kx,2) == 0 ) then
       do li_i = ai_kx + 1, ai_nx
          apr_tx ( li_i ) = apr_taux ( li_i - ai_kx/2 ) 
          
       end do
     else
        
        do li_i = ai_kx + 1, ai_nx
           apr_tx ( li_i ) = &
                0.5*( apr_taux ( li_i - (ai_kx-1)/2 ) + &
                apr_taux ( li_i -1 - (ai_kx-1)/2 ) )
           
        end do
     
     end if
     apr_Bcoef = 0.0_8
     do li_i = 1, ai_nx
        do li_j = 1, ai_ny
           apr_Bcoef ( li_i, li_j ) = apr_g ( li_i, li_j )
        end do
     end do
     !  *** construct b-coefficients of interpolant
     !
     apr_ty = 0.0_f64
     
     if ( mod(ai_ky,2) == 0 ) then
        do li_i = ai_ky + 1, ai_ny
           apr_ty ( li_i ) = apr_tauy ( li_i - ai_ky/2 ) 
           
        end do
     else
        
        do li_i = ai_ky + 1, ai_ny
           apr_ty ( li_i ) = &
                0.5*( apr_tauy ( li_i - (ai_ky-1)/2 ) + &
                apr_tauy ( li_i -1 - (ai_ky-1)/2 ) )
           
        end do
        
     end if
     apr_ty ( 1 : ai_ky ) = apr_tauy ( 1 )
     apr_ty ( ai_ny + 1 : ai_ny + ai_ky ) = apr_tauy ( ai_ny )
     
     apr_ty_bis = apr_tauy(1:ai_ny)
     
     lpr_work5_ptr => lpr_work5
     call spli2d ( &
          apr_taux,&
          apr_Bcoef,&
          apr_tx, &
          ai_nx,&
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
     
     apr_ty_bis_ptr => apr_ty_bis
     call spli2d ( &
          apr_ty_bis_ptr,&
          lpr_work5_ptr,&
          apr_ty,&
          ai_ny, &
          ai_ky, &
          ai_nx, &
          lpr_work4, &
          lpr_work32,&
          apr_bcoef, &
          li_iflag )
     
   end subroutine spli2d_custom
   



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

   subroutine spli2d_custom_der1 ( &
     ai_nx,&
     ai_nx_der,&
     ai_kx,&
     apr_taux,&
     apr_taux_der,&
     ai_ny,&
     ai_ky,&
     apr_tauy,&
     apr_g,&
     apr_g_der1,&
     apr_Bcoef,&
     apr_tx,&
     apr_ty )
     implicit none
     ! INPUT
     sll_int32  :: ai_nx, ai_kx, ai_ny, ai_ky
     sll_int32  :: ai_nx_der
     sll_real64, dimension(:),pointer :: apr_taux !!ai_nx
     sll_real64, dimension(:),pointer :: apr_tauy !! ai_ny
     sll_int32,  dimension(:),pointer :: apr_taux_der !!ai_nx_der
     sll_real64, dimension(:,:),pointer :: apr_g    ! ai_nx,ai_ny
     sll_real64, dimension(:,:),pointer :: apr_g_der1 ! ai_nx_der,ai_ny
 
     ! OUTPUT
     sll_real64, dimension(:,:),pointer::apr_Bcoef!ai_nx + ai_nx_der,ai_ny
     sll_real64, dimension( : ),pointer:: apr_tx ! ai_nx + ai_kx + ai_nx_der
     sll_real64, dimension( : ),pointer:: apr_ty ! ai_ny + ai_ky
    ! LOCAL VARIABLES		
    sll_real64, dimension ( ai_nx + ai_nx_der , ai_ny) :: lpr_work1
    sll_real64, dimension ( ai_nx + ai_nx_der ) :: lpr_work2
    sll_real64, dimension ( (ai_nx + ai_nx_der)* (ai_ny) ) :: lpr_work3
    sll_real64, dimension ( (ai_nx+ai_nx_der) *( 2*ai_kx-1) ) :: lpr_work31
    sll_real64, dimension (( 2*ai_ky-1) * (ai_ny) ) :: lpr_work32
    sll_real64, dimension ( ai_ny) :: lpr_work4
    sll_real64, dimension (1:ai_ny,1:ai_nx +ai_nx_der),target :: lpr_work5 !  ai_ny , ai_nx 
    sll_real64, dimension (:,:),pointer :: lpr_work5_ptr
    sll_int32  :: li_i, li_iflag
    sll_int32 :: ierr
    
    
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
    
    apr_ty = 0.0_f64
    
    if ( mod(ai_ky,2) == 0 ) then
       do li_i = ai_ky + 1, ai_ny
          apr_ty ( li_i ) = apr_tauy ( li_i - ai_ky/2 ) 
          
       end do
    else
       
       do li_i = ai_ky + 1, ai_ny
          apr_ty ( li_i ) = &
                0.5*( apr_tauy ( li_i - (ai_ky-1)/2 ) + &
                apr_tauy ( li_i -1 - (ai_ky-1)/2 ) )
          
       end do
       
    end if
    apr_ty ( 1 : ai_ky ) = apr_tauy ( 1 )
    apr_ty ( ai_ny + 1 : ai_ny + ai_ky ) = apr_tauy ( ai_ny )
    
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
     

     call spli2d ( &
          apr_tauy,&
          lpr_work5_ptr,&
          apr_ty,&
          ai_ny, &
          ai_ky, &
          ai_nx+ai_nx_der, &
          lpr_work4, &
          lpr_work32,&
          apr_bcoef, &
          li_iflag )

   end subroutine spli2d_custom_der1
   

   subroutine spli2d_custom_der2 ( &
     ai_nx,&
     ai_kx,&
     ai_ny,&
     ai_ny_der,&
     ai_ky,&
     apr_tauy,&
     apr_tauy_der,&
     apr_g,&
     apr_g_der2,&
     apr_Bcoef,&
     apr_tx,&
     apr_ty )
     implicit none
     ! INPUT
     sll_int32  :: ai_nx, ai_kx, ai_ny, ai_ky
     sll_int32  :: ai_ny_der
     sll_real64, dimension(:),pointer :: apr_taux !!ai_nx
     sll_real64, dimension(:),pointer :: apr_tauy !! ai_ny
     sll_int32,  dimension(:),pointer :: apr_tauy_der !!ai_ny_der
     sll_real64, dimension(:,:),pointer :: apr_g ! ai_nx,ai_ny
     sll_real64, dimension(:,:),pointer :: apr_g_der2 !ai_ny_der,ai_nx
     ! OUTPUT
     sll_real64, dimension(:,:),pointer::apr_Bcoef!ai_nx ,ai_ny+ ai_ny_der 
     sll_real64, dimension( : ),pointer:: apr_tx ! ai_nx + ai_kx 
     sll_real64, dimension( : ),pointer:: apr_ty ! ai_ny + ai_ky + ai_ny_der
     ! LOCAL VARIABLES		
     sll_real64, dimension ( ai_nx , ai_ny + ai_ny_der) :: lpr_work1
    sll_real64, dimension ( ai_nx ) :: lpr_work2
    sll_real64, dimension ( (ai_nx )* (ai_ny+ai_ny_der) ) :: lpr_work3
    sll_real64, dimension ( (ai_nx) *( 2*ai_kx-1) ) :: lpr_work31
    sll_real64, dimension (( 2*ai_ky-1) * (ai_ny+ai_ny_der) ) :: lpr_work32
    sll_real64, dimension ( ai_ny +ai_ny_der) :: lpr_work4
    sll_real64, dimension (1:ai_ny,1:ai_nx),target :: lpr_work5 !  ai_ny , ai_nx 
    sll_real64, dimension (:,:),pointer :: lpr_work5_ptr
    sll_int32  :: li_i,li_iflag
    sll_int32 :: ierr
    
    
    lpr_work1(:,:) = 0.0
    lpr_work5(:,:) = 0.0
    
    ! *** set up knots
    !     interpolate between knots
    
    apr_tx ( 1 : ai_kx ) = apr_taux ( 1 )
    apr_tx ( ai_nx + 1 : ai_nx + ai_kx ) = apr_taux ( ai_nx )
    
    if ( mod(ai_kx,2) == 0 ) then
       do li_i = ai_kx + 1, ai_nx
          apr_tx ( li_i ) = apr_taux ( li_i - ai_kx/2 ) 
          
       end do
     else
        
        do li_i = ai_kx + 1, ai_nx
           apr_tx ( li_i ) = &
                0.5*( apr_taux ( li_i - (ai_kx-1)/2 ) + &
                apr_taux ( li_i -1 - (ai_kx-1)/2 ) )
           
        end do
     
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
    
    lpr_work5_ptr = lpr_work5
    call spli2d( &
         apr_taux,&
         apr_g,&
         apr_tx, &
         ai_nx,&
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
          ai_nx, &
          lpr_work4, &
          lpr_work32,&
          apr_bcoef, &
          li_iflag )

   end subroutine spli2d_custom_der2

   subroutine spli2d_perdir (&
        ar_L,&
        ai_nx,&
        ai_kx,&
        apr_taux,&
        ai_ny,&
        ai_ky,&
        apr_tauy,&
        apr_g,&
        apr_Bcoef,&
        apr_tx,&
        apr_ty )
     ! CALLED WHEN WE WANT TO INTERPOL WITH A PERIODIC FIRST PARAM WITH A PERIOD = ar_L
     implicit none
     ! INPUT
     sll_real64 :: ar_L 
     sll_int32  :: ai_nx, ai_kx, ai_ny, ai_ky
     sll_real64, dimension ( :),pointer :: apr_taux ! ai_nx- 1
     sll_real64, dimension ( :),pointer :: apr_tauy ! ai_ny		
     sll_real64, dimension ( :,:) :: apr_g !ai_nx - 1, ai_ny
     ! OUTPUT
     sll_real64, dimension (:,:),pointer :: apr_Bcoef !  ai_nx , ai_ny	
     sll_real64, dimension (:),pointer :: apr_tx !  ai_nx + ai_kx
     sll_real64, dimension (:),pointer :: apr_ty ! ai_ny + ai_ky
     ! LOCAL VARIABLES		
     sll_real64, dimension (1:ai_nx),target :: lpr_taux !  ai_nx
     sll_real64, dimension (:),pointer :: lpr_taux_ptr
     sll_real64, dimension (1:ai_nx,1:ai_ny),target :: lpr_g !  ai_nx ,ai_ny
     sll_real64, dimension (:,:),pointer :: lpr_g_ptr
     sll_int32 :: ierr

     if ( ar_L == 0.0_8 ) then
        print*,'Error spli2d_per : called with a period = 0 '
        stop
     end if
     
    
     lpr_taux ( 1 : ai_nx - 1 ) = apr_taux ( 1 : ai_nx-1)
     lpr_taux ( ai_nx ) = apr_taux ( 1 ) + ar_L

     lpr_g ( 1 : ai_nx - 1 , 1 : ai_ny ) = apr_g ( 1 : ai_nx - 1 , 1 : ai_ny )
     lpr_g ( ai_nx , 1 : ai_ny ) = apr_g ( 1 , 1 : ai_ny )

     lpr_taux_ptr => lpr_taux
     lpr_g_ptr => lpr_g
     
     call spli2d_custom ( &
          ai_nx, &
          ai_kx, &
          lpr_taux_ptr,&
          ai_ny,&
          ai_ky, &
          apr_tauy, &
          lpr_g_ptr,&
          apr_Bcoef,&
          apr_tx,&
          apr_ty )

     
   end subroutine spli2d_perdir
   
   subroutine spli2d_dirper (&
        ai_nx,&
        ai_kx,&
        apr_taux,&
        ar_L, &
        ai_ny,&
        ai_ky, &
        apr_tauy,&
        apr_g,&
        apr_Bcoef,&
        apr_tx,&
        apr_ty )
     ! CALLED WHEN WE WANT TO INTERPOL WITH A PERIODIC second PARAM WITH A PERIOD = ar_L
     implicit none
     ! INPUT
     sll_real64 :: ar_L
     sll_int32  :: ai_nx, ai_kx, ai_ny, ai_ky
     sll_real64, dimension ( :),pointer :: apr_taux ! ai_nx
     sll_real64, dimension (:),pointer :: apr_tauy !  ai_ny -1
     sll_real64, dimension ( :,:) :: apr_g ! ai_nx , ai_ny-1
     ! OUTPUT
     sll_real64, dimension (:,:),pointer :: apr_Bcoef !  ai_nx , ai_ny
     sll_real64, dimension ( :),pointer :: apr_tx ! ai_nx + ai_kx	
     sll_real64, dimension (:),pointer :: apr_ty ! ai_ny + ai_ky 
     ! LOCAL VARIABLES
     sll_real64, dimension (1:ai_ny),target :: lpr_tauy ! ai_ny	
     sll_real64, dimension (1:ai_nx,1:ai_ny),target :: lpr_g  !  ai_nx ,ai_ny
     sll_real64, dimension (:),pointer :: lpr_tauy_ptr ! ai_ny	
     sll_real64, dimension (:,:),pointer :: lpr_g_ptr
     sll_int32 :: ierr
     
     
     if ( ar_L == 0.0_8 ) then
        print*,'Error spli2d_per : called with a period = 0 '
        stop
     end if
     
     
     lpr_tauy ( 1 : ai_ny - 1 ) = apr_tauy ( 1 : ai_ny - 1 )
     lpr_tauy ( ai_ny ) = apr_tauy ( 1 ) + ar_L
     
     lpr_g ( 1 : ai_nx , 1 : ai_ny -1 ) = apr_g ( 1 : ai_nx , 1 : ai_ny -1)
     lpr_g (1: ai_nx , ai_ny ) = apr_g ( 1 : ai_nx, 1 )
     
     lpr_tauy_ptr => lpr_tauy
     lpr_g_ptr => lpr_g
     call spli2d_custom (&
          ai_nx,&
          ai_kx,&
          apr_taux,&
          ai_ny, &
          ai_ky,&
          lpr_tauy_ptr, &
          lpr_g_ptr, &
          apr_Bcoef,&
          apr_tx,&
          apr_ty )

  
   end subroutine spli2d_dirper
   

   subroutine spli2d_perper(&
     ar_Lx,&
     ai_nx,&
     ai_kx,&
     apr_taux,&
     ar_Ly,&
     ai_ny,&
     ai_ky,&
     apr_tauy,&
     apr_g,&
     apr_Bcoef,&
     apr_tx,&
     apr_ty )
  ! CALLED WHEN WE WANT TO INTERPOL WITH A PERIODIC FIRST PARAM 
     !WITH A PERIOD = ar_L
     implicit none
     ! INPUT
     sll_real64 :: ar_Lx
     sll_real64 :: ar_Ly
     sll_int32  :: ai_nx, ai_kx, ai_ny, ai_ky
     sll_real64, dimension (:),pointer :: apr_taux ! ai_nx -1
     sll_real64, dimension (:),pointer :: apr_tauy ! ai_ny - 1
     sll_real64, dimension (:,:),pointer :: apr_g !  ai_nx  - 1, ai_ny - 1
     ! OUTPUT
     sll_real64, dimension ( :,:),pointer :: apr_Bcoef ! ai_nx , ai_ny
     sll_real64, dimension (:),pointer :: apr_tx !  ai_nx + ai_kx
     sll_real64, dimension ( :),pointer :: apr_ty ! ai_ny + ai_ky
     ! LOCAL VARIABLES
     sll_real64, dimension (1:ai_nx),target :: lpr_taux ! tmp_ty
     sll_real64, dimension (1:ai_ny),target :: lpr_tauy
     sll_real64, dimension(1:ai_nx,1:ai_ny),target :: lpr_g !  ( ai_nx ,ai_ny)
     sll_real64, dimension (:),pointer :: lpr_taux_ptr ! tmp_ty
     sll_real64, dimension (:),pointer :: lpr_tauy_ptr
     sll_real64, dimension(:,:),pointer :: lpr_g_ptr
     sll_int32 :: ierr
     
     if ( ar_Lx == 0.0_8 ) then
        print*,'Error spli2d_perper : called with a period = 0 '
        stop
     end if
     if ( ar_Ly == 0.0_8 ) then
        print*,'Error spli2d_perper : called with a period = 0 '
        stop
     end if
  
     lpr_taux ( 1 : ai_nx - 1 ) = apr_taux ( 1 : ai_nx - 1 )
     lpr_taux ( ai_nx ) = apr_taux ( 1 ) + ar_Lx
     
     
     lpr_tauy ( 1 : ai_ny - 1 ) = apr_tauy ( 1 : ai_ny - 1 )
     lpr_tauy ( ai_ny ) = apr_tauy ( 1 ) + ar_Ly
     
     
     lpr_g ( 1 : ai_nx - 1 , 1 : ai_ny - 1 ) = &
          apr_g ( 1 : ai_nx - 1 , 1 : ai_ny -1 )
     lpr_g ( ai_nx , 1 : ai_ny -1 ) = apr_g ( 1 , 1 : ai_ny -1 )
     lpr_g ( 1 : ai_nx -1 , ai_ny ) = apr_g ( 1 : ai_nx -1, 1 )
     lpr_g ( ai_nx , ai_ny ) = apr_g ( 1 , 1 )

     lpr_taux_ptr => lpr_taux
     lpr_tauy_ptr => lpr_tauy
     lpr_g_ptr => lpr_g
     call spli2d_custom ( &
          ai_nx,&
          ai_kx,&
          lpr_taux_ptr,&
          ai_ny,&
          ai_ky,&
          lpr_tauy_ptr,&
          lpr_g_ptr,&
          apr_Bcoef,&
          apr_tx,&
          apr_ty )
     
   end subroutine spli2d_perper
   




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


end module sll_module_arbitrary_degree_spline_interpolator_2d
