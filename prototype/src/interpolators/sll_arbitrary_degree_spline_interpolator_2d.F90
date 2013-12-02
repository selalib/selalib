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

module sll_arbitrary_degree_spline_interpolator_2d_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h" 
  use sll_timer
#ifdef STDF95
use sll_boundary_condition_descriptors
!use sll_constants
#else
use sll_module_interpolators_2d_base
#endif
  implicit none

  ! in what follows, the direction '1' is in the contiguous memory direction.
#ifdef STDF95
  type                                    :: arb_deg_2d_interpolator           
#else
  type, extends(sll_interpolator_2d_base) :: arb_deg_2d_interpolator           
#endif
     sll_int32  :: num_pts1
     sll_int32  :: num_pts2
     sll_real64 :: eta1_min
     sll_real64 :: eta1_max
     sll_real64 :: eta2_min
     sll_real64 :: eta2_max
     sll_int32  :: bc_left
     sll_int32  :: bc_right
     sll_int32  :: bc_bottom
     sll_int32  :: bc_top
     sll_int32  :: spline_degree1
     sll_int32  :: spline_degree2
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
     ! table contains the coeff spline of the function in boundary 
     ! in the case of dirichlet boundary condition non homogene 
     sll_real64, dimension(:),pointer :: slope_left
     sll_real64, dimension(:),pointer :: slope_right
     sll_real64, dimension(:),pointer :: slope_bottom
     sll_real64, dimension(:),pointer :: slope_top

#ifndef STDF95
   contains
    procedure, pass(interpolator) :: initialize=>initialize_ad2d_interpolator
    procedure, pass(interpolator) :: set_coefficients => set_coefficients_ad2d
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
#endif
  end type arb_deg_2d_interpolator

  interface delete
     module procedure delete_arbitrary_degree_2d_interpolator
  end interface delete

contains

  subroutine delete_arbitrary_degree_2d_interpolator( interpolator )
#ifdef STDF95
    type (arb_deg_2d_interpolator), intent(inout) :: interpolator
#else
    class(arb_deg_2d_interpolator), intent(inout) :: interpolator
#endif
    sll_int32 :: ierr
    SLL_DEALLOCATE(interpolator%knots1,ierr)
    SLL_DEALLOCATE(interpolator%knots2,ierr)
    SLL_DEALLOCATE(interpolator%t1,ierr)
    SLL_DEALLOCATE(interpolator%t2,ierr)
    SLL_DEALLOCATE(interpolator%coeff_splines,ierr)
    SLL_DEALLOCATE(interpolator%slope_left,ierr)
    SLL_DEALLOCATE(interpolator%slope_right,ierr)
    SLL_DEALLOCATE(interpolator%slope_bottom,ierr)
    SLL_DEALLOCATE(interpolator%slope_top,ierr)
  end subroutine delete_arbitrary_degree_2d_interpolator

  function new_arbitrary_degree_spline_interpolator_2d( &
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
    spline_degree2,&
    slope_left,&
    slope_right,&
    slope_bottom,&
    slope_top) result( res )

    type(arb_deg_2d_interpolator), pointer :: res
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
    sll_real64, dimension(:),optional :: slope_left
    sll_real64, dimension(:),optional :: slope_right
    sll_real64, dimension(:),optional :: slope_bottom
    sll_real64, dimension(:),optional :: slope_top
    sll_int32 :: ierr

    SLL_ALLOCATE(res,ierr)
    call initialize_ad2d_interpolator( &
         res, &
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
         spline_degree2,&
         slope_left,&
         slope_right,&
         slope_bottom,&
         slope_top)
  end function new_arbitrary_degree_spline_interpolator_2d

#ifdef STDF95
  subroutine arbitrary_degree_spline_interp2d_initialize( &
#else
  subroutine initialize_ad2d_interpolator( &
#endif
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
    spline_degree2,&
    slope_left,&
    slope_right,&
    slope_bottom,&
    slope_top)
    use sll_arbitrary_degree_spline_interpolator_1d_module

#ifdef STDF95
    type (arb_deg_2d_interpolator), intent(inout) :: interpolator
#else
    class(arb_deg_2d_interpolator), intent(inout) :: interpolator
#endif
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
    ! In the case of Dirichlet boundary conditions we can have 
    ! non homogene and homogene case
    ! slope_  represente the value of a function in the nodes of boundary
    ! if the user put anything we consider that is equal to 0 
    sll_real64, dimension(:),optional :: slope_left
    sll_real64, dimension(:),optional :: slope_right
    sll_real64, dimension(:),optional :: slope_bottom
    sll_real64, dimension(:),optional :: slope_top
    type(arb_deg_1d_interpolator),pointer :: interp1d_slope_left
    type(arb_deg_1d_interpolator),pointer :: interp1d_slope_right
    type(arb_deg_1d_interpolator),pointer :: interp1d_slope_bottom
    type(arb_deg_1d_interpolator),pointer :: interp1d_slope_top
    sll_int32 :: ierr
    sll_int32 :: tmp1
    sll_int32 :: tmp2
    sll_int64 :: bc_selector
    sll_int32 :: sz_slope_left,sz_slope_right,sz_slope_bottom,sz_slope_top
    ! only for troubleshooting
!!$    type(sll_time_mark) :: tm
!!$    sll_real64 :: time


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

   
    interpolator%spline_degree1 = spline_degree1
    interpolator%spline_degree2 = spline_degree2
    interpolator%eta1_min = eta1_min
    interpolator%eta1_max = eta1_max
    interpolator%eta2_min = eta2_min
    interpolator%eta2_max = eta2_max
    interpolator%bc_left  = bc_left
    interpolator%bc_right = bc_right
    interpolator%bc_left  = bc_bottom
    interpolator%bc_right = bc_top
    interpolator%bc_selector = bc_selector
    interpolator%num_pts1 = num_pts1
    interpolator%num_pts2 = num_pts2

    SLL_ALLOCATE(interpolator%slope_left  (num_pts2),ierr)
    SLL_ALLOCATE(interpolator%slope_right (num_pts2),ierr)
    SLL_ALLOCATE(interpolator%slope_bottom(num_pts1),ierr)
    SLL_ALLOCATE(interpolator%slope_top   (num_pts1),ierr)

    ! tmp1 and tmp2 is the maximun (not absolue) for the size of coefficients
    select case (bc_selector)
    case (0) ! 1. periodic-periodic
       SLL_ALLOCATE( interpolator%knots1(2*spline_degree1+2),ierr )
       SLL_ALLOCATE( interpolator%knots2(2*spline_degree2+2),ierr )
       !  tmp1 et tmp2 sont des limites suffisantes mais pas absolu 
       tmp1 = num_pts1 + 4*spline_degree1!*num_pts1 !+ 2*spline_degree1
       tmp2 = num_pts2 + 4*spline_degree2!*num_pts2 !+ 2*spline_degree2
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

    case (9) ! 2. dirichlet-left, dirichlet-right, periodic
       SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
       SLL_ALLOCATE( interpolator%knots2(2*spline_degree2+2),ierr )
       tmp1 = num_pts1+ 4*spline_degree1!*num_pts1! + spline_degree1 !- 1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2! + 2*spline_degree2
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)
       if (present(slope_left)) then 
          sz_slope_left = size(slope_left)
          if ( sz_slope_left .ne. interpolator%num_pts2 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_left must have the size of numbers of pts in direction 2 '
             stop
          end if
          
          call interp1d_slope_left%initialize(&
               interpolator%num_pts2, &
               interpolator%eta2_min, &
               interpolator%eta2_max, &
               interpolator%bc_bottom, &
               interpolator%bc_top, &
               interpolator%spline_degree2 )
          
          call interp1d_slope_left%compute_interpolants( &
               slope_left(1:sz_slope_left))
          
          interpolator%slope_left(1:sz_slope_left) = interp1d_slope_left%coeff_splines(1:sz_slope_left)
          call delete(interp1d_slope_left)
       else
          interpolator%slope_left(:) = 0.0_f64
       end if

       if (present(slope_right)) then 
          sz_slope_right = size(slope_right)
          if ( sz_slope_right .ne. interpolator%num_pts2 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_right must have the size of numbers of pts in direction 2 '
             stop
          end if
          
          call interp1d_slope_right%initialize(&
               interpolator%num_pts2, &
               interpolator%eta2_min, &
               interpolator%eta2_max, &
               interpolator%bc_bottom, &
               interpolator%bc_top, &
               interpolator%spline_degree2 )
          
          call interp1d_slope_right%compute_interpolants( &
               slope_right(1:sz_slope_right))
          
          interpolator%slope_right(1:sz_slope_right) = interp1d_slope_right%coeff_splines(1:sz_slope_right)
          call delete(interp1d_slope_right)
       else
          interpolator%slope_right(:) = 0.0_f64
       end if

    case (576) ! 3. periodic, dirichlet-bottom, dirichlet-top
       SLL_ALLOCATE( interpolator%knots1(2*spline_degree1+2),ierr )
       SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )
       tmp1 = num_pts1+ 4*spline_degree1!*num_pts1! + 2*spline_degree1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2 + spline_degree2 !- 1
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)

       if (present(slope_bottom)) then 
          sz_slope_bottom = size(slope_bottom)
          if ( sz_slope_bottom .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_bottom must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          call interp1d_slope_bottom%initialize(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_slope_bottom%compute_interpolants( &
               slope_bottom(1:sz_slope_bottom))
          
          interpolator%slope_bottom(1:sz_slope_bottom) = interp1d_slope_bottom%coeff_splines(1:sz_slope_bottom)
          call delete(interp1d_slope_bottom)
       else
          interpolator%slope_bottom(:) = 0.0_f64
       end if

       if (present(slope_top)) then 
          sz_slope_top = size(slope_top)
          if ( sz_slope_top .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_top must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          call  interp1d_slope_top%initialize(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1 )
          
          call interp1d_slope_top%compute_interpolants(&
               slope_top(1:sz_slope_top))
          
          interpolator%slope_top(1:sz_slope_top) = interp1d_slope_top%coeff_splines(1:sz_slope_top)
          call delete(interp1d_slope_top)
       else
          interpolator%slope_top(:) = 0.0_f64
       end if

    case (585) ! 4. dirichlet in all sides
       SLL_ALLOCATE( interpolator%knots1(num_pts1+2*spline_degree1),ierr )
       SLL_ALLOCATE( interpolator%knots2(num_pts2+2*spline_degree2),ierr )
       tmp1 = num_pts1+ 4*spline_degree1!*num_pts1! + spline_degree1 !- 1
       tmp2 = num_pts2+ 4*spline_degree2!*num_pts2! + spline_degree2 !- 1
       SLL_ALLOCATE( interpolator%coeff_splines(tmp1,tmp2),ierr)


       if (present(slope_left)) then 
          sz_slope_left = size(slope_left)
          if ( sz_slope_left .ne. interpolator%num_pts2 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_left must have the size of numbers of pts in direction 2 '
             stop
          end if
          
          call interp1d_slope_left%initialize(&
               interpolator%num_pts2, &
               interpolator%eta2_min, &
               interpolator%eta2_max, &
               interpolator%bc_bottom, &
               interpolator%bc_top, &
               interpolator%spline_degree2,&
               slope_left(1),&
               slope_left(sz_slope_left))
          
          call interp1d_slope_left%compute_interpolants( &
               slope_left(1:sz_slope_left))
          
          interpolator%slope_left(1:sz_slope_left) = interp1d_slope_left%coeff_splines(1:sz_slope_left)
          call delete(interp1d_slope_left)
       else
          interpolator%slope_left(:) = 0.0_f64
       end if
       
       if (present(slope_right)) then 
          sz_slope_right = size(slope_right)
          if ( sz_slope_right .ne. interpolator%num_pts2 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_right must have the size of numbers of pts in direction 2 '
             stop
          end if
          
          call interp1d_slope_right%initialize(&
               interpolator%num_pts2, &
               interpolator%eta2_min, &
               interpolator%eta2_max, &
               interpolator%bc_bottom, &
               interpolator%bc_top, &
               interpolator%spline_degree2,&
               slope_right(1),&
               slope_right(sz_slope_right))
          
          
          call interp1d_slope_right%compute_interpolants( &
               slope_right(1:sz_slope_right))
          
          interpolator%slope_right(1:sz_slope_right) = interp1d_slope_right%coeff_splines(1:sz_slope_right)
          call delete(interp1d_slope_right)
       else
          interpolator%slope_right(:) = 0.0_f64
       end if
       
       
       if (present(slope_bottom)) then 
          sz_slope_bottom = size(slope_bottom)
          if ( sz_slope_bottom .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_bottom must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          call interp1d_slope_bottom%initialize(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1,&
               slope_bottom(1),&
               slope_bottom(sz_slope_bottom))
         
          
          call interp1d_slope_bottom%compute_interpolants( &
               slope_bottom(1:sz_slope_bottom))
          
          interpolator%slope_bottom(1:sz_slope_bottom) = interp1d_slope_bottom%coeff_splines(1:sz_slope_bottom)
          call delete(interp1d_slope_bottom)
       else
          interpolator%slope_bottom(:) = 0.0_f64
       end if
       
       if (present(slope_top)) then 
          sz_slope_top = size(slope_top)
          if ( sz_slope_top .ne. interpolator%num_pts1 ) then 
             print*, ' problem in the initialization of arb_deg_spline 2d'
             print*, ' slope_top must have the size of numbers of pts in direction 1 '
             stop
          end if
          
          call interp1d_slope_top%initialize(&
               interpolator%num_pts1, &
               interpolator%eta1_min, &
               interpolator%eta1_max, &
               interpolator%bc_left, &
               interpolator%bc_right, &
               interpolator%spline_degree1,&
               slope_top(1),&
               slope_top(sz_slope_top))
          
          
          call interp1d_slope_top%compute_interpolants( &
               slope_top(1:sz_slope_top))
          
          interpolator%slope_top(1:sz_slope_top) = interp1d_slope_top%coeff_splines(1:sz_slope_top)
          call delete(interp1d_slope_top)
       else
          interpolator%slope_top(:) = 0.0_f64
       end if
       


    case default
       print *, 'initialize_ad2d_interpolator: BC combination not implemented.'
    end select

    interpolator%coeff_splines(:,:) = 0.0_f64
    SLL_ALLOCATE( interpolator%t1(tmp1),ierr)!num_pts1*num_pts1),ierr)
    !+ 2*(spline_degree1 + 1)), ierr)
    SLL_ALLOCATE( interpolator%t2(tmp2),ierr)!num_pts2*num_pts2),ierr) 
    !+ 2*(spline_degree2 + 1)), ierr)

    interpolator%t1(:) = 0.0_f64
    interpolator%t2(:) = 0.0_f64
    !print*,'SIZE',  num_pts1 + 2*(spline_degree1 + 1)
    
  end subroutine !initialize_ad2d_interpolator

!!$  subroutine compute_interpolants_ad2d( &
!!$    interpolator, &
!!$    data_array )
!!$
!!$    class(arb_deg_2d_interpolator), intent(inout) :: interpolator
!!$    sll_real64, dimension(:,:), intent(in)         :: data_array
!!$
!!$    print *, 'compute interpolants not implemented'
!!$
!!$  end subroutine compute_interpolants_ad2d

#ifdef STDF95
  subroutine arbitrary_degree_spline_interp2_set_coefficients( &
#else
  subroutine set_coefficients_ad2d( &
#endif
   interpolator, &
   coeffs_1d, &
   coeffs_2d )

#ifdef STDF95
   type (arb_deg_2d_interpolator), intent(inout)  :: interpolator
#else
   class(arb_deg_2d_interpolator), intent(inout)  :: interpolator
#endif
   sll_real64, dimension(:), intent(in), optional :: coeffs_1d
   sll_real64, dimension(:,:), intent(in), optional :: coeffs_2d
   sll_int32 :: sp_deg1
   sll_int32 :: sp_deg2
   sll_int32 :: num_cells1
   sll_int32 :: num_cells2
   sll_int32 :: tmp1, tmp2
   sll_int32 :: i, j
   sll_real64 :: eta1_min, eta1_max
   sll_real64 :: eta2_min, eta2_max
   sll_real64 :: delta1
   sll_real64 :: delta2
   sll_int32  ::  nb_spline_eta1
   sll_int32  ::  nb_spline_eta2
   sll_real64 :: eta1
   sll_real64 :: eta2

   if(present(coeffs_2d)) then
      print *, 'set_coefficients_ad2d(), ERROR: this function has not ', &
           'implmented the option to set the coeffcients provided as a 2d ', &
           'array.'
      stop
   end if

   sp_deg1    = interpolator%spline_degree1
   sp_deg2    = interpolator%spline_degree2
   num_cells1 = interpolator%num_pts1 - 1
   num_cells2 = interpolator%num_pts2 - 1
   eta1_min = interpolator%eta1_min
   eta2_min = interpolator%eta2_min
   eta1_max = interpolator%eta1_max
   eta2_max = interpolator%eta2_max
   delta1 = (eta1_max - eta1_min)/num_cells1
   delta2 = (eta2_max - eta2_min)/num_cells2

   tmp1 = (sp_deg1 + 1)/2
   tmp2 = (sp_deg2 + 1)/2
   !print*, tmp1,tmp2
   ! The interpretation and further filling of the spline coefficients array
   ! depends on the boundary conditions.
   select case (interpolator%bc_selector)
   case(0) ! periodic-periodic
      
      interpolator%size_coeffs1=  num_cells1 + sp_deg1
      interpolator%size_coeffs2=  num_cells2 + sp_deg2
      interpolator%size_t1 = 2*sp_deg1 + num_cells1 +1 
      interpolator%size_t2 = 2*sp_deg2 + num_cells2 +1
      ! allocation and definition of knots
      do i = -sp_deg1, num_cells1 + sp_deg1
         interpolator%t1( i + sp_deg1 + 1 ) = eta1_min + i*delta1
      end do
      
      do i = -sp_deg2, num_cells2 + sp_deg2
         interpolator%t2( i + sp_deg2 + 1 ) = eta2_min + i*delta2
      end do

      do i = 1,num_cells1
         do j = 1,num_cells2
            interpolator%coeff_splines(i,j) = &
                 coeffs_1d( i + num_cells1 *(j-1) )
         end do
      end do

      do j = 1, sp_deg2
         do i = 1,num_cells1
            
            interpolator%coeff_splines(i ,num_cells2 + j ) = &
                 coeffs_1d(i+num_cells1*(j-1))
         end do
      end do
      do i = 1, sp_deg1
         do j = 1,num_cells2
            
            interpolator%coeff_splines(num_cells1 + i ,j) = &
                 coeffs_1d(i+num_cells1 *(j-1) )
            !nb_spline_eta1 - (tmp1 -i)  + nb_spline_eta1 *(j-1) )
         end do
      end do

      do i= 1,sp_deg1
         do j=1,sp_deg2
            
            interpolator%coeff_splines(num_cells1 +  i ,num_cells2 + j) = &
                 interpolator%coeff_splines(i,j)!(sp_deg1-(i-1)),(sp_deg2-(j-1)))
         end do
      end do

   case (9) ! 2. dirichlet-left, dirichlet-right, periodic
      interpolator%size_coeffs1=  num_cells1 + sp_deg1
      interpolator%size_coeffs2=  num_cells2 + sp_deg2
      interpolator%size_t1 = 2*sp_deg1 + num_cells1 + 1
      interpolator%size_t2 = 2*sp_deg2 + num_cells2 + 1
      nb_spline_eta1 = num_cells1 + sp_deg1 - 2
      nb_spline_eta2 = num_cells2
      ! allocation and definition of knots
      do i = - sp_deg2, num_cells2 + sp_deg2
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
       
       
       do i = 1 ,nb_spline_eta1
          do j = 1,nb_spline_eta2
             interpolator%coeff_splines(i+1,j) = &
                  coeffs_1d(i+nb_spline_eta1*(j-1))
         end do
      end do
     
      
      do j = 1, sp_deg2
         do i = 1,nb_spline_eta1

            interpolator%coeff_splines(i + 1 ,nb_spline_eta2 + j ) = &
                 coeffs_1d(i+nb_spline_eta1*(j-1))
         end do
      end do

      ! achtung ! normaly interpolator%slope_left(:) and interpolator%slope_right(:)
      interpolator%coeff_splines(1,:) = 0.0_8
      interpolator%coeff_splines(nb_spline_eta1+2,:) = 0.0_8
      
   case(576)!3. periodic, dirichlet-bottom, dirichlet-top
      interpolator%size_coeffs1=  num_cells1 + sp_deg1 
      interpolator%size_coeffs2=  num_cells2 + sp_deg2
      interpolator%size_t1 = 2.0_f64*sp_deg1 + num_cells1 + 1
      interpolator%size_t2 = 2.0_f64*sp_deg2 + num_cells2 + 1
      nb_spline_eta1 = num_cells1
      nb_spline_eta2 = num_cells2 + sp_deg2 - 2
      
      ! allocation and definition of knots
      do i = - sp_deg1, nb_spline_eta1 + sp_deg1

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
      do i = num_cells2 + sp_deg2 + 2, num_cells2 + 1 + 2*sp_deg2
         interpolator%t2(i) = eta2
      enddo
       

      
      do i = 1 , nb_spline_eta1
         do j = 1,nb_spline_eta2
            
            interpolator%coeff_splines(i ,j+1) = &
                 coeffs_1d(i+nb_spline_eta1 *(j-1) )
         end do
      end do
      
      do i = 1, sp_deg1
         do j = 1,nb_spline_eta2
            
            interpolator%coeff_splines(nb_spline_eta1 + i ,j+1) = &
                 coeffs_1d(i+nb_spline_eta1 *(j-1) )
            !nb_spline_eta1 - (tmp1 -i)  + nb_spline_eta1 *(j-1) )
         end do
      end do
      ! achtung ! normaly interpolator%slope_bottom(:) and interpolator%slope_top(:)
      interpolator%coeff_splines(:,1) = 0.0_8
      interpolator%coeff_splines(:,nb_spline_eta2+2) = 0.0_8
      
   case(585) ! 4. dirichlet in all sides
      interpolator%size_coeffs1=  num_cells1 + sp_deg1
      interpolator%size_coeffs2=  num_cells2 + sp_deg2
      interpolator%size_t1 = 2.0_f64*sp_deg1 + num_cells1 + 1
      interpolator%size_t2 = 2.0_f64*sp_deg2 + num_cells2 + 1
      nb_spline_eta1 = num_cells1 + sp_deg1 - 2
      nb_spline_eta2 = num_cells2 + sp_deg2 - 2
      
      ! allocation and definition of knots

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
      
      ! achtung ! normaly interpolator%slope_left(:) and interpolator%slope_right(:)
      ! achtung ! normaly interpolator%slope_bottom(:) and interpolator%slope_top(:)
      interpolator%coeff_splines(:,:) = 0.0_8
      ! allocation coefficient spline
      do i = 1,nb_spline_eta1
         do j = 1,nb_spline_eta2
            
            interpolator%coeff_splines(i+1,j+1) = &
                 coeffs_1d( i + nb_spline_eta1 *(j-1))
         end do
      end do
      
   case default
      print *, 'arbitrary_degree_spline_2d() error: set_spline_coefficients ',&
           'not recognized.'
      stop
   end select
 end subroutine !set_coefficients_ad2d

#ifdef STDF95
  subroutine arbitrary_degree_spline_interp2d_compute_interpolants( &
#else
  subroutine compute_interpolants_ad2d( &
#endif
    interpolator, &
    data_array, &
    eta1_coords, &
    size_eta1_coords, &
    eta2_coords, &
    size_eta2_coords )

#ifdef STDF95
    type (arb_deg_2d_interpolator), intent(inout)  :: interpolator
#else
    class(arb_deg_2d_interpolator), intent(inout)  :: interpolator
#endif
    sll_real64, dimension(:,:), intent(in)         :: data_array
    sll_real64, dimension(:), intent(in),optional  :: eta1_coords
    sll_real64, dimension(:), intent(in),optional  :: eta2_coords
    sll_int32, intent(in),optional                 :: size_eta1_coords
    sll_int32, intent(in),optional                 :: size_eta2_coords
    sll_real64, dimension(:),pointer               :: point_location_eta1
    sll_real64, dimension(:),pointer               :: point_location_eta2
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
       
       SLL_ALLOCATE(point_location_eta1(sz1),ierr)
       SLL_ALLOCATE(point_location_eta2(sz2),ierr)
       point_location_eta1 = eta1_coords
       point_location_eta2 = eta2_coords

    else ! size depends on BC combination, filled out at initialization.

       select case (interpolator%bc_selector)
       case (0) ! 1. periodic-periodic
          sz1 = interpolator%num_pts1!-1
          sz2 = interpolator%num_pts2!-1
          
       case (9) ! 2. dirichlet-left, dirichlet-right, periodic
          sz1 = interpolator%num_pts1
          sz2 = interpolator%num_pts2!-1
          
       case (576) ! 3. periodic, dirichlet-bottom, dirichlet-top
          sz1 = interpolator%num_pts1!-1
          sz2 = interpolator%num_pts2
       
       case (585) ! 4. dirichlet in all sides
          sz1 = interpolator%num_pts1
          sz2 = interpolator%num_pts2
   
       case default
          print *, 'compute_interpolants_ad2d():BC combination not implemented.'
       end select

       delta_eta1 = (interpolator%eta1_max - interpolator%eta1_min)&
            /(interpolator%num_pts1 -1)
       delta_eta2 = (interpolator%eta2_max - interpolator%eta2_min)&
            /(interpolator%num_pts2 -1)
       SLL_ALLOCATE(point_location_eta1(sz1),ierr)
       SLL_ALLOCATE(point_location_eta2(sz2),ierr)
       
       do i = 1,sz1
          point_location_eta1(i) = interpolator%eta1_min + delta_eta1*(i-1)
       end do
       do i = 1,sz2
          point_location_eta2(i) = interpolator%eta2_min + delta_eta2*(i-1)
       end do
    end if
    
    !print*, 'point location1',point_location_eta1
    !print*, 'point location2',point_location_eta2
    SLL_ASSERT(sz1 .le. interpolator%num_pts1* interpolator%num_pts1)
    SLL_ASSERT(sz2 .le. interpolator%num_pts2* interpolator%num_pts2)
    SLL_ASSERT(size(data_array,1) .ge. sz1)
    SLL_ASSERT(size(data_array,2) .ge. sz2)
    SLL_ASSERT(size(point_location_eta1)  .ge. sz1)
    SLL_ASSERT(size(point_location_eta2)  .ge. sz2)
    
    order1  = interpolator%spline_degree1 + 1
    order2  = interpolator%spline_degree2 + 1
    period1 = interpolator%eta1_max - interpolator%eta1_min
    period2 = interpolator%eta2_max - interpolator%eta2_min
    
   ! print*, 'pointlocation',point_location_eta2
    select case (interpolator%bc_selector)
    case (0) ! periodic-periodic
       interpolator%size_coeffs1 = sz1!+1
       interpolator%size_coeffs2 = sz2!+1
       interpolator%size_t1 = order1 + sz1 !+ 1
       interpolator%size_t2 = order2 + sz2 !+ 1 

       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2

       call spli2d_perper( &
            period1, sz1, order1, point_location_eta1(1:sz1-1), & !+1
            period2, sz2, order2, point_location_eta2(1:sz2-1), & !+1
            data_array(1:sz1-1,1:sz2-1), interpolator%coeff_splines(1:sz1,1:sz2),&!(1:sz1+1,1:sz2+1),&
            interpolator%t1(1:order1 + sz1 ), &!+ 1), &
            interpolator%t2(1:order2 + sz2 ))!+ 1) )
   
       
    case (9) ! 2. dirichlet-left, dirichlet-right, periodic
       interpolator%size_coeffs1 = sz1
       interpolator%size_coeffs2 = sz2!+1
       interpolator%size_t1 = order1 + sz1
       interpolator%size_t2 = order2 + sz2 !+ 1
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       call spli2d_dirper( sz1, order1, point_location_eta1, &
            period2, sz2, order2, point_location_eta2(1:sz2-1), & !+1
            data_array(1:sz1,1:sz2-1), interpolator%coeff_splines(1:sz1,1:sz2),&!+1
            interpolator%t1(1:sz1+order1), &
            interpolator%t2(1:sz2+order2) ) !+1

       ! boundary condition non homogene
       interpolator%coeff_splines(1,1:sz2)   = interpolator%slope_left(1:sz2)
       interpolator%coeff_splines(sz1,1:sz2) = interpolator%slope_right(1:sz2)
  
    case(576) !  3. periodic, dirichlet-bottom, dirichlet-top
       interpolator%size_coeffs1 = sz1!+1
       interpolator%size_coeffs2 = sz2
       interpolator%size_t1 = order1 + sz1 !+ 1
       interpolator%size_t2 = order2 + sz2 
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       call spli2d_perdir( period1, sz1, order1, point_location_eta1(1:sz1-1), & !+ 1
            sz2, order2, point_location_eta2, &
            data_array(1:sz1-1,1:sz2), interpolator%coeff_splines(1:sz1,1:sz2),& !+ 1
            interpolator%t1(1:sz1+order1), & ! + 1
            interpolator%t2(1:sz2+order2) )

       ! boundary condition non homogene
       interpolator%coeff_splines(1:sz1,1)   = interpolator%slope_bottom(1:sz1)
       interpolator%coeff_splines(1:sz1,sz2) = interpolator%slope_top(1:sz1)
       
    case (585) ! 4. dirichlet in all sides
       !print*, 'her'
       interpolator%size_coeffs1 = sz1
       interpolator%size_coeffs2 = sz2
       interpolator%size_t1 = order1 + sz1 
       interpolator%size_t2 = order2 + sz2 
 
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       call spli2d_custom( sz1, order1, point_location_eta1, &
            sz2, order2, point_location_eta2, &
            data_array(1:sz1,1:sz2), interpolator%coeff_splines(1:sz1,1:sz2),&
            interpolator%t1(1:sz1+order1), &
            interpolator%t2(1:sz2+order2) )

       ! boundary condition non homogene
       interpolator%coeff_splines(1,1:sz2)   = interpolator%slope_left(1:sz2)
       interpolator%coeff_splines(sz1,1:sz2) = interpolator%slope_right(1:sz2)
       ! boundary condition non homogene
       interpolator%coeff_splines(1:sz1,1)   = interpolator%slope_bottom(1:sz1)
       interpolator%coeff_splines(1:sz1,sz2) = interpolator%slope_top(1:sz1)

    end select

    SLL_DEALLOCATE(point_location_eta2,ierr)
    SLL_DEALLOCATE(point_location_eta1,ierr)
  end subroutine !compute_interpolants_ad2d

#ifdef STDF95
  function arbitrary_degree_spline_interp2d_interpolate_value( &
#else
  function interpolate_value_ad2d( &
#endif
    interpolator, &
    eta1, &
    eta2 ) result(val)

#ifdef STDF95
    type (arb_deg_2d_interpolator), intent(in)  :: interpolator
#else
    class(arb_deg_2d_interpolator), intent(in)  :: interpolator
#endif
    sll_real64, intent(in)         :: eta1
    sll_real64, intent(in)         :: eta2
    sll_real64                     :: val
    sll_int32 :: size_coeffs1
    sll_int32 :: size_coeffs2
    sll_real64 :: bvalue2d
    sll_real64 :: res1,res2
 

    size_coeffs1 = interpolator%size_coeffs1
    size_coeffs2 = interpolator%size_coeffs2

    res1 = eta1
    res2 = eta2


    select case (interpolator%bc_selector)
    case (0) ! periodic-periodic

       if( res1 < interpolator%eta1_min ) then
          res1 = res1+interpolator%eta1_max-interpolator%eta1_min
       else if( res1 >  interpolator%eta1_max ) then
          res1 = res1+interpolator%eta1_min-interpolator%eta1_max
       end if
        if( res2 < interpolator%eta2_min ) then
          res2 = res2+interpolator%eta2_max-interpolator%eta2_min
       else if( res2 >  interpolator%eta2_max ) then
          res2 = res2+interpolator%eta2_min-interpolator%eta2_max
       end if
    case (9) ! 2. dirichlet-left, dirichlet-right, periodic

       if( res2 < interpolator%eta2_min ) then
          res2 = res2+interpolator%eta2_max-interpolator%eta2_min
       else if( res2 >  interpolator%eta2_max ) then
          res2 = res2+interpolator%eta2_min-interpolator%eta2_max
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
          res1 = res1+interpolator%eta1_max-interpolator%eta1_min
       else if( res1 >  interpolator%eta1_max ) then
          res1 = res1+interpolator%eta1_min-interpolator%eta1_max
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
       
    case (585) ! dirichlet-dirichlet 
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

    end select

    val = bvalue2d( &
         res1, &
         res2, &
         size_coeffs1, &
         interpolator%spline_degree1+1, &
         size_coeffs2, &
         interpolator%spline_degree2+1, &
         interpolator%coeff_splines(1:size_coeffs1,1:size_coeffs2), &
         interpolator%t1(1:interpolator%size_t1), &
         interpolator%t2(1:interpolator%size_t2))

  end function interpolate_value_ad2d


#ifdef STDF95
  function arbitrary_degree_spline_interp2d_interpolate_derivative1( &
#else
  function interpolate_derivative1_ad2d( &
#endif
    interpolator, &
    eta1, &
    eta2 ) result(val)

#ifdef STDF95
    type (arb_deg_2d_interpolator), intent(in)  :: interpolator
#else
    class(arb_deg_2d_interpolator), intent(in)  :: interpolator
#endif
    sll_real64, intent(in)         :: eta1
    sll_real64, intent(in)         :: eta2
    sll_real64                     :: val
    sll_int32 :: size_coeffs1
    sll_int32 :: size_coeffs2
    sll_real64 :: dvalue2d
    sll_real64 :: res1,res2
    
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
          res1 = res1+interpolator%eta1_max-interpolator%eta1_min
       else if( res1 >  interpolator%eta1_max ) then
          res1 = res1+interpolator%eta1_min-interpolator%eta1_max
       end if
        if( res2 < interpolator%eta2_min ) then
          res2 = res2+interpolator%eta2_max-interpolator%eta2_min
       else if( res2 >  interpolator%eta2_max ) then
          res2 = res2+interpolator%eta2_min-interpolator%eta2_max
       end if
    case (9) ! 2. dirichlet-left, dirichlet-right, periodic

       if( res2 < interpolator%eta2_min ) then
          res2 = res2+interpolator%eta2_max-interpolator%eta2_min
       else if( res2 >  interpolator%eta2_max ) then
          res2 = res2+interpolator%eta2_min-interpolator%eta2_max
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
!!$       if ( res1 .ge. interpolator%eta1_max ) then 
!!$          res1 = res1 -(interpolator%eta1_max-interpolator%eta1_min)
!!$       end if
       if( res1 < interpolator%eta1_min ) then
          res1 = res1+interpolator%eta1_max-interpolator%eta1_min
       else if( res1 >  interpolator%eta1_max ) then
          res1 = res1+interpolator%eta1_min-interpolator%eta1_max
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
       
    case (585) ! dirichlet-dirichlet 
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
       end select
    
    val = dvalue2d( &
         res1, &
         res2, &
         size_coeffs1, &
         interpolator%spline_degree1+1, &
         size_coeffs2, &
         interpolator%spline_degree2+1, &
         interpolator%coeff_splines(1:size_coeffs1,1:size_coeffs2), &
         interpolator%t1(1:interpolator%size_t1), &
         interpolator%t2(1:interpolator%size_t2),&
         1,0)
    
  end function interpolate_derivative1_ad2d
  

#ifdef STDF95
  function arbitrary_degree_spline_interp2d_interpolate_derivative2( &
#else
  function interpolate_derivative2_ad2d( &
#endif
    interpolator, &
    eta1, &
    eta2 ) result(val)

#ifdef STDF95
    type (arb_deg_2d_interpolator), intent(in)  :: interpolator
#else
    class(arb_deg_2d_interpolator), intent(in)  :: interpolator
#endif
    sll_real64, intent(in)         :: eta1
    sll_real64, intent(in)         :: eta2
    sll_real64                     :: val
    sll_int32 :: size_coeffs1
    sll_int32 :: size_coeffs2
    sll_real64 :: dvalue2d
    sll_real64 :: res1,res2

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
          res1 = res1+interpolator%eta1_max-interpolator%eta1_min
       else if( res1 >  interpolator%eta1_max ) then
          res1 = res1+interpolator%eta1_min-interpolator%eta1_max
       end if
        if( res2 < interpolator%eta2_min ) then
          res2 = res2+interpolator%eta2_max-interpolator%eta2_min
       else if( res2 >  interpolator%eta2_max ) then
          res2 = res2+interpolator%eta2_min-interpolator%eta2_max
       end if
    case (9) ! 2. dirichlet-left, dirichlet-right, periodic

        if( res2 < interpolator%eta2_min ) then
          res2 = res2+interpolator%eta2_max-interpolator%eta2_min
       else if( res2 >  interpolator%eta2_max ) then
          res2 = res2+interpolator%eta2_min-interpolator%eta2_max
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
!!$       if ( res1 .ge. interpolator%eta1_max ) then 
!!$          res1 = res1 -(interpolator%eta1_max-interpolator%eta1_min)
!!$       end if
       if( res1 < interpolator%eta1_min ) then
          res1 = res1+interpolator%eta1_max-interpolator%eta1_min
       else if( res1 >  interpolator%eta1_max ) then
          res1 = res1+interpolator%eta1_min-interpolator%eta1_max
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
       
    case (585) ! dirichlet-dirichlet 
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
    end select
    
    val = dvalue2d( &
         res1, &
         res2, &
         size_coeffs1, &
         interpolator%spline_degree1+1, &
         size_coeffs2, &
         interpolator%spline_degree2+1, &
         interpolator%coeff_splines(1:size_coeffs1,1:size_coeffs2), &
         interpolator%t1(1:interpolator%size_t1), &
         interpolator%t2(1:interpolator%size_t2),&
         0,1)

  end function !interpolate_derivative2_ad2d

#ifdef STDF95
  function arbitrary_degree_spline_interp2d_interpolate_array( &
#else
  function interpolate_array_ad2d( &
#endif
    this, &
    num_points1, &
    num_points2, &
    data_in, &
    eta1, &
    eta2 ) result(res)

#ifdef STDF95
    type (arb_deg_2d_interpolator), intent(in)  :: this
#else
    class(arb_deg_2d_interpolator), intent(in)  :: this
#endif
    sll_real64,  dimension(:,:), intent(in)         :: eta1
    sll_real64,  dimension(:,:), intent(in)         :: eta2
    sll_real64, dimension(:,:), intent(in)         :: data_in
    sll_int32, intent(in)         :: num_points1
    sll_int32, intent(in)         :: num_points2

    sll_real64, dimension(num_points1,num_points2) :: res
 
    print *, 'interpolate_array_ad2d: not implemented'
    res = -1000000._f64
  end function !interpolate_array_ad2d
  
#ifdef STDF95
  function arbitrary_degree_spline_interp2d_interpolate_2d_array_disp( &
#else
  function interpolate_2d_array_disp_ad2d( &
#endif
       this,        &
       num_points1, &
       num_points2, &
       data_in,     &
       alpha1,      &
       alpha2) result(res)
      
#ifdef STDF95
    type (arb_deg_2d_interpolator), intent(in)    :: this
#else
    class(arb_deg_2d_interpolator), intent(in)    :: this
#endif
    sll_int32, intent(in)                          :: num_points1  
    sll_int32, intent(in)                          :: num_points2 
    sll_real64, dimension(:,:), intent(in)         :: data_in
    sll_real64, dimension(:,:), intent(in)         :: alpha1
    sll_real64, dimension(:,:), intent(in)         :: alpha2  
    sll_real64, dimension(num_points1,num_points2) :: res
    
    print *, 'interpolate_2d_array_disp_ad2d: not implemented.'
    res = -1000000._f64
  end function !interpolate_2d_array_disp_ad2d
    
   
#ifdef STDF95
  function arbitrary_degree_spline_interp2d_get_coefficients(interpolator)
    type (arb_deg_2d_interpolator), intent(in)    :: interpolator
    sll_real64, dimension(:,:), pointer           :: arbitrary_degree_spline_interp2d_get_coefficients     

    arbitrary_degree_spline_interp2d_get_coefficients => interpolator%coeff_splines
  end function arbitrary_degree_spline_interp2d_get_coefficients
#else 
  function get_coefficients_ad2d(interpolator)
    class(arb_deg_2d_interpolator), intent(in)    :: interpolator
    sll_real64, dimension(:,:), pointer           :: get_coefficients_ad2d     

    get_coefficients_ad2d => interpolator%coeff_splines
  end function get_coefficients_ad2d
#endif
  
end module sll_arbitrary_degree_spline_interpolator_2d_module
