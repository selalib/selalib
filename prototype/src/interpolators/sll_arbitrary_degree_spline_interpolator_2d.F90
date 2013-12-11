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
#ifdef STDF95
use sll_boundary_condition_descriptors
!use sll_constants
#else
use sll_module_interpolators_2d_base
#endif
use sll_utilities
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
     logical    :: coefficients_set = .false.
#ifndef STDF95
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
  end subroutine delete_arbitrary_degree_2d_interpolator



  ! -----------------------------------------------
  ! This subroutine allocate the type of interpolator
  !    the  arbitrary_spline_interp2d
  ! -----------------------------------------------
  function new_arbitrary_degree_spline_interp2d(&
       num_pts1,  & ! number of points in direction 1
       num_pts2,  & ! number of points in direction 2
       eta1_min,  & ! minimun of direction 1
       eta1_max,  & ! maximun of direction 1
       eta2_min,  & ! minimun of direction 2
       eta2_max,  & ! maximun of direction 2
       bc_left,   & ! boundary condition on the left of direction 1
       bc_right,  & ! boundary condition on the rigth of direction 1
       bc_bottom, & ! boundary condition on the bottom of direction 2 
       bc_top,    & ! boundary condition on the top of direction 1
       spline_degree1, & ! degree spline in direction 1
       spline_degree2 ) result(arbitrary_degree_spline_interp2d)  ! degree spline in direction 2
       
    ! INPUT VARIABLES

    sll_int32 , intent(in) :: num_pts1
    sll_int32 , intent(in) :: num_pts2
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_real64, intent(in) :: eta2_min
    sll_real64, intent(in) :: eta2_max
    sll_int32 , intent(in) :: bc_left
    sll_int32 , intent(in) :: bc_right
    sll_int32 , intent(in) :: bc_bottom
    sll_int32 , intent(in) :: bc_top
    sll_int32 , intent(in) :: spline_degree1
    sll_int32 , intent(in) :: spline_degree2

    ! LOCAL VARIABLES
#ifdef STDF95
    type (arb_deg_2d_interpolator), pointer :: arbitrary_degree_spline_interp2d
#else
    class(arb_deg_2d_interpolator), pointer :: arbitrary_degree_spline_interp2d
#endif
    sll_int32 :: ierr
    
    SLL_ALLOCATE(arbitrary_degree_spline_interp2d, ierr)

    call initialize_ad2d_interpolator( &
         arbitrary_degree_spline_interp2d, &
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
         spline_degree2 )
  end function new_arbitrary_degree_spline_interp2d


  ! -----------------------------------------------
  ! This subroutine allocate the type of interpolator
  !    the  arbitrary_spline_interp2d
  ! -----------------------------------------------
#ifdef STDF95
  subroutine arbitrary_degree_spline_interp2d_initialize( &
#else
       subroutine initialize_ad2d_interpolator( &
#endif
       interpolator, & ! the interpolator as a pointer
       num_pts1,  & ! number of points in direction 1
       num_pts2,  & ! number of points in direction 2
       eta1_min,  & ! minimun of direction 1
       eta1_max,  & ! maximun of direction 1
       eta2_min,  & ! minimun of direction 2
       eta2_max,  & ! maximun of direction 2
       bc_left,   & ! boundary condition on the left of direction 1
       bc_right,  & ! boundary condition on the rigth of direction 1
       bc_bottom, & ! boundary condition on the bottom of direction 2 
       bc_top,    & ! boundary condition on the top of direction 1
       spline_degree1, & ! degree spline in direction 1
       spline_degree2 )  ! degree spline in direction 2
    

#ifdef STDF95
    type (arb_deg_2d_interpolator):: interpolator
#else
    class(arb_deg_2d_interpolator):: interpolator
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
    interpolator%bc_left  = bc_bottom
    interpolator%bc_right = bc_top
    interpolator%bc_selector = bc_selector
    interpolator%num_pts1 = num_pts1
    interpolator%num_pts2 = num_pts2

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

    case default
       print *, 'initialize_ad2d_interpolator: BC combination not implemented.'
    end select

    ! knots and coeff splines allocations 
    interpolator%coeff_splines(:,:) = 0.0_f64
    ! the minimun is to be of class C^0 everywhere on the knots
    ! i.e. each knot have multiplicity (spline_degree1+1) 
    ! so the maximun number of knots is num_pts1*(spline_degree1+1)
    SLL_ALLOCATE( interpolator%t1(num_pts1*(spline_degree1+1)),ierr)
    SLL_ALLOCATE( interpolator%t2(num_pts2*(spline_degree2+1)),ierr) 

    interpolator%t1(:) = 0.0_f64
    interpolator%t2(:) = 0.0_f64

  end subroutine !initialize_ad2d_interpolator



  ! -------------------------------------------------------------
  !  subroutine initializing the coefficients of splines
  !  in the cas of linearization of them i.e. if we have 
  !  a table in 1d corresponding of coefficients 2d
  ! 
  ! -------------------------------------------------------------
#ifdef STDF95
  subroutine arbitrary_degree_spline_interp2_set_coefficients( &
#else
  subroutine set_coefficients_ad2d( &
#endif
   interpolator, &
   coeffs_1d, &
   coeffs_2d,&
   coeff2d_size1,&
   coeff2d_size2,&
   knots1,&
   size_knots1,&
   knots2,&
   size_knots2)

#ifdef STDF95
   type (arb_deg_2d_interpolator), intent(inout)  :: interpolator
#else
   class(arb_deg_2d_interpolator), intent(inout)  :: interpolator
#endif
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
         
         interpolator%size_coeffs1 =  num_cells1 + sp_deg1
         interpolator%size_coeffs2 =  num_cells2 + sp_deg2
         interpolator%size_t1      =  2*sp_deg1 + num_cells1 +1 
         interpolator%size_t2      =  2*sp_deg2 + num_cells2 +1
         
         if ( size( coeffs_1d,1) .ne. num_cells1*num_cells2) then
            print*, 'Problem in set_coefficients in arbitrary_degree_spline_2d'
            print*, ' Problem with the size coeffs_1d must have the size equal to '
            print*, ' num_cells1*num_cells2=', num_cells1*num_cells2
            stop
         end if
         ! ------------------------------------------------------------
         ! allocation and definition of knots
         ! ------------------------------------------------------------
         
         do i = -sp_deg1, num_cells1 + sp_deg1
            interpolator%t1( i + sp_deg1 + 1 ) = eta1_min + i*delta1
         end do
         
         do i = -sp_deg2, num_cells2 + sp_deg2
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
            end do
         end do
         
         do i= 1,sp_deg1
            do j=1,sp_deg2
               
               interpolator%coeff_splines(num_cells1 +  i ,num_cells2 + j) = &
                    interpolator%coeff_splines(i,j)
            end do
         end do
      ! ------------------------------------------------------------
      case (9) ! 2. dirichlet-left, dirichlet-right, periodic
         
         
         
         interpolator%size_coeffs1 =  num_cells1 + sp_deg1
         interpolator%size_coeffs2 =  num_cells2 + sp_deg2
         interpolator%size_t1      =  2*sp_deg1 + num_cells1 + 1
         interpolator%size_t2      =  2*sp_deg2 + num_cells2 + 1
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
         
         ! ------------------------------------------------------------
         ! reorganization of spline coefficients 1D in coefficients 2D 
         ! ------------------------------------------------------------
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
         
         interpolator%coeff_splines(1,:) = 0.0_8
         interpolator%coeff_splines(nb_spline_eta1+2,:) = 0.0_8
         ! ------------------------------------------------------------
      case(576)!3. periodic, dirichlet-bottom, dirichlet-top
       
         
         interpolator%size_coeffs1 =  num_cells1 + sp_deg1 
         interpolator%size_coeffs2 =  num_cells2 + sp_deg2
         interpolator%size_t1      = 2.0_f64*sp_deg1 + num_cells1 + 1
         interpolator%size_t2      = 2.0_f64*sp_deg2 + num_cells2 + 1
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
       
       
       ! ------------------------------------------------------------
       ! reorganization of spline coefficients 1D in coefficients 2D 
       ! -----------------------------------------------------------
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
             
          end do
       end do
         
       interpolator%coeff_splines(:,1) = 0.0_8
       interpolator%coeff_splines(:,nb_spline_eta2+2) = 0.0_8
       ! ------------------------------------------------------------
       
      case(585) ! 4. dirichlet in all sides
         interpolator%size_coeffs1=  num_cells1 + sp_deg1
         interpolator%size_coeffs2=  num_cells2 + sp_deg2
         interpolator%size_t1 = 2.0_f64*sp_deg1 + num_cells1 + 1
         interpolator%size_t2 = 2.0_f64*sp_deg2 + num_cells2 + 1
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
         
         interpolator%coeff_splines(:,:) = 0.0_8
         ! allocation coefficient spline
         do i = 1,nb_spline_eta1
            do j = 1,nb_spline_eta2
               
               interpolator%coeff_splines(i+1,j+1) = &
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
       point_location_eta1(1:sz1) = eta1_coords(1:sz1)
       point_location_eta2(1:sz2) = eta2_coords(1:sz2)

    else ! size depends on BC combination, filled out at initialization.

       select case (interpolator%bc_selector)
       case (0) ! 1. periodic-periodic
          sz1 = interpolator%num_pts1-1
          sz2 = interpolator%num_pts2-1
          
       case (9) ! 2. dirichlet-left, dirichlet-right, periodic
          sz1 = interpolator%num_pts1
          sz2 = interpolator%num_pts2-1
          
       case (576) ! 3. periodic, dirichlet-bottom, dirichlet-top
          sz1 = interpolator%num_pts1-1
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
    
    select case (interpolator%bc_selector)
    case (0) ! periodic-periodic
       interpolator%size_coeffs1 = sz1+1
       interpolator%size_coeffs2 = sz2+1
       interpolator%size_t1 = order1 + sz1 + 1
       interpolator%size_t2 = order2 + sz2 + 1 

       
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       
       call spli2d_perper( &
            period1, sz1+1, order1, point_location_eta1, &
            period2, sz2+1, order2, point_location_eta2, &
            data_array(1:sz1,1:sz2), interpolator%coeff_splines(1:sz1+1,1:sz2+1),&
            interpolator%t1(1:order1 + sz1 + 1), &
            interpolator%t2(1:order2 + sz2 + 1) )
       
       
    case (9) ! 2. dirichlet-left, dirichlet-right, periodic
       interpolator%size_coeffs1 = sz1
       interpolator%size_coeffs2 = sz2+1
       interpolator%size_t1 = order1 + sz1
       interpolator%size_t2 = order2 + sz2 + 1
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       call spli2d_dirper( sz1, order1, point_location_eta1, &
            period2, sz2+1, order2, point_location_eta2, &
            data_array(1:sz1,1:sz2), interpolator%coeff_splines(1:sz1,1:sz2+1),&
            interpolator%t1(1:sz1+order1), &
            interpolator%t2(1:sz2+order2+1) )
  
       
    case(576) !  3. periodic, dirichlet-bottom, dirichlet-top
       interpolator%size_coeffs1 = sz1+1
       interpolator%size_coeffs2 = sz2
       interpolator%size_t1 = order1 + sz1 + 1
       interpolator%size_t2 = order2 + sz2 
       !  data_array must have the same dimension than 
       !  size(  point_location_eta1 ) x  size(  point_location_eta2 )
       !  i.e  data_array must have the dimension sz1 x sz2
       call spli2d_perdir( period1, sz1+1, order1, point_location_eta1, &
            sz2, order2, point_location_eta2, &
            data_array(1:sz1,1:sz2), interpolator%coeff_splines(1:sz1+1,1:sz2),&
            interpolator%t1(1:sz1+order1+1), &
            interpolator%t2(1:sz2+order2) )
       
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

    end select
    interpolator%coefficients_set = .true.
    SLL_DEALLOCATE(point_location_eta2,ierr)
    SLL_DEALLOCATE(point_location_eta1,ierr)
  end subroutine !compute_interpolants_ad2d

  function coefficients_are_set_ad2d( interpolator ) result(res)
    class(arb_deg_2d_interpolator), intent(in)  :: interpolator
    logical :: res
    res = interpolator%coefficients_set
  end function coefficients_are_set_ad2d


  !  ----------------------------------------------------------
  !  Interpolation on the points eta1 and eta2 
  !  ---------------------------------------------------------
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
       if ( res1 > interpolator%eta1_max ) then 
          res1 = res1 -(interpolator%eta1_max-interpolator%eta1_min)
       end if
       if ( res2 > interpolator%eta2_max ) then 
          res2 = res2 -(interpolator%eta2_max-interpolator%eta2_min)
       end if

       if ( res1 < interpolator%eta1_min ) then 
          res1 = res1 + (interpolator%eta1_max-interpolator%eta1_min)
       end if
       if ( res2 < interpolator%eta2_min ) then 
          res2 = res2 +(interpolator%eta2_max-interpolator%eta2_min)
       end if
       

    case (9) ! 2. dirichlet-left, dirichlet-right, periodic
       if ( res2 > interpolator%eta2_max ) then 
          res2 = res2 - (interpolator%eta2_max-interpolator%eta2_min)
       end if
       if ( res2 < interpolator%eta2_min ) then 
          res2 = res2 +(interpolator%eta2_max-interpolator%eta2_min)
       end if
       SLL_ASSERT( res1 >= interpolator%eta1_min )
       SLL_ASSERT( res1 <= interpolator%eta1_max )
  
    case(576) !  3. periodic, dirichlet-bottom, dirichlet-top
       if ( res1 > interpolator%eta1_max ) then 
          res1 = res1 -(interpolator%eta1_max-interpolator%eta1_min)
       end if
       if ( res1 < interpolator%eta1_min ) then 
          res1 = res1 + (interpolator%eta1_max-interpolator%eta1_min)
       end if
       SLL_ASSERT( res2 >= interpolator%eta2_min )
       SLL_ASSERT( res2 <= interpolator%eta2_max )
       
    case (585) ! dirichlet-dirichlet 
       SLL_ASSERT( res1 >= interpolator%eta1_min )
       SLL_ASSERT( res1 <= interpolator%eta1_max )
       SLL_ASSERT( res2 >= interpolator%eta2_min )
       SLL_ASSERT( res2 <= interpolator%eta2_max )
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
       if ( res1 > interpolator%eta1_max ) then 
          res1 = res1 -(interpolator%eta1_max-interpolator%eta1_min)
       end if
       if ( res2 > interpolator%eta2_max ) then 
          res2 = res2 -(interpolator%eta2_max-interpolator%eta2_min)
       end if
       
       if ( res1 < interpolator%eta1_min ) then 
          res1 = res1 + (interpolator%eta1_max-interpolator%eta1_min)
       end if
       if ( res2 < interpolator%eta2_min ) then 
          res2 = res2 +(interpolator%eta2_max-interpolator%eta2_min)
       end if

    case (9) ! 2. dirichlet-left, dirichlet-right, periodic
       if ( res2 > interpolator%eta2_max ) then 
          res2 = res2 - (interpolator%eta2_max-interpolator%eta2_min)
       end if
       if ( res2 < interpolator%eta2_min ) then 
          res2 = res2 +(interpolator%eta2_max-interpolator%eta2_min)
       end if
       SLL_ASSERT( res1 >= interpolator%eta1_min )
       SLL_ASSERT( res1 <= interpolator%eta1_max )
       
    case(576) !  3. periodic, dirichlet-bottom, dirichlet-top
       if ( res1 > interpolator%eta1_max ) then 
          res1 = res1 -(interpolator%eta1_max-interpolator%eta1_min)
       end if
       if ( res1 < interpolator%eta1_min ) then 
          res1 = res1 + (interpolator%eta1_max-interpolator%eta1_min)
       end if
       SLL_ASSERT( res2 >= interpolator%eta2_min )
       SLL_ASSERT( res2 <= interpolator%eta2_max )
       
    case (585) ! dirichlet-dirichlet 
       SLL_ASSERT( res1 >= interpolator%eta1_min )
       SLL_ASSERT( res1 <= interpolator%eta1_max )
       SLL_ASSERT( res2 >= interpolator%eta2_min )
       SLL_ASSERT( res2 <= interpolator%eta2_max )
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
       if ( res1 > interpolator%eta1_max ) then 
          res1 = res1 -(interpolator%eta1_max-interpolator%eta1_min)
       end if
       if ( res2 > interpolator%eta2_max ) then 
          res2 = res2 -(interpolator%eta2_max-interpolator%eta2_min)
       end if
       if ( res1 < interpolator%eta1_min ) then 
          res1 = res1 + (interpolator%eta1_max-interpolator%eta1_min)
       end if
       if ( res2 < interpolator%eta2_min ) then 
          res2 = res2 +(interpolator%eta2_max-interpolator%eta2_min)
       end if
    case (9) ! 2. dirichlet-left, dirichlet-right, periodic
       if ( res2 > interpolator%eta2_max ) then 
          res2 = res2 - (interpolator%eta2_max-interpolator%eta2_min)
       end if
       if ( res2 < interpolator%eta2_min ) then 
          res2 = res2 +(interpolator%eta2_max-interpolator%eta2_min)
       end if
       SLL_ASSERT( res1 >= interpolator%eta1_min )
       SLL_ASSERT( res1 <= interpolator%eta1_max )
       
    case(576) !  3. periodic, dirichlet-bottom, dirichlet-top
       if ( res1 > interpolator%eta1_max ) then 
          res1 = res1 -(interpolator%eta1_max-interpolator%eta1_min)
       end if
       
       if ( res1 < interpolator%eta1_min ) then 
          res1 = res1 + (interpolator%eta1_max-interpolator%eta1_min)
       end if
       SLL_ASSERT( res2 >= interpolator%eta2_min )
       SLL_ASSERT( res2 <= interpolator%eta2_max )
       
    case (585) ! dirichlet-dirichlet 
       SLL_ASSERT( res1 >= interpolator%eta1_min )
       SLL_ASSERT( res1 <= interpolator%eta1_max )
       SLL_ASSERT( res2 >= interpolator%eta2_min )
       SLL_ASSERT( res2 <= interpolator%eta2_max )
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
    
  end function interpolate_derivative2_ad2d !interpolate_derivative2_ad2d

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
