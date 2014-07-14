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

module sll_arbitrary_degree_spline_interpolator_1d_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h" 
  use sll_module_interpolators_1d_base
  use sll_module_deboor_splines_1d
  

  implicit none
  
  ! in what follows, the direction '1' is in the contiguous memory direction.
  type, extends(sll_interpolator_1d_base) :: sll_arb_deg_1d_interpolator           
     sll_int32  :: num_pts
     sll_real64 :: eta_min
     sll_real64 :: eta_max
     sll_int32  :: bc_left
     sll_int32  :: bc_right
     sll_int32  :: spline_degree
     sll_real64, dimension(:), pointer :: knots
     ! some knot-like arrays needed by the spli1d_per routine
     sll_real64, dimension(:), pointer :: t
     sll_int32  :: size_t 
     sll_int64  :: bc_selector ! this is set in initialization
     sll_real64, dimension(:), pointer :: coeff_splines
     sll_int32  :: size_coeffs
     sll_real64 :: slope_left
     sll_real64 :: slope_right
     sll_real64 :: value_left
     sll_real64 :: value_right
     logical    :: compute_slope_left = .TRUE.
     logical    :: compute_slope_right= .TRUE.
     logical    :: compute_value_left = .TRUE.
     logical    :: compute_value_right= .TRUE.

   contains
    procedure, pass(interpolator) :: initialize=>initialize_ad1d_interpolator
    procedure, pass :: set_coefficients => set_coefficients_ad1d 
! better: pre-compute-interpolation-information or something...
    procedure :: compute_interpolants => compute_interpolants_ad1d
   ! procedure,  pass(interpolator) :: compute_spline_coefficients => &
    !     compute_spline_coefficients_ad2d
    !procedure, pass:: compute_spline_coefficients =>compute_spline_coefficients_ad2d
    procedure :: interpolate_value => interpolate_value_ad1d
    procedure :: interpolate_array_values => interpolate_values_ad1d
    procedure :: interpolate_pointer_values => interpolate_pointer_values_ad1d
    procedure :: interpolate_derivative_eta1 => interpolate_derivative_ad1d
    procedure :: interpolate_array_derivatives => interpolate_derivatives_ad1d
    procedure :: interpolate_pointer_derivatives =>interpolate_pointer_derivatives_ad1d
    procedure, pass:: interpolate_array => interpolate_array_ad1d
    procedure, pass:: interpolate_array_disp => interpolate_1d_array_disp_ad1d 
    procedure, pass:: get_coefficients => get_coefficients_ad1d 
    procedure, pass:: reconstruct_array
 end type sll_arb_deg_1d_interpolator

  interface delete
     module procedure delete_arbitrary_degree_1d_interpolator
  end interface delete


contains

  !> @brief delete interpolator arbitrary degree splines.
  !> @details   
  !> 
  !> The parameters are
  !> @param interpolator the type sll_arb_deg_1d_interpolator

  subroutine delete_arbitrary_degree_1d_interpolator( interpolator )
    class(sll_arb_deg_1d_interpolator), intent(inout) :: interpolator
    sll_int32 :: ierr
    SLL_DEALLOCATE(interpolator%knots,ierr)
    SLL_DEALLOCATE(interpolator%t,ierr)
    SLL_DEALLOCATE(interpolator%coeff_splines,ierr)
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

  function new_arbitrary_degree_1d_interpolator(&
       num_pts, &
       eta_min, &
       eta_max, &
       bc_left, &
       bc_right, &
       spline_degree) result(interpolator)
    
    class(sll_arb_deg_1d_interpolator),pointer :: interpolator
    sll_int32, intent(in) :: num_pts
    sll_real64, intent(in) :: eta_min
    sll_real64, intent(in) :: eta_max
    sll_int32, intent(in) :: bc_left
    sll_int32, intent(in) :: bc_right
    sll_int32, intent(in) :: spline_degree
    sll_int32 :: ierr


    SLL_ALLOCATE(interpolator,ierr)
    
    call initialize_ad1d_interpolator( &
         interpolator, &
         num_pts, &
         eta_min, &
         eta_max, &
         bc_left, &
         bc_right, &
         spline_degree)
    
  end function new_arbitrary_degree_1d_interpolator

  !> @brief Initialization of interpolator arbitrary degree splines 1d.
  !> @details To have the interpolator arbitrary degree splines 1d
  !> 
  !> The parameters are
  !> @params interpolator the type sll_arb_deg_1d_interpolator 
  !> @param[in] num_pts the number of points
  !> @param[in] eta_min the minimun
  !> @param[in] eta_max the maximun
  !> @param[in] bc_left  the boundary condition at left
  !> @param[in] bc_right the boundary condition at right
  !> @param[in] spline_degree the degree of B-spline
  !> @return the type arb_deg_1d_interpolator

  subroutine initialize_ad1d_interpolator( &
       interpolator, &
       num_pts, &
       eta_min, &
       eta_max, &
       bc_left, &
       bc_right, &
       spline_degree)

    class(sll_arb_deg_1d_interpolator), intent(inout) :: interpolator
    sll_int32, intent(in) :: num_pts
    sll_real64, intent(in) :: eta_min
    sll_real64, intent(in) :: eta_max
    sll_int32, intent(in) :: bc_left
    sll_int32, intent(in) :: bc_right
    sll_int32, intent(in) :: spline_degree
    sll_int32 :: ierr
    sll_int32 :: tmp
    sll_int64 :: bc_selector
    
    ! do some argument checking...
    if(((bc_left  == SLL_PERIODIC).and.(bc_right.ne. SLL_PERIODIC))) then
       print *, 'initialize_arbitrary_degree_1d_interpolator, ERROR: ', &
            'if one boundary condition is specified as periodic, then ', &
            'both must be. Error in first direction.'
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
    
    interpolator%spline_degree = spline_degree
    interpolator%eta_min       = eta_min
    interpolator%eta_max       = eta_max
    interpolator%bc_left       = bc_left
    interpolator%bc_right      = bc_right
    interpolator%bc_selector   = bc_selector
    interpolator%num_pts       = num_pts
    
    select case (bc_selector)
    case (0) ! 1. periodic
       SLL_ALLOCATE( interpolator%knots(2*spline_degree+2),ierr )
       tmp = num_pts*num_pts
       SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)
       
    case (9) ! 2. dirichlet-left, dirichlet-right
       SLL_ALLOCATE( interpolator%knots(num_pts+2*spline_degree),ierr )
       tmp = num_pts*num_pts
       SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)
       interpolator%value_left = 0.0_f64
       interpolator%value_right = 0.0_f64
    case(10) ! Neumann - Dirichlet
       SLL_ALLOCATE( interpolator%knots(num_pts+2*spline_degree),ierr )
       tmp = num_pts*num_pts
       SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)
       interpolator%slope_left  = 0.0_f64
       interpolator%value_right = 0.0_f64

    case(12) ! Hermite- Dirichlet
       SLL_ALLOCATE( interpolator%knots(num_pts+2*spline_degree),ierr )
       tmp = num_pts*num_pts
       SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)
       interpolator%slope_left = 0.0_f64
       interpolator%value_right = 0.0_f64

    case(17) ! Dirichlet-Neumann
       SLL_ALLOCATE( interpolator%knots(num_pts+2*spline_degree),ierr )
       tmp = num_pts*num_pts
       SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)
       interpolator%slope_right = 0.0_f64
       interpolator%value_left = 0.0_f64

    case(18) ! Neumann - Neumann 
       SLL_ALLOCATE( interpolator%knots(num_pts+2*spline_degree),ierr )
       tmp = num_pts*num_pts
       SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)
       interpolator%slope_right = 0.0_f64
       interpolator%slope_left = 0.0_f64

    case(20) ! Hermite - Neumann
       SLL_ALLOCATE( interpolator%knots(num_pts+2*spline_degree),ierr )
       tmp = num_pts*num_pts
       SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)
       interpolator%slope_right = 0.0_f64
       interpolator%slope_left = 0.0_f64

    case(33) ! Dirichlet - Hermite
       SLL_ALLOCATE( interpolator%knots(num_pts+2*spline_degree),ierr )
       tmp = num_pts*num_pts
       SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)
       interpolator%slope_right = 0.0_f64
       interpolator%value_left = 0.0_f64

    case(34) ! Neumann- Hermite
       SLL_ALLOCATE( interpolator%knots(num_pts+2*spline_degree),ierr )
       tmp = num_pts*num_pts
       SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)
       interpolator%slope_right = 0.0_f64
       interpolator%slope_left = 0.0_f64

    case(36) ! Hermitte - Hermite
       SLL_ALLOCATE( interpolator%knots(num_pts+2*spline_degree),ierr )
       tmp = num_pts*num_pts
       SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)
       interpolator%slope_right = 0.0_f64
       interpolator%slope_left = 0.0_f64
       
    case default
       print *, 'initialize_ad1d_interpolator: BC combination not implemented.'
    end select
    
    interpolator%coeff_splines(:) = 0.0_f64
    SLL_ALLOCATE( interpolator%t(num_pts*num_pts),ierr)
    interpolator%t(:) = 0.0_f64
  end subroutine initialize_ad1d_interpolator
  
  !> @brief initializing the coefficients of splines.
  !> @details  initializing the coefficients of splines
  !>  fot the arbitrary degree splines interpolator 1d 
  !> The parameters are
  !> @param interpolator the type sll_arb_deg_1d_interpolator
  !> @param[in],optional, coeffs the 1d arrays corresponding of the splines coefficients
  !> @return the type arb_deg_1d_interpolator

  subroutine set_coefficients_ad1d( &
   interpolator, &
   coeffs)

   class(sll_arb_deg_1d_interpolator), intent(inout)  :: interpolator
   sll_real64, dimension(:), intent(in), optional :: coeffs
   sll_int32 :: sp_deg
   sll_int32 :: num_cells
   sll_int32 :: tmp
   sll_int32 :: i
   sll_real64 :: eta_min, eta_max
   sll_real64 :: delta
   sll_int32  ::  nb_spline_eta
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
      if ( size(coeffs) .ne.  num_cells + 1 ) then
         print*, 'problem in set_coeff_1d_arb_deg_spline '
         print*, 'size coeffs must be equal to ',num_cells + 1
         print*, 'and not =', size(coeffs)
         stop
      endif
      ! allocation and definition of knots
      do i = -sp_deg, num_cells + sp_deg
         interpolator%t( i + sp_deg + 1 ) = eta_min + i*delta
      end do
      

      do i = 1,num_cells
         interpolator%coeff_splines(i) = coeffs( i )
      end do

      do i = 1, sp_deg
         
         interpolator%coeff_splines(num_cells + i ) = coeffs(i )
      end do

      do i= 1,sp_deg
         
         interpolator%coeff_splines(num_cells +  i ) = &
              interpolator%coeff_splines(sp_deg-(i-1))
      end do
      
   case (9) ! 2. dirichlet-left, dirichlet-right
      interpolator%size_coeffs=  num_cells + sp_deg
      interpolator%size_t = 2*sp_deg + num_cells + 1
      nb_spline_eta = num_cells + sp_deg - 2
      if ( size(coeffs) .ne.  nb_spline_eta ) then
         print*, 'problem in set_coeff_1d_arb_deg_spline '
         print*, 'size coeffs must be equal to ',nb_spline_eta
         print*, 'and not =', size(coeffs)
         stop
      endif
      ! allocation and definition of knots
      
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

      if ( size(coeffs) .ne.  nb_spline_eta ) then
         print*, 'problem in set_coeff_1d_arb_deg_spline '
         print*, 'size coeffs must be equal to ',nb_spline_eta
         print*, 'and not =', size(coeffs)
         stop
      endif

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
      if ( size(coeffs) .ne.  nb_spline_eta ) then
         print*, 'problem in set_coeff_1d_arb_deg_spline '
         print*, 'size coeffs must be equal to ',nb_spline_eta
         print*, 'and not =', size(coeffs)
         stop
      endif
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
      if ( size(coeffs) .ne.  nb_spline_eta ) then
         print*, 'problem in set_coeff_1d_arb_deg_spline '
         print*, 'size coeffs must be equal to ',nb_spline_eta
         print*, 'and not =', size(coeffs)
         stop
      endif

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
   case(18) ! Neumann - Neumann 
      interpolator%size_coeffs= num_cells + sp_deg
      interpolator%size_t     = 2*sp_deg + num_cells + 1
      nb_spline_eta = num_cells + sp_deg 
      if ( size(coeffs) .ne.  nb_spline_eta ) then
         print*, 'problem in set_coeff_1d_arb_deg_spline '
         print*, 'size coeffs must be equal to ',nb_spline_eta
         print*, 'and not =', size(coeffs)
         stop
      endif
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
         interpolator%coeff_splines(i ) =  coeffs(i)
      end do
   case(20) ! Hermite - Neumann
      interpolator%size_coeffs= num_cells + sp_deg
      interpolator%size_t     = 2*sp_deg + num_cells + 1
      nb_spline_eta = num_cells + sp_deg 
      if ( size(coeffs) .ne.  nb_spline_eta ) then
         print*, 'problem in set_coeff_1d_arb_deg_spline '
         print*, 'size coeffs must be equal to ',nb_spline_eta
         print*, 'and not =', size(coeffs)
         stop
      endif
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
         interpolator%coeff_splines(i) =  coeffs(i)
      end do
   case(33) ! Dirichlet - Hermite
      interpolator%size_coeffs= num_cells + sp_deg
      interpolator%size_t     = 2*sp_deg + num_cells + 1
      nb_spline_eta = num_cells + sp_deg - 1
      if ( size(coeffs) .ne.  nb_spline_eta ) then
         print*, 'problem in set_coeff_1d_arb_deg_spline '
         print*, 'size coeffs must be equal to ',nb_spline_eta
         print*, 'and not =', size(coeffs)
         stop
      endif
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
   case(34) ! Neumann- Hermite
      interpolator%size_coeffs= num_cells + sp_deg
      interpolator%size_t     = 2*sp_deg + num_cells + 1
      nb_spline_eta = num_cells + sp_deg 
      if ( size(coeffs) .ne.  nb_spline_eta ) then
         print*, 'problem in set_coeff_1d_arb_deg_spline '
         print*, 'size coeffs must be equal to ',nb_spline_eta
         print*, 'and not =', size(coeffs)
         stop
      endif
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
         interpolator%coeff_splines(i) =  coeffs(i)
      end do
   case(36)! Hermite - Hermite

      interpolator%size_coeffs= num_cells + sp_deg
      interpolator%size_t     = 2*sp_deg + num_cells + 1
      nb_spline_eta = num_cells + sp_deg
      if ( size(coeffs) .ne.  nb_spline_eta ) then
         print*, 'problem in set_coeff_1d_arb_deg_spline '
         print*, 'size coeffs must be equal to ',nb_spline_eta
         print*, 'and not =', size(coeffs)
         stop
      endif
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
         interpolator%coeff_splines(i ) =  coeffs(i)
      end do
   case default
      print *, 'arbitrary_degree_spline_1d() error: set_spline_coefficients ',&
           'not recognized.'
      stop
   end select
 end subroutine set_coefficients_ad1d

 !> @brief Initialization of the boundary for interpolator arbitrary degree splines 1d.
 !> @details Initialization of the boundary
 !>  interpolator sll_arbitrary degree splines 1d
 !> The parameters are
 !> @params interpolator the type sll_arb_deg_1d_interpolator
 !> @param[in],optional,value_left a array contains the value in the left
 !> @param[in],optional,value_right a array contains the value in the right
 !> @param[in],optional,slope_left a array contains the value in the left for derivative
 !> @param[in],optional,slope_right, array contains the value in the right for derivative
 !> @return the type sll_arb_deg_1d_interpolator
 subroutine set_values_at_boundary1d(&
      interpolator,&
      value_left,&
      value_right,&
      slope_left,&
      slope_right)
   
   class(sll_arb_deg_1d_interpolator), intent(inout) :: interpolator
   sll_real64, intent(in), optional     :: value_left
   sll_real64, intent(in), optional     :: value_right
   sll_real64, intent(in), optional     :: slope_left
   sll_real64, intent(in), optional     :: slope_right
   sll_int32 :: bc_left
   sll_int32 :: bc_right
   sll_int64 :: bc_selector
   
   bc_left = interpolator%bc_left
   bc_right= interpolator%bc_right
   bc_selector = interpolator%bc_selector
   
   ! modify this such as slope
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
  !> @details computing the coefficients spline with a given 
  !>  data_array 1D coorespondind at the values of a function 
  !>  on eta_coords of size size_eta_coords
  !>  if the eta_coords and eta_coords is not given 
  !>  we consider that the values of the function is on the points in the mesh_1d
  !> 
  !> The parameters are
  !> @param interpolator the type sll_arb_deg_1d_interpolator
  !> @param[in] data_array the 1d arrays corresponding at the values of a function
  !> @param[in],optional, eta_coords the 1d arrays  
  !> @param[in],optional, size_eta_coords the size of eta_coords
  !> @return the type sll_arb_deg_1d_interpolator
 subroutine compute_interpolants_ad1d( &
      interpolator,data_array, &
      eta_coords, &
      size_eta_coords)
   
   class(sll_arb_deg_1d_interpolator), intent(inout)  :: interpolator
   sll_real64, dimension(:), intent(in)           :: data_array
   sll_real64, dimension(:),pointer               :: data_array_derivative
   sll_real64, dimension(:), intent(in),optional  :: eta_coords
   sll_int32, intent(in),optional                 :: size_eta_coords
   sll_real64, dimension(:),pointer               :: point_locate_eta
   sll_int32, dimension(:),pointer               :: point_locate_eta_derivative
   sll_real64 :: delta_eta
   sll_int32  :: sz,sz_deriv
   sll_real64 :: period
   sll_int32  :: order
   sll_int32  :: ierr
   sll_int32  :: i
   logical    :: user_coords
   
   if(present(eta_coords) .and. (.not. present(size_eta_coords))) then
      print *, 'compute_interpolants_ad1d(), ERROR: if eta_coords is ', &
           'passed, its size must be specified as well through ', &
           'size_eta_coords.'
      stop
   end if
   


   if( present(eta_coords) ) then
       user_coords = .true.
    else
       user_coords = .false.
    end if


    if (user_coords .eqv. .true.) then
       sz = size(data_array)!size_eta_coords
       
       SLL_ALLOCATE(point_locate_eta(sz),ierr)
       point_locate_eta = eta_coords

    else ! size depends on BC combination
       sz = interpolator%num_pts

       delta_eta = (interpolator%eta_max - interpolator%eta_min)&
            /(interpolator%num_pts -1)
       SLL_ALLOCATE(point_locate_eta(sz),ierr)
       
       do i = 1,sz
          point_locate_eta(i) = interpolator%eta_min + delta_eta*(i-1)
       end do
    end if

    
    if (interpolator%compute_slope_left .eqv. .true.) then
       interpolator%slope_left = forward_fd_5pt( data_array,point_locate_eta)
    end if
    if (interpolator%compute_slope_right .eqv. .true.) then
       interpolator%slope_right = backward_fd_5pt( data_array,point_locate_eta,sz)
    end if

    if (interpolator%compute_value_left .eqv. .true.) then
       interpolator%value_left = data_array(1)
    end if
    if (interpolator%compute_value_right .eqv. .true.) then
       interpolator%value_right = data_array(sz)
    end if



    SLL_ASSERT(sz .le. interpolator%num_pts* interpolator%num_pts)
    SLL_ASSERT(size(data_array) .le. sz)
    SLL_ASSERT(size(point_locate_eta)  .ge. sz)


    order  = interpolator%spline_degree + 1
    period = interpolator%eta_max - interpolator%eta_min

    select case (interpolator%bc_selector)
    case (0) ! periodic
       interpolator%size_coeffs = sz !+ 1
       interpolator%size_t = order + sz !+ 1 
       call spli1d_per( & 
            period, sz, order, point_locate_eta, &
            data_array, interpolator%coeff_splines(1:sz),&!+1),&
            interpolator%t(1:order + sz ))!+ 1))
       
       
       
    case (9) ! 2. dirichlet-left, dirichlet-right
       interpolator%size_coeffs = sz
       interpolator%size_t = order + sz 

       call spli1d_dir( sz, order, point_locate_eta, &
            data_array, interpolator%coeff_splines(1:sz),&
            interpolator%t(1:sz+order) )
 
       ! test dirichlet non homogene
       interpolator%coeff_splines(1)  = interpolator%value_left
       interpolator%coeff_splines(sz) = interpolator%value_right

    case(10) ! Neumann - Dirichlet
       ! -----------------------------------
       !!! It is only for cubic spline !!!!
       ! -----------------------------------
       sz_deriv = 2
       interpolator%size_coeffs = sz + sz_deriv
       interpolator%size_t = order + sz + sz_deriv 
       
       SLL_ALLOCATE(point_locate_eta_derivative(sz_deriv),ierr)
       SLL_ALLOCATE(data_array_derivative(sz_deriv),ierr)

       point_locate_eta_derivative(1) = 1
       data_array_derivative(1)       = 0.0_f64
       point_locate_eta_derivative(2) = sz
       data_array_derivative(2)       = interpolator%slope_right

       call spli1d_der(sz,sz_deriv,order,&
            point_locate_eta,&
            data_array,&
            point_locate_eta_derivative,&
            data_array_derivative,&
            interpolator%coeff_splines(1:sz+sz_deriv),&
            interpolator%t(1:sz+order+sz_deriv))
       
       ! test dirichlet non homogene
       interpolator%coeff_splines(sz+sz_deriv) = interpolator%value_right

       SLL_DEALLOCATE(point_locate_eta_derivative,ierr)
       SLL_DEALLOCATE(data_array_derivative,ierr)

    case(12) ! Hermite - Dirichlet

       ! -----------------------------------
       !!! It is only for cubic spline !!!!
       ! -----------------------------------
       sz_deriv = 2
       interpolator%size_coeffs = sz + sz_deriv
       interpolator%size_t = order + sz + sz_deriv 
       
       SLL_ALLOCATE(point_locate_eta_derivative(sz_deriv),ierr)
       SLL_ALLOCATE(data_array_derivative(sz_deriv),ierr)

       point_locate_eta_derivative(1) = 1
       data_array_derivative(1)       = interpolator%slope_left
       point_locate_eta_derivative(2) = sz
       data_array_derivative(2)       = interpolator%slope_right

       call spli1d_der(sz,sz_deriv,order,&
            point_locate_eta,&
            data_array,&
            point_locate_eta_derivative,&
            data_array_derivative,&
            interpolator%coeff_splines(1:sz+sz_deriv),&
            interpolator%t(1:sz+order+sz_deriv))
       
       ! test dirichlet non homogene
       interpolator%coeff_splines(sz+sz_deriv) = interpolator%value_right
       SLL_DEALLOCATE(point_locate_eta_derivative,ierr)
       SLL_DEALLOCATE(data_array_derivative,ierr)

    case(17) ! Dirichlet - Neumann

        ! -----------------------------------
       !!! It is only for cubic spline !!!!
       ! -----------------------------------
       sz_deriv = 2
       interpolator%size_coeffs = sz + sz_deriv
       interpolator%size_t = order + sz + sz_deriv 

       SLL_ALLOCATE(point_locate_eta_derivative(sz_deriv),ierr)
       SLL_ALLOCATE(data_array_derivative(sz_deriv),ierr)

       point_locate_eta_derivative(1) = 1
       data_array_derivative(1)       = interpolator%slope_left
       point_locate_eta_derivative(2) = sz
       data_array_derivative(2)       = 0.0_f64
       
       call spli1d_der(sz,sz_deriv,order,&
            point_locate_eta,&
            data_array,&
            point_locate_eta_derivative,&
            data_array_derivative,&
            interpolator%coeff_splines(1:sz+sz_deriv),&
            interpolator%t(1:sz+order+sz_deriv))
       
       ! test dirichlet non homogene
       interpolator%coeff_splines(1) = interpolator%value_left
       SLL_DEALLOCATE(point_locate_eta_derivative,ierr)
       SLL_DEALLOCATE(data_array_derivative,ierr)
    case(18) ! Neumann - Neumann

        ! -----------------------------------
       !!! It is only for cubic spline !!!!
       ! -----------------------------------
       sz_deriv = 2
       interpolator%size_coeffs = sz + sz_deriv
       interpolator%size_t = order + sz + sz_deriv 

       SLL_ALLOCATE(point_locate_eta_derivative(sz_deriv),ierr)
       SLL_ALLOCATE(data_array_derivative(sz_deriv),ierr)

       point_locate_eta_derivative(1) = 1
       data_array_derivative(1)       = 0.0_f64
       point_locate_eta_derivative(2) = sz
       data_array_derivative(2)       = 0.0_f64
       call spli1d_der(sz,sz_deriv,order,&
            point_locate_eta,&
            data_array,&
            point_locate_eta_derivative,&
            data_array_derivative,&
            interpolator%coeff_splines(1:sz+sz_deriv),&
            interpolator%t(1:sz+order+sz_deriv))

       SLL_DEALLOCATE(point_locate_eta_derivative,ierr)
       SLL_DEALLOCATE(data_array_derivative,ierr)
 
    case(20) ! Hermite - Neumann

        ! -----------------------------------
       !!! It is only for cubic spline !!!!
       ! -----------------------------------
       sz_deriv = 2
       interpolator%size_coeffs = sz + sz_deriv
       interpolator%size_t = order + sz + sz_deriv 

       SLL_ALLOCATE(point_locate_eta_derivative(sz_deriv),ierr)
       SLL_ALLOCATE(data_array_derivative(sz_deriv),ierr)

       point_locate_eta_derivative(1) = 1
       data_array_derivative(1)       = interpolator%slope_left
       point_locate_eta_derivative(2) = sz
       data_array_derivative(2)       = 0.0_f64

       call spli1d_der(sz,sz_deriv,order,&
            point_locate_eta,&
            data_array,&
            point_locate_eta_derivative,&
            data_array_derivative,&
            interpolator%coeff_splines(1:sz+sz_deriv),&
            interpolator%t(1:sz+order+sz_deriv))
       SLL_DEALLOCATE(point_locate_eta_derivative,ierr)
       SLL_DEALLOCATE(data_array_derivative,ierr)

    case(33) ! Dirichlet - Hermite

        ! -----------------------------------
       !!! It is only for cubic spline !!!!
       ! -----------------------------------
       sz_deriv = 2
       interpolator%size_coeffs = sz + sz_deriv
       interpolator%size_t = order + sz + sz_deriv 

       SLL_ALLOCATE(point_locate_eta_derivative(sz_deriv),ierr)
       SLL_ALLOCATE(data_array_derivative(sz_deriv),ierr)

       point_locate_eta_derivative(1) = 1
       data_array_derivative(1)       = interpolator%slope_left
       point_locate_eta_derivative(2) = sz
       data_array_derivative(2)       = interpolator%slope_right

       call spli1d_der(sz,sz_deriv,order,&
            point_locate_eta,&
            data_array,&
            point_locate_eta_derivative,&
            data_array_derivative,&
            interpolator%coeff_splines(1:sz+sz_deriv),&
            interpolator%t(1:sz+order+sz_deriv))
       
       ! test dirichlet non homogene
       interpolator%coeff_splines(1) = interpolator%value_left
       SLL_DEALLOCATE(point_locate_eta_derivative,ierr)
       SLL_DEALLOCATE(data_array_derivative,ierr)
    case(34) ! Neumann - Hermite 
       
        ! -----------------------------------
       !!! It is only for cubic spline !!!!
       ! -----------------------------------

       sz_deriv = 2
       interpolator%size_coeffs = sz + sz_deriv
       interpolator%size_t = order + sz + sz_deriv 

       SLL_ALLOCATE(point_locate_eta_derivative(sz_deriv),ierr)
       SLL_ALLOCATE(data_array_derivative(sz_deriv),ierr)

       point_locate_eta_derivative(1) = 1
       data_array_derivative(1)       = 0.0_f64
       point_locate_eta_derivative(2) = sz
       data_array_derivative(2)       = interpolator%slope_right
       
       call spli1d_der(sz,sz_deriv,order,&
            point_locate_eta,&
            data_array,&
            point_locate_eta_derivative,&
            data_array_derivative,&
            interpolator%coeff_splines(1:sz+sz_deriv),&
            interpolator%t(1:sz+order+sz_deriv))

       SLL_DEALLOCATE(point_locate_eta_derivative,ierr)
       SLL_DEALLOCATE(data_array_derivative,ierr)
       
    case(36) ! Hermite - Hermite

        ! -----------------------------------
       !!! It is only for cubic spline !!!!
       ! -----------------------------------
       sz_deriv = 2
       interpolator%size_coeffs = sz + sz_deriv
       interpolator%size_t = order + sz + sz_deriv 

        SLL_ALLOCATE(point_locate_eta_derivative(sz_deriv),ierr)
       SLL_ALLOCATE(data_array_derivative(sz_deriv),ierr)

       point_locate_eta_derivative(1) = 1
       data_array_derivative(1)       = interpolator%slope_left
       point_locate_eta_derivative(2) = sz
       data_array_derivative(2)       = interpolator%slope_right
       
       call spli1d_der(sz,sz_deriv,order,&
            point_locate_eta,&
            data_array,&
            point_locate_eta_derivative,&
            data_array_derivative,&
            interpolator%coeff_splines(1:sz+sz_deriv),&
            interpolator%t(1:sz+order+sz_deriv))

       SLL_DEALLOCATE(point_locate_eta_derivative,ierr)
       SLL_DEALLOCATE(data_array_derivative,ierr)
 
    end select
    SLL_DEALLOCATE(point_locate_eta,ierr)

  end subroutine compute_interpolants_ad1d
  

  !> @brief Interpolation on the points eta using 
  !> the arbitrary degree splines interpolator 1d 
  !> @details computing the values with the interpolator arbitrary degree splines 1d
  !>  on the points eta of arbitrary degree splines 1d
  !> 
  !> The parameters are
  !> @param interpolator the type sll_arb_deg_1d_interpolator
  !> @param[in] eta the point 
  !> @return val the values on the points eta
  function interpolate_value_ad1d( &
       interpolator, &
       eta1) result(val)
    

    class(sll_arb_deg_1d_interpolator), intent(inout)  :: interpolator
    sll_real64, intent(in)         :: eta1
    sll_real64                     :: val
    sll_int32 :: size_coeffs
    !sll_real64 :: bvalue
    sll_real64 :: res
    sll_real64,dimension(:),pointer :: knot_tmp
    sll_real64,dimension(:),pointer :: coef_tmp
    size_coeffs = interpolator%size_coeffs

    res = eta1
    select case (interpolator%bc_selector)
    case (0) ! periodic

       if( res < interpolator%eta_min ) then
          res = res+interpolator%eta_max-interpolator%eta_min
       else if( res >  interpolator%eta_max ) then
          res = res+interpolator%eta_min-interpolator%eta_max
       end if
    case (9) ! 2. dirichlet-left, dirichlet-right
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then 
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then 
          print*, 'problem  x < eta_min'
          stop
       end if
    case(10) ! Neumann - Dirichlet
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then 
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then 
          print*, 'problem  x < eta_min'
          stop
       end if
    case(12) ! Hermite - Dirichlet
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then 
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then 
          print*, 'problem  x < eta_min'
          stop
       end if
    case(17) ! Dirichlet - Neumann
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then 
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then 
          print*, 'problem  x < eta_min'
          stop
       end if
    case(18) ! Neumann - Neumann
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then 
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then 
          print*, 'problem  x < eta_min'
          stop
       end if
    case(20) ! Hermite - Neumann
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then 
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then 
          print*, 'problem  x < eta_min'
          stop
       end if
    case(33) ! Dirichlet - Hermite
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then 
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then 
          print*, 'problem  x < eta_min'
          stop
       end if
    case(34) ! Neumann - Hermite
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then 
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then 
          print*, 'problem  x < eta_min'
          stop
       end if
    case(36) ! Hermite - Hermite
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then 
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then 
          print*, 'problem  x < eta_min'
          stop
       end if
    
    end select
       
    knot_tmp => interpolator%t(1:interpolator%size_t)
    coef_tmp => interpolator%coeff_splines(1:size_coeffs)
    val = bvalue( &
         knot_tmp,&
         coef_tmp,&
         size_coeffs, &
         interpolator%spline_degree+1, &
         res,0)
  end function interpolate_value_ad1d


  !> @brief First derivative interpolation on the point eta 
  !> @details computing the values of the first derivative
  !> with the interpolator arbitrary degree splines 1d
  !> on the points eta of arbitrary degree splines 1d
  !> 
  !> The parameters are
  !> @param interpolator the type arb_deg_1d_interpolator
  !> @param[in] eta the point 
  !> @return val the values on the point eta of the first derivative

  function interpolate_derivative_ad1d( &
    interpolator, &
    eta1 ) result(val)

    class(sll_arb_deg_1d_interpolator), intent(inout)  :: interpolator
    sll_real64, intent(in)         :: eta1
    sll_real64                     :: val
    sll_int32 :: size_coeffs
   ! sll_real64 :: dvalue1d
    sll_real64 :: res
    sll_real64, dimension(:),pointer :: knot_tmp 
    sll_real64, dimension(:),pointer :: coef_tmp

    SLL_ASSERT( eta1 .ge. interpolator%eta_min )
    SLL_ASSERT( eta1 .le. interpolator%eta_max )

    size_coeffs = interpolator%size_coeffs

    res = eta1
    
    select case (interpolator%bc_selector)
    case (0) ! periodic
!
       if( res < interpolator%eta_min ) then
          res = res+interpolator%eta_max-interpolator%eta_min
       else if( res >  interpolator%eta_max ) then
          res = res+interpolator%eta_min-interpolator%eta_max
       end if
    case (9) ! 2. dirichlet-left, dirichlet-right
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then 
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then 
          print*, 'problem  x < eta_min'
          stop
       end if
    case(10) ! Neumann - Dirichlet
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then 
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then 
          print*, 'problem  x < eta_min'
          stop
       end if
    case(12) ! Hermite - Dirichlet
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then 
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then 
          print*, 'problem  x < eta_min'
          stop
       end if
    case(17) ! Dirichlet - Neumann
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then 
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then 
          print*, 'problem  x < eta_min'
          stop
       end if
    case(18) ! Neumann - Neumann
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then 
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then 
          print*, 'problem  x < eta_min'
          stop
       end if
    case(20) ! Hermite - Neumann
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then 
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then 
          print*, 'problem  x < eta_min'
          stop
       end if
    case(33) ! Dirichlet - Hermite
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then 
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then 
          print*, 'problem  x < eta_min'
          stop
       end if
    case(34) ! Neumann - Hermite
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then 
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then 
          print*, 'problem  x < eta_min'
          stop
       end if
    case(36) ! Hermite - Hermite
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then 
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then 
          print*, 'problem  x < eta_min'
          stop
       end if
       
    end select
    
    knot_tmp =>  interpolator%t(1:interpolator%size_t)
    coef_tmp => interpolator%coeff_splines(1:size_coeffs)
    val = dvalue1d( &
         res, &
         size_coeffs, &
         interpolator%spline_degree+1, &
         coef_tmp, &
         knot_tmp,&
         1)
  end function interpolate_derivative_ad1d
  
  function interpolate_array_ad1d( &
       this, num_points, data, coordinates)&
    result(data_out)
    class(sll_arb_deg_1d_interpolator),  intent(in)       :: this
    !class(sll_cubic_spline_1D),  intent(in)      :: this
    sll_int32,  intent(in)                 :: num_points
    sll_real64, dimension(:), intent(in)   :: coordinates
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_points)      :: data_out
    
    print *, 'interpolate_array_ad1d: not implemented'
    data_out = -1000000._f64*data*coordinates*this%spline_degree
  end function interpolate_array_ad1d
  
  function interpolate_1d_array_disp_ad1d( &
       this,        &
       num_points, &
       data,     &
       alpha) result(res)
      
    class(sll_arb_deg_1d_interpolator), intent(in)    :: this
    sll_int32, intent(in)                          :: num_points 
    sll_real64, dimension(:), intent(in)         :: data
    sll_real64, intent(in)         :: alpha  
    sll_real64, dimension(num_points) :: res
    
    print *, 'interpolate_1d_array_disp_ad1d: not implemented.'
    res = -1000000._f64*alpha*data*this%spline_degree
  end function interpolate_1d_array_disp_ad1d
    
    
  function get_coefficients_ad1d(interpolator)
    class(sll_arb_deg_1d_interpolator), intent(in)    :: interpolator
    sll_real64, dimension(:), pointer            :: get_coefficients_ad1d     
    
    get_coefficients_ad1d => interpolator%coeff_splines
  end function get_coefficients_ad1d

  subroutine interpolate_values_ad1d( &
       interpolator, &
       num_pts, &
       vals_to_interpolate, &
       output_array )
    
    class(sll_arb_deg_1d_interpolator),  intent(in) :: interpolator
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(:), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(:), intent(out)  :: output_array
    
    print*, 'interpolate_values_ad1d NOT iMPLEMENTED YET'
    output_array = -1000000._f64*num_pts&
         *vals_to_interpolate*interpolator%spline_degree
  end subroutine interpolate_values_ad1d

  subroutine interpolate_pointer_values_ad1d( &
       interpolator, &
       num_pts, &
       vals_to_interpolate, &
       output )
    
    class(sll_arb_deg_1d_interpolator),  intent(in) :: interpolator
    sll_int32,  intent(in)            :: num_pts
    sll_real64, dimension(:), pointer :: vals_to_interpolate
    sll_real64, dimension(:), pointer :: output
    
    print*, 'interpolate_pointer_values_ad1d NOT iMPLEMENTED YET'
    output = -1000000._f64*num_pts&
         *vals_to_interpolate*interpolator%spline_degree
  end subroutine interpolate_pointer_values_ad1d
  
  subroutine interpolate_derivatives_ad1d( &
       interpolator, &
       num_pts, &
       vals_to_interpolate, &
       output_array )
    
    class(sll_arb_deg_1d_interpolator),  intent(in) :: interpolator
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(:), intent(in) :: vals_to_interpolate
    sll_real64, dimension(:), intent(out) :: output_array
    
    print*, 'interpolate_derivatives_ad1d NOT iMPLEMENTED YET'
    output_array = -1000000._f64*num_pts &
         *vals_to_interpolate*interpolator%spline_degree
  end subroutine interpolate_derivatives_ad1d
  
  subroutine interpolate_pointer_derivatives_ad1d( &
       interpolator, &
       num_pts, &
       vals_to_interpolate, &
       output )
#ifdef STDF95
    type(sll_arb_deg_1d_interpolator),  intent(in) :: interpolator
#else
    class(sll_arb_deg_1d_interpolator),  intent(in) :: interpolator
#endif
    sll_int32,  intent(in)              :: num_pts
    sll_real64, dimension(:), pointer   :: vals_to_interpolate
    sll_real64, dimension(:), pointer   :: output
    output = -1000000.0_f64*num_pts&
         *vals_to_interpolate*interpolator%spline_degree
    print*, 'interpolate_pointer_derivatives_ad1d NOT iMPLEMENTED YET'
  end subroutine interpolate_pointer_derivatives_ad1d
  
  function reconstruct_array(this, num_points, data) result(res)
    ! dummy procedure

    class(sll_arb_deg_1d_interpolator),  intent(in) :: this
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
  function forward_fd_5pt( data,eta) result(res)
    sll_real64, dimension(:), intent(in) :: data
    sll_real64, dimension(:), intent(in) :: eta
    sll_real64 :: res
    
    res = (-(25.0_f64/12.0_f64)*data(1)*(eta(2) - eta(1)) &
                      + 4.0_f64*data(2)*(eta(3) - eta(2)) &
                      - 3.0_f64*data(3)*(eta(4) - eta(3)) &
            + (4.0_f64/3.0_f64)*data(4)*(eta(5) - eta(4)) &
                      -0.25_f64*data(5)*(eta(6) - eta(5)))
  end function forward_fd_5pt

  function backward_fd_5pt( data,eta,li)result(res)
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
end module sll_arbitrary_degree_spline_interpolator_1d_module
