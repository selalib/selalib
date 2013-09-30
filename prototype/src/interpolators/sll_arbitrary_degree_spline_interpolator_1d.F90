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
  implicit none
  
  ! in what follows, the direction '1' is in the contiguous memory direction.
  type, extends(sll_interpolator_1d_base) :: arb_deg_1d_interpolator           
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
     sll_int32                  :: size_coeffs

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
    procedure :: interpolate_pointer_derivatives => interpolate_pointer_derivatives_ad1d
    procedure, pass:: interpolate_array => interpolate_array_ad1d
    procedure, pass:: interpolate_array_disp => interpolate_1d_array_disp_ad1d 
    procedure, pass:: get_coefficients => get_coefficients_ad1d 
    procedure, pass:: reconstruct_array
 end type arb_deg_1d_interpolator

  interface delete
     module procedure delete_arbitrary_degree_1d_interpolator
  end interface delete

contains

  subroutine delete_arbitrary_degree_1d_interpolator( interpolator )
    class(arb_deg_1d_interpolator), intent(inout) :: interpolator
    sll_int32 :: ierr
    SLL_DEALLOCATE(interpolator%knots,ierr)
    SLL_DEALLOCATE(interpolator%t,ierr)
    SLL_DEALLOCATE(interpolator%coeff_splines,ierr)
  end subroutine delete_arbitrary_degree_1d_interpolator

  subroutine initialize_ad1d_interpolator( &
    interpolator, &
    num_pts, &
    eta_min, &
    eta_max, &
    bc_left, &
    bc_right, &
    spline_degree )

    class(arb_deg_1d_interpolator), intent(inout) :: interpolator
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
    interpolator%eta_min = eta_min
    interpolator%eta_max = eta_max
    interpolator%bc_left  = bc_left
    interpolator%bc_right = bc_right
    interpolator%bc_selector = bc_selector
    interpolator%num_pts = num_pts

    select case (bc_selector)
    case (0) ! 1. periodic
       SLL_ALLOCATE( interpolator%knots(2*spline_degree+2),ierr )
       tmp = num_pts*num_pts
       SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)

    case (9) ! 2. dirichlet-left, dirichlet-right
       SLL_ALLOCATE( interpolator%knots(num_pts+2*spline_degree),ierr )
       tmp = num_pts*num_pts
       SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)

       
    case default
       print *, 'initialize_ad1d_interpolator: BC combination not implemented.'
    end select

    interpolator%coeff_splines(:) = 0.0_f64
    SLL_ALLOCATE( interpolator%t(num_pts*num_pts),ierr)
    interpolator%t(:) = 0.0_f64
  end subroutine initialize_ad1d_interpolator
  
  subroutine set_coefficients_ad1d( &
   interpolator, &
   coeffs)

   class(arb_deg_1d_interpolator), intent(inout)  :: interpolator
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
   eta_min = interpolator%eta_min
   eta_max = interpolator%eta_max
   delta = (eta_max - eta_min)/num_cells

   tmp = (sp_deg + 1)/2

   ! The interpretation and further filling of the spline coefficients array
   ! depends on the boundary conditions.
   select case (interpolator%bc_selector)
   case(0) ! periodic
      
      interpolator%size_coeffs =  num_cells + sp_deg
      interpolator%size_t = 2*sp_deg + num_cells +1
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
       
      
      do i = 1 ,nb_spline_eta
         interpolator%coeff_splines(i+1) = coeffs(i)
      end do
      
      do i = 1,nb_spline_eta
         
         interpolator%coeff_splines(i + 1 ) =  coeffs(i)
      end do
      
      interpolator%coeff_splines(1) = 0.0_8
      interpolator%coeff_splines(nb_spline_eta+2) = 0.0_8
      
      
   case default
      print *, 'arbitrary_degree_spline_1d() error: set_spline_coefficients ',&
           'not recognized.'
      stop
   end select
 end subroutine set_coefficients_ad1d

 subroutine compute_interpolants_ad1d( &
      interpolator,data_array, &
      eta_coords, &
      size_eta_coords)
   
   class(arb_deg_1d_interpolator), intent(inout)  :: interpolator
   sll_real64, dimension(:), intent(in)         :: data_array
   sll_real64, dimension(:), intent(in),optional  :: eta_coords
   sll_int32, intent(in),optional                 :: size_eta_coords
   sll_real64, dimension(:),pointer               :: point_locate_eta
   sll_real64 :: delta_eta
   sll_int32  :: sz
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
       
       select case (interpolator%bc_selector)
       case (0) ! 1. periodic
          sz = interpolator%num_pts-1
          
       case (9) ! 2. dirichlet-left, dirichlet-right
          sz = interpolator%num_pts
          
          
       case default
          print *, 'compute_interpolants_ad1d():BC combination not implemented.'
       end select

       delta_eta = (interpolator%eta_max - interpolator%eta_min)&
            /(interpolator%num_pts -1)
       SLL_ALLOCATE(point_locate_eta(sz),ierr)
       
       do i = 1,sz
          point_locate_eta(i) = interpolator%eta_min + delta_eta*(i-1)
       end do
    end if
    
    SLL_ASSERT(sz .le. interpolator%num_pts* interpolator%num_pts)
    SLL_ASSERT(size(data_array) .ge. sz)
    SLL_ASSERT(size(point_locate_eta)  .ge. sz)
    
    order  = interpolator%spline_degree + 1
    period = interpolator%eta_max - interpolator%eta_min
    
   ! print*, 'pointlocate',point_locate_eta1
    select case (interpolator%bc_selector)
    case (0) ! periodic
       interpolator%size_coeffs = sz+1
       interpolator%size_t = order + sz + 1 
       call spli1d_per( & ! a implementer
            period, sz+1, order, point_locate_eta, &
            data_array, interpolator%coeff_splines(1:sz+1),&
            interpolator%t(1:order + sz + 1))
       
       
     !  print*, 'moyenne', sum( interpolator%coeff_splines(1:sz+1))
       
    case (9) ! 2. dirichlet-left, dirichlet-right
       interpolator%size_coeffs = sz
       interpolator%size_t = order + sz ! a implementer
       !print*, 'data',data_array
       !print*, 'de',point_locate_eta
       call spli1d_dir( sz, order, point_locate_eta, &
            data_array, interpolator%coeff_splines(1:sz),&
            interpolator%t(1:sz+order) )
 
  
    end select
  end subroutine compute_interpolants_ad1d


  function interpolate_value_ad1d( &
       interpolator, &
       eta1) result(val)

    class(arb_deg_1d_interpolator), intent(inout)  :: interpolator
    sll_real64, intent(in)         :: eta1
    sll_real64                     :: val
    sll_int32 :: size_coeffs
    sll_real64 :: bvalue
    sll_real64 :: res

    size_coeffs = interpolator%size_coeffs

    res = eta1
    select case (interpolator%bc_selector)
    case (0) ! periodic
       if ( res .ge. interpolator%eta_max ) then 
          res = res -(interpolator%eta_max-interpolator%eta_min)
       end if
    case (9) ! 2. dirichlet-left, dirichlet-right
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
  
    end select
       
    !print*, 'coucou'
    val = bvalue( &
         interpolator%t(1:interpolator%size_t),&
         interpolator%coeff_splines(1:size_coeffs),&
         size_coeffs, &
         interpolator%spline_degree+1, &
         res,0)
  end function interpolate_value_ad1d


  function interpolate_derivative_ad1d( &
    interpolator, &
    eta1 ) result(val)

    class(arb_deg_1d_interpolator), intent(inout)  :: interpolator
    sll_real64, intent(in)         :: eta1
    sll_real64                     :: val
    sll_int32 :: size_coeffs
    sll_real64 :: dvalue1d
    sll_real64 :: res
    
    SLL_ASSERT( eta1 .ge. interpolator%eta_min )
    SLL_ASSERT( eta1 .le. interpolator%eta_max )

    size_coeffs = interpolator%size_coeffs

    res = eta1
    
    select case (interpolator%bc_selector)
    case (0) ! periodic
       if ( res .ge. interpolator%eta_max ) then 
          res = res -(interpolator%eta_max-interpolator%eta_min)
       end if
    case (9) ! 2. dirichlet-left, dirichlet-right
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       
    end select
    
    val = dvalue1d( &
         res, &
         size_coeffs, &
         interpolator%spline_degree+1, &
         interpolator%coeff_splines(1:size_coeffs), &
         interpolator%t(1:interpolator%size_t),&
         1)
    
  end function interpolate_derivative_ad1d
  
  function interpolate_array_ad1d( &
       this, num_points, data, coordinates)&
    result(data_out)
    class(arb_deg_1d_interpolator),  intent(in)       :: this
    !class(sll_cubic_spline_1D),  intent(in)      :: this
    sll_int32,  intent(in)                 :: num_points
    sll_real64, dimension(:), intent(in)   :: coordinates
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_points)      :: data_out
    
    print *, 'interpolate_array_ad1d: not implemented'
  end function interpolate_array_ad1d
  
  function interpolate_1d_array_disp_ad1d( &
       this,        &
       num_points, &
       data,     &
       alpha) result(res)
      
    class(arb_deg_1d_interpolator), intent(in)    :: this
    sll_int32, intent(in)                          :: num_points 
    sll_real64, dimension(:), intent(in)         :: data
    sll_real64, intent(in)         :: alpha  
    sll_real64, dimension(num_points) :: res
    
    print *, 'interpolate_1d_array_disp_ad1d: not implemented.'
  end function interpolate_1d_array_disp_ad1d
    
    
  function get_coefficients_ad1d(interpolator)
    class(arb_deg_1d_interpolator), intent(in)    :: interpolator
    sll_real64, dimension(:), pointer            :: get_coefficients_ad1d     
    
    get_coefficients_ad1d => interpolator%coeff_splines
  end function get_coefficients_ad1d

  subroutine interpolate_values_ad1d( &
       interpolator, &
       num_pts, &
       vals_to_interpolate, &
       output_array )
    
    class(arb_deg_1d_interpolator),  intent(in) :: interpolator
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(:), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(:), intent(out)  :: output_array
    
    print*, 'interpolate_values_ad1d NOT iMPLEMENTED YET'
  end subroutine interpolate_values_ad1d

  subroutine interpolate_pointer_values_ad1d( &
       interpolator, &
       num_pts, &
       vals_to_interpolate, &
       output )
    
    class(arb_deg_1d_interpolator),  intent(in) :: interpolator
    sll_int32,  intent(in)            :: num_pts
    sll_real64, dimension(:), pointer :: vals_to_interpolate
    sll_real64, dimension(:), pointer :: output
    
    print*, 'interpolate_pointer_values_ad1d NOT iMPLEMENTED YET'
  end subroutine interpolate_pointer_values_ad1d
  
  subroutine interpolate_derivatives_ad1d( &
       interpolator, &
       num_pts, &
       vals_to_interpolate, &
       output_array )
    
    class(arb_deg_1d_interpolator),  intent(in) :: interpolator
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(:), intent(in) :: vals_to_interpolate
    sll_real64, dimension(:), intent(out) :: output_array
    
    print*, 'interpolate_derivatives_ad1d NOT iMPLEMENTED YET'
  end subroutine interpolate_derivatives_ad1d
  
  subroutine interpolate_pointer_derivatives_ad1d( &
       interpolator, &
       num_pts, &
       vals_to_interpolate, &
       output )
#ifdef STDF95
    type(arb_deg_1d_interpolator),  intent(in) :: interpolator
#else
    class(arb_deg_1d_interpolator),  intent(in) :: interpolator
#endif
    sll_int32,  intent(in)              :: num_pts
    sll_real64, dimension(:), pointer   :: vals_to_interpolate
    sll_real64, dimension(:), pointer   :: output
    
    print*, 'interpolate_pointer_derivatives_ad1d NOT iMPLEMENTED YET'
  end subroutine interpolate_pointer_derivatives_ad1d
  
  function reconstruct_array(this, num_points, data) result(res)
    ! dummy procedure

    class(arb_deg_1d_interpolator),  intent(in) :: this
    sll_int32, intent(in)                :: num_points! size of output array
    sll_real64, dimension(:), intent(in) :: data   ! data to be interpolated 
    sll_real64, dimension(num_points)    :: res
    res(:) = 0.0_f64
  end function reconstruct_array
end module sll_arbitrary_degree_spline_interpolator_1d_module
