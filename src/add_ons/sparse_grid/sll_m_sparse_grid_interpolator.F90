!> @ingroup sparse_grid
!> @brief
!> Dimension-independent functions for sparse grid with polynomial basis functions
!> @details
!> Implements the sll_m_sparse_grid_interpolator interface
module sll_m_sparse_grid_interpolator
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_fftw.h"

  use iso_c_binding, only: &
    c_double_complex, &
    c_size_t

  use sll_m_fft, only : &
    sll_t_fft, &
    sll_s_fft_init_c2c_1d, &
    sll_s_fft_exec_c2c_1d, &
    sll_s_fft_free, &
    sll_f_fft_allocate_aligned_complex, &
    sll_p_fft_forward, &
    sll_p_fft_backward, &
    sll_p_fft_measure

  use sll_m_interpolators_1d_base, only: &
    sll_c_interpolator_1d

  use sll_m_periodic_interpolator_1d, only: &
    sll_t_periodic_interpolator_1d

  implicit none

  public :: &
    sll_t_sparse_grid_interpolator

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> class to hold values for hierarchical fft computations
  type :: fft_hierarchical
     type( sll_t_fft ) :: fw
     type( sll_t_fft ) :: bw
     complex(c_double_complex), dimension(:), pointer :: in
     complex(c_double_complex), dimension(:), pointer :: out
     integer(c_size_t) :: sz_in
     integer(c_size_t) :: sz_out
  end type fft_hierarchical

  !> Data type for sparse grid node
  type sparsegrid_node
     sll_real64, dimension(:), allocatable   :: coordinate !< \a coordinate holds the coordinates of the node
     sll_int32, dimension(:), allocatable     :: parent !< \a parent holds the neighbors of the node on the level of the node
     ! type(sparsegrid_node_ptr),dimension(:), pointer :: parents
     sll_int32, dimension(:), allocatable     :: children !< \a children holds the children of the node (child 2i-1 (2i): left (right) neighbor along dimension i
     sll_int32, dimension(:), allocatable     :: function_type !< \a function_type: integer defining the basis function along each dimension
     sll_int32, dimension(:), allocatable     :: level !< \a level defines the refinement level along each dimension
     sll_int32, dimension(:), allocatable     :: index_on_level !< \a index_on_level defines the number of the node within the level hierarchy
     !sll_real64, dimension(:), allocatable    :: values
  end type sparsegrid_node


  type :: interpolator_base_ptr
     class(sll_c_interpolator_1d), pointer :: ptr
  end type interpolator_base_ptr



!> Class defining the sparse grid data structure
  type :: sll_t_sparse_grid_interpolator
     sll_int32                          :: dim !< \a dim defines the dimension of the sparse grid 
     sll_real64,dimension(:), pointer    :: eta_min !< \a eta_min defines the lower bound of the computational domain
     sll_real64, dimension(:), pointer            :: eta_max !< \a eta_max defines the upper bound of the computational domain
     sll_real64,dimension(:), pointer             :: length !< \a length defines the length of the computational domain
     sll_real64                         :: volume !< \a volumne defines the volume of the computational domain
     sll_int32                           :: max_level !< \a max_level is the upper bound of the l_infinity norm of the level vector
     sll_int32, dimension(:), pointer     :: levels !< \a levels is the upper bound on the level in each direction
     sll_int32                           :: order !< \a order it the (maximal) order of basis functions
     sll_int32                           :: interpolation !< \a specifies type of interpolation
     sll_int32                           :: no_basis_functions !< \a no_basis_functions is the number of different basis functions
     sll_int32                           :: size_basis !< \a size_basis is the number of grid points
     sll_int32                           :: modified !< \a modified specifies the number of boundary points in sparse grid (0 for usual sparse grid, 1 for as many boundary as middle points)
     sll_int32                           :: boundary !< \a boundary specifies if the boundary is periodic (0) or usual sparse grid (1)

     class(sparsegrid_node), dimension(:), pointer     :: hierarchy !< \a hierarchy is an array containing the nodes

     class(fft_hierarchical), dimension(:), pointer :: fft_object !< \a fft_object is the object for hierarchical fft
     
     !sll_int32, dimension(:,:,:,:), pointer  :: index
     sll_real64,dimension(:), pointer :: stripe,stripe_out !< \a stripe and \a stripe are internal arrays to handle the 1d interpolations
     sll_real64,dimension(:), pointer :: hs_weights !< \a hs_weights contains the weights for the computation of the hierarchical surplus
     sll_int32, dimension(:), pointer ::  hs_weights_index !< \a hs_weights_index is the index vector for \a hs_weights
     type(sll_t_periodic_interpolator_1d),dimension(:,:), pointer :: interp_per !< \a interp_per is the periodic interpolator object along the stripes
     !real,dimension(:), target      :: interp_per_x
     !type(per_1d_interpolator),dimension(:,:), pointer :: interp_v
    ! type(odd_degree_spline_1d_interpolator),dimension(:,:), pointer :: interp_v
     !type(lagrange_1d_interpolator),dimension(:,:), pointer :: interpl_v
     type(interpolator_base_ptr), dimension(:,:), pointer  :: interp !< \a interp is the interpolator object for the 1d interpolations along the stripes
     !type(sll_c_interpolator_1d), dimension(:,:), pointer  :: interp
     sll_int32, dimension(:), pointer :: level_mapping !< \a level_mapping is an index pointing the the start of each level

   contains
     procedure :: compute_hierarchical_surplus !< Compute the hierarchical surplus (to the order of the sparse grid)
     procedure :: compute_linear_hierarchical_surplus !< Compute the linear hierarchical surplus
     procedure :: compute_dehierarchical !< Compute values from (order) hierarchical surplus
     procedure :: dehierarchical
     procedure :: dehierarchical_part
     procedure :: hierarchical_part
     procedure :: interpolate_disp !< Interpolation function for interpolation at (constantly) displaced grid points; displacement only in dimension dim
     procedure :: integrate_trapezoidal
     procedure :: integrate_trapezoidal2
     procedure :: extract_periodic
     procedure :: hierarchical_stripe
     procedure :: dehierarchical_stripe
     procedure :: insert_periodic
     procedure :: interpolate_disp_1d_periodic_self
     procedure, nopass :: basis_function
     procedure, nopass :: basis_function_derivative
     procedure :: dehierarchical_part_order
     procedure :: interpolate_disp_1d_periodic_for_neighbor
     procedure :: interpolate_disp_1d_periodic
     procedure :: displace1d
     procedure :: tonodal1d
     procedure :: tonodal1d_comp
     procedure :: todehi1d
     procedure :: tohira1d
     procedure :: tohierarchical1d
     procedure :: tohierarchical1d_comp
     procedure :: initialize_sg
     procedure :: hierarchical_part_order
     procedure :: dehierarchical_order
     procedure :: free => free_sparse_grid

  end type sll_t_sparse_grid_interpolator

contains


!------------------------------------------------------------------------------!
!!!! Helper functions for initialization

 subroutine initialize_sg( &
    interpolator, &
    levels, &
    order, &
    interpolation, &
    interpolation_type, &
    eta_min, &
    eta_max)
    class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator

    sll_real64, dimension(:), intent(in)          :: eta_min !< \a eta_min defines the lower bound of the domain
    sll_real64, dimension(:),  intent(in)         :: eta_max !< \a eta_max defines the upper bound of the domain
    sll_int32, dimension(:), intent(in)           :: levels !< \a levels defines the maximum level in the sparse grid
    sll_int32, intent(in)                         :: order !< \a order of the sparse grid functions
    sll_int32, intent(in)                         :: interpolation !< \a Order of the interpolator
    sll_int32, intent(in)                         :: interpolation_type !< Choose spline (\a interpolation_type = 0) or Lagrange (\a interpolation_type = 1) interpolation for the 1D interpolators if not traditional sparse grid interpolation is used.

    sll_int32                                     :: i,j
    sll_int32                                     :: ierr

    SLL_ALLOCATE(interpolator%eta_min(interpolator%dim), ierr);
    SLL_ALLOCATE(interpolator%eta_max(interpolator%dim), ierr);
    SLL_ALLOCATE(interpolator%length(interpolator%dim), ierr);
    SLL_ALLOCATE(interpolator%levels(interpolator%dim), ierr);

    !SLL_ALLOCATE( interpolator, ierr )
    interpolator%volume = 1.0_f64;
    do j=1,interpolator%dim
       interpolator%eta_min(j) = eta_min(j);
       interpolator%eta_max(j) = eta_max(j);
       interpolator%length(j)  = eta_max(j) - eta_min(j);
       interpolator%volume = interpolator%volume * interpolator%length(j);
    end do
    interpolator%levels     = levels;
    SLL_ALLOCATE(interpolator%level_mapping(0:interpolator%max_level+1),ierr);
    interpolator%order = order
    interpolator%interpolation = interpolation;
    if(order == 1) then
       interpolator%no_basis_functions = 2
    elseif(order == 2) then
       interpolator%no_basis_functions = 3
    elseif(order == 3) then
       interpolator%no_basis_functions = 7
    elseif(order == 4) then
       interpolator%no_basis_functions = 15
    end if

    print*, "Size of the basis: ", interpolator%size_basis
    SLL_ALLOCATE(interpolator%hierarchy(interpolator%size_basis),ierr)
    do j=1,interpolator%size_basis
       SLL_ALLOCATE(interpolator%hierarchy(j)%coordinate(interpolator%dim),ierr);
       SLL_ALLOCATE(interpolator%hierarchy(j)%parent(2*interpolator%dim),ierr);
       SLL_ALLOCATE(interpolator%hierarchy(j)%children(2*interpolator%dim),ierr);
       SLL_ALLOCATE(interpolator%hierarchy(j)%level(interpolator%dim),ierr);
       SLL_ALLOCATE(interpolator%hierarchy(j)%index_on_level(interpolator%dim),ierr);
       SLL_ALLOCATE(interpolator%hierarchy(j)%function_type(interpolator%dim),ierr);
    end do

    SLL_ALLOCATE(interpolator%stripe(2**interpolator%max_level+1),ierr)
    SLL_ALLOCATE(interpolator%stripe_out(2**interpolator%max_level+1),ierr)

    do j=1,interpolator%size_basis
       interpolator%hierarchy(j)%children = -1
    end do


    if(interpolator%order>1) then
       SLL_ALLOCATE(interpolator%hs_weights_index(interpolator%order+1),ierr)
       interpolator%hs_weights_index(1) = 1
       interpolator%hs_weights_index(2) = 1
       interpolator%hs_weights_index(3) = 3
       if(interpolator%order==2) then
          SLL_ALLOCATE(interpolator%hs_weights(2),ierr)
       elseif(interpolator%order==3) then
          interpolator%hs_weights_index(4) = 7
          SLL_ALLOCATE(interpolator%hs_weights(6),ierr)
          interpolator%hs_weights(3) = -0.125_f64
          interpolator%hs_weights(4) = -interpolator%hs_weights(3)
          interpolator%hs_weights(5) = interpolator%hs_weights(4)
          interpolator%hs_weights(6) = interpolator%hs_weights(3)
       elseif(interpolator%order==4) then
          interpolator%hs_weights_index(4) = 7
          interpolator%hs_weights_index(5) = 15
          SLL_ALLOCATE(interpolator%hs_weights(14),ierr)
          interpolator%hs_weights(3) = -0.125_f64
          interpolator%hs_weights(4) = -interpolator%hs_weights(3)
          interpolator%hs_weights(5) = interpolator%hs_weights(4)
          interpolator%hs_weights(6) = interpolator%hs_weights(3)
          interpolator%hs_weights(7) = -1.0_f64/16.0_f64
          interpolator%hs_weights(8) = 5.0_f64/112.0_f64
          interpolator%hs_weights(9) = interpolator%hs_weights(7)
          interpolator%hs_weights(10) = 7.0_f64/80.0_f64
          interpolator%hs_weights(11) = 7.0_f64/80.0_f64
          interpolator%hs_weights(12) = interpolator%hs_weights(7)
          interpolator%hs_weights(13) = interpolator%hs_weights(8)
          interpolator%hs_weights(14) = interpolator%hs_weights(7)
       end if
       interpolator%hs_weights(1) = -0.25_f64
       interpolator%hs_weights(2) = interpolator%hs_weights(1)
    end if

    SLL_ALLOCATE(interpolator%interp_per(interpolator%dim,interpolator%max_level),ierr);
    !SLL_ALLOCATE(interpolator%interp_v(interpolator%dim/2,interpolator%max_level),ierr);
    !SLL_ALLOCATE(interpolator%interpl_v(interpolator%dim/2,interpolator%max_level),ierr);
    SLL_ALLOCATE(interpolator%interp(interpolator%dim, interpolator%max_level), ierr);

    do i=1,interpolator%max_level
       do j=1,interpolator%dim
          if (interpolation_type == 0) then
             call interpolator%interp_per(j,i)%initialize( 2**i + 1, &
                  interpolator%eta_min(j), &
                  interpolator%eta_max(j),&
                  1, interpolation);
            ! call interpolator%interp_v(j,i)%initialize( 2**i + 1, &
            !      interpolator%eta_min(j+1), &
            !      interpolator%eta_max(j+1),&
            !      interpolation-1);
 
          else
             call interpolator%interp_per(j,i)%initialize( 2**i +1,&
                  interpolator%eta_min(j), &
                  interpolator%eta_max(j),&
                  2, interpolation);
             !call interpolator%interpl_v(j,i)%initialize( 2**i +1,&
             !     interpolator%eta_min(j+interpolator%dim/2), &
             !     interpolator%eta_max(j+interpolator%dim/2),&
             !     HERMITE_LAGRANGE,interpolation/2);
          end if

 !         interpolator%interp(j,i)=interpolator%interp_x(j,i);
 !         interpolator%interp(j+interpolator%dim/2,i) = interpolator%interp_v(j,i);
       end do
    end do

    ! For FFT
    ALLOCATE(interpolator%fft_object(interpolator%max_level+1))
    call fft_initialize(interpolator%fft_object, interpolator%max_level)

  end subroutine initialize_sg



!!!! End helper functions for initialization
!------------------------------------------------------------------------------!






!------------------------------------------------------------------------------!
!!!! Functions to build hierarchical surplus !!!!

! Public function to be called
 subroutine compute_hierarchical_surplus( interpolator, data_array )
    class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
    sll_real64, dimension(:), intent(inout) :: data_array
    sll_int32                              :: i

    call hierarchical(interpolator,data_array);

    do i=2,interpolator%order
       call hierarchical_order(interpolator,data_array,i);
    end do

  end subroutine compute_hierarchical_surplus

 subroutine compute_linear_hierarchical_surplus( interpolator, data_array )
    class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
    sll_real64, dimension(:), intent(inout) :: data_array
 
    call hierarchical(interpolator,data_array);

  end subroutine compute_linear_hierarchical_surplus

! Public function to be called
 subroutine compute_dehierarchical( interpolator, data_array )
    class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
    sll_real64, dimension(:), intent(inout) :: data_array
    sll_int32                              :: i


    do i=interpolator%order,2,-1
       call dehierarchical_order(interpolator,data_array,i);
    end do

    call dehierarchical(interpolator,data_array);

  end subroutine compute_dehierarchical


!!!! End functions to build hierarchical surplus !!!!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!!!! Helper functions for interpolation routines !!!!

!!!! End helper functions for interpolation routines !!!!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!!!! Evaluation of basis functions !!!!

  subroutine basis_function(x,fx,type) 
    sll_real64,  intent(in)                 :: x
    sll_real64,  intent(inout)              :: fx
    sll_int32, intent(in)                   :: type

    select case (type)
    case(-1)
       !print*, -1
       if ((x>=-1.0_f64) .AND. (x<0.0_f64)) then
          fx = 1.0_f64+x
       elseif  ((x>=0.0_f64) .AND. (x<1.0_f64)) then
          fx = 1.0_f64-x
       else
          fx = 0.0_f64
       end if
    case(0)
       !print*, 0
       fx = 1.0_f64; ! CHANGE_CONSTANT
       !if ((x>=-1.0_f64) .AND. (x<0.0_f64)) then
        !  fx = -x
       !elseif  ((x>=0.0_f64) .AND. (x<1.0_f64)) then
       !   fx = x
       !else
       !   fx = 0.0_f64
       !end if
    case(1, 2)
       if ((x>-1.0_f64) .AND. (x<1.0_f64)) then
          fx = (1.0_f64-x)*(x+1.0_f64)
       else
          fx = 0.0_f64
       end if
    case(3, 5)
       if ((x>-1.0_f64) .AND. (x<1.0_f64)) then
          fx = (1.0_f64-x)*(x+1.0_f64)*(1.0_f64-x/3.0_f64)
       else
          fx = 0.0_f64
       end if
    case(4, 6)
       if ((x>-1.0_f64) .AND. (x<1.0_f64)) then
          fx = (1.0_f64-x)*(x+1.0_f64)*(1.0_f64+x/3.0_f64)
       else
          fx = 0.0_f64
       end if
    case(7, 11)
       if ((x>-1.0_f64) .AND. (x<1.0_f64)) then
          fx = (1.0_f64-x)*(1.0_f64+x)*(1.0_f64-x/3.0_f64)*(1.0_f64-x/7.0_f64)
       else
          fx = 0.0_f64
       end if
    case(8, 12)
       if ((x>-1.0_f64) .AND. (x<1.0_f64)) then
          fx = (1.0_f64-x)*(1.0_f64+x)*(1.0_f64+x/3.0_f64)*(1.0_f64-x/5.0_f64)
       else
          fx = 0.0_f64
       end if
    case(9, 13)
       if ((x>-1.0_f64) .AND. (x<1.0_f64)) then
          fx = (1.0_f64-x)*(1.0_f64+x)*(1.0_f64-x/3.0_f64)*(1.0_f64+x/5.0_f64)
       else
          fx = 0.0_f64
       end if
    case(10,  14)
       if ((x>-1.0_f64) .AND. (x<1.0_f64)) then
          fx = (1.0_f64-x)*(1.0_f64+x)*(1.0_f64+x/3.0_f64)*(1.0_f64+x/7.0_f64)
       else
          fx = 0.0_f64
       end if
    end select
  end subroutine basis_function

  subroutine basis_function_derivative(x,fx,type) 
    sll_real64,  intent(in)                 :: x
    sll_real64,  intent(inout)              :: fx
    sll_int32, intent(in)                   :: type

    select case (type)
    case(0)       
       fx = 0.0_f64! CHANGE_CONSTANT
       !if ((x>=-1.0_f64) .AND. (x<0.0_f64)) then
       !   fx = -1.0_f64
       !elseif  ((x>=0.0_f64) .AND. (x<1.0_f64)) then
       !   fx = 1.0_f64
       !else
       !   fx = 0.0_f64
       !end if
    case(1)
       if ((x>=-1.0_f64) .AND. (x<1.0_f64)) then
          fx = -2.0_f64*x
       else
          fx = 0.0_f64
       end if
    case(2)
       if ((x>=-1.0_f64) .AND. (x<1.0_f64)) then
          fx = -1.0_f64/3.0_f64-2.0_f64*x+x*x
       else
          fx = 0.0_f64
       end if
    case(3)
       if ((x>=-1.0_f64) .AND. (x<1.0_f64)) then
          fx = 1.0_f64/3.0_f64-2.0_f64*x-x*x
       else
          fx = 0.0_f64
       end if
    end select
 
  end subroutine basis_function_derivative


!!!! End evaluation of basis functions !!!!
!------------------------------------------------------------------------------!



!------------------------------------------------------------------------------!
!!!! Helper functions on 1D stripes (extract, insert, displace, interpolate) !!!


! Displace functions for 1D
subroutine displace_on_stripe_periodic_for_neighbor(interpolator,displacement,dim, max_level,max_level_neighbor)
    class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64, intent(in) :: displacement
  sll_int32 , intent(in) :: max_level,max_level_neighbor,dim
  sll_int32 :: cell,size_fraction, size_neighbor
  
  call  displace_on_stripe_periodic(interpolator,displacement,dim, max_level)

  size_neighbor = 2**max_level_neighbor
  size_fraction = 2**max_level/size_neighbor

  do cell=1,size_neighbor
     interpolator%stripe_out(cell) = interpolator%stripe_out(1+(cell-1)*size_fraction)
  end do
  
end subroutine displace_on_stripe_periodic_for_neighbor

! end functions for neighbors


subroutine displace_on_stripe_periodic(interpolator,displacement,dim, max_level)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64, intent(in) :: displacement
  sll_int32 , intent(in) :: max_level,dim
  sll_int32 :: size
  
  size = 2**max_level+1;
  if(max_level == 0) then
     interpolator%stripe_out(1) = interpolator%stripe(1)
  else
     interpolator%stripe(size) = interpolator%stripe(1);
     call interpolator%interp_per(dim,max_level)%interpolate_array_disp(size, interpolator%stripe(1:size), displacement, interpolator%stripe_out(1:size))
  end if

end subroutine displace_on_stripe_periodic


! Interpolator functions in 1D

! Interpolate along dimension dim only (periodic boundary conditions), interpolation at a displaced function value
subroutine interpolate_disp_1d_periodic(interpolator,displacement,dim,max_level,index,data_in,data_out,hiera)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(in) :: data_in
  sll_real64, dimension(:), intent(out) :: data_out
  sll_real64, intent(in) :: displacement
  sll_int32, intent(in) :: dim,max_level,index
  logical, intent(in) :: hiera

  call extract_periodic(interpolator,dim,max_level,&
       index,data_in,interpolator%stripe)

  call dehierarchical_stripe(interpolator,&
       interpolator%stripe,max_level)

  call displace_on_stripe_periodic(interpolator,displacement,dim,max_level)

  if (hiera) then
     call hierarchical_stripe (interpolator,&
          interpolator%stripe_out, max_level);
  end if

  call insert_periodic(interpolator,dim,max_level,&
       index,interpolator%stripe_out,data_out)

end subroutine interpolate_disp_1d_periodic

recursive subroutine interpolate_disp_recursive(interpolator,no_dims,dim,node_index, displacement,data_in, data_out,hiera)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64, intent(in) :: displacement
  sll_real64, dimension(:), intent(in) :: data_in
  sll_real64, dimension(:), intent(out) :: data_out
  sll_int32, intent(in) :: dim, no_dims
  sll_int32, intent(in) :: node_index
  logical, intent(in) :: hiera
  sll_int32 :: j, child_index,max_level


  max_level = interpolator%max_level!interpolator%levels(dim)
  do j=1,no_dims
     if (j .NE. dim) then
        max_level = max_level - interpolator%hierarchy(node_index)%level(j);
     end if
  end do
  max_level = min(max_level,interpolator%levels(dim))

  call interpolate_disp_1d_periodic(interpolator,displacement,dim,&
       max_level,node_index,data_in,data_out,hiera)
  
  do j=1,2*no_dims
     if((j+1)/2 .NE. dim) then
        child_index = interpolator%hierarchy(node_index)%children(j);
        if (child_index > 0) then
           call interpolate_disp_recursive(interpolator,no_dims,dim,&
                child_index,displacement,data_in,data_out,hiera)
        end if
     end if
  end do

end subroutine interpolate_disp_recursive

! Interpolation function for interpolation at (constantly) displaced grid points; displacement only in dimension dim
subroutine interpolate_disp(interpolator,dim,displacement,data_in, data_out,hiera)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(inout) :: data_in
  sll_real64, dimension(:), intent(out) :: data_out
  sll_int32, intent(in) :: dim
  sll_real64, intent(in) ::displacement
  sll_int32 :: j,jj
  sll_int32, dimension(:), allocatable :: dorder
  logical, intent(in) :: hiera

  SLL_ALLOCATE(dorder(interpolator%dim),j);
  dorder(1) = dim;
  jj = 2;
  do j=1,interpolator%dim
     if(j .NE. dim) then
        dorder(jj) = j;
        jj = jj+1;
     end if
  end do  

  call interpolate_disp_recursive(interpolator,interpolator%dim,dim,1,&
       displacement, data_in, data_out,hiera);

  if (hiera .EQV. .FALSE.) then
     ! Dehierarchization along dimension dorder(1) only
     do j=interpolator%order,2,-1
        call dehierarchical_part_order&
             (interpolator,data_out,&
             interpolator%dim,2,dorder,j)
     end do

     call dehierarchical_part(interpolator,data_out,&
          interpolator%dim,2,dorder)
  end if


!!$  SLL_ALLOCATE(dorder(interpolator%dim),j);
!!$  dorder = (1:interpolator%dim);
!!$  dorder(1) = dim;
!!$  jj = 2;
!!$  do j=1,interpolator%dim
!!$     if(j .NE. dim) then
!!$        dorder(jj) = j;
!!$        jj = jj+1;
!!$     end if
!!$  end do
!!$
!!$  call interpolate_disp_recursive_self(interpolator,interpolator%dim,dim,1,displacement, data_in, data_out)
!!$
!!$
!!$  ! Dehierarchization along dimension dorder(1) only
!!$  do j=interpolator%order,2,-1
!!$     call dehierarchical_part_order&
!!$          (interpolator%sparse_grid_interpolator,data_in,&
!!$          interpolator%dim,2,dorder,j)
!!$  end do
!!$
!!$  call dehierarchical_part(interpolator%sparse_grid_interpolator,data_in,&
!!$       interpolator%dim,2,dorder)


end subroutine Interpolate_disp


!!!!!
! Interpolates along dimensions "dim" on the stripe defined by "index". Then result is only computed for every "size_fraction" element and inserted into the global vector for the stripe defined by "index_neighbor".
!!!!!!!!
subroutine interpolate_disp_1d_periodic_for_neighbor(interpolator,displacement,factor,dim,max_level,max_level_neighbor,index,index_neighbor,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(in) :: data_in
  sll_real64, dimension(:), intent(out) :: data_out
  sll_real64, intent(in) :: displacement, factor
  sll_int32, intent(in) :: dim,max_level,index, index_neighbor,max_level_neighbor


  call extract_periodic(interpolator,dim,max_level,&
       index,data_in,interpolator%stripe)
  !call dehi_periodic(interpolator,&
  !     interpolator%stripe,max_level)

  call displace_on_stripe_periodic_for_neighbor(interpolator,displacement,dim,&
       max_level,max_level_neighbor)

  !if(index_neighbor == 4) then
  !   print*, interpolator%stripe_out(1:2**max_level_neighbor)
  !end if


  !call hira_periodic (interpolator,&
  !     interpolator%stripe_out, max_level)
  call insert_periodic_additive(interpolator,factor,dim,max_level_neighbor,&
       index_neighbor,interpolator%stripe_out,data_out)

end subroutine interpolate_disp_1d_periodic_for_neighbor


subroutine interpolate_disp_1d_periodic_self(interpolator,displacement,dim,max_level,index,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(in) :: data_in
  sll_real64, dimension(:), intent(out) :: data_out
  sll_real64, intent(in) :: displacement
  sll_int32, intent(in) :: dim,max_level,index

  call extract_periodic(interpolator,dim,max_level,&
       index,data_in,interpolator%stripe)
 
  call displace_on_stripe_periodic(interpolator,displacement,dim,max_level)

  !call hira_periodic (interpolator,&
  !     interpolator%stripe_out, max_level)
  call insert_periodic(interpolator,dim,max_level,&
       index,interpolator%stripe_out,data_out)

end subroutine interpolate_disp_1d_periodic_self

! Extract functions (extract 1D stripe from dD)




subroutine extract_periodic(sparsegrid,dim,max_level,index,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(in) :: sparsegrid
  sll_int32, intent(in) :: dim,max_level,index
  sll_int32             :: n_points,index_running
  sll_real64,dimension(:),intent(in) :: data_in
  sll_real64,dimension(:),intent(out) :: data_out

  n_points = 2**(max_level)
  data_out(1) = data_in(index)
  if (max_level>0) then
     index_running = sparsegrid%hierarchy(index)%children(dim*2)
     data_out(n_points/2+1) = data_in(index_running)
     if (max_level>1) then   
        call extract_recursive(sparsegrid,n_points/2+1,n_points/4,&
             index_running,dim,data_in,data_out)
     end if
  end if


end subroutine extract_periodic

recursive subroutine extract_recursive(sparsegrid,index_stripe,stride,index_sg,dim,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(in) :: sparsegrid
  sll_int32, intent(in) :: index_stripe,stride,index_sg,dim
  sll_real64,dimension(:),intent(in) :: data_in
  sll_real64,dimension(:),intent(inout) :: data_out
 !print*, 'expr', index_sg,index_stripe
 
  data_out(index_stripe-stride) = &
       data_in(sparsegrid%hierarchy(index_sg)%children(dim*2-1))
  data_out(index_stripe+stride) = &
       data_in(sparsegrid%hierarchy(index_sg)%children(dim*2))
  if (stride>1) then
     call extract_recursive(sparsegrid,index_stripe-stride,stride/2,&
          sparsegrid%hierarchy(index_sg)%children(dim*2-1),dim,&
          data_in,data_out)
     call extract_recursive(sparsegrid,index_stripe+stride,stride/2,&
          sparsegrid%hierarchy(index_sg)%children(dim*2),dim,&
          data_in,data_out)
  end if
end subroutine extract_recursive


! Insert functions (write 1D stripe back to dD)



subroutine insert_periodic(sparsegrid,dim,max_level,index,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(in) :: sparsegrid
  sll_int32, intent(in) :: dim,max_level,index
  sll_int32             :: n_points,index_running
  sll_real64,dimension(:),intent(in) :: data_in
  sll_real64,dimension(:),intent(inout) :: data_out

  
  n_points = 2**(max_level)
  data_out(index) = data_in(1)
  if (max_level>0) then
     index_running = sparsegrid%hierarchy(index)%children(dim*2)
     
     data_out(index_running) = data_in(n_points/2+1)
     if (max_level>1) then
        
        call insert_recursive(sparsegrid,n_points/2+1,n_points/4,index_running,&
             dim,data_in,data_out)
     end if
  end if


end subroutine insert_periodic

recursive subroutine insert_recursive(sparsegrid,index_stripe,stride,index_sg,dim,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(in) :: sparsegrid
  sll_int32, intent(in) :: index_stripe,stride,index_sg,dim
  sll_real64,dimension(:),intent(in) :: data_in
  sll_real64,dimension(:),intent(inout) :: data_out

  data_out(sparsegrid%hierarchy(index_sg)%children(dim*2-1)) =&
       (data_in(index_stripe-stride)) 
  data_out(sparsegrid%hierarchy(index_sg)%children(dim*2)) =&
       (data_in(index_stripe+stride))
  if (stride>1) then
     call insert_recursive(sparsegrid,index_stripe-stride,stride/2,&
          sparsegrid%hierarchy(index_sg)%children(dim*2-1),dim,&
          data_in,data_out)
     call insert_recursive(sparsegrid,index_stripe+stride,stride/2,&
          sparsegrid%hierarchy(index_sg)%children(dim*2),dim,&
          data_in,data_out)
  end if
end subroutine insert_recursive


subroutine insert_periodic_additive(sparsegrid,factor,dim,max_level,index,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(in) :: sparsegrid
  sll_real64, intent(in) :: factor
  sll_int32, intent(in) :: dim,max_level,index
  sll_int32             :: n_points,index_running
  sll_real64,dimension(:),intent(in) :: data_in
  sll_real64,dimension(:),intent(inout) :: data_out

  
  n_points = 2**(max_level)
  data_out(index) =  data_out(index)+ data_in(1)*factor
  if (max_level>0) then
     index_running = sparsegrid%hierarchy(index)%children(dim*2)
     
     data_out(index_running) = data_out(index_running)+data_in(n_points/2+1)*factor
     if (max_level>1) then
        
        call insert_additive_recursive(sparsegrid,factor,n_points/2+1,n_points/4,index_running,&
             dim,data_in,data_out)
     end if
  end if


end subroutine insert_periodic_additive

recursive subroutine insert_additive_recursive(sparsegrid,factor,index_stripe,stride,index_sg,dim,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(in) :: sparsegrid
  sll_real64, intent(in) :: factor
  sll_int32, intent(in) :: index_stripe,stride,index_sg,dim
  sll_real64,dimension(:),intent(in) :: data_in
  sll_real64,dimension(:),intent(inout) :: data_out

  data_out(sparsegrid%hierarchy(index_sg)%children(dim*2-1)) =&
       data_out(sparsegrid%hierarchy(index_sg)%children(dim*2-1))+(data_in(index_stripe-stride))*factor 
  !TODO: Should this be there or not?
  !if(sparsegrid%hierarchy(index_sg)%children(dim*2) .NE. sparsegrid%hierarchy(index_sg)%children(dim*2-1)) then
     data_out(sparsegrid%hierarchy(index_sg)%children(dim*2)) =&
          data_out(sparsegrid%hierarchy(index_sg)%children(dim*2))+(data_in(index_stripe+stride))*factor
  !end if
  if (stride>1) then
     call insert_additive_recursive(sparsegrid,factor,index_stripe-stride,stride/2,&
          sparsegrid%hierarchy(index_sg)%children(dim*2-1),dim,&
          data_in,data_out)
     call insert_additive_recursive(sparsegrid,factor,index_stripe+stride,stride/2,&
          sparsegrid%hierarchy(index_sg)%children(dim*2),dim,&
          data_in,data_out)
  end if
end subroutine insert_additive_recursive

!!!! End helper function on 1D stripes !!!!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!!!! General sparse grid helper functions !!!!

! Various dehierarchical routines (for computing hierarchical surplus)

! Hierarchization of the sparse grid linear hierarchical surplus (data on out) (along all dimensions)
subroutine hierarchical(interpolator,data)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(inout)   :: data
  !sll_real64, dimension(:), intent(out)  :: data_out
  sll_int32                                :: counter
  sll_real64 :: factor
  
 
  do counter=interpolator%size_basis,1,-1
     !write(16,*), 'Now:', l
     !data_out(counter) = 0.0_f64!data_array(counter,1)
     factor = 1.0_f64
     call dehierarchical_d_dimension&
          (interpolator,data(counter),&
          data,interpolator%hierarchy(counter)%level,factor,counter,0)
  end do
end subroutine hierarchical


!  Hierarchization of the higher order sparse grid hierarchical surplus (data on output) (along all dimensions)
subroutine hierarchical_order(interpolator,data,order)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(inout)   :: data
  sll_int32, intent(in) :: order
  sll_int32             :: counter,start_level,order_level_size,j
  sll_int32, dimension(:), allocatable :: k
  sll_real64 :: factor

  start_level = order-1
  order_level_size = 2**(order-1)
  SLL_ALLOCATE(k(interpolator%dim),j);

  do counter=interpolator%size_basis,1,-1
     factor = 1.0_f64
     
     do j=1,interpolator%dim
        k(j) = modulo(interpolator%hierarchy(counter)%function_type(j)-&
             interpolator%hs_weights_index(order), &
             order_level_size)+interpolator%hs_weights_index(order);
     end do
     call dehierarchical_order_d_dimension&
          (interpolator,data(counter),&
          data,start_level,&
          interpolator%hierarchy(counter)%level,&
          k,&
          factor,counter,0)
  end do

end subroutine hierarchical_order



! Dehierarchization of the sparse grid linear hierarchical surplus (data) (along all dimensions)
subroutine dehierarchical(interpolator,data)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(inout)   :: data
  !sll_real64, dimension(:), intent(out)  :: data_out
  sll_int32                                :: counter
  sll_real64 :: factor
  
 
  do counter=1,interpolator%size_basis
     !write(16,*), 'Now:', l
     !data_out(counter) = 0.0_f64!data_array(counter,1)
     factor = -1.0_f64
     call dehierarchical_d_dimension&
          (interpolator,data(counter),&
          data,interpolator%hierarchy(counter)%level,factor,counter,0)
  end do
end subroutine dehierarchical

!  Dehierarchization of the higher order sparse grid hierarchical surplus (data) (along all dimensions)
subroutine dehierarchical_order(interpolator,data,order)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(inout)   :: data
  sll_int32, intent(in) :: order
  sll_int32             :: counter,start_level,order_level_size,j
  sll_int32, dimension(:), allocatable :: k
  sll_real64 :: factor

  start_level = order-1
  order_level_size = 2**(order-1)
  SLL_ALLOCATE(k(interpolator%dim),j);

  do counter=1, interpolator%size_basis
     factor = -1.0_f64
     
     do j=1,interpolator%dim
        k(j) = modulo(interpolator%hierarchy(counter)%function_type(j)-&
             interpolator%hs_weights_index(order), &
             order_level_size)+interpolator%hs_weights_index(order);
     end do
     call dehierarchical_order_d_dimension&
          (interpolator,data(counter),&
          data,start_level,&
          interpolator%hierarchy(counter)%level,&
          k,&
          factor,counter,0)
  end do

end subroutine dehierarchical_order

! Recursive worker function for dehierarchical
recursive subroutine dehierarchical_d_dimension(interpolator,surplus,data_array,level,factor,index,d)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64, intent(inout) :: surplus
  sll_real64, dimension(:), intent(in) :: data_array
  sll_int32, dimension(:), intent(in) :: level
  sll_real64, intent(in) :: factor
  sll_int32, intent(in) :: index
  sll_int32, intent(in) :: d
  sll_int32 :: j

 if (index .NE. -1) then
   
    if (d>0) then
       surplus = surplus + data_array(index)*factor
    end if
     do j=d+1,interpolator%dim
        if (level(j)>0) then
           call dehierarchical_d_dimension(interpolator,surplus,data_array,&
                level,&
                -0.5_f64*factor,interpolator%hierarchy(index)%parent(2*j-1),j)
           call dehierarchical_d_dimension(interpolator,surplus,data_array,&
                level,&
                -0.5_f64*factor,interpolator%hierarchy(index)%parent(2*j),j)
        end if
     end do
  end if

end subroutine dehierarchical_d_dimension

! Recursive worker function for dehierarchical order
recursive subroutine dehierarchical_order_d_dimension(interpolator,surplus,data_array,start_level,level,k,factor,index,d)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64, intent(inout) :: surplus
  sll_real64, dimension(:), intent(in) :: data_array
  sll_int32, intent(in) :: start_level
  !sll_real64,dimension(:), intent(in) :: weights
  sll_int32, dimension(:), intent(in) :: level,k
  sll_real64, intent(in) :: factor
  sll_int32, intent(in) :: index
  sll_int32, intent(in) :: d
  sll_int32 :: j,father

  if (d>0) then
     surplus = surplus + data_array(index)*factor
  end if
  do j=d+1,interpolator%dim
     if (level(j)>start_level) then
        father = max(interpolator%hierarchy(index)%parent(2*j-1),&
                           interpolator%hierarchy(index)%parent(2*j))
        if (father .NE. -1) then
           call dehierarchical_order_d_dimension(&
                interpolator,surplus,data_array,start_level,level,k, &
                interpolator%hs_weights(k(j))*factor,father,j)
        end if
     end if
  end do

end subroutine dehierarchical_order_d_dimension


!!! (De)hierarchization functions on a 1D stripe

subroutine hierarchical_stripe(sparsegrid,data,max_level)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: sparsegrid
  sll_real64, dimension(:), intent(inout) :: data
  sll_int32, intent(in) :: max_level
  sll_int32 :: index, stride, index_run,j, upper,od,level,weights_index,weights_number,factor

  stride = 1;
  index = 2;
  upper = 2**max_level;
  do level = max_level,1,-1
     index_run = index
     do j=1,2**(level-1)
        data(index_run) = data(index_run) - &
             data(index_run-stride)*0.5_f64-&
             data(modulo(index_run+stride-1,upper)+1)*0.5_f64
        index_run = index_run + 2*stride
     end do
     index = index+stride
     stride = stride*2
  end do

 do od=2,sparsegrid%order
     weights_index = sparsegrid%hs_weights_index(od);
     weights_number = sparsegrid%hs_weights_index(od+1)-&
          sparsegrid%hs_weights_index(od)
     stride = 1;
     index = 2;
     do level = max_level,od,-1
        index_run = index
        factor = 1
        do j=1,2**(level-1)
           data(index_run) = data(index_run) + &
                data(index_run+factor*stride)*&
                sparsegrid%hs_weights(weights_index+modulo(j-1,weights_number))
           index_run = index_run + 2*stride
           factor = -factor
        end do
        index = index+stride
        stride = stride*2
     end do
  end do


end subroutine hierarchical_stripe


subroutine dehierarchical_stripe(sparsegrid,data, max_level)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: sparsegrid
  sll_real64, dimension(:), intent(inout) :: data
  sll_int32, intent(in) :: max_level
  sll_int32 :: index_stripe, stride,od,weights_index,weights_number,index,factor,level,j,index_run

  do od=sparsegrid%order,2,-1
     weights_index = sparsegrid%hs_weights_index(od);
     weights_number = sparsegrid%hs_weights_index(od+1)-sparsegrid%hs_weights_index(od)
     stride = 2**(max_level-od);
     index = stride+1;
     do level = od,max_level
        index_run = index
        factor = 1
        do j=1,2**(level-1)
           data(index_run) = data(index_run) - &
                data(index_run+factor*stride)*&
                sparsegrid%hs_weights(weights_index+modulo(j-1,weights_number))
           index_run = index_run + 2*stride
           factor = -factor
        end do
        stride = stride/2
        index = index-stride
     end do
  end do


  if(max_level>0) then
     index_stripe = 2**(max_level-1)
     
     data(index_stripe+1) = data(index_stripe+1)+data(1)
     if(max_level>1) then
        index_stripe = index_stripe/2
        
        call dehierarchical_stripe_recursive(index_stripe+1,index_stripe,index_stripe*4,data)
        
        call  dehierarchical_stripe_recursive(index_stripe*3+1,index_stripe,index_stripe*4,data)
     end if
  end if
end subroutine dehierarchical_stripe

recursive subroutine dehierarchical_stripe_order_recursive(index,stride,data_out)
  sll_int32, intent(in) :: index,stride
  sll_real64,dimension(:),intent(inout) :: data_out

  data_out(index-stride) = data_out(index-stride) + 0.25_f64*data_out(index)
  data_out(index+stride) = data_out(index+stride) + 0.25_f64*data_out(index)
  
  if(stride>1) then
     call dehierarchical_stripe_order_recursive(index-stride,stride/2,data_out)
     call dehierarchical_stripe_order_recursive(index+stride,stride/2,data_out)
  end if

end subroutine dehierarchical_stripe_order_recursive

recursive subroutine dehierarchical_stripe_recursive(index,stride,upper,data_out)
  sll_int32, intent(in) :: index,stride,upper
  sll_real64,dimension(:),intent(inout) :: data_out

  data_out(index) = data_out(index) + 0.5_f64*data_out(index-stride)&
       + 0.5_f64*data_out(modulo(index+stride-1,upper)+1);

  if (stride>1) then 
     call dehierarchical_stripe_recursive(index-stride/2,stride/2,upper,&
          data_out)
     call dehierarchical_stripe_recursive(index+stride/2,stride/2,upper,&
          data_out)
  end if
end subroutine dehierarchical_stripe_recursive


!!!! (De)hierarchization functions involving only parts of the dimensions

! Recursive worker function for dehierarchical_part
recursive subroutine dehierarchical_part_d_dimension(interpolator,surplus,data_array,level,factor,index,dmax,dmin,dim,dorder)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64, intent(inout) :: surplus
  sll_real64, dimension(:), intent(in) :: data_array
  sll_int32, dimension(:), intent(in) :: level
  sll_real64, intent(in) :: factor
  sll_int32, intent(in) :: index
  sll_int32, intent(in) :: dmax,dmin,dim
  sll_int32, dimension(:), intent(in) :: dorder
  sll_int32 :: j,jj

  if (index .NE. -1) then
     if (dim>dmin-1) then
        surplus = surplus + data_array(index)*factor
     end if
     do j=dim+1,dmax
        jj = dorder(j)
        if (level(jj)>0) then
           call dehierarchical_part_d_dimension(interpolator,surplus,data_array,&
                level,&
                -0.5_f64*factor,interpolator%hierarchy(index)%parent(2*jj-1),&
                dmax,dmin,j,dorder)
           call dehierarchical_part_d_dimension(interpolator,surplus,data_array,&
                level,&
                -0.5_f64*factor,interpolator%hierarchy(index)%parent(2*jj),&
                dmax,dmin,j,dorder)
        end if
     end do
  end if

end subroutine dehierarchical_part_d_dimension


! Recursive worker function for dehierarchical_part_order
recursive subroutine dehierarchical_part_order_d_dimension(interpolator,surplus,data_array,start_level,level,k,factor,index,dmax,dmin,d,dorder)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64, intent(inout) :: surplus
  sll_real64, dimension(:), intent(in) :: data_array
  sll_int32, intent(in) :: start_level
  sll_int32, dimension(:), intent(in) :: level,k, dorder
  sll_real64, intent(in) :: factor
  sll_int32, intent(in) :: index
  sll_int32, intent(in) :: dmax, dmin, d
  sll_int32 :: j,father,jj

  if (d>dmin-1) then
     surplus = surplus + data_array(index)*factor
  end if
  do j=d+1,dmax
     jj = dorder(j)
     if (level(jj)>start_level) then
        father = max(interpolator%hierarchy(index)%parent(2*jj-1),&
             interpolator%hierarchy(index)%parent(2*jj))
        if (father .NE. -1) then
           call dehierarchical_part_order_d_dimension(&
                interpolator,surplus,data_array,start_level,level,k, &
                interpolator%hs_weights(k(jj))*factor,father,dmax,dmin,j,dorder)
        end if
     end if
  end do

end subroutine dehierarchical_part_order_d_dimension

! Dehierarchization from linear hierarchical surplus only applied to parts of the dimensions. 
! Dimensions dehiearchized: dorder(dmin:dmax)
subroutine dehierarchical_part(interpolator,data,dmax,dmin,dorder)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(inout)   :: data
  sll_int32, intent(in) :: dmax,dmin
  sll_int32, dimension(:), intent(in) :: dorder
  !sll_real64, dimension(:), intent(out)  :: data_out
  sll_int32                                :: counter
  sll_real64 :: factor

  do counter=1,interpolator%size_basis
     factor = -1.0_f64
     call dehierarchical_part_d_dimension&
          (interpolator,data(counter),&
          data,interpolator%hierarchy(counter)%level,&
          factor,counter,dmax,dmin,dmin-1,dorder)
  end do
end subroutine dehierarchical_part


! Dehierarchization from higher order hierarchical surplus only applied to parts of the dimensions. 
! Dimensions dehiearchized: dorder(dmin:dmax)
subroutine dehierarchical_part_order(interpolator,data,dmax,dmin,dorder,order)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(inout)   :: data
  sll_int32, intent(in) :: dmax,dmin, order
  sll_int32, dimension(:), intent(in) :: dorder
  !sll_real64, dimension(:), intent(out)  :: data_out
  sll_int32                                :: counter,j,start_level,order_level_size
  sll_int32, dimension(:), allocatable :: k
  sll_real64 :: factor

  SLL_ALLOCATE(k(interpolator%dim),j);  

  start_level = order-1
  order_level_size = 2**start_level

  do counter=1, interpolator%size_basis
     factor = -1.0_f64
     
     do j=1,interpolator%dim
        k(j) = modulo(interpolator%hierarchy(counter)%function_type(j)-&
             interpolator%hs_weights_index(order), &
             order_level_size)+interpolator%hs_weights_index(order);
     end do
     call dehierarchical_part_order_d_dimension&
          (interpolator,data(counter),&
          data,start_level,interpolator%hierarchy(counter)%level,k,&
          factor,counter,dmax,dmin,dmin-1,dorder)
  end do

end subroutine dehierarchical_part_order

! Hierarchization to linear hierarchical surplus only applied to parts of the dimensions. 
! Dimensions hiearchized: dorder(dmin:dmax)
subroutine hierarchical_part(interpolator,data,dmax,dmin,dorder)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(inout)   :: data
  sll_int32, intent(in) :: dmax,dmin
  sll_int32, dimension(:), intent(in) :: dorder
  !sll_real64, dimension(:), intent(out)  :: data_out
  sll_int32                                :: counter
  sll_real64 :: factor
  
  do counter=interpolator%size_basis,1,-1
     factor = 1.0_f64
     call dehierarchical_part_d_dimension&
          (interpolator,data(counter),&
          data,interpolator%hierarchy(counter)%level,&
          factor,counter,dmax,dmin,dmin-1,dorder)
  end do

end subroutine hierarchical_part

! Hierarchization to higher order hierarchical surplus only applied to parts of the dimensions. 
! Dimensions hiearchized: dorder(dmin:dmax)
subroutine hierarchical_part_order(interpolator,data,dmax,dmin,dorder,order)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(inout)   :: data
  sll_int32, intent(in) :: dmax,dmin,order
  sll_int32, dimension(:), intent(in) :: dorder
  !sll_real64, dimension(:), intent(out)  :: data_out
  sll_int32                                :: counter,j
  sll_int32 :: start_level, order_level_size
  sll_int32, dimension(:), allocatable :: k
  sll_real64 :: factor

  SLL_ALLOCATE(k(interpolator%dim),j);

  start_level = order-1
  order_level_size = 2**start_level

  
  do counter = interpolator%size_basis,1, -1
     do j=1,interpolator%dim
        k(j) = modulo(interpolator%hierarchy(counter)%function_type(j)-&
             interpolator%hs_weights_index(order), &
             order_level_size)+interpolator%hs_weights_index(order);
     end do
     factor = 1.0_f64
     call dehierarchical_part_order_d_dimension&
          (interpolator,data(counter),&
          data,start_level,&
          interpolator%hierarchy(counter)%level,k,factor,counter,dmax,dmin,&
          dmin-1,dorder)
  end do
end subroutine hierarchical_part_order

!!!! End general sparse grid helper functions !!!!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!!!! Trapezoidal integrator on sparse grid
subroutine integrate_trapezoidal( interpolator,data_in,val)
  class(sll_t_sparse_grid_interpolator), intent(in) :: interpolator
  sll_int32 :: i,j
  sll_real64,intent(inout) :: val
  sll_real64, dimension(:), intent(in)   :: data_in
  sll_real64 :: phix1

  val = 0.0_f64
  do i=0, interpolator%max_level
     if (interpolator%boundary == 0) then
        phix1 = 1.0_f64/(real((2**(i)),f64))
     else
        phix1 = 1.0_f64/(real(max(2**(i),2),f64))
     end if
     do j=interpolator%level_mapping(i),interpolator%level_mapping(i+1)-1
        val = val + data_in(j)*phix1
     end do
  end do
  val = val*interpolator%volume
end subroutine integrate_trapezoidal

!!!! Trapezoidal integrator on sparse grid working with the semi-hierarchical surplus (hierarchical in (d-1) dimensions, nodal along one dimension). Note this functions should usually not be used since it is only working if the maximum number of levels along the nodal dimension is the maximum total levels
subroutine integrate_trapezoidal2( interpolator,dorder,data_in,val)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64,intent(inout) :: val
  sll_real64, dimension(:), intent(inout)   :: data_in
  sll_int32, dimension(:), intent(in) :: dorder

  call hierarchical_part(interpolator,data_in,interpolator%dim,2,dorder)
  val = sum(data_in)*interpolator%volume/(real((2**(interpolator%levels(dorder(1)))),f64));
  call dehierarchical_part(interpolator,data_in,interpolator%dim,2,dorder);
  
end subroutine integrate_trapezoidal2



!------------------------------------------------------------------------------!
!!!!!! Various SGFFT helper functions. Need a clean up. !!!!

subroutine extract_real_to_comp(sparsegrid,dim,max_level,index,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(in) :: sparsegrid
  sll_int32, intent(in) :: dim,max_level,index
  sll_int32             :: n_points,index_running
  sll_real64,dimension(:),intent(in) :: data_in
  sll_comp64,dimension(:),intent(out) :: data_out

  
  n_points = 2**(max_level)
  data_out(1) = cmplx(data_in(index), 0.0_f64, kind=f64)
  if (max_level>0) then
     index_running = sparsegrid%hierarchy(index)%children(dim*2)
     
     data_out(n_points/2+1) = cmplx(data_in(index_running), 0.0_f64, kind=f64)
     if (max_level>1) then
        
        call extract_recursive_real_to_comp(sparsegrid,n_points/2+1,n_points/4,index_running,&
             dim,data_in,data_out)
     end if
  end if


end subroutine extract_real_to_comp


recursive subroutine extract_recursive_real_to_comp(sparsegrid,index_stripe,stride,index_sg,dim,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(in) :: sparsegrid
  sll_int32, intent(in) :: index_stripe,stride,index_sg,dim
  sll_real64,dimension(:),intent(in) :: data_in
  sll_comp64,dimension(:),intent(inout) :: data_out

  data_out(index_stripe-stride) = &
       cmplx(data_in(sparsegrid%hierarchy(index_sg)%children(dim*2-1)), 0.0_f64, kind=f64)
  data_out(index_stripe+stride) = &
       cmplx(data_in(sparsegrid%hierarchy(index_sg)%children(dim*2)), 0.0_f64, kind=f64)
  if (stride>1) then
     call extract_recursive_real_to_comp(sparsegrid,index_stripe-stride,stride/2,&
          sparsegrid%hierarchy(index_sg)%children(dim*2-1),dim,&
          data_in,data_out)
     call extract_recursive_real_to_comp(sparsegrid,index_stripe+stride,stride/2,&
          sparsegrid%hierarchy(index_sg)%children(dim*2),dim,&
          data_in,data_out)
  end if
end subroutine extract_recursive_real_to_comp



subroutine extract_comp(sparsegrid,dim,max_level,index,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(in) :: sparsegrid
  sll_int32, intent(in) :: dim,max_level,index
  sll_int32             :: n_points,index_running
  sll_comp64,dimension(:),intent(in) :: data_in
  sll_comp64,dimension(:),intent(out) :: data_out

  
  n_points = 2**(max_level)
  data_out(1) = data_in(index)
  if (max_level>0) then
     index_running = sparsegrid%hierarchy(index)%children(dim*2)
     
     data_out(n_points/2+1) = data_in(index_running)
     if (max_level>1) then
        
        call extract_recursive_comp(sparsegrid,n_points/2+1,n_points/4,index_running,&
             dim,data_in,data_out)
     end if
  end if


end subroutine extract_comp


recursive subroutine extract_recursive_comp(sparsegrid,index_stripe,stride,index_sg,dim,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(in) :: sparsegrid
  sll_int32, intent(in) :: index_stripe,stride,index_sg,dim
  sll_comp64,dimension(:),intent(in) :: data_in
  sll_comp64,dimension(:),intent(inout) :: data_out

  data_out(index_stripe-stride) = &
       data_in(sparsegrid%hierarchy(index_sg)%children(dim*2-1))
  data_out(index_stripe+stride) = &
       data_in(sparsegrid%hierarchy(index_sg)%children(dim*2))
  if (stride>1) then
     call extract_recursive_comp(sparsegrid,index_stripe-stride,stride/2,&
          sparsegrid%hierarchy(index_sg)%children(dim*2-1),dim,&
          data_in,data_out)
     call extract_recursive_comp(sparsegrid,index_stripe+stride,stride/2,&
          sparsegrid%hierarchy(index_sg)%children(dim*2),dim,&
          data_in,data_out)
  end if
end subroutine extract_recursive_comp


subroutine extract_fourier(sparsegrid,dim,max_level,index,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(in) :: sparsegrid
  sll_int32, intent(in) :: dim,max_level,index
  sll_int32             :: n_points,index_running
  sll_comp64,dimension(:),intent(in) :: data_in
  sll_comp64,dimension(:),intent(out) :: data_out

  
  n_points = 2**(max_level)
  data_out(1) = data_in(index)
  if (max_level>0) then
     index_running = sparsegrid%hierarchy(index)%children(dim*2)   
     data_out(2) = data_in(index_running)
     if (max_level>1) then
        call extract_recursive_fourier(sparsegrid,index_running,0,&
             2,max_level,dim,data_in,data_out)
     end if
  end if
end subroutine extract_fourier

recursive subroutine extract_recursive_fourier(sparsegrid,index_sg,ind,level,max_level,dim,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(in) :: sparsegrid
  sll_int32, intent(in) :: level,max_level,index_sg,dim,ind
  sll_comp64,dimension(:),intent(in) :: data_in
  sll_comp64,dimension(:),intent(inout) :: data_out

  data_out(2**(level-1)+1+2*ind) = &
       data_in(sparsegrid%hierarchy(index_sg)%children(dim*2-1))
  data_out(2**(level-1)+1+2*ind+1) = &
       data_in(sparsegrid%hierarchy(index_sg)%children(dim*2))
  if (level<max_level) then
     call extract_recursive_fourier(sparsegrid,&
          sparsegrid%hierarchy(index_sg)%children(dim*2-1),2*ind,&
          level+1,max_level,dim,data_in,data_out)
     call extract_recursive_fourier(sparsegrid,&
          sparsegrid%hierarchy(index_sg)%children(dim*2),2*ind+1,&
          level+1,max_level,dim,data_in,data_out)
  end if
end subroutine extract_recursive_fourier

subroutine insert_fourier(sparsegrid,dim,max_level,index,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(in) :: sparsegrid
  sll_int32, intent(in) :: dim,max_level,index
  sll_int32             :: n_points,index_running
  sll_comp64,dimension(:),intent(in) :: data_in
  sll_comp64,dimension(:),intent(out) :: data_out

  
  n_points = 2**(max_level)
  data_out(index) = data_in(1)
  if (max_level>0) then
     index_running = sparsegrid%hierarchy(index)%children(dim*2)   
     data_out(index_running) = data_in(2)
     if (max_level>1) then
        call insert_recursive_fourier(sparsegrid,index_running,0,&
             2,max_level,dim,data_in,data_out)
     end if
  end if
end subroutine insert_fourier

recursive subroutine insert_recursive_fourier(sparsegrid,index_sg,ind,level,max_level,dim,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(in) :: sparsegrid
  sll_int32, intent(in) :: level,max_level,index_sg,dim,ind
  sll_comp64,dimension(:),intent(in) :: data_in
  sll_comp64,dimension(:),intent(inout) :: data_out

  data_out(sparsegrid%hierarchy(index_sg)%children(dim*2-1)) = &
       data_in(2**(level-1)+1+2*ind)
  data_out(sparsegrid%hierarchy(index_sg)%children(dim*2)) = &
       data_in(2**(level-1)+1+2*ind+1)
  if (level<max_level) then
     call insert_recursive_fourier(sparsegrid,&
          sparsegrid%hierarchy(index_sg)%children(dim*2-1),2*ind,&
          level+1,max_level,dim,data_in,data_out)
     call insert_recursive_fourier(sparsegrid,&
          sparsegrid%hierarchy(index_sg)%children(dim*2),2*ind+1,&
          level+1,max_level,dim,data_in,data_out)
  end if
end subroutine insert_recursive_fourier

!PN DEFINED BUT NOT USED
!subroutine sort_hiera_1d(max_level,data_in,data_out)
!  sll_int32, intent(in) :: max_level
!  sll_real64,dimension(:),intent(in) :: data_in
!  sll_real64,dimension(:),intent(inout) :: data_out
!  sll_int32 :: size,l,j,initial,stride,counter
!
!  size = 2**max_level;
!
!  counter = 1
!  data_out(1) = data_in(1);
!  if(max_level>1) then
!     data_out(2) = data_in(size/2+1)
!     initial = size/4
!     stride = size/2
!     do l=2,max_level
!        do j=0,2**(l-1)
!           data_out(counter) = data_in(initial+1+j*stride)
!           counter = counter +1
!        end do
!        initial = initial/2
!        stride = stride/2
!     end do
!  end if
!
!end subroutine sort_hiera_1d

subroutine insert_comp_to_real(sparsegrid,dim,max_level,index,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(in) :: sparsegrid
  sll_int32, intent(in) :: dim,max_level,index
  sll_int32             :: n_points,index_running
  sll_comp64,dimension(:),intent(in) :: data_in
  sll_real64,dimension(:),intent(out) :: data_out

  
  n_points = 2**(max_level)
  data_out(index) = real(data_in(1))
  if (max_level>0) then
     index_running = sparsegrid%hierarchy(index)%children(dim*2)
     
     data_out(index_running) = real(data_in(n_points/2+1))
     if (max_level>1) then
        
        call insert_recursive_comp_to_real(sparsegrid,n_points/2+1,n_points/4,index_running,&
             dim,data_in,data_out)
     end if
  end if


end subroutine insert_comp_to_real

recursive subroutine insert_recursive_comp_to_real(sparsegrid,index_stripe,stride,index_sg,dim,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(in) :: sparsegrid
  sll_int32, intent(in) :: index_stripe,stride,index_sg,dim
  sll_comp64,dimension(:),intent(in) :: data_in
  sll_real64,dimension(:),intent(inout) :: data_out

  data_out(sparsegrid%hierarchy(index_sg)%children(dim*2-1)) =&
       real(data_in(index_stripe-stride)) 
  data_out(sparsegrid%hierarchy(index_sg)%children(dim*2)) =&
       real(data_in(index_stripe+stride))
  if (stride>1) then
     call insert_recursive_comp_to_real(sparsegrid,index_stripe-stride,stride/2,&
          sparsegrid%hierarchy(index_sg)%children(dim*2-1),dim,&
          data_in,data_out)
     call insert_recursive_comp_to_real(sparsegrid,index_stripe+stride,stride/2,&
          sparsegrid%hierarchy(index_sg)%children(dim*2),dim,&
          data_in,data_out)
  end if
end subroutine insert_recursive_comp_to_real


subroutine insert_comp(sparsegrid,dim,max_level,index,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(in) :: sparsegrid
  sll_int32, intent(in) :: dim,max_level,index
  sll_int32             :: n_points,index_running
  sll_comp64,dimension(:),intent(in) :: data_in
  sll_comp64,dimension(:),intent(out) :: data_out

  
  n_points = 2**(max_level)
  data_out(index) = data_in(1)
  if (max_level>0) then
     index_running = sparsegrid%hierarchy(index)%children(dim*2)
     
     data_out(index_running) = data_in(n_points/2+1)
     if (max_level>1) then
        
        call insert_recursive_comp(sparsegrid,n_points/2+1,n_points/4,index_running,&
             dim,data_in,data_out)
     end if
  end if


end subroutine insert_comp

recursive subroutine insert_recursive_comp(sparsegrid,index_stripe,stride,index_sg,dim,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(in) :: sparsegrid
  sll_int32, intent(in) :: index_stripe,stride,index_sg,dim
  sll_comp64,dimension(:),intent(in) :: data_in
  sll_comp64,dimension(:),intent(inout) :: data_out

  data_out(sparsegrid%hierarchy(index_sg)%children(dim*2-1)) =&
       data_in(index_stripe-stride) 
  data_out(sparsegrid%hierarchy(index_sg)%children(dim*2)) =&
       data_in(index_stripe+stride)
  if (stride>1) then
     call insert_recursive_comp(sparsegrid,index_stripe-stride,stride/2,&
          sparsegrid%hierarchy(index_sg)%children(dim*2-1),dim,&
          data_in,data_out)
     call insert_recursive_comp(sparsegrid,index_stripe+stride,stride/2,&
          sparsegrid%hierarchy(index_sg)%children(dim*2),dim,&
          data_in,data_out)
  end if
end subroutine insert_recursive_comp



subroutine hira(data,max_level)
  sll_comp64, dimension(:), intent(inout) :: data
  sll_int32, intent(in) :: max_level
  sll_int32 ::l,i

  do l=max_level,1,-1
     do i=2**l,2**(l-1)+1,-1
        data(2**l+1-i) =  data(2**l+1-i)+data(i)
     end do
  end do

end subroutine hira

subroutine dehi(data,max_level)
  sll_comp64, dimension(:), intent(inout) :: data
  sll_int32, intent(in) :: max_level
  sll_int32 ::l,i

  do l = 1,max_level
     do i = 2**(l-1)+1,2**l
        data(2**l+1-i) =  data(2**l+1-i)-data(i)
     end do
  end do

end subroutine dehi

subroutine fft_on_stripe(interpolator,level) 
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_int32, intent(in) ::level
  sll_int32 :: i,no_points

  no_points = 2**level

  call sll_s_fft_exec_c2c_1d( interpolator%fft_object(level+1)%fw,&
       interpolator%fft_object(level+1)%in,interpolator%fft_object(level+1)%out )
  ! Scale the Fourier coefficients
  do i=1,no_points
     interpolator%fft_object(level+1)%out(i) = &
          interpolator%fft_object(level+1)%out(i)/no_points
  end do

end subroutine fft_on_stripe

subroutine ifft_on_stripe(interpolator,level)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_int32, intent(in) :: level
  
  call sll_s_fft_exec_c2c_1d(interpolator%fft_object(level+1)%bw,&
       interpolator%fft_object(level+1)%out,interpolator%fft_object(level+1)%in)

end subroutine ifft_on_stripe


subroutine fft_initialize(fft_object,levels)
  class(fft_hierarchical), dimension(:),pointer :: fft_object
  sll_int32,intent(in) :: levels 
  sll_int32 :: l,size
#ifndef FFTW_F2003
  sll_int32 :: ierr
#endif

  size = 1
  do l=1,levels+1
     
     fft_object(l)%sz_in = int(size,C_SIZE_T)
     fft_object(l)%sz_out = int(size,C_SIZE_T)
     fft_object(l)%out => sll_f_fft_allocate_aligned_complex(size)
     fft_object(l)%in => sll_f_fft_allocate_aligned_complex(size)
     call sll_s_fft_init_c2c_1d( fft_object(l)%fw, size, fft_object(l)%in, &
          fft_object(l)%out, sll_p_fft_backward, normalized = .false., &
          aligned = .true., optimization = sll_p_fft_measure )
     
     call sll_s_fft_init_c2c_1d( fft_object(l)%bw, size, fft_object(l)%out, &
          fft_object(l)%in, sll_p_fft_forward, normalized = .false., &
          aligned = .true., optimization = sll_p_fft_measure )
     
     size = size*2
  end do
     

end subroutine fft_initialize



subroutine fft_finalize(fft_object,levels)
  class(fft_hierarchical), dimension(:),pointer :: fft_object
  sll_int32, intent(in) :: levels
  sll_int32 :: l
 
  do l=1,levels+1
     call sll_s_fft_free(fft_object(l)%fw)
     call sll_s_fft_free(fft_object(l)%bw)
  end do

end subroutine fft_finalize

!> Finalize sparse grid
subroutine free_sparse_grid(interpolator)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator

  call fft_finalize(interpolator%fft_object, interpolator%max_level)
  deallocate(interpolator%hierarchy)
  deallocate(interpolator%stripe)
  deallocate(interpolator%stripe_out)
  deallocate(interpolator%hs_weights)
  deallocate(interpolator%hs_weights_index)
  deallocate(interpolator%interp)
  deallocate(interpolator%level_mapping)
  deallocate(interpolator%eta_min)
  deallocate(interpolator%eta_max)
  deallocate(interpolator%length)

end subroutine free_sparse_grid


!> Compute Fourier coefficients on sparse grid along dimension \a dim. 
subroutine ToHierarchical1D(interpolator,dim,max_level,index,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_real64, dimension(:), intent(in) :: data_in
  sll_comp64, dimension(:), intent(out) :: data_out
  sll_int32, intent(in) :: dim,max_level,index

  call extract_real_to_comp(interpolator,dim,max_level,&
       index,data_in,interpolator%fft_object(max_level+1)%in)
  call fft_on_stripe(interpolator,max_level)

  call fft_to_centered(2**max_level,interpolator%fft_object(max_level+1)%out,&
       interpolator%fft_object(max_level+1)%in)

  call hira(interpolator%fft_object(max_level+1)%in,max_level)

  call insert_fourier(interpolator,dim,max_level,&
       index,interpolator%fft_object(max_level+1)%in,data_out)
end subroutine ToHierarchical1D

!> Complex version of \a ToHierarchical1d_comp
subroutine ToHierarchical1D_comp(interpolator,dim,max_level,index,data)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_comp64, dimension(:), intent(inout) :: data
  sll_int32, intent(in) :: dim,max_level,index

  call extract_comp(interpolator,dim,max_level,&
       index,data,interpolator%fft_object(max_level+1)%in)
  call fft_on_stripe(interpolator,max_level)

  call fft_to_centered(2**max_level,interpolator%fft_object(max_level+1)%out,&
       interpolator%fft_object(max_level+1)%in)
  call hira(interpolator%fft_object(max_level+1)%in,max_level)
  call insert_fourier(interpolator,dim,max_level,&
       index,interpolator%fft_object(max_level+1)%in,data)
end subroutine ToHierarchical1D_comp

!>
subroutine ToHira1D(interpolator,dim,max_level,index,data)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_comp64, dimension(:), intent(inout) :: data
  sll_int32, intent(in) :: dim,max_level,index

  call extract_fourier(interpolator,dim,max_level,index,&
       data,interpolator%fft_object(max_level+1)%in)
  call hira(interpolator%fft_object(max_level+1)%in,max_level)
  call insert_fourier(interpolator,dim,max_level,index,&
       interpolator%fft_object(max_level+1)%in,data)

end subroutine ToHira1D

!>
subroutine ToNodal1D(interpolator,dim,max_level,index,data_in,data_out)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_comp64, dimension(:), intent(in) :: data_in
  sll_real64, dimension(:), intent(out) :: data_out
  sll_int32, intent(in) :: dim,max_level,index

  call extract_fourier(interpolator,dim,&
       max_level,index,data_in,interpolator%fft_object(max_level+1)%in) 
  call dehi(interpolator%fft_object(max_level+1)%in,max_level)
  call fft_to_inorder(2**max_level,interpolator%fft_object(max_level+1)%in,&
       interpolator%fft_object(max_level+1)%out)
  !call dehi(interpolator%fft_object(max_level+1)%out,max_level)
  call ifft_on_stripe(interpolator,max_level)
  call insert_comp_to_real(interpolator,dim,max_level,index,&
       interpolator%fft_object(max_level+1)%in,data_out)

end subroutine ToNodal1D

!>
subroutine ToNodal1D_comp(interpolator,dim,max_level,index,data)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_comp64, dimension(:), intent(inout) :: data
  sll_int32, intent(in) :: dim,max_level,index

  call extract_fourier(interpolator,dim,&
       max_level,index,data,interpolator%fft_object(max_level+1)%in)  
  call dehi(interpolator%fft_object(max_level+1)%in,max_level)
  !print*, 'dehi', dim
  !print*, interpolator%fft_object(max_level+1)%in
  call fft_to_inorder(2**max_level,interpolator%fft_object(max_level+1)%in,&
       interpolator%fft_object(max_level+1)%out)
  !call dehi(interpolator%fft_object(max_level+1)%out,max_level)
  call ifft_on_stripe(interpolator,max_level)
  call insert_comp(interpolator,dim,max_level,index,&
       interpolator%fft_object(max_level+1)%in,data)

end subroutine ToNodal1D_comp

subroutine ToDehi1D(interpolator,dim,max_level,index,data_array)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_comp64, dimension(:), intent(inout) :: data_array
  sll_int32, intent(in) :: dim,max_level,index

  call extract_fourier(interpolator,dim,max_level,index,&
       data_array,interpolator%fft_object(max_level+1)%out)
  call dehi(interpolator%fft_object(max_level+1)%out,max_level)
  call insert_fourier(interpolator,dim,max_level,index,&
       interpolator%fft_object(max_level+1)%out,data_array)

end subroutine ToDehi1D



subroutine FFT_TO_CENTERED(length,data_in,data_out)
  sll_int32, intent(in) :: length
  sll_comp64, dimension(:), intent(in) :: data_in
  sll_comp64, dimension(:), intent(out) :: data_out
  sll_int32 :: k

  data_out(1) = data_in(1)
  data_out(length) = data_in(length/2+1)
  
  do k=1,length/2-1
     data_out(2*k) = data_in(k+1)
     data_out(2*k+1) = data_in(length-k+1)
  end do

end subroutine FFT_TO_CENTERED

subroutine FFT_TO_INORDER(length,data_in,data_out)
  sll_int32, intent(in) :: length
  sll_comp64, dimension(:), intent(in) :: data_in
  sll_comp64, dimension(:), intent(out) :: data_out
  sll_int32 :: k

  data_out(1) = data_in(1)
  data_out(length/2+1) = data_in(length)
  
  do k=1,length/2-1
     data_out(k+1) = data_in(2*k)
     data_out(length-k+1) = data_in(2*k+1)
  end do

end subroutine FFT_TO_INORDER

subroutine displace_fourier_coeffs(d_scale,size,data)!displacement)
  sll_comp64, dimension(:), intent(inout) :: data
  !sll_real64, intent(in) :: displacement
  sll_real64, intent(in) :: d_scale
  sll_int32, intent(in) :: size
  sll_int32 :: k
  sll_real64 :: sin_val,cos_val

  do k=1,size/2
     sin_val = sin(-k*d_scale)
     cos_val = cos(-k*d_scale)
     data(2*k) = cmplx(real(data(2*k))*cos_val - aimag(data(2*k))*sin_val, &
          real(data(2*k))*sin_val+aimag(data(2*k))*cos_val, kind=f64)
     if(k<size/2) then
        sin_val = -sin_val!sin(k*d_scale)
        !cos_val = cos(k*d_scale)
        data(2*k+1) = cmplx(real(data(2*k+1))*cos_val-aimag(data(2*k+1))*sin_val, &
             real(data(2*k+1))*sin_val+aimag(data(2*k+1))*cos_val, kind=f64)
     end if
  end do

end subroutine displace_fourier_coeffs



subroutine Displace1D(interpolator,dim,max_level,index,displacement,data)
  class(sll_t_sparse_grid_interpolator), intent(inout) :: interpolator
  sll_comp64, dimension(:), intent(inout) :: data
  sll_int32, intent(in) :: dim,max_level,index
  sll_real64 ,intent(in) :: displacement

  call extract_fourier(interpolator,dim,max_level,&
       index,data,interpolator%fft_object(max_level+1)%in)
  call displace_fourier_coeffs(displacement,2**max_level,&
       interpolator%fft_object(max_level+1)%in)
  call insert_fourier(interpolator,dim,max_level,&
       index,interpolator%fft_object(max_level+1)%in,data)
end subroutine Displace1D




!!!! End SGFFT helper functions !!!!
!------------------------------------------------------------------------------!



end module sll_m_sparse_grid_interpolator
