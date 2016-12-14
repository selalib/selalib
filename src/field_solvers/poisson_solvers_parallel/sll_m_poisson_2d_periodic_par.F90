!> @ingroup poisson_solvers
!> @brief 
!> Selalib periodic 2D poisson solver for cartesian coordinates.
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!> @details
!> While the solution region is 2D, the following solver is meant to be used 
!> with 4D arrays with size 1 in the last two dimensions. The reason for this
!> is that the input data, the rho array, is the result of a reduction
!> process in the last 2 dimensions. The computation of the potential
!> which may be an in-place operation, needs to be treated as a 4D array
!> in order to remap it appropriately with the 4D distribution function.
!> There might be ways around this, like using 'reshape'...
module sll_m_poisson_2d_periodic_par

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_t_collective_t, &
    sll_f_get_collective_size

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_fft, only: &
    sll_s_fft_exec_c2c_1d, &
    sll_s_fft_exec_c2c_2d, &
    sll_p_fft_backward, &
    sll_s_fft_free, &
    sll_p_fft_forward, &
    sll_s_fft_init_c2c_1d, &
    sll_t_fft

  use sll_m_remapper, only: &
    sll_o_apply_remap_2d, &
    sll_o_compute_local_sizes, &
    sll_o_get_layout_collective, &
    sll_o_initialize_layout_with_distributed_array, &
    sll_t_layout_2d, &
    sll_o_local_to_global, &
    sll_f_new_layout_2d, &
    sll_o_new_remap_plan, &
    sll_t_remap_plan_2d_comp64, &
    sll_o_delete

  use sll_m_utilities, only: &
    sll_f_is_power_of_two

  implicit none

  public :: &
    sll_t_poisson_2d_periodic_par, &
    sll_f_poisson_2d_periodic_par_new, &
    sll_f_poisson_2d_periodic_par_new_alt, &
    sll_s_poisson_2d_periodic_par_solve, &
    sll_s_poisson_2d_periodic_par_solve_alt, &
    sll_s_poisson_2d_periodic_par_free

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Structure to store data from Poisson solver. This
  !> solver is parallel on structured cartesian mesh. Numerical method
  !> uses FFT transforms.
  type sll_t_poisson_2d_periodic_par
     sll_int32                           :: ncx    !< number of cells  
     sll_int32                           :: ncy    !< number of cells
     sll_real64                          :: Lx     !< domain length 
     sll_real64                          :: Ly     !< domain length
     type(sll_t_fft)        :: px     !< fft plan in x
     type(sll_t_fft)         :: py     !< fft plan in y
     type(sll_t_fft)         :: px_inv !< inverse fft plan in x
     type(sll_t_fft)         :: py_inv !< inverse fft plan in y
     type(sll_t_layout_2d), pointer           :: layout_seq_x1 !< layout sequential in x
     type(sll_t_layout_2d), pointer           :: layout_seq_x2 !< layout sequential in y
     sll_int32                           :: seq_x1_local_sz_x1 !< local size 
     sll_int32                           :: seq_x1_local_sz_x2 !< local size 
     sll_int32                           :: seq_x2_local_sz_x1 !< local size 
     sll_int32                           :: seq_x2_local_sz_x2 !< local size 
     sll_comp64, dimension(:,:), pointer :: fft_x_array !< 2d array for fft in x
     sll_comp64, dimension(:,:), pointer :: fft_y_array !< array for fft in y
     sll_comp64, dimension(:),   pointer :: fft_x_1d_array !< 1d array for sliced fft in x
     sll_comp64, dimension(:),   pointer :: fft_y_1d_array !< 1d array for sliced fft in y
     type(sll_t_remap_plan_2d_comp64), pointer :: rmp_xy !< remap to transpose from x to y
     type(sll_t_remap_plan_2d_comp64), pointer :: rmp_yx !< remap to transpose from y to x
  end type sll_t_poisson_2d_periodic_par

contains

  !> Presently, this function receives the geometric information as 
  !> individual arguments. We should consider passing the 'simple geometry'
  !> object that we have for the cartesian cases.
  !> @return
  function sll_f_poisson_2d_periodic_par_new( &
    start_layout, &   
    ncx, &            
    ncy, &            
    Lx, &    
    Ly ) result(plan)

    type (sll_t_poisson_2d_periodic_par), pointer :: plan
    type(sll_t_layout_2d), pointer         :: start_layout !< First layout
    sll_int32                        :: ncx          !< number of cells in x
    sll_int32                        :: ncy          !< number of cells in y
    sll_real64                       :: Lx           !< length x
    sll_real64                       :: Ly           !< length y
    sll_int64                        :: colsz ! collective size
    type(sll_t_collective_t), pointer  :: collective
    ! number of processors
    sll_int32                        :: nprocx1
    sll_int32                        :: nprocx2
    sll_int32                        :: ierr 
    sll_int32                        :: loc_sz_x1
    sll_int32                        :: loc_sz_x2
    !sll_int32                        :: seq_x1_local_sz_x1
    !sll_int32                        :: seq_x1_local_sz_x2
    !sll_int32                        :: seq_x2_local_sz_x1
    !sll_int32                        :: seq_x2_local_sz_x2

    ! The collective to be used is the one that comes with the given layout.
    collective => sll_o_get_layout_collective( start_layout )
    colsz      = int(sll_f_get_collective_size( collective ),i64)

    if ( (.not.sll_f_is_power_of_two(int(ncx,i64))) .and. &
         (.not.sll_f_is_power_of_two(int(ncy,i64))) ) then     
       print *, 'This test needs to run with numbers of cells which are',  &
                'powers of 2.'
       print *, 'Exiting...'
       stop
    end if

    SLL_ALLOCATE(plan, ierr)

    ! We use the number of cells since due to periodicity, the last point is
    ! not considered. 
    plan%ncx = ncx
    plan%ncy = ncy
    plan%Lx  = Lx
    plan%Ly  = Ly

    ! Layout and local sizes for FFTs in x-direction
    plan%layout_seq_x1 => start_layout
    call sll_o_compute_local_sizes( &
         plan%layout_seq_x1, &
         loc_sz_x1, &
         loc_sz_x2 )

    plan%seq_x1_local_sz_x1 = loc_sz_x1 
    plan%seq_x1_local_sz_x2 = loc_sz_x2

    SLL_ALLOCATE( plan%fft_x_array(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE( plan%fft_x_1d_array(loc_sz_x1),ierr)

    ! For FFTs (in x-direction)
    call sll_s_fft_init_c2c_1d( &
         plan%px, &
         ncx, &
         plan%fft_x_1d_array, &
         plan%fft_x_1d_array, &
         sll_p_fft_forward)!+FFT_NORMALIZE )

    call sll_s_fft_init_c2c_1d( &
         plan%px_inv, &
         ncx, &
         plan%fft_x_1d_array, &
         plan%fft_x_1d_array, &
         sll_p_fft_backward) !+FFT_NORMALIZE )

    ! Layout and local sizes for FFTs in y-direction (x2)
    plan%layout_seq_x2 => sll_f_new_layout_2d( collective )
    nprocx1 = int(colsz,kind=4)
    nprocx2 = 1

    call sll_o_initialize_layout_with_distributed_array( &
         ncx+1, &
         ncy+1, &
         nprocx1, &
         nprocx2, &
         plan%layout_seq_x2 )

    call sll_o_compute_local_sizes( &
         plan%layout_seq_x2, &
         loc_sz_x1, &
         loc_sz_x2 )

    plan%seq_x2_local_sz_x1 = loc_sz_x1
    plan%seq_x2_local_sz_x2 = loc_sz_x2
    SLL_ALLOCATE( plan%fft_y_array(loc_sz_x1,loc_sz_x2), ierr )
    SLL_ALLOCATE( plan%fft_y_1d_array(loc_sz_x2), ierr )

    ! For FFTs (in y-direction)

    call sll_s_fft_init_c2c_1d( &
         plan%py, &
         ncy, &
         plan%fft_y_1d_array, &
         plan%fft_y_1d_array, &
         sll_p_fft_forward)! + FFT_NORMALIZE )

    call sll_s_fft_init_c2c_1d( &
         plan%py_inv, &
         ncy, &
         plan%fft_y_1d_array, &
         plan%fft_y_1d_array, &
         sll_p_fft_backward)! + FFT_NORMALIZE )

    plan%rmp_xy => &
     sll_o_new_remap_plan(plan%layout_seq_x1, plan%layout_seq_x2, plan%fft_x_array)
    plan%rmp_yx => &
     sll_o_new_remap_plan(plan%layout_seq_x2, plan%layout_seq_x1, plan%fft_y_array)
  end function sll_f_poisson_2d_periodic_par_new


  !> Presently, this function receives the geometric information as 
  !> individual arguments. We should consider passing the 'simple geometry'
  !> object that we have for the cartesian cases. The 'alt' version does not
  !> consider the last point in the problem, the arrays involved and the layout
  !> that represents them should also not include this last point
  !> @return
  function sll_f_poisson_2d_periodic_par_new_alt( &
    start_layout, &   
    ncx, &            
    ncy, &            
    Lx, &    
    Ly ) result(plan)

    type (sll_t_poisson_2d_periodic_par), pointer :: plan
    type(sll_t_layout_2d), pointer         :: start_layout !< First layout
    sll_int32                        :: ncx          !< number of cells in x
    sll_int32                        :: ncy          !< number of cells in y
    sll_real64                       :: Lx           !< length x
    sll_real64                       :: Ly           !< length y
    sll_int64                        :: colsz ! collective size
    type(sll_t_collective_t), pointer  :: collective
    ! number of processors
    sll_int32                        :: nprocx1
    sll_int32                        :: nprocx2
    sll_int32                        :: ierr 
    sll_int32                        :: loc_sz_x1
    sll_int32                        :: loc_sz_x2
    !sll_int32                        :: seq_x1_local_sz_x1
    !sll_int32                        :: seq_x1_local_sz_x2
    !sll_int32                        :: seq_x2_local_sz_x1
    !sll_int32                        :: seq_x2_local_sz_x2

    ! The collective to be used is the one that comes with the given layout.
    collective => sll_o_get_layout_collective( start_layout )
    colsz      = int(sll_f_get_collective_size( collective ),i64)

    if ( (.not.sll_f_is_power_of_two(int(ncx,i64))) .and. &
         (.not.sll_f_is_power_of_two(int(ncy,i64))) ) then     
       print *, 'This test needs to run with numbers of cells which are',  &
                'powers of 2.'
       print *, 'Exiting...'
       stop
    end if

    SLL_ALLOCATE(plan, ierr)

    ! We use the number of cells since due to periodicity, the last point is
    ! not considered. It should not even exist...
    plan%ncx = ncx
    plan%ncy = ncy
    plan%Lx  = Lx
    plan%Ly  = Ly

    ! Layout and local sizes for FFTs in x-direction
    plan%layout_seq_x1 => start_layout
    call sll_o_compute_local_sizes( &
         plan%layout_seq_x1, &
         loc_sz_x1, &
         loc_sz_x2 )

    plan%seq_x1_local_sz_x1 = loc_sz_x1  ! ncx
    plan%seq_x1_local_sz_x2 = loc_sz_x2

    SLL_ALLOCATE( plan%fft_x_array(loc_sz_x1,loc_sz_x2),ierr)

    ! For FFTs (in x-direction)
    call sll_s_fft_init_c2c_1d( &
         plan%px, &
         ncx, &
         plan%fft_x_1d_array, &
         plan%fft_x_1d_array, &
         sll_p_fft_forward)!+FFT_NORMALIZE )

    call sll_s_fft_init_c2c_1d( &
         plan%px_inv, &
         ncx, &
         plan%fft_x_1d_array, &
         plan%fft_x_1d_array, &
         sll_p_fft_backward) !+FFT_NORMALIZE )

    ! Layout and local sizes for FFTs in y-direction (x2)
    plan%layout_seq_x2 => sll_f_new_layout_2d( collective )
    nprocx1 = int(colsz,kind=4)
    nprocx2 = 1

    call sll_o_initialize_layout_with_distributed_array( &
         ncx, &
         ncy, &
         nprocx1, &
         nprocx2, &
         plan%layout_seq_x2 )

    call sll_o_compute_local_sizes( &
         plan%layout_seq_x2, &
         loc_sz_x1, &
         loc_sz_x2 )

    plan%seq_x2_local_sz_x1 = loc_sz_x1
    plan%seq_x2_local_sz_x2 = ncy ! loc_sz_x2
    SLL_ALLOCATE( plan%fft_y_array(loc_sz_x1,loc_sz_x2), ierr )

    ! For FFTs (in y-direction)

    call sll_s_fft_init_c2c_1d( &
         plan%py, &
         ncy, &
         plan%fft_y_1d_array, &
         plan%fft_y_1d_array, &
         sll_p_fft_forward)! + FFT_NORMALIZE )

    call sll_s_fft_init_c2c_1d( &
         plan%py_inv, &
         ncy, &
         plan%fft_y_1d_array, &
         plan%fft_y_1d_array, &
         sll_p_fft_backward)! + FFT_NORMALIZE )

    plan%rmp_xy => &
     sll_o_new_remap_plan(plan%layout_seq_x1, plan%layout_seq_x2, plan%fft_x_array)
    plan%rmp_yx => &
     sll_o_new_remap_plan(plan%layout_seq_x2, plan%layout_seq_x1, plan%fft_y_array)
  end function sll_f_poisson_2d_periodic_par_new_alt



  !> Note that the equation that is solved is: \f$ \Delta \phi = \rho \f$
  !> Thus the user is responsible for giving the proper sign to the source term.
  subroutine sll_s_poisson_2d_periodic_par_solve(plan, rho, phi)
    type (sll_t_poisson_2d_periodic_par), pointer :: plan !< self object
    sll_real64, dimension(:,:)        :: rho      !< charge density
    sll_real64, dimension(:,:)        :: phi      !< electric potential
    sll_int32                         :: ncx      !< global size
    sll_int32                         :: ncy      !< global size
    sll_int32                         :: npx_loc
    sll_int32                         :: npy_loc
    sll_int32                         :: i, j
    ! Reciprocals of domain lengths.
    sll_real64                        :: r_Lx, r_Ly
    sll_real64                        :: kx, ky
    !sll_comp64                        :: val
    sll_real64                        :: normalization
    type(sll_t_layout_2d), pointer          :: layout_x
    type(sll_t_layout_2d), pointer          :: layout_y
    sll_int32, dimension(1:2)         :: global
    sll_int32                         :: gi, gj

    ncx  = plan%ncx
    ncy  = plan%ncy
    r_Lx = 1.0_f64/plan%Lx
    r_Ly = 1.0_f64/plan%Ly

    !print*, 1.0_f64/plan%Lx, 1.0_f64/plan%Ly
    ! Get layouts to compute FFTs (in each direction)
    layout_x => plan%layout_seq_x1
    layout_y => plan%layout_seq_x2
    !call verify_argument_sizes_par(layout_x, rho, phi)

    ! FFTs in x-direction
    npx_loc = plan%seq_x1_local_sz_x1 
    npy_loc = plan%seq_x1_local_sz_x2 

    ! The input is handled internally as a complex array
    plan%fft_x_array = cmplx(rho, 0.0_f64, kind=f64)

    do j=1,npy_loc
       call sll_s_fft_exec_c2c_1d(plan%px, plan%fft_x_array(:,j), plan%fft_x_array(:,j))
    end do
    ! FFTs in y-direction
    npx_loc = plan%seq_x2_local_sz_x1
    npy_loc = plan%seq_x2_local_sz_x2

    call sll_o_apply_remap_2d( plan%rmp_xy, plan%fft_x_array, plan%fft_y_array )

    do i=1,npx_loc
       call sll_s_fft_exec_c2c_1d(plan%py, plan%fft_y_array(i,:), plan%fft_y_array(i,:)) 
    end do

    ! This should be inside the FFT plan...
    normalization = 1.0_f64/(ncx*ncy)

    ! Apply the kernel 
    do j=1, npy_loc-1 ! last point was not transformed
       do i=1, npx_loc-1
          ! Make sure that the first mode is set to zero so that we get an
          ! answer with zero mean value. This step assumes that the (1,1) point
          ! will always be at the border of any splitting of the domains. This 
          ! seems safe in this case.
          
          global = sll_o_local_to_global( layout_y, (/i, j/))
          gi = global(1)
          gj = global(2)
          
          if( (gi == 1) .and. (gj == 1) ) then
             plan%fft_y_array(1,1) = (0.0_f64,0.0_f64)
          else
             kx  = real(gi-1,f64)
             ky  = real(gj-1,f64)
             ! Crucial step: While the results of the FFT are basically a
             ! list of the modes in the 0..n-1 sense, the algorithm that we
             ! are using uses a decomposition with negative frequencies,
             ! i.e.: -n/2..n/2-1. Therefore a step is needed in which we are
             ! effectively reinterpreting the FFT results in this light.
             ! More specifically, we take the second half of the transformed
             ! array and apply a modified factor, proper of a negative 
             ! frequency.
             if( kx .ge. ncx/2 ) then
                kx = kx - ncx
             end if

             if( ky .ge. ncy/2 ) then
                ky = ky - ncy
             end if

              plan%fft_y_array(i,j) = -plan%fft_y_array(i,j)*normalization / &
                  ( ( (kx*r_Lx)**2 + (ky*r_Ly)**2)*4.0_f64*sll_p_pi**2)
          end if
       enddo
    enddo

    ! Inverse FFTs in y-direction
    do i=1,npx_loc
       call sll_s_fft_exec_c2c_1d(plan%py_inv, plan%fft_y_array(i,:), plan%fft_y_array(i,:))
    end do

    ! Force the periodicity condition in the y-direction. CAN'T USE THE FFT
    ! INTERFACE SINCE THIS POINT FALLS OUTSIDE OF THE POINTS IN THE ARRAY
    ! TOUCHED BY THE FFT. This is another reason to permit not including the
    ! last point in the periodic cases...

    plan%fft_y_array(:,npy_loc) = plan%fft_y_array(:,1)

    ! Prepare to take inverse FFTs in x-direction
    call sll_o_apply_remap_2d( plan%rmp_yx, plan%fft_y_array, plan%fft_x_array )

    npx_loc = plan%seq_x1_local_sz_x1 
    npy_loc = plan%seq_x1_local_sz_x2 

    do j=1,npy_loc
       call sll_s_fft_exec_c2c_1d(plan%px_inv, plan%fft_x_array(:,j), plan%fft_x_array(:,j))
    end do

    ! Also ensure the periodicity in x
    plan%fft_x_array(npx_loc,:) = plan%fft_x_array(1,:)

    !print*, 'before copying the last line, '
    !print*, 'npx_loc = ', npx_loc, 'npy_loc = ', npy_loc

    phi(1:npx_loc,1:npy_loc) = real(plan%fft_x_array(1:npx_loc,1:npy_loc),f64)

  end subroutine sll_s_poisson_2d_periodic_par_solve

  !> Note that the equation that is solved is: \f$ \Delta \phi = \rho \f$
  !> Thus the user is responsible for giving the proper sign to the source term.
  !> The 'alt' version of this function considers only a domain that does not
  !> include the last, periodic, point
  subroutine sll_s_poisson_2d_periodic_par_solve_alt(plan, rho, phi)
    type (sll_t_poisson_2d_periodic_par), pointer :: plan !< self object
    sll_real64, dimension(:,:)        :: rho      !< charge density
    sll_real64, dimension(:,:)        :: phi      !< electric potential
    sll_int32                         :: ncx      !< global size
    sll_int32                         :: ncy      !< global size
    sll_int32                         :: npx_loc
    sll_int32                         :: npy_loc
    sll_int32                         :: i, j
    ! Reciprocals of domain lengths.
    sll_real64                        :: r_Lx, r_Ly
    sll_real64                        :: kx, ky
    !sll_comp64                        :: val
    sll_real64                        :: normalization
    type(sll_t_layout_2d), pointer          :: layout_x
    type(sll_t_layout_2d), pointer          :: layout_y
    sll_int32, dimension(1:2)         :: global
    sll_int32                         :: gi, gj

    ncx  = plan%ncx
    ncy  = plan%ncy
    r_Lx = 1.0_f64/plan%Lx
    r_Ly = 1.0_f64/plan%Ly

    !print*, 1.0_f64/plan%Lx, 1.0_f64/plan%Ly
    ! Get layouts to compute FFTs (in each direction)
    layout_x => plan%layout_seq_x1
    layout_y => plan%layout_seq_x2
    !call verify_argument_sizes_par(layout_x, rho, phi)

    ! FFTs in x-direction
    npx_loc = plan%seq_x1_local_sz_x1 
    npy_loc = plan%seq_x1_local_sz_x2 

    ! The input is handled internally as a complex array
    plan%fft_x_array = cmplx(rho, 0.0_f64, kind=f64)

    call sll_s_fft_exec_c2c_2d(plan%px, plan%fft_x_array, plan%fft_x_array)
    ! FFTs in y-direction
    npx_loc = plan%seq_x2_local_sz_x1
    npy_loc = plan%seq_x2_local_sz_x2

    call sll_o_apply_remap_2d( plan%rmp_xy, plan%fft_x_array, plan%fft_y_array )

    call sll_s_fft_exec_c2c_2d(plan%py, plan%fft_y_array, plan%fft_y_array) 

    ! This should be inside the FFT plan...
    normalization = 1.0_f64/(ncx*ncy)

    ! Apply the kernel 
    do j=1, npy_loc
       do i=1, npx_loc
          ! Make sure that the first mode is set to zero so that we get an
          ! answer with zero mean value. This step assumes that the (1,1) point
          ! will always be at the border of any splitting of the domains. This 
          ! seems safe in this case.
          
          global = sll_o_local_to_global( layout_y, (/i, j/))
          gi = global(1)
          gj = global(2)
          
          if( (gi == 1) .and. (gj == 1) ) then
             plan%fft_y_array(1,1) = (0.0_f64,0.0_f64)
          else
             kx  = real(gi-1,f64)
             ky  = real(gj-1,f64)
             ! Crucial step: While the results of the FFT are basically a
             ! list of the modes in the 0..n-1 sense, the algorithm that we
             ! are using uses a decomposition with negative frequencies,
             ! i.e.: -n/2..n/2-1. Therefore a step is needed in which we are
             ! effectively reinterpreting the FFT results in this light.
             ! More specifically, we take the second half of the transformed
             ! array and apply a modified factor, proper of a negative 
             ! frequency.
             if( kx .ge. ncx/2 ) then
                kx = kx - ncx
             end if

             if( ky .ge. ncy/2 ) then
                ky = ky - ncy
             end if

              plan%fft_y_array(i,j) = -plan%fft_y_array(i,j)*normalization / &
                  ( ( (kx*r_Lx)**2 + (ky*r_Ly)**2)*4.0_f64*sll_p_pi**2)
          end if
       enddo
    enddo

    ! Inverse FFTs in y-direction
    call sll_s_fft_exec_c2c_2d(plan%py_inv, plan%fft_y_array, plan%fft_y_array) 

    ! Prepare to take inverse FFTs in x-direction
    call sll_o_apply_remap_2d( plan%rmp_yx, plan%fft_y_array, plan%fft_x_array )

    npx_loc = plan%seq_x1_local_sz_x1 
    npy_loc = plan%seq_x1_local_sz_x2 

    call sll_s_fft_exec_c2c_2d(plan%px_inv, plan%fft_x_array, plan%fft_x_array)

    phi(1:npx_loc,1:npy_loc) = real(plan%fft_x_array(1:npx_loc,1:npy_loc),f64)
  end subroutine sll_s_poisson_2d_periodic_par_solve_alt




!> Delete the Poisson solver object
  subroutine sll_s_poisson_2d_periodic_par_free(plan)
    type (sll_t_poisson_2d_periodic_par), pointer :: plan
    sll_int32                                              :: ierr

    if( .not. associated(plan) ) then
       print *, 'ERROR, delete_poisson_3d_periodic_plan_par(): ', &
            'passed plan is not associated.'
       STOP
    end if

    call sll_s_fft_free(plan%px)
    call sll_s_fft_free(plan%py)
    call sll_s_fft_free(plan%px_inv)
    call sll_s_fft_free(plan%py_inv)

!    call delete( plan%layout_x ) ! can't delete this, the plan does not own it
    call sll_o_delete( plan%layout_seq_x1 )
    call sll_o_delete( plan%layout_seq_x2 )
    SLL_DEALLOCATE_ARRAY(plan%fft_x_array, ierr)
    SLL_DEALLOCATE_ARRAY(plan%fft_y_array, ierr)
    call sll_o_delete( plan%rmp_xy )
    call sll_o_delete( plan%rmp_yx )
    SLL_DEALLOCATE(plan, ierr)
  end subroutine sll_s_poisson_2d_periodic_par_free

!> Check that arrays match layout properties
  subroutine verify_argument_sizes_par(layout, rho, phi)
    type(sll_t_layout_2d), pointer       :: layout !< layout for remap
    sll_real64, dimension(:,:)     :: rho    !< charge density
    sll_real64, dimension(:,:)     :: phi    !< electric potential
    sll_int32                      :: nx
    sll_int32                      :: ny
    !sll_int32                      :: i

    ! Note that this checks for strict sizes, not an array being bigger
    ! than a certain size, but exactly a desired size... This may be a bit
    ! too stringent.
    call sll_o_compute_local_sizes( layout, nx, ny )
    ! Verify the first direction
    if ( nx /= size(rho,1) ) then
       print*, 'ERROR: sll_s_poisson_2d_periodic_par_solve()', &
            'size of rho does not match expected size. ', &
            'Expected size according to layout = ', nx, 'Received size = ',&
            size(rho,1)
       print *, 'Exiting...'
       stop
    end if
    if ( nx /= size(phi,1) ) then
       print*, 'ERROR: sll_s_poisson_2d_periodic_par_solve()', &
            'size of phi does not match expected size. ', &
            'Expected size according to layout = ', nx, 'Received size = ',&
            size(phi,1)
       print *, 'Exiting...'
       stop
    end if
    ! Verify the second direction
    if ( ny /= size(rho,2) ) then
       print*, 'ERROR: sll_s_poisson_2d_periodic_par_solve()', &
            'size of rho does not match expected size. ', &
            'Expected size according to layout = ', ny, 'Received size = ',&
            size(rho,2)
       print *, 'Exiting...'
       stop
    end if
    if ( ny /= size(phi,2) ) then
       print*, 'ERROR: sll_s_poisson_2d_periodic_par_solve()', &
            'size of phi does not match expected size. ', &
            'Expected size according to layout = ', ny, 'Received size = ',&
            size(phi,2)
       print *, 'Exiting...'
       stop
    end if
  end subroutine verify_argument_sizes_par


end module sll_m_poisson_2d_periodic_par
