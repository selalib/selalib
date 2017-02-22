!> @ingroup poisson_solvers
!> @brief 
!> Periodic 3D poisson solver (parallel version)
!> @details
!> depends on sll_m_remapper
module sll_m_poisson_3d_periodic_par

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_t_collective_t, &
    sll_f_get_collective_size

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_fft, only: &
    sll_s_fft_exec_c2c_1d, &
    sll_p_fft_backward, &
    sll_s_fft_free, &
    sll_p_fft_forward, &
    sll_s_fft_init_c2c_1d, &
    sll_t_fft

  use sll_m_remapper, only: &
    sll_o_apply_remap_3d, &
    sll_o_compute_local_sizes, &
    sll_o_get_layout_collective, &
    sll_o_initialize_layout_with_distributed_array, &
    sll_t_layout_3d, &
    sll_o_local_to_global, &
    sll_f_new_layout_3d, &
    sll_o_new_remap_plan, &
    sll_t_remap_plan_3d_comp64, &
    sll_t_remap_plan_3d_real64, &
    sll_o_delete, &
    sll_s_factorize_in_two_powers_of_two

  use sll_m_utilities, only: &
    sll_f_is_even, &
    sll_f_is_power_of_two

  implicit none

  public :: &
    sll_s_poisson_3d_periodic_par_compute_e_from_phi, &
    sll_s_poisson_3d_periodic_par_compute_e_from_phi_layoutseq3, &
    sll_s_poisson_3d_periodic_par_free, &
    sll_s_poisson_3d_periodic_par_init, &
    sll_t_poisson_3d_periodic_par, &
    sll_s_poisson_3d_periodic_par_solve

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Structure to solve Poisson equation on 3d mesh with periodic boundary
  !> conditions. Solver is parallel and numerical method is based on fft 
  !> transform.  Number of cells, which in this periodic case is equal to 
  !> the number of points.
  type sll_t_poisson_3d_periodic_par
     sll_int32                             :: ncx       !< number of cells in x
     sll_int32                             :: ncy       !< number of cells in y
     sll_int32                             :: ncz       !< number of cells in z
     sll_real64                            :: Lx        !< x domain length
     sll_real64                            :: Ly        !< y domain length
     sll_real64                            :: Lz        !< z domain length
     sll_real64                            :: dx        !< x domain grid spacing
     sll_real64                            :: dy        !< y domain grid spacing
     sll_real64                            :: dz        !< z domain grid spacing
     type(sll_t_fft)           :: px        !< fft plan in x
     type(sll_t_fft)           :: py        !< fft plan in y
     type(sll_t_fft)           :: pz        !< fft plan in z
     type(sll_t_fft)           :: px_inv    !< inverse fft in x
     type(sll_t_fft)           :: py_inv    !< inverse fft in y
     type(sll_t_fft)           :: pz_inv    !< inverse fft in z
     type(sll_t_layout_3d),  pointer             :: layout_split !< split layout for remap comming from program
     type(sll_t_layout_3d),  pointer             :: layout_x  !< x layout for remap
     type(sll_t_layout_3d),  pointer             :: layout_y  !< y layout for remap
     type(sll_t_layout_3d),  pointer             :: layout_z  !< z layout for remap

     sll_int32, dimension(3,3)             :: loc_sizes !< local sizes
     sll_comp64, dimension(:,:,:), pointer :: array_x   !< x array component
     sll_comp64, dimension(:,:,:), pointer :: array_y   !< y array component
     sll_comp64, dimension(:,:,:), pointer :: array_z   !< z array component
     sll_comp64, allocatable               :: array1d_x(:)
     sll_comp64, allocatable               :: array1d_y(:)
     sll_comp64, allocatable               :: array1d_z(:)
     sll_real64, allocatable               :: phi_x1(:,:,:)
     sll_real64, allocatable               :: phi_x2(:,:,:)
     sll_real64, allocatable               :: phi_x3(:,:,:)
     sll_real64, allocatable               :: ex_x1(:,:,:)
     sll_real64, allocatable               :: ey_x2(:,:,:)
     sll_real64, allocatable               :: ez_x3(:,:,:)

     type(sll_t_remap_plan_3D_comp64), pointer   :: rmp3_xy   !< transpose from x to y
     type(sll_t_remap_plan_3D_comp64), pointer   :: rmp3_yz   !< transpose from y to z
     type(sll_t_remap_plan_3D_comp64), pointer   :: rmp3_zy   !< transpose from z to y
     type(sll_t_remap_plan_3D_comp64), pointer   :: rmp3_yx   !< transpose from y to x
     type(sll_t_remap_plan_3D_real64), pointer   :: rmp_x1_to_split
     type(sll_t_remap_plan_3D_real64), pointer   :: rmp_split_to_x1
     type(sll_t_remap_plan_3D_real64), pointer   :: rmp_x2_to_split
     type(sll_t_remap_plan_3D_real64), pointer   :: rmp_x3_to_split
     type(sll_t_remap_plan_3D_real64), pointer   :: rmp_x1_to_x2
     type(sll_t_remap_plan_3D_real64), pointer   :: rmp_x2_to_x3
     type(sll_t_remap_plan_3D_real64), pointer   :: rmp_x3_to_x1
  end type sll_t_poisson_3d_periodic_par

contains


  !> Allocate the structure for the 3d parallel Poisson solver
  !> @return
  subroutine sll_s_poisson_3d_periodic_par_init( &
    start_layout, &
    ncx, &
    ncy, &
    ncz, &
    Lx, &
    Ly, &
    Lz, & 
    plan)

    type(sll_t_layout_3d),  pointer,     intent( in    )  :: start_layout !< intiial layout
    sll_int32,                           intent( in    )  :: ncx !< number of cells in x
    sll_int32,                           intent( in    )  :: ncy !< number of cells in y
    sll_int32,                           intent( in    )  :: ncz !< number of cells in z
    sll_real64,                          intent( in    )  :: Lx  !< x length
    sll_real64,                          intent( in    )  :: Ly  !< y length
    sll_real64,                          intent( in    )  :: Lz  !< z length
    type(sll_t_poisson_3d_periodic_par), intent(   out )  :: plan !< Poisson solver object
  
    
    sll_comp64                                            :: x(ncx)   !< 1d array in x
    sll_comp64                                            :: y(ncy)   !< 1d array in y
    sll_comp64                                            :: z(ncz)   !< 1d array in z
    sll_int64                                             :: colsz ! collective size
    type(sll_t_collective_t), pointer                     :: collective
    ! npx, npy, npz are the numbers of processors in directions x, y, z
    sll_int32                                             :: npx, npy, npz
    sll_int32                                             :: ierr
    sll_int32, dimension(3,3)                             :: loc_sizes
    sll_real64, allocatable                               :: tmp(:,:,:)

    collective => sll_o_get_layout_collective(start_layout)
    colsz      = int(sll_f_get_collective_size(collective),i64)

    if ( int(colsz,i32) > min(ncx*ncy,ncy*ncz,ncz*ncx) ) then
       print *, 'This test needs to run in a number of processes which',  &
                ' is less than or equal', min(ncx,ncy,ncz), ' in order to ', &
                'be able properly split the arrays.'
       print *, 'Exiting...'
       stop
    end if

    plan%ncx = ncx
    plan%ncy = ncy
    plan%ncz = ncz
    plan%Lx  = Lx
    plan%Ly  = Ly
    plan%Lz  = Lz

    plan%dx = Lx/ncx
    plan%dy = Ly/ncy
    plan%dz = Lz/ncz

    ! For FFTs (in each direction)
    call sll_s_fft_init_c2c_1d( plan%px, ncx, x, x, sll_p_fft_forward )
    call sll_s_fft_init_c2c_1d( plan%py, ncy, y, y, sll_p_fft_forward )
    call sll_s_fft_init_c2c_1d( plan%pz, ncz, z, z, sll_p_fft_forward )

    ! For inverse FFTs (in each direction)
    call sll_s_fft_init_c2c_1d( plan%px_inv, ncx, x, x, sll_p_fft_backward )
    call sll_s_fft_init_c2c_1d( plan%py_inv, ncy, y, y, sll_p_fft_backward )
    call sll_s_fft_init_c2c_1d( plan%pz_inv, ncz, z, z, sll_p_fft_backward )

    ! Layout and local sizes for FFTs in x-direction
    plan%layout_split => start_layout

    ! Layout and local sizes for FFTs in x-direction
    plan%layout_x => sll_f_new_layout_3d( collective )
    npx = 1
    call sll_s_factorize_in_two_powers_of_two( int(colsz,i32), npy, npz)
    
    call sll_o_initialize_layout_with_distributed_array( &
         ncx, &
         ncy, &
         ncz, &
         npx, &
         npy, &
         npz, &
         plan%layout_x )
    call sll_o_compute_local_sizes( &
         plan%layout_x, &
         loc_sizes(1,1), &
         loc_sizes(1,2), &
         loc_sizes(1,3)  )

    ! Layout and local sizes for FFTs in y-direction
    plan%layout_y => sll_f_new_layout_3d( collective )
    npx = npy
    npy = 1
    call sll_o_initialize_layout_with_distributed_array( &
         ncx, &
         ncy, &
         ncz, &
         npx, &
         npy, &
         npz, &
         plan%layout_y )

    call sll_o_compute_local_sizes( &
         plan%layout_y,  &
         loc_sizes(2,1), &
         loc_sizes(2,2), &
         loc_sizes(2,3) )

    ! Layout and local sizes for FFTs in z-direction
    plan%layout_z => sll_f_new_layout_3d( collective )
    ! npx remains the same. Exchange npy and npz.
    npy = npz
    npz = 1
    call sll_o_initialize_layout_with_distributed_array( &
         ncx, &
         ncy, &
         ncz, &
         npx, &
         npy, &
         npz, &
         plan%layout_z )

    call sll_o_compute_local_sizes( &
         plan%layout_z, &
         loc_sizes(3,1), &
         loc_sizes(3,2), &
         loc_sizes(3,3) )

    plan%loc_sizes = loc_sizes

    SLL_ALLOCATE(plan%array_x(loc_sizes(1,1),loc_sizes(1,2),loc_sizes(1,3)),ierr)
    SLL_ALLOCATE(plan%array_y(loc_sizes(2,1),loc_sizes(2,2),loc_sizes(2,3)),ierr)
    SLL_ALLOCATE(plan%array_z(loc_sizes(3,1),loc_sizes(3,2),loc_sizes(3,3)),ierr)

    allocate(plan%array1d_x(loc_sizes(1,1)))
    allocate(plan%array1d_y(loc_sizes(2,2)))
    allocate(plan%array1d_z(loc_sizes(3,3)))

    plan%rmp3_xy => sll_o_new_remap_plan(plan%layout_x, plan%layout_y, plan%array_x)
    plan%rmp3_yz => sll_o_new_remap_plan(plan%layout_y, plan%layout_z, plan%array_y)
    plan%rmp3_zy => sll_o_new_remap_plan(plan%layout_z, plan%layout_y, plan%array_z)
    plan%rmp3_yx => sll_o_new_remap_plan(plan%layout_y, plan%layout_x, plan%array_y)

    SLL_ALLOCATE(plan%phi_x1(loc_sizes(1,1), loc_sizes(1,2), loc_sizes(1,3)), ierr)
    SLL_ALLOCATE(plan%phi_x2(loc_sizes(2,1), loc_sizes(2,2), loc_sizes(2,3)), ierr)
    SLL_ALLOCATE(plan%phi_x3(loc_sizes(3,1), loc_sizes(3,2), loc_sizes(3,3)), ierr)
    SLL_ALLOCATE(plan%ex_x1(loc_sizes(1,1), loc_sizes(1,2), loc_sizes(1,3)), ierr)
    SLL_ALLOCATE(plan%ey_x2(loc_sizes(2,1), loc_sizes(2,2), loc_sizes(2,3)), ierr)
    SLL_ALLOCATE(plan%ez_x3(loc_sizes(3,1), loc_sizes(3,2), loc_sizes(3,3)), ierr)

    plan%rmp_x1_to_split => sll_o_new_remap_plan(plan%layout_x, plan%layout_split, plan%ex_x1)
    plan%rmp_x2_to_split => sll_o_new_remap_plan(plan%layout_y, plan%layout_split, plan%ey_x2)
    plan%rmp_x3_to_split => sll_o_new_remap_plan(plan%layout_z, plan%layout_split, plan%ez_x3)
    plan%rmp_x1_to_x2 => sll_o_new_remap_plan(plan%layout_x, plan%layout_y, plan%ex_x1)
    plan%rmp_x2_to_x3 => sll_o_new_remap_plan(plan%layout_y, plan%layout_z, plan%ey_x2)
    plan%rmp_x3_to_x1 => sll_o_new_remap_plan(plan%layout_z, plan%layout_x, plan%ez_x3)

    call sll_o_compute_local_sizes( &
         plan%layout_split, &
         loc_sizes(1,1), &
         loc_sizes(1,2), &
         loc_sizes(1,3) )

    SLL_ALLOCATE(tmp(loc_sizes(1,1),loc_sizes(1,2),loc_sizes(1,3)),ierr)
    plan%rmp_split_to_x1 => sll_o_new_remap_plan(plan%layout_split, plan%layout_x, tmp)

  end subroutine sll_s_poisson_3d_periodic_par_init


  !> Compute the 3d potential from the Poisson equation with periodic
  !> boundary conditions.
  subroutine sll_s_poisson_3d_periodic_par_solve(plan, rho, phi)
    type (sll_t_poisson_3d_periodic_par)         :: plan !< Solver structure
    sll_real64, dimension(:,:,:)                 :: rho  !< Charge density
    sll_real64, dimension(:,:,:)                 :: phi  !< Electric potential
    sll_int32                                    :: nx, ny, nz
    ! nx, ny, nz are the numbers of points - 1 in directions x, y ,z
    sll_int32                                    :: nx_loc, ny_loc, nz_loc
    sll_int32                                    :: i, j, k
    sll_real64                                   :: Lx, Ly, Lz
    sll_real64                                   :: ind_x, ind_y, ind_z
    type(sll_t_layout_3d), pointer               :: layout_x
    type(sll_t_layout_3d), pointer               :: layout_y
    type(sll_t_layout_3d), pointer               :: layout_z
    sll_int32, dimension(1:3)                    :: global
    sll_int32                                    :: gi, gj, gk

    ! Get geometry information
    nx = plan%ncx
    ny = plan%ncy
    nz = plan%ncz
    Lx = plan%Lx
    Ly = plan%Ly
    Lz = plan%Lz

    ! Get layouts to compute FFTs (in each direction)
    layout_x => plan%layout_x
    layout_y => plan%layout_y
    layout_z => plan%layout_z

    call sll_o_apply_remap_3d( plan%rmp_split_to_x1, rho, plan%ex_x1)

    !call verify_argument_sizes_par(layout_x, rho, phi)

    ! FFTs in x-direction
    nx_loc = plan%loc_sizes(1,1) 
    ny_loc = plan%loc_sizes(1,2) 
    nz_loc = plan%loc_sizes(1,3)
    !plan%array_x = cmplx(plan%ex_x1, 0_f64, kind=f64)
    do k=1,nz_loc
       do j=1,ny_loc
          plan%array1d_x = cmplx(plan%ex_x1(:,j,k), 0_f64, kind=f64)
          call sll_s_fft_exec_c2c_1d(plan%px, plan%array1d_x, plan%array1d_x)
          !call sll_s_fft_apply_plan_c2c_1d(plan%px, plan%array_x(:,j,k), plan%array_x(:,j,k))
          plan%array_x(:,j,k) = plan%array1d_x
       enddo
    enddo

    ! FFTs in y-direction
    nx_loc = plan%loc_sizes(2,1) 
    ny_loc = plan%loc_sizes(2,2) 
    nz_loc = plan%loc_sizes(2,3)
    call sll_o_apply_remap_3d( plan%rmp3_xy, plan%array_x, plan%array_y ) 
    do k=1,nz_loc
       do i=1,nx_loc
          plan%array1d_y = plan%array_y(i,:,k)
          call sll_s_fft_exec_c2c_1d(plan%py, plan%array1d_y, plan%array1d_y)
          plan%array_y(i,:,k) = plan%array1d_y
       enddo
    enddo

    ! FFTs in z-direction
    nx_loc = plan%loc_sizes(3,1) 
    ny_loc = plan%loc_sizes(3,2) 
    nz_loc = plan%loc_sizes(3,3)
    call sll_o_apply_remap_3d( plan%rmp3_yz, plan%array_y, plan%array_z ) 
    do j=1,ny_loc
       do i=1,nx_loc
          plan%array1d_z = plan%array_z(i,j,:)
          call sll_s_fft_exec_c2c_1d(plan%pz, plan%array1d_z, plan%array1d_z)
          plan%array_z(i,j,:) = plan%array1d_z
       enddo
    enddo

    ! move this normalization elsewhere where the modes are modified
    plan%array_z = plan%array_z/(nx*ny*nz)

    ! Compute hat_phi, phi = inv_fft(hat_phi)
    do k=1,nz_loc
       do j=1,ny_loc
          do i=1,nx_loc
             global = sll_o_local_to_global( layout_z, (/i, j, k/))
             gi = global(1)
             gj = global(2)
             gk = global(3)
             if( (gi==1) .and. (gj==1) .and. (gk==1) ) then
                plan%array_z(1,1,1) = (0.0_f64,0.0_f64)
             else
                if (gi<=nx/2) then
                   ind_x = real(gi-1,f64)
                else
                   ind_x = real(nx-(gi-1),f64)
                endif
                if (gj<=ny/2) then
                   ind_y = real(gj-1,f64)
                else
                   ind_y = real(ny-(gj-1),f64)
                endif
                if (gk<=nz/2) then
                   ind_z = real(gk-1,f64)
                else
                   ind_z = real(nz-(gk-1),f64)
                endif
                plan%array_z(i,j,k) = plan%array_z(i,j,k)/(4*sll_p_pi**2 * &
                            ((ind_x/Lx)**2 + (ind_y/Ly)**2+(ind_z/Lz)**2))
             endif
          enddo
       enddo
    enddo

    ! Inverse FFTs in z-direction
    do j=1,ny_loc
       do i=1,nx_loc
          plan%array1d_z = plan%array_z(i,j,:)
          call sll_s_fft_exec_c2c_1d( &
               plan%pz_inv, &
               plan%array1d_z, &
               plan%array1d_z)
          plan%array_z(i,j,:) = plan%array1d_z
       enddo
    enddo

    ! Inverse FFTs in y-direction
    nx_loc = plan%loc_sizes(2,1) 
    ny_loc = plan%loc_sizes(2,2) 
    nz_loc = plan%loc_sizes(2,3)
    call sll_o_apply_remap_3d( plan%rmp3_zy, plan%array_z, plan%array_y )
    do k=1,nz_loc
       do i=1,nx_loc
          plan%array1d_y = plan%array_y(i,:,k)
          call sll_s_fft_exec_c2c_1d( &
               plan%py_inv, &
               plan%array1d_y, &
               plan%array1d_y )
          plan%array_y(i,:,k) = plan%array1d_y
       enddo
    enddo

    ! Inverse FFTs in x-direction
    nx_loc = plan%loc_sizes(1,1) 
    ny_loc = plan%loc_sizes(1,2) 
    nz_loc = plan%loc_sizes(1,3)
    call sll_o_apply_remap_3d( plan%rmp3_yx, plan%array_y, plan%array_x ) 
    do k=1,nz_loc
       do j=1,ny_loc
          plan%array1d_x = plan%array_x(:,j,k)
          call sll_s_fft_exec_c2c_1d( &
               plan%px_inv, &
               plan%array1d_x, &
               plan%array1d_x )
          plan%array_x(:,j,k) = plan%array1d_x
       enddo
    enddo

    phi = real(plan%array_x, f64)

  end subroutine sll_s_poisson_3d_periodic_par_solve

  subroutine sll_s_poisson_3d_periodic_par_compute_e_from_phi(plan, phi, ex, ey, ez)
    type (sll_t_poisson_3d_periodic_par),        intent(inout) :: plan !< Solver structure
    sll_real64,                                  intent(in   ) :: phi(:,:,:)  !< Electric potential
    sll_real64,                                  intent(  out) :: ex(:,:,:)
    sll_real64,                                  intent(  out) :: ey(:,:,:)
    sll_real64,                                  intent(  out) :: ez(:,:,:)
   

    ! compute the values of the electric field. rho is configured for
    ! sequential operations in x1, thus we start by computing the E_x
    ! component.
    ! The following call is inefficient and unnecessary. The local sizes for
    ! the arrays should be kept around as parameters basically and not on
    ! variables whose content could be anything... This will have to do for
    ! now.
    call compute_electric_field_x1_3d( plan, phi, plan%ex_x1)
    call sll_o_apply_remap_3d( plan%rmp_x1_to_split, plan%ex_x1, ex)

    call sll_o_apply_remap_3d( plan%rmp_x1_to_x2, phi, plan%phi_x2)
    call compute_electric_field_x2_3d( plan, plan%phi_x2, plan%ey_x2)
    call sll_o_apply_remap_3d( plan%rmp_x2_to_split, plan%ey_x2, ey)


    call sll_o_apply_remap_3d( plan%rmp_x2_to_x3, plan%phi_x2, plan%phi_x3)
    call compute_electric_field_x3_3d( plan, plan%phi_x3, plan%ez_x3)
    call sll_o_apply_remap_3d( plan%rmp_x3_to_split, plan%ez_x3, ez)

  end subroutine sll_s_poisson_3d_periodic_par_compute_e_from_phi


  subroutine sll_s_poisson_3d_periodic_par_compute_e_from_phi_layoutseq3(plan, phi, ex, ey, ez)
    type (sll_t_poisson_3d_periodic_par),        intent(inout) :: plan !< Solver structure
    sll_real64,                                  intent(in   ) :: phi(:,:,:)  !< Electric potential
    sll_real64,                                  intent(  out) :: ex(:,:,:)
    sll_real64,                                  intent(  out) :: ey(:,:,:)
    sll_real64,                                  intent(  out) :: ez(:,:,:)
   

    ! compute the values of the electric field. rho is configured for
    ! sequential operations in x1, thus we start by computing the E_x
    ! component.
    ! The following call is inefficient and unnecessary. The local sizes for
    ! the arrays should be kept around as parameters basically and not on
    ! variables whose content could be anything... This will have to do for
    ! now.
    call compute_electric_field_x3_3d( plan, phi, plan%ez_x3)
    call sll_o_apply_remap_3d( plan%rmp_x3_to_split, plan%ez_x3, ez)

    call sll_o_apply_remap_3d( plan%rmp_x3_to_x1, phi, plan%phi_x1)
    call compute_electric_field_x1_3d( plan, plan%phi_x1, plan%ex_x1)
    call sll_o_apply_remap_3d( plan%rmp_x1_to_split, plan%ex_x1, ex)

    call sll_o_apply_remap_3d( plan%rmp_x1_to_x2, plan%phi_x1, plan%phi_x2)
    call compute_electric_field_x2_3d( plan, plan%phi_x2, plan%ey_x2)
    call sll_o_apply_remap_3d( plan%rmp_x2_to_split, plan%ey_x2, ey)


  end subroutine sll_s_poisson_3d_periodic_par_compute_e_from_phi_layoutseq3


  !> Delete the solver structure
  subroutine sll_s_poisson_3d_periodic_par_free(plan)
    type (sll_t_poisson_3d_periodic_par)     :: plan
    sll_int32                                    :: ierr


    call sll_s_fft_free(plan%px)
    call sll_s_fft_free(plan%py)
    call sll_s_fft_free(plan%pz)

    call sll_s_fft_free(plan%px_inv)
    call sll_s_fft_free(plan%py_inv)
    call sll_s_fft_free(plan%pz_inv)

    call sll_o_delete( plan%layout_x )
    call sll_o_delete( plan%layout_y )
    call sll_o_delete( plan%layout_z )

    SLL_DEALLOCATE_ARRAY(plan%array_x, ierr)
    SLL_DEALLOCATE_ARRAY(plan%array_y, ierr)
    SLL_DEALLOCATE_ARRAY(plan%array_z, ierr)
    
  end subroutine sll_s_poisson_3d_periodic_par_free

  !> Check sizes of arrays in input
  subroutine verify_argument_sizes_par(layout, rho, phi)
    type(sll_t_layout_3d), pointer   :: layout
    sll_real64, dimension(:,:,:) :: rho
    sll_real64, dimension(:,:,:) :: phi
    sll_int32,  dimension(3)     :: n ! nx_loc, ny_loc, nz_loc
    sll_int32                    :: i

    call sll_o_compute_local_sizes( layout, n(1), n(2), n(3) )

    do i=1,3
       if ( (n(i)/=size(rho,i)) .or. (n(i)/=size(phi,i))  ) then
          print*, 'Input sizes passed to sll_s_poisson_3d_periodic_par_solve do not match'
          if (i==1) then
             print*, 'Input sizes passed to "sll_s_poisson_3d_periodic_par_solve" ', &
                  'do not match in direction x'
          elseif (i==2) then
             print*, 'Input sizes passed to "sll_s_poisson_3d_periodic_par_solve" ', &
                  'do not match in direction y'
          else
             print*, 'Input sizes passed to "sll_s_poisson_3d_periodic_par_solve" ', &
                  'do not match in direction z'
          endif
          print *, 'Exiting...'
          stop
       endif
    enddo

  end subroutine verify_argument_sizes_par


!------------------------------------------------------------------------------!

 subroutine compute_electric_field_x1_3d( &
      plan, &
      phi_x1, &
      efield_x1)
   type(sll_t_poisson_3d_periodic_par), intent(inout) :: plan !< Solver structure
   sll_real64,                         intent(in)    :: phi_x1(:,:,:)
   sll_real64,                         intent(out)   :: efield_x1(:,:,:)

   ! local variables
   sll_int32                                 :: num_pts_x1
   sll_int32                                 :: i
   sll_int32                                 :: j
   sll_int32                                 :: k
   sll_real64                                :: r_delta
   sll_real64                                :: ex
   ! FIXME: arg checking

    num_pts_x1 = plan%loc_sizes(1,1)

    r_delta = 1.0_f64/plan%dx

    ! Compute the electric field values on the left and right edges.
    do k=1,plan%loc_sizes(1,3)
       do j=1,plan%loc_sizes(1,2)
          ! i=1 plane:
          ex = r_delta*(-1.5_f64*phi_x1(1,j,k) + &
                         2.0_f64*phi_x1(2,j,k) - &
                         0.5_f64*phi_x1(3,j,k) )
          efield_x1(1,j,k) = -ex
          ! i=num_pts_x1 plane:
          ex = r_delta*(0.5_f64*phi_x1(num_pts_x1-2,j,k) - &
                        2.0_f64*phi_x1(num_pts_x1-1,j,k) + &
                        1.5_f64*phi_x1(num_pts_x1  ,j,k) )
          efield_x1(num_pts_x1,j,k) = -ex
       end do
    end do

    ! Electric field in interior points
    do k=1,plan%loc_sizes(1,3)
       do j=1,plan%loc_sizes(1,2)
          do i=2, num_pts_x1-1
             ex = r_delta*0.5_f64*(phi_x1(i+1,j,k) - phi_x1(i-1,j,k))
             efield_x1(i,j,k) = -ex
          end do
       end do
    end do

  end subroutine compute_electric_field_x1_3d


  ! This function only sets the Ey component of the electric field.
  subroutine compute_electric_field_x2_3d( & 
      plan, &
      phi_x2, &
      efield_x2)
   type(sll_t_poisson_3d_periodic_par), intent(inout) :: plan !< Solver structure
   sll_real64,                         intent(in)    :: phi_x2(:,:,:)
   sll_real64,                         intent(out)   :: efield_x2(:,:,:)

   ! local variables
   sll_int32                                 :: num_pts_x2
   sll_int32                                 :: i
   sll_int32                                 :: j
   sll_int32                                 :: k
   sll_real64                                :: r_delta
   sll_real64                                :: ey

    ! FIXME: arg checking
    num_pts_x2 = plan%loc_sizes(2,2)
    r_delta = 1.0_f64/plan%dy

    ! Compute the electric field values on the bottom and top edges.
    do k=1,plan%loc_sizes(2,3)
       do i=1,plan%loc_sizes(2,1)
          ! bottom:
          ey = r_delta*(-1.5_f64*phi_x2(i,1,k) + 2.0_f64*phi_x2(i,2,k) - &
                         0.5_f64*phi_x2(i,3,k))
          efield_x2(i,1,k) = -ey
          ! top:
          ey = r_delta*(0.5_f64*phi_x2(i,num_pts_x2-2,k) - &
                        2.0_f64*phi_x2(i,num_pts_x2-1,k) + &
                        1.5_f64*phi_x2(i,num_pts_x2  ,k))
          efield_x2(i,num_pts_x2,k) = -ey
       end do
    end do

    ! Electric field in interior points
    do k=1,plan%loc_sizes(2,3)
       do j=2,num_pts_x2-1
          do i=1,plan%loc_sizes(2,1)
             ey = r_delta*0.5_f64*(phi_x2(i,j+1,k) - phi_x2(i,j-1,k))
             efield_x2(i,j,k) = -ey
          end do
       end do
    end do

  end subroutine compute_electric_field_x2_3d

  ! This function only sets the Ey component of the electric field.
  subroutine compute_electric_field_x3_3d( & 
      plan, &
      phi_x3, &
      efield_x3)
   type(sll_t_poisson_3d_periodic_par), intent(inout) :: plan !< Solver structure
   sll_real64,                         intent(in)    :: phi_x3(:,:,:)
   sll_real64,                         intent(out)   :: efield_x3(:,:,:)

   ! local variables
   sll_int32                                 :: num_pts_x1
   sll_int32                                 :: num_pts_x2
   sll_int32                                 :: num_pts_x3
   sll_int32                                 :: i
   sll_int32                                 :: j
   sll_int32                                 :: k
   sll_real64                                :: r_delta
   sll_real64                                :: ez

   ! FIXME: arg checking
   num_pts_x1 = plan%loc_sizes(3,1)
   num_pts_x2 = plan%loc_sizes(3,2)
   num_pts_x3 = plan%loc_sizes(3,3)


   r_delta = 1.0_f64/plan%dz

    ! Compute the electric field values on the end faces.
    do j=1,num_pts_x2
       do i=1,num_pts_x1
          ez = r_delta*(-1.5_f64*phi_x3(i,j,1) + 2.0_f64*phi_x3(i,j,2) - &
                         0.5_f64*phi_x3(i,j,3))
          efield_x3(i,j,1) = -ez
          ! top:
          ez = r_delta*(0.5_f64*phi_x3(i,j,num_pts_x3-2) - &
                        2.0_f64*phi_x3(i,j,num_pts_x3-1) + &
                        1.5_f64*phi_x3(i,j,num_pts_x3  ))
          efield_x3(i,j,num_pts_x3) = -ez
       end do
    end do

    ! Electric field in interior points
    do k=2,num_pts_x3-1
       do j=1,num_pts_x2
          do i=1, num_pts_x1
             ez = r_delta*0.5_f64*(phi_x3(i,j,k+1) - phi_x3(i,j,k-1))
             efield_x3(i,j,k) = -ez
          end do
       end do
    end do

  end subroutine compute_electric_field_x3_3d



end module sll_m_poisson_3d_periodic_par
