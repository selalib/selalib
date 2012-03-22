
!**************************************************************************
!
! Selalib 2012     
! Module: sll_poisson_3d_periodic.F90
!
!> @brief 
!> Selalib periodic 3D poisson solver
!> Start date: Feb. 08, 2012
!> Last modification: March 22, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!**************************************************************************

module sll_poisson_3d_periodic_par

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
#include "sll_assert.h"
#include "sll_remap.h"

  use sll_fft
  use numeric_constants
  use sll_collective

  implicit none

  type poisson_3d_periodic_plan_par
     sll_int32                   :: nx ! Number of points-1 in x-direction
     sll_int32                   :: ny ! Number of points-1 in y-direction
     sll_int32                   :: nz ! Number of points-1 in z-direction
     sll_real64                  :: Lx
     sll_real64                  :: Ly
     sll_real64                  :: Lz
     type(sll_fft_plan), pointer :: px
     type(sll_fft_plan), pointer :: py
     type(sll_fft_plan), pointer :: pz
     type(sll_fft_plan), pointer :: px_inv
     type(sll_fft_plan), pointer :: py_inv
     type(sll_fft_plan), pointer :: pz_inv
     type(layout_3D_t),  pointer :: layout_x
     type(layout_3D_t),  pointer :: layout_y
     type(layout_3D_t),  pointer :: layout_z
     sll_int32, dimension(3,3)   :: loc_sizes
  end type poisson_3d_periodic_plan_par

contains


  function new_poisson_3d_periodic_plan_par(start_layout, nx, ny, nz, Lx, &
                                                       Ly, Lz) result(plan)

    type(layout_3D_t),  pointer                  :: start_layout
    sll_int32                                    :: nx, ny, nz
    ! nx, ny, nz are the numbers of points - 1 in directions x, y, z
    sll_comp64,                    dimension(nx) :: x
    sll_comp64,                    dimension(ny) :: y
    sll_comp64,                    dimension(nz) :: z
    sll_real64                                   :: Lx, Ly, Lz
    sll_int64                                    :: colsz ! collective size
    sll_int32                                    :: npx, npy, npz
    ! npx, npy, npz are the numbers of processors in directions x, y, z
    sll_int32                                    :: e
    sll_int32                                    :: ierr
    type (poisson_3d_periodic_plan_par), pointer :: plan
    sll_int32, dimension(3,3)                    :: loc_sizes

    if ( (.not.is_power_of_two(int(nx,i64))) .and. (.not.is_power_of_two( & 
              int(ny,i64))) .and. (.not.is_power_of_two(int(nz,i64)))) then     
       print *, 'This test needs to run in numbers of points which are',  &
                'powers of 2.' 
       stop
    end if

    colsz  = sll_get_collective_size(sll_world_collective)
    if ( colsz > min(nx,ny,nz) ) then     
       print *, 'This test needs to run in a number of processes which',  &
                'is less than', min(nx,ny,nz)
       stop
    end if
    if ( (.not.is_power_of_two(int(nx,i64))) .and. (.not.is_power_of_two( &
               int(ny,i64))) .and.(.not.is_power_of_two(int(nz,i64)))) then     
       print *, 'This test needs to run in numbers of points which are',  &
                'powers of 2.'
       stop
    end if

    SLL_ALLOCATE(plan, ierr)

    ! Geometry informations
    plan%nx = nx
    plan%ny = ny
    plan%nz = nz
    plan%Lx = Lx
    plan%Ly = Ly
    plan%Lz = Lz

    ! For FFTs (in each direction)
    plan%px => new_plan_c2c_1d( nx, x, x, FFT_FORWARD )
    plan%py => new_plan_c2c_1d( ny, y, y, FFT_FORWARD )
    plan%pz => new_plan_c2c_1d( nz, z, z, FFT_FORWARD )

    ! For inverse FFTs (in each direction)
    plan%px_inv => new_plan_c2c_1d( nx, x, x, FFT_INVERSE )
    plan%py_inv => new_plan_c2c_1d( ny, y, y, FFT_INVERSE )
    plan%pz_inv => new_plan_c2c_1d( nz, z, z, FFT_INVERSE )

    ! Layout and local sizes for FFTs in x-direction
    plan%layout_x => start_layout
    call compute_local_sizes( plan%layout_x, loc_sizes(1,1), &
                              loc_sizes(1,2), loc_sizes(1,3) )

    ! Layout and local sizes for FFTs in y-direction
    e = int(log(real(colsz))/log(2.))
    plan%layout_y => new_layout_3D( sll_world_collective )
    npx = 2**(e/2)
    npy = 1
    npz = int(colsz)/npx
    call initialize_layout_with_distributed_3D_array( nx, ny, &
                             nz, npx, npy, npz, plan%layout_y )
    call compute_local_sizes( plan%layout_y, loc_sizes(2,1),  &
                               loc_sizes(2,2), loc_sizes(2,3) )

    ! Layout and local sizes for FFTs in z-direction
    plan%layout_z => new_layout_3D( sll_world_collective )
    npy = npz
    npz = 1
    call initialize_layout_with_distributed_3D_array( nx, ny, &
                             nz, npx, npy, npz, plan%layout_z )
    call compute_local_sizes( plan%layout_z, loc_sizes(3,1),  &
                               loc_sizes(3,2), loc_sizes(3,3) )

    plan%loc_sizes = loc_sizes

  end function new_poisson_3d_periodic_plan_par


  subroutine solve_poisson_3d_periodic_par(plan, rho, phi)

    type (poisson_3d_periodic_plan_par), pointer :: plan
    sll_real64, dimension(:,:,:)                 :: rho
    sll_real64, dimension(:,:,:)                 :: phi
    sll_comp64, dimension(:,:,:), allocatable    :: hat_rho, tmp, hat_phi
    sll_int32                                    :: nx, ny, nz
    ! nx, ny, nz are the numbers of points - 1 in directions x, y ,z
    sll_int32                                    :: nx_loc, ny_loc, nz_loc
    sll_int32                                    :: i, j, k
    sll_int32                                    :: ierr
    sll_real64                                   :: Lx, Ly, Lz
    sll_real64                                   :: ind_x, ind_y, ind_z
    sll_int32                                    :: myrank
    sll_int64                                    :: colsz!collective size
    type(layout_3D_t), pointer                   :: layout_x
    type(layout_3D_t), pointer                   :: layout_y
    type(layout_3D_t), pointer                   :: layout_z
    type(remap_plan_3D_t), pointer               :: rmp3
    sll_int32, dimension(1:3)                    :: global
    sll_int32                                    :: gi, gj, gk

    ! Get geometry informations
    nx = plan%nx
    ny = plan%ny
    nz = plan%nz
    Lx = plan%Lx
    Ly = plan%Ly
    Lz = plan%Lz

    ! Get layouts to compute FFTs (in each direction) and poisson solver 
    !kernel
    layout_x => plan%layout_x
    layout_y => plan%layout_y
    layout_z => plan%layout_z

    ! FFTs in x-direction
    nx_loc = plan%loc_sizes(1,1) 
    ny_loc = plan%loc_sizes(1,2) 
    nz_loc = plan%loc_sizes(1,3)
    call verify_argument_sizes_par(layout_x, rho, phi)
    SLL_ALLOCATE(hat_rho(nx_loc,ny_loc,nz_loc), ierr)
    hat_rho = cmplx(rho, 0_f64, kind=f64)
    do k=1,nz_loc
       do j=1,ny_loc
          call apply_fft_c2c_1d( plan%px, hat_rho(:,j,k), hat_rho(:,j,k) )
       enddo
    enddo

    ! FFTs in y-direction
    nx_loc = plan%loc_sizes(2,1) 
    ny_loc = plan%loc_sizes(2,2) 
    nz_loc = plan%loc_sizes(2,3)
    SLL_ALLOCATE(tmp(nx_loc,ny_loc,nz_loc), ierr)
    rmp3 => NEW_REMAPPER_PLAN_3D( layout_x, layout_y, hat_rho )
    call apply_remap_3D( rmp3, hat_rho, tmp ) 
    do k=1,nz_loc
       do i=1,nx_loc
          call apply_fft_c2c_1d( plan%py, tmp(i,:,k), tmp(i,:,k) )
       enddo
    enddo
    SLL_DEALLOCATE_ARRAY(hat_rho, ierr)

    ! FFTs in z-direction
    nx_loc = plan%loc_sizes(3,1) 
    ny_loc = plan%loc_sizes(3,2) 
    nz_loc = plan%loc_sizes(3,3)
    SLL_ALLOCATE(hat_rho(nx_loc,ny_loc,nz_loc), ierr)
    rmp3 => NEW_REMAPPER_PLAN_3D( layout_y, layout_z, tmp )
    call apply_remap_3D( rmp3, tmp, hat_rho ) 
    do j=1,ny_loc
       do i=1,nx_loc
          call apply_fft_c2c_1d( plan%pz, hat_rho(i,j,:), hat_rho(i,j,:) )
       enddo
    enddo
    SLL_DEALLOCATE_ARRAY(tmp, ierr)
    hat_rho = hat_rho/(nx*ny*nz)

    ! Compute hat_phi, phi = inv_fft(hat_phi)
    SLL_ALLOCATE(hat_phi(nx_loc,ny_loc,nz_loc), ierr)
    do k=1,nz_loc
       do j=1,ny_loc
          do i=1,nx_loc
             global = local_to_global_3D( layout_z, (/i, j, k/))
             gi = global(1)
             gj = global(2)
             gk = global(3)
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
             if ( (ind_x==0) .and. (ind_y==0) .and. (ind_z==0) ) then
                hat_phi(i,j,k) = 0.d0
             else
                hat_phi(i,j,k) = hat_rho(i,j,k)/(4*sll_pi**2*((ind_x/Lx)**2 &
                                 + (ind_y/Ly)**2+(ind_z/Lz)**2))
             endif
          enddo
       enddo
    enddo
    SLL_DEALLOCATE_ARRAY(hat_rho, ierr)

    ! Inverse FFTs in z-direction
    do j=1,ny_loc
       do i=1,nx_loc
          call apply_fft_c2c_1d( plan%pz_inv, hat_phi(i,j,:), hat_phi(i,j,:) )
       enddo
    enddo

    ! Inverse FFTs in y-direction
    nx_loc = plan%loc_sizes(2,1) 
    ny_loc = plan%loc_sizes(2,2) 
    nz_loc = plan%loc_sizes(2,3)
    SLL_ALLOCATE(tmp(nx_loc,ny_loc,nz_loc), ierr)
    rmp3 => NEW_REMAPPER_PLAN_3D( layout_z, layout_y, tmp )
    call apply_remap_3D( rmp3, hat_phi, tmp )
    do k=1,nz_loc
       do i=1,nx_loc
          call apply_fft_c2c_1d( plan%py_inv, tmp(i,:,k), tmp(i,:,k) )
       enddo
    enddo
    SLL_DEALLOCATE_ARRAY(hat_phi, ierr)

    ! Inverse FFTs in x-direction
    nx_loc = plan%loc_sizes(1,1) 
    ny_loc = plan%loc_sizes(1,2) 
    nz_loc = plan%loc_sizes(1,3)
    SLL_ALLOCATE(hat_phi(nx_loc,ny_loc,nz_loc), ierr)
    rmp3 => NEW_REMAPPER_PLAN_3D( layout_y, layout_x, hat_phi )
    call apply_remap_3D( rmp3, tmp, hat_phi ) 
    do k=1,nz_loc
       do j=1,ny_loc
          call apply_fft_c2c_1d( plan%px_inv, hat_phi(:,j,k), hat_phi(:,j,k) )
       enddo
    enddo

    phi = real(hat_phi, f64)

    SLL_DEALLOCATE_ARRAY(hat_phi, ierr)
    SLL_DEALLOCATE_ARRAY(tmp, ierr)

  end subroutine solve_poisson_3d_periodic_par


  subroutine delete_poisson_3d_periodic_plan_par(plan)

    type (poisson_3d_periodic_plan_par), pointer :: plan
    sll_int32                                    :: ierr

    ! Fixme: some error checking, whether the poisson pointer is associated
    ! for instance
    SLL_ASSERT( associated(plan) )

    call delete_fft_plan1d(plan%px)
    call delete_fft_plan1d(plan%py)
    call delete_fft_plan1d(plan%pz)

    call delete_fft_plan1d(plan%px_inv)
    call delete_fft_plan1d(plan%py_inv)
    call delete_fft_plan1d(plan%pz_inv)

    call delete_layout_3D( plan%layout_x )
    call delete_layout_3D( plan%layout_y )
    call delete_layout_3D( plan%layout_z )

    SLL_DEALLOCATE(plan, ierr)

  end subroutine delete_poisson_3d_periodic_plan_par


  subroutine compute_local_sizes( layout, loc_sz_i, loc_sz_j, loc_sz_k )
    type(layout_3D_t), pointer :: layout
    sll_int32, intent(out) :: loc_sz_i
    sll_int32, intent(out) :: loc_sz_j
    sll_int32, intent(out) :: loc_sz_k
    sll_int32 :: i_min
    sll_int32 :: i_max
    sll_int32 :: j_min
    sll_int32 :: j_max
    sll_int32 :: k_min
    sll_int32 :: k_max
    sll_int32 :: my_rank
    if( .not. associated(layout) ) then
       print *, 'not-associated layout passed to new_distributed_mesh_3D'
       print *, 'Exiting...'
       STOP
    end if
    my_rank = sll_get_collective_rank(get_layout_3D_collective(layout))
    i_min = get_layout_3D_i_min( layout, my_rank )
    i_max = get_layout_3D_i_max( layout, my_rank )
    j_min = get_layout_3D_j_min( layout, my_rank )
    j_max = get_layout_3D_j_max( layout, my_rank )
    k_min = get_layout_3D_k_min( layout, my_rank )
    k_max = get_layout_3D_k_max( layout, my_rank )
    loc_sz_i = i_max - i_min + 1
    loc_sz_j = j_max - j_min + 1
    loc_sz_k = k_max - k_min + 1
  end subroutine compute_local_sizes

  subroutine verify_argument_sizes_par(layout, rho, phi)


    type(layout_3D_t), pointer   :: layout
    sll_real64, dimension(:,:,:) :: rho
    sll_real64, dimension(:,:,:) :: phi
    sll_int32,  dimension(3)     :: n ! nx_loc, ny_loc, nz_loc
    sll_int32                    :: i

    call compute_local_sizes( layout, n(1), n(2), n(3) )

    do i=1,3
       if ( (n(i)/=size(rho,i)) .or. (n(i)/=size(phi,i))  ) then
          print*, 'Input sizes passed to solve_poisson_3d_periodic_par do not match'
          if (i==1) then
             print*, 'Input sizes passed to "solve_poisson_3d_periodic_par" ', &
                  'do not match in direction x'
          elseif (i==2) then
             print*, 'Input sizes passed to "solve_poisson_3d_periodic_par" ', &
                  'do not match in direction y'
          else
             print*, 'Input sizes passed to "solve_poisson_3d_periodic_par" ', &
                  'do not match in direction z'
          endif
          print *, 'Exiting...'
          stop
       endif
    enddo

  end subroutine verify_argument_sizes_par


end module sll_poisson_3d_periodic_par
