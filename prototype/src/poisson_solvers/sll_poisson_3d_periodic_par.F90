
!***************************************************************************
!
! Selalib 2012     
! Module: sll_poisson_3d_periodic.F90
!
!> @brief 
!> Selalib 3D poisson solver
!> Start date: Feb. 08, 2012
!> Last modification: Feb. 23, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!***************************************************************************

module sll_poisson_3d_periodic_par

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
#include "sll_assert.h"
#include "sll_remap.h"

  use sll_fft
  use numeric_constants
  use sll_collective
  use sll_poisson_3d_periodic_util

  implicit none

contains

  subroutine solve_poisson_3d_periodic_par(plan, rho, phi)

    type (poisson_3d_periodic_plan), pointer  :: plan
    sll_real64, dimension(:,:,:)              :: rho, phi
    sll_comp64, dimension(:,:,:), allocatable :: hat_rho, tmp, hat_phi
    sll_int64                                 :: nx, ny, nz
    sll_int32                                 :: npx, npy, npz
    sll_int32                                 :: e, ex, ey, ez
    sll_int64                                 :: i, j, k
    sll_int32                                 :: ierr
    sll_real64                                :: Lx, Ly, Lz
    sll_real64                                :: ind_x, ind_y, ind_z
    sll_int32                                 :: myrank
    sll_int64                                 :: colsz ! collective size
    type(layout_3D_t), pointer                :: layout1
    type(layout_3D_t), pointer                :: layout2
    type(remap_plan_3D_t), pointer            :: rmp3
    sll_int32, dimension(1:3)                 :: global
    sll_int32                                 :: gi, gj

    nx = plan%nx
    ny = plan%ny
    nz = plan%nz

    if ( (.not.is_power_of_two(nx)) .and. (.not.is_power_of_two(ny)) &
         .and.(.not.is_power_of_two(nz))) then     
       print *, 'This test needs to run in numbers of points which are powers of 2.'
       stop
    end if

    Lx = plan%Lx
    Ly = plan%Ly
    Lz = plan%Lz

    colsz  = sll_get_collective_size(sll_world_collective)
    myrank = sll_get_collective_rank(sll_world_collective)

    if ( colsz > min(nx,ny,nz) ) then     
       print *, 'This test needs to run in a number of processes which is less than', min(nx,ny,nz)
       stop
    end if

    if ( (.not.is_power_of_two(nx)) .and. (.not.is_power_of_two(ny)) &
         .and.(.not.is_power_of_two(nz))) then     
       print *, 'This test needs to run in numbers of points which are powers of 2.'
       stop
    end if

    ! FFTs in z-direction

    e = int(log(real(colsz))/log(2.))
    npx = 2**(e/2)
    npy = int(colsz)/npx
    npz = 1
    layout1  => new_layout_3D( sll_world_collective ) 
    call initialize_layout_with_distributed_3D_array( int(nx), int(ny), int(nz), npx, npy, npz, layout1 )

    SLL_ALLOCATE(hat_rho(nx/npx,ny/npy,nz/npz), ierr)
    do j=1,ny/npy
       do i=1,nx/npx
       global = local_to_global_3D( layout1, (/int(i), int(j), 1/))
       gi =global(1)
       gj =global(2)
       hat_rho(i,j,:) = cmplx(rho(gi,gj,:), 0_f64, kind=f64)
          call apply_fft_c2c_1d( plan%pz, hat_rho(i,j,:), hat_rho(i,j,:) )
       enddo
    enddo       

    ! FFTs in y-direction

    npz = npy
    npy = 1
    layout2  => new_layout_3D( sll_world_collective )
    call initialize_layout_with_distributed_3D_array( int(nx), int(ny), int(nz), npx, npy, npz, layout2 )

    SLL_ALLOCATE(tmp(nx/npx,ny/npy,nz/npz), ierr)
    rmp3 => NEW_REMAPPER_PLAN_3D( layout1, layout2, hat_rho )
    call apply_remap_3D( rmp3, hat_rho, tmp) 

    do k=1,nz/npz
       do i=1,nx/npx
          call apply_fft_c2c_1d( plan%py, tmp(i,:,k), tmp(i,:,k) )
       enddo
    enddo

    SLL_DEALLOCATE_ARRAY(hat_rho, ierr)

    ! FFTs in x-direction

    npy = npx
    npx = 1
    layout1  => new_layout_3D( sll_world_collective )
    call initialize_layout_with_distributed_3D_array( int(nx), int(ny), int(nz), npx, npy, npz, layout1 )

    SLL_ALLOCATE(hat_rho(nx/npx,ny/npy,nz/npz), ierr)
    rmp3 => NEW_REMAPPER_PLAN_3D( layout2, layout1, tmp )
    call apply_remap_3D( rmp3, tmp, hat_rho) 

    do k=1,nz/npz
       do j=1,ny/npy
          call apply_fft_c2c_1d( plan%px, hat_rho(:,j,k), hat_rho(:,j,k) )
       enddo
    enddo

    SLL_DEALLOCATE_ARRAY(tmp, ierr)
    hat_rho = hat_rho/(nx*ny*nz)

    ! Compute hat_phi, phi = inv_fft(hat_phi)

    ex = e/3
    ey = (e-ex)/2
    ez = e - (ex+ey)
    npx = 2**ex
    npy = 2**ey
    npz = 2**ez
    layout2  => new_layout_3D( sll_world_collective ) 
    call initialize_layout_with_distributed_3D_array( int(nx), int(ny), int(nz), npx, npy, npz, layout2 )

    SLL_ALLOCATE(tmp(nx/npx,ny/npy,nz/npz), ierr)
    rmp3 => NEW_REMAPPER_PLAN_3D( layout1, layout2, hat_rho )
    call apply_remap_3D( rmp3, hat_rho, tmp) 

    SLL_ALLOCATE(hat_phi(nx/npx,ny/npy,nz/npz), ierr)
    do k=1,nz/npz
       do j=1,ny/npy
          do i=1,nx/npx
             if (i<=nx/2) then
                ind_x = real(i-1,f64)
             else
                ind_x = real(nx-(i-1),f64)
             endif
             if (j<=ny/2) then
                ind_y = real(j-1,f64)
             else
                ind_y = real(ny-(j-1),f64)
             endif
             if (k<=nz/2) then
                ind_z = real(k-1,f64)
             else
                ind_z = real(nz-(k-1),f64)
             endif                
             if ( (ind_x==0) .and. (ind_y==0) .and. (ind_z==0) ) then
                hat_phi(i,j,k) = 0.d0
             else
                hat_phi(i,j,k) = tmp(i,j,k) / (4*sll_pi**2*((ind_x/Lx)**2+(ind_y/Ly)**2+(ind_z/Lz)**2))
                ! tmp is hat_rho remapped
             endif
          enddo
       enddo
    enddo

    SLL_DEALLOCATE_ARRAY(tmp, ierr)

    ! Inverse FFTs in z-direction

    npx = 2**(e/2)
    npy = int(colsz)/npx
    npz = 1
    layout1  => new_layout_3D( sll_world_collective ) 
    call initialize_layout_with_distributed_3D_array( int(nx), int(ny), int(nz), npx, npy, npz, layout1 )

    SLL_ALLOCATE(tmp(nx/npx,ny/npy,nz/npz), ierr)
    rmp3 => NEW_REMAPPER_PLAN_3D( layout2, layout1, hat_phi )
    call apply_remap_3D( rmp3, hat_phi, tmp ) 

    do j=1,ny/npy
       do i=1,nx/npx
          call apply_fft_c2c_1d( plan%pz_inv, tmp(i,j,:), tmp(i,j,:) )
       enddo
    enddo

    SLL_DEALLOCATE_ARRAY(hat_phi, ierr)

    ! Inverse FFTs in y-direction

    npz = npy
    npy = 1
    layout2  => new_layout_3D( sll_world_collective )
    call initialize_layout_with_distributed_3D_array( int(nx), int(ny), int(nz), npx, npy, npz, layout2 )

    SLL_ALLOCATE(hat_phi(nx/npx,ny/npy,nz/npz), ierr)
    rmp3 => NEW_REMAPPER_PLAN_3D( layout1, layout2, tmp )
    call apply_remap_3D( rmp3, tmp, hat_phi)

    do k=1,nz/npz
       do i=1,nx/npx
          call apply_fft_c2c_1d( plan%py_inv, hat_phi(i,:,k), hat_phi(i,:,k) )
       enddo
    enddo

    SLL_DEALLOCATE_ARRAY(tmp, ierr)

    ! Inverse FFTs in x-direction

    npy = npx
    npx = 1
    layout1  => new_layout_3D( sll_world_collective )
    call initialize_layout_with_distributed_3D_array( int(nx), int(ny), int(nz), npx, npy, npz, layout1 )

    SLL_ALLOCATE(tmp(nx/npx,ny/npy,nz/npz), ierr)
    rmp3 => NEW_REMAPPER_PLAN_3D( layout2, layout1, hat_phi )
    call apply_remap_3D( rmp3, hat_phi, tmp ) 

    do k=1,nz/npz
       do j=1,ny/npy
          call apply_fft_c2c_1d( plan%px_inv, tmp(:,j,k), tmp(:,j,k) )
       enddo
    enddo

    phi = real(tmp, f64)

    SLL_DEALLOCATE_ARRAY(tmp, ierr)
    call delete_layout_3D( layout1 )
    call delete_layout_3D( layout2 )

  end subroutine solve_poisson_3d_periodic_par

end module sll_poisson_3d_periodic_par
