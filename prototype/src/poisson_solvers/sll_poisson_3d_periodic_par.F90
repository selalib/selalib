
!***************************************************************************
!
! Selalib 2012     
! Module: sll_poisson_3d_periodic.F90
!
!> @brief 
!> Selalib periodic 3D poisson solver
!> Start date: Feb. 08, 2012
!> Last modification: March 09, 2012
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

    type (poisson_3d_periodic_plan_par), pointer :: plan
    sll_real64, dimension(:,:,:)                 :: rho
    sll_real64, dimension(:,:,:), allocatable    :: phi
    sll_comp64, dimension(:,:,:), allocatable    :: hat_rho, tmp, hat_phi
    sll_int64                                    :: nx, ny, nz
    sll_int64                                    :: nx_loc, ny_loc, nz_loc
    sll_int32                                    :: e, ex, ey, ez
    sll_int64                                    :: i, j, k
    sll_int32                                    :: ierr
    sll_real64                                   :: Lx, Ly, Lz
    sll_real64                                   :: ind_x, ind_y, ind_z
    sll_int32                                    :: myrank
    sll_int64                                    :: colsz ! collective size
    type(layout_3D_t), pointer                   :: layout_x, layout_y, layout_z, layout_kernel
    type(remap_plan_3D_t), pointer               :: rmp3
    sll_int32, dimension(1:3)                    :: global
    sll_int32                                    :: gi, gj, gk
    sll_int32, dimension(4,3)                    :: loc_sizes ! local sizes in the 4 layouts

    ! Get geometry informations
    nx = plan%nx
    ny = plan%ny
    nz = plan%nz
    Lx = plan%Lx
    Ly = plan%Ly
    Lz = plan%Lz

    ! Get layouts to compute FFTs (in each direction) and poisson solver kernel
    layout_x => plan%layout_x
    layout_y => plan%layout_y
    layout_z => plan%layout_z
    layout_kernel => plan%layout_kernel

    ! Get loc_sizes in the 4 layouts
    loc_sizes = plan%loc_sizes

    if ( (.not.is_power_of_two(nx)) .and. (.not.is_power_of_two(ny)) &
         .and.(.not.is_power_of_two(nz))) then     
       print *, 'This test needs to run in numbers of points which are powers of 2.'
       stop
    end if

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

    ! FFTs in x-direction

    nx_loc = loc_sizes(1,1)
    ny_loc = loc_sizes(1,2)
    nz_loc = loc_sizes(1,3)
    SLL_ALLOCATE(hat_rho(nx_loc,ny_loc,nz_loc), ierr)

    do k=1,nz_loc
       do j=1,ny_loc
          global = local_to_global_3D( layout_x, (/1, int(j), int(k)/))
          gj = global(2)
          gk = global(3)
          hat_rho(:,j,k) = cmplx(rho(:,gj,gk), 0_f64, kind=f64)
          call apply_fft_c2c_1d( plan%px, hat_rho(:,j,k), hat_rho(:,j,k) )
       enddo
    enddo

    ! FFTs in y-direction

    nx_loc = loc_sizes(2,1)
    ny_loc = loc_sizes(2,2)
    nz_loc = loc_sizes(2,3)
    SLL_ALLOCATE(tmp(nx_loc,ny_loc,nz_loc), ierr)

    rmp3 => NEW_REMAPPER_PLAN_3D( layout_x, layout_y, hat_rho )
    call apply_remap_3D( rmp3, hat_rho, tmp) 
    do k=1,nz_loc
       do i=1,nx_loc
          call apply_fft_c2c_1d( plan%py, tmp(i,:,k), tmp(i,:,k) )
       enddo
    enddo

    SLL_DEALLOCATE_ARRAY(hat_rho, ierr)

    ! FFTs in z-direction

    nx_loc = loc_sizes(3,1)
    ny_loc = loc_sizes(3,2)
    nz_loc = loc_sizes(3,3)
    SLL_ALLOCATE(hat_rho(nx_loc,ny_loc,nz_loc), ierr)

    rmp3 => NEW_REMAPPER_PLAN_3D( layout_y, layout_z, tmp )

    call apply_remap_3D( rmp3, tmp, hat_rho) 
    do j=1,ny_loc
       do i=1,nx_loc
          call apply_fft_c2c_1d( plan%pz, hat_rho(i,j,:), hat_rho(i,j,:) )
       enddo
    enddo

    SLL_DEALLOCATE_ARRAY(tmp, ierr)
    hat_rho = hat_rho/(nx*ny*nz)

    ! Compute hat_phi, phi = inv_fft(hat_phi)

    nx_loc = loc_sizes(4,1)
    ny_loc = loc_sizes(4,2)
    nz_loc = loc_sizes(4,3)
    SLL_ALLOCATE(tmp(nx_loc,ny_loc,nz_loc), ierr)

    rmp3 => NEW_REMAPPER_PLAN_3D( layout_z, layout_kernel, hat_rho )
    call apply_remap_3D( rmp3, hat_rho, tmp) 

    SLL_ALLOCATE(hat_phi(nx_loc,ny_loc,nz_loc), ierr)
    do k=1,nz_loc
       do j=1,ny_loc
          do i=1,nx_loc
             global = local_to_global_3D( layout_kernel, (/int(i), int(j), int(k)/))
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
                hat_phi(i,j,k) = tmp(i,j,k)/(4*sll_pi**2*((ind_x/Lx)**2+(ind_y/Ly)**2+(ind_z/Lz)**2))
                ! tmp is hat_rho remapped
             endif
          enddo
       enddo
    enddo

    SLL_DEALLOCATE_ARRAY(tmp, ierr)

    ! Inverse FFTs in x-direction

    nx_loc = loc_sizes(1,1)
    ny_loc = loc_sizes(1,2)
    nz_loc = loc_sizes(1,3)
    SLL_ALLOCATE(tmp(nx_loc,ny_loc,nz_loc), ierr)

    rmp3 => NEW_REMAPPER_PLAN_3D( layout_kernel, layout_x, hat_phi )
    call apply_remap_3D( rmp3, hat_phi, tmp )
    do k=1,nz_loc
       do j=1,ny_loc
          call apply_fft_c2c_1d( plan%px_inv, tmp(:,j,k), tmp(:,j,k) )
       enddo
    enddo
    SLL_DEALLOCATE_ARRAY(hat_phi, ierr)

    ! Inverse FFTs in y-direction

    nx_loc = loc_sizes(2,1)
    ny_loc = loc_sizes(2,2)
    nz_loc = loc_sizes(2,3)
    SLL_ALLOCATE(hat_phi(nx_loc,ny_loc,nz_loc), ierr)

    rmp3 => NEW_REMAPPER_PLAN_3D( layout_x, layout_y, tmp )
    call apply_remap_3D( rmp3, tmp, hat_phi)
    do k=1,nz_loc
       do i=1,nx_loc
          call apply_fft_c2c_1d( plan%py_inv, hat_phi(i,:,k), hat_phi(i,:,k) )
       enddo
    enddo
    SLL_DEALLOCATE_ARRAY(tmp, ierr)

    ! Inverse FFTs in z-direction

    nx_loc = loc_sizes(3,1)
    ny_loc = loc_sizes(3,2)
    nz_loc = loc_sizes(3,3)
    SLL_ALLOCATE(tmp(nx_loc,ny_loc,nz_loc), ierr)

    rmp3 => NEW_REMAPPER_PLAN_3D( layout_y, layout_z, hat_phi )
    call apply_remap_3D( rmp3, hat_phi, tmp ) 
    do j=1,ny_loc
       do i=1,nx_loc
          call apply_fft_c2c_1d( plan%pz_inv, tmp(i,j,:), tmp(i,j,:) )
       enddo
    enddo

    SLL_ALLOCATE(phi(nx_loc,ny_loc,nz_loc), ierr)
    phi = real(tmp, f64)

    SLL_DEALLOCATE_ARRAY(tmp, ierr)

  end subroutine solve_poisson_3d_periodic_par

end module sll_poisson_3d_periodic_par
