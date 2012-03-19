
!***************************************************************************
!
! Selalib 2012     
! Module: sll_poisson_3d_periodic.F90
!
!> @brief 
!> Selalib 3D poisson solver
!> Start date: Feb. 08, 2012
!> Last modification: March 19, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!***************************************************************************

module sll_poisson_3d_periodic_seq

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
#include "sll_assert.h"

  use sll_fft
  use numeric_constants
  use sll_poisson_3d_periodic_util

  implicit none

contains


  subroutine solve_poisson_3d_periodic_seq(plan, rho, phi)

    type (poisson_3d_periodic_plan_seq), pointer :: plan
    sll_real64, dimension(:,:,:)                 :: rho, phi
    sll_comp64, dimension(:,:,:), allocatable    :: hat_rho, hat_phi
    sll_int32                                    :: nx, ny, nz
    sll_int32                                    :: i, j, k, ierr
    sll_real64                                   :: Lx, Ly, Lz
    sll_real64                                   :: ind_x, ind_y, ind_z

    call if_sizes_do_not_match(plan, rho, phi)

    nx = plan%nx
    ny = plan%ny
    nz = plan%nz
    Lx = plan%Lx
    Ly = plan%Ly
    Lz = plan%Lz

    SLL_ALLOCATE(hat_rho(nx,ny,nz), ierr)
    SLL_ALLOCATE(hat_phi(nx,ny,nz), ierr)

    ! FFTs in x-direction
    hat_rho = cmplx(rho, 0_f64, kind=f64)
    do k=1,nz
       do j=1,ny
          call apply_fft_c2c_1d( plan%px, hat_rho(:,j,k), hat_rho(:,j,k) )
       enddo
    enddo       

    ! FFTs in y-direction
    do k=1,nz
       do i=1,nx
          call apply_fft_c2c_1d( plan%py, hat_rho(i,:,k), hat_rho(i,:,k) )
       enddo
    enddo

    ! FFTs in z-direction
    do j=1,ny
       do i=1,nx
          call apply_fft_c2c_1d( plan%pz, hat_rho(i,j,:), hat_rho(i,j,:) )
       enddo
    enddo

    hat_rho = hat_rho/(nx*ny*nz)

    ! Compute hat_phi, phi = inv_fft(hat_phi)
    do k=1,nz
       do j=1,ny
          do i=1,nx
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
                hat_phi(i,j,k) = hat_rho(i,j,k) / (4*sll_pi**2*((ind_x/Lx)**2+(ind_y/Ly)**2+(ind_z/Lz)**2))
             endif
          enddo
       enddo
    enddo

    ! Inverse FFTs in x-direction
    do k=1,nz
       do j=1,ny
          call apply_fft_c2c_1d( plan%px_inv, hat_phi(:,j,k), hat_phi(:,j,k) )
       enddo
    enddo

    ! Inverse FFTs in y-direction
    do k=1,nz
       do i=1,nx
          call apply_fft_c2c_1d( plan%py_inv, hat_phi(i,:,k), hat_phi(i,:,k) )
       enddo
    enddo

    ! Inverse FFTs in z-direction
    do j=1,ny
       do i=1,nx
          call apply_fft_c2c_1d( plan%pz_inv, hat_phi(i,j,:), hat_phi(i,j,:) )
       enddo
    enddo

    phi = real(hat_phi, f64)

    SLL_DEALLOCATE_ARRAY(hat_rho, ierr)
    SLL_DEALLOCATE_ARRAY(hat_phi, ierr)

  end subroutine solve_poisson_3d_periodic_seq


  subroutine if_sizes_do_not_match(plan, rho, phi)

    type (poisson_3d_periodic_plan_seq), pointer :: plan
    sll_real64, dimension(:,:,:)                 :: rho
    sll_real64, dimension(:,:,:)                 :: phi
    sll_int32,  dimension(3)                     :: n ! nx_loc, ny_loc, nz_loc
    sll_int32                                    :: i

    n(1) = plan%nx
    n(2) = plan%ny
    n(3) = plan%nz

    do i=1,3
       if ( (n(i)/=size(rho,i)) .or. (n(i)/=size(phi,i))  ) then
          print*, 'Input sizes passed to solve_poisson_3d_periodic_par do not match'
          print *, 'Exiting...'
          stop
       endif
    enddo

  end subroutine if_sizes_do_not_match


end module sll_poisson_3d_periodic_seq
