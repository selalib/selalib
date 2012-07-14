
!*************************************************************************
!
! Selalib 2012     
! Module: sll_poisson_3d_periodic.F90
!
!> @brief 
!> Selalib 3D poisson solver
!> Start date: Feb. 08, 2012
!> Last modification: April 10, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!*************************************************************************

module sll_poisson_3d_periodic_seq

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
#include "sll_assert.h"

  use sll_fft
  use numeric_constants

  implicit none

  type poisson_3d_periodic_plan_seq
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
  end type poisson_3d_periodic_plan_seq

contains


  function new_poisson_3d_periodic_plan_seq(nx ,ny ,nz, Lx, Ly, Lz) &
                                                         result(plan)

    sll_int32                                    :: nx, ny, nz
    sll_comp64,                    dimension(nx) :: x
    sll_comp64,                    dimension(ny) :: y
    sll_comp64,                    dimension(nz) :: z
    ! nx, ny, nz are the numbers of points - 1 in directions x, y, z
    sll_int32                                    :: ierr
    sll_real64                                   :: Lx, Ly, Lz
    type (poisson_3d_periodic_plan_seq), pointer :: plan

    SLL_ALLOCATE(plan, ierr)

    ! Geometry informations
    plan%nx = nx
    plan%ny = ny
    plan%nz = nz
    plan%Lx = Lx
    plan%Ly = Ly
    plan%Lz = Lz

    ! For FFTs (in each direction)
    plan%px => fft_new_plan( nx, x, x, FFT_FORWARD )
    plan%py => fft_new_plan( ny, y, y, FFT_FORWARD )
    plan%pz => fft_new_plan( nz, z, z, FFT_FORWARD )

    ! For inverse FFTs (in each direction)
    plan%px_inv => fft_new_plan( nx, x, x, FFT_INVERSE )
    plan%py_inv => fft_new_plan( ny, y, y, FFT_INVERSE )
    plan%pz_inv => fft_new_plan( nz, z, z, FFT_INVERSE )

  end function new_poisson_3d_periodic_plan_seq


  subroutine solve_poisson_3d_periodic_seq(plan, rho, phi)

    type (poisson_3d_periodic_plan_seq), pointer   :: plan
    sll_real64, dimension(:,:,:)                   :: rho, phi
    sll_comp64, dimension(plan%nx,plan%ny,plan%nz) :: hat_rho, hat_phi
    sll_int32                                      :: nx, ny, nz
    ! nx, ny, nz are the numbers of points - 1 in directions x, y, z
    sll_int32                                      :: i, j, k
    sll_real64                                     :: Lx, Ly, Lz
    sll_real64                                     :: ind_x, ind_y, ind_z

    ! Checking input arguments consistency
    if ( rho(1,1,1) /= 0.d0 ) then     
       print *, '3D periodic poisson cannot be solved without rho(1,1,1)=0'
       print *, 'Exiting...'
       stop
    endif
    call verify_argument_sizes_seq(plan, rho, phi)

    nx = plan%nx
    ny = plan%ny
    nz = plan%nz
    Lx = plan%Lx
    Ly = plan%Ly
    Lz = plan%Lz

    ! FFTs in x-direction
    hat_rho = cmplx(rho, 0_f64, kind=f64)
    do k=1,nz
       do j=1,ny
          call fft_apply_plan( plan%px, hat_rho(:,j,k), hat_rho(:,j,k) )
       enddo
    enddo

    ! FFTs in y-direction
    do k=1,nz
       do i=1,nx
          call fft_apply_plan( plan%py, hat_rho(i,:,k), hat_rho(i,:,k) )
       enddo
    enddo

    ! FFTs in z-direction
    do j=1,ny
       do i=1,nx
          call fft_apply_plan( plan%pz, hat_rho(i,j,:), hat_rho(i,j,:) )
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
                hat_phi(i,j,k) = hat_rho(i,j,k)/(4*sll_pi**2*((ind_x/Lx)**2 &
                                 + (ind_y/Ly)**2+(ind_z/Lz)**2))
             endif
          enddo
       enddo
    enddo

    ! Inverse FFTs in z-direction
    do j=1,ny
       do i=1,nx
          call fft_apply_plan( plan%pz_inv, hat_phi(i,j,:), hat_phi(i,j,:) )
       enddo
    enddo

    ! Inverse FFTs in y-direction
    do k=1,nz
       do i=1,nx
          call fft_apply_plan( plan%py_inv, hat_phi(i,:,k), hat_phi(i,:,k) )
       enddo
    enddo

    ! Inverse FFTs in x-direction
    do k=1,nz
       do j=1,ny
          call fft_apply_plan( plan%px_inv, hat_phi(:,j,k), hat_phi(:,j,k) )
       enddo
    enddo

    phi = real(hat_phi, f64)

  end subroutine solve_poisson_3d_periodic_seq


  subroutine delete_poisson_3d_periodic_plan_seq(plan)

    type (poisson_3d_periodic_plan_seq), pointer :: plan
    sll_int32                                    :: ierr

    ! Fixme: some error checking, whether the poisson pointer is associated
    ! for instance
    SLL_ASSERT( associated(plan) )

    call fft_delete_plan(plan%px)
    call fft_delete_plan(plan%py)
    call fft_delete_plan(plan%pz)

    call fft_delete_plan(plan%px_inv)
    call fft_delete_plan(plan%py_inv)
    call fft_delete_plan(plan%pz_inv)

    SLL_DEALLOCATE(plan, ierr)

  end subroutine delete_poisson_3d_periodic_plan_seq


  subroutine verify_argument_sizes_seq(plan, rho, phi)

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

  end subroutine verify_argument_sizes_seq


end module sll_poisson_3d_periodic_seq
