
!***************************************************************************
!
! Selalib 2012     
! Module: sll_poisson_3d_periodic.F90
!
!> @brief 
!> Selalib 3D poisson solver
!> Start date: Feb. 08, 2012
!> Last modification: Feb. 14, 2012
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

  implicit none

  type poisson_3d_periodic_plan
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
  end type poisson_3d_periodic_plan

contains


  function new_poisson_3d_periodic_plan(array, Lx, Ly, Lz)

    sll_comp64, dimension(:,:,:)             :: array
    sll_int32                                :: nx, ny, nz
    sll_int32                                 :: i, j, k, ierr
    sll_real64                                :: Lx, Ly, Lz
    type (poisson_3d_periodic_plan), pointer :: new_poisson_3d_periodic_plan

    nx = size(array,1)
    ny = size(array,2)
    nz = size(array,3)

    SLL_ALLOCATE(new_poisson_3d_periodic_plan, ierr)

    new_poisson_3d_periodic_plan%nx = nx
    new_poisson_3d_periodic_plan%ny = ny
    new_poisson_3d_periodic_plan%nz = nz
    new_poisson_3d_periodic_plan%Lx = Lx
    new_poisson_3d_periodic_plan%Ly = Ly
    new_poisson_3d_periodic_plan%Lz = Lz

    new_poisson_3d_periodic_plan%px => new_plan_c2c_1d( nx, array(:,1,1), array(:,1,1), FFT_FORWARD )
    new_poisson_3d_periodic_plan%py => new_plan_c2c_1d( ny, array(1,:,1), array(1,:,1), FFT_FORWARD )
    new_poisson_3d_periodic_plan%pz => new_plan_c2c_1d( nz, array(1,1,:), array(1,1,:), FFT_FORWARD )
    new_poisson_3d_periodic_plan%px_inv => new_plan_c2c_1d( nx, array(:,1,1), array(:,1,1), FFT_INVERSE)
    new_poisson_3d_periodic_plan%py_inv => new_plan_c2c_1d( ny, array(1,:,1), array(1,:,1), FFT_INVERSE )
    new_poisson_3d_periodic_plan%pz_inv => new_plan_c2c_1d( nz, array(1,1,:), array(1,1,:), FFT_INVERSE )

  end function new_poisson_3d_periodic_plan


  subroutine solve_poisson_3d_periodic(plan, rho, phi)

    type (poisson_3d_periodic_plan), pointer  :: plan
    sll_real64, dimension(:,:,:)              :: rho, phi
    sll_comp64, dimension(:,:,:), allocatable :: hat_rho, hat_phi
    sll_int32                                 :: nx, ny, nz
    sll_int32                                 :: i, j, k, ierr
    sll_real64                                :: Lx, Ly, Lz
    sll_real64                                :: ind_x, ind_y, ind_z

    nx = plan%nx
    ny = plan%ny
    nz = plan%nz
    Lx = plan%Lx
    Ly = plan%Ly
    Lz = plan%Lz

    SLL_ALLOCATE(hat_rho(nx,ny,nz), ierr)
    SLL_ALLOCATE(hat_phi(nx,ny,nz), ierr)

    ! FFTs in z-direction
    hat_rho = cmplx(rho, 0_f64, kind=f64)
    do j=1,ny
       do i=1,nx
          call apply_fft_c2c_1d( plan%pz, hat_rho(i,j,:), hat_rho(i,j,:) )
       enddo
    enddo       

    ! FFTs in y-direction
    do k=1,nz
       do i=1,nx
          call apply_fft_c2c_1d( plan%py, hat_rho(i,:,k), hat_rho(i,:,k) )
       enddo
    enddo

    ! FFTs in x-direction
    do k=1,nz
       do j=1,ny
          call apply_fft_c2c_1d( plan%px, hat_rho(:,j,k), hat_rho(:,j,k) )
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

    ! Inverse FFTs in z-direction
    do j=1,ny
       do i=1,nx
          call apply_fft_c2c_1d( plan%pz_inv, hat_phi(i,j,:), hat_phi(i,j,:) )
       enddo
    enddo

    ! Inverse FFTs in y-direction
    do k=1,nz
       do i=1,nx
          call apply_fft_c2c_1d( plan%py_inv, hat_phi(i,:,k), hat_phi(i,:,k) )
       enddo
    enddo

    ! Inverse FFTs in x-direction
    do k=1,nz
       do j=1,ny
          call apply_fft_c2c_1d( plan%px_inv, hat_phi(:,j,k), hat_phi(:,j,k) )
       enddo
    enddo

    phi = real(hat_phi, f64) !Inverse FFTs are not normalized

    SLL_DEALLOCATE_ARRAY(hat_rho, ierr)
    SLL_DEALLOCATE_ARRAY(hat_phi, ierr)

  end subroutine solve_poisson_3d_periodic


  subroutine delete_poisson_3d_periodic_plan(plan)

    type (poisson_3d_periodic_plan), pointer :: plan
    sll_int32                                :: ierr

    ! Fixme: some error checking, whether the poisson pointer is associated
    ! for instance
    SLL_ASSERT( associated(plan) )

    call delete(plan%px)
    call delete(plan%py)
    call delete(plan%pz)

    call delete(plan%px_inv)
    call delete(plan%py_inv)
    call delete(plan%pz_inv)

    SLL_DEALLOCATE(plan, ierr)

  end subroutine delete_poisson_3d_periodic_plan


end module sll_poisson_3d_periodic_seq
