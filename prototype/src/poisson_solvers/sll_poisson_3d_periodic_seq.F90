
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
     type(sll_fft_plan), pointer :: px
     type(sll_fft_plan), pointer :: py
     type(sll_fft_plan), pointer :: pz
     type(sll_fft_plan), pointer :: px_inv
     type(sll_fft_plan), pointer :: py_inv
     type(sll_fft_plan), pointer :: pz_inv
  end type poisson_3d_periodic_plan

contains


  function new_poisson_3d_periodic_plan(array)

    sll_comp64, dimension(:,:,:)             :: array
    sll_int32                                :: nx, ny, nz
    type (poisson_3d_periodic_plan), pointer :: new_poisson_3d_periodic_plan
    sll_int32                                :: ierr

    nx = size(array,1)
    ny = size(array,2)
    nz = size(array,3)

    SLL_ALLOCATE(new_poisson_3d_periodic_plan, ierr)

    new_poisson_3d_periodic_plan%nx = nx
    new_poisson_3d_periodic_plan%ny = ny
    new_poisson_3d_periodic_plan%nz = nz

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

    nx = plan%nx
    ny = plan%ny
    nz = plan%nz

    SLL_ALLOCATE(hat_rho(nx,ny,nz), ierr)
    SLL_ALLOCATE(hat_phi(nx,ny,nz), ierr)

    ! FFTs in z-direction
 !call random_number(rho)
    hat_rho = cmplx(rho, 0_f64, kind=f64)
!print*, 'test1:',hat_rho(1,1,:)
    do j=1,ny
       do i=1,nx
          call apply_fft_c2c_1d( plan%pz, hat_rho(i,j,:), hat_rho(i,j,:) )
       enddo
    enddo
!print*, 'test2',hat_rho(1,1,:)
!stop
    do k=1,nz-1,2
       hat_rho(:,:,k) = -hat_rho(:,:,k)
    enddo

    ! FFTs in y-direction
    do k=1,nz
       do i=1,nx
          call apply_fft_c2c_1d( plan%py, hat_rho(i,:,k), hat_rho(i,:,k) )
       enddo
    enddo

    do j=1,ny-1,2
       hat_rho(:,j,:) = -hat_rho(:,j,:)
    enddo

    ! FFTs in x-direction
    do k=1,nz
       do j=1,ny
          call apply_fft_c2c_1d( plan%px, hat_rho(:,j,k), hat_rho(:,j,k) )
       enddo
    enddo
    do i=1,nx-1,2
       hat_rho(i,:,:) = -hat_rho(i,:,:)
    enddo

    hat_rho = hat_rho/(nx*ny*nz)

    ! Compute hat_phi, phi = inv_fft(hat_phi)
    do k=-nz/2,nz/2-1
       do j=-ny/2,ny/2-1
          do i=-nx/2,nx/2-1
             if ( k==0 ) then
                hat_phi(i+nx/2+1,j+ny/2+1,k+nz/2+1) = 0.d0
             else
                hat_phi(i+nx/2+1,j+ny/2+1,k+nz/2+1) = hat_rho(i+nx/2+1,j+ny/2+1,k+nz/2+1) / ( &
                  4*sll_pi**2*( (real(i,f64)/nx)**2 + (real(j,f64)/ny)**2 + (real(k,f64)/nz)**2 ) )
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
    do k=1,nz-1,2
       hat_phi(:,:,k) = -hat_phi(:,:,k)
    enddo

    ! Inverse FFTs in y-direction
    do k=1,nz
       do i=1,nx
          call apply_fft_c2c_1d( plan%py_inv, hat_phi(i,:,k), hat_phi(i,:,k) )
       enddo
    enddo
    do j=1,ny-1,2
       hat_phi(:,j,:) = -hat_phi(:,j,:)
    enddo

    ! Inverse FFTs in x-direction
    do k=1,nz
       do j=1,ny
          call apply_fft_c2c_1d( plan%px_inv, hat_phi(:,j,k), hat_phi(:,j,k) )
       enddo
    enddo
    do i=1,nx-1,2
       hat_phi(i,:,:) = -hat_phi(i,:,:)
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
