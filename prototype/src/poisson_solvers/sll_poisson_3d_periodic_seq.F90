
!*******************************************************************
!
! Selalib 2012     
! Module: sll_poisson_3d_periodic.F90
!
!> @brief 
!> Selalib 3D poisson solver
!> Start date: Feb. 08, 2012
!> Last modification: Feb. 09, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!*******************************************************************

module sll_poisson_3d_periodic_seq

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
#include "sll_assert.h"
  
  use sll_fft
  use numeric_constants

  implicit none

  type poisson_3d_periodic_plan
     sll_int32 :: npx
     sll_int32 :: npy
     sll_int32 :: npz
  end type poisson_3d_periodic_plan

contains


  function new_poisson_3d_periodic_plan()

    type (poisson_3d_periodic_plan), pointer :: new_poisson_3d_periodic_plan

  end function new_poisson_3d_periodic_plan


  subroutine initialize_poisson_3d_periodic_plan(npx, npy, npz, plan)

    sll_int32                                :: npx, npy, npz
    type (poisson_3d_periodic_plan), pointer :: plan

    plan%npx = npx
    plan%npy = npy
    plan%npz = npz

  end subroutine initialize_poisson_3d_periodic_plan


  subroutine solve_poisson_3d_periodic(plan, rho, phi)

    type (poisson_3d_periodic_plan), pointer  :: plan
    sll_real64, dimension(:,:,:)              :: rho, phi
    sll_comp64, dimension(:,:,:), allocatable :: hat_rho, hat_phi
    sll_int32                                 :: mx, my, mz
    sll_int32                                 :: i, j, k, ierr
    type(sll_fft_plan), pointer               :: p

    mx = plan%npx - 1
    my = plan%npy - 1
    mz = plan%npz - 1

    SLL_ALLOCATE(hat_rho(mx,my,mz), ierr)
    SLL_ALLOCATE(hat_phi(mx,my,mz), ierr)

    ! FFTs in x-direction
    hat_rho = cmplx(rho, 0_f64, kind=f64)
    p => new_plan_c2c_1d( mx, hat_rho(:,1,1), hat_rho(:,1,1), FFT_FORWARD )
    do k=1,mz
       do j=1,my
          call apply_fft_c2c_1d( p, hat_rho(:,j,k), hat_rho(:,j,k) )
       enddo
    enddo
    call delete(p)

    ! FFTs in y-direction
    p => new_plan_c2c_1d( my, hat_rho(1,:,1), hat_rho(1,:,1), FFT_FORWARD )
    do k=1,mz
       do i=1,mx
          call apply_fft_c2c_1d( p, hat_rho(i,:,k), hat_rho(i,:,k) )
       enddo
    enddo
    call delete(p)

    ! FFTs in z-direction
    p => new_plan_c2c_1d( mz, hat_rho(1,1,:), hat_rho(1,1,:), FFT_FORWARD )
    do j=1,mx
       do i=1,mx
          call apply_fft_c2c_1d( p, hat_rho(i,j,:), hat_rho(i,j,:) )
       enddo
    enddo
    call delete(p)

    do k=1,mz
       do j=1,my
          do i=1,mx
             hat_phi(i,j,k) = hat_rho(i,j,k) / ( 4*sll_pi * &
             ( ((i-1)/mx)**2 + ((j-1)/my)**2 + ((k-1)/mz)**2 ) )
          enddo
       enddo
    enddo

    ! Inverse FFTs in z-direction
    p => new_plan_c2c_1d( mz, hat_phi(1,1,:), hat_phi(1,1,:), FFT_INVERSE )
    do j=1,mx
       do i=1,mx
          call apply_fft_c2c_1d( p, hat_phi(i,j,:), hat_phi(i,j,:) )
       enddo
    enddo
    call delete(p)

    ! Inverse FFTs in y-direction
    p => new_plan_c2c_1d( my, hat_phi(1,:,1), hat_phi(1,:,1), FFT_INVERSE )
    do k=1,mz
       do i=1,mx
          call apply_fft_c2c_1d( p, hat_phi(i,:,k), hat_phi(i,:,k) )
       enddo
    enddo
    call delete(p)

    ! Inverse FFTs in x-direction
    p => new_plan_c2c_1d( mx, hat_phi(:,1,1), hat_phi(:,1,1), FFT_INVERSE )
    do k=1,mz
       do j=1,my
          call apply_fft_c2c_1d( p, hat_phi(:,j,k), hat_phi(:,j,k) )
       enddo
    enddo
    call delete(p)

    phi = real(hat_phi, f64)

    SLL_DEALLOCATE_ARRAY(hat_rho, ierr)
    SLL_DEALLOCATE_ARRAY(hat_phi, ierr)

  end subroutine solve_poisson_3d_periodic


  subroutine delete_poisson_3d_periodic_plan(poisson)

    type (poisson_3d_periodic_plan), pointer :: poisson

    ! Fixme: some error checking, whether the poisson pointer is associated
    ! for instance
    SLL_ASSERT( associated(poisson) )

  end subroutine delete_poisson_3d_periodic_plan


end module sll_poisson_3d_periodic_seq
