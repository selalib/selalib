!***************************************************************************
!
! Selalib 2012     
! Module: sll_qns_2d_with_finite_diff.F90
!
!> @brief 
!> Selalib 2D (r, theta) quasi-neutral solver with finite differences
!> Some arrays are here in 3D for remap utilities
!> Start date: March 13, 2012
!> Last modification: March 22, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!***************************************************************************

module sll_qns_2d_with_finite_diff_par

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
#include "sll_assert.h"
#include "sll_remap.h"

  use numeric_constants
  use sll_fft
  use sll_tridiagonal
  use sll_collective

  implicit none


  type qns_2d_with_finite_diff_plan_par
     character(len=100)                      :: bc ! Boundary_conditions
     sll_int32                               :: nr ! Number of points in r-direction
     sll_int32                               :: ntheta ! Number of points in theta-direction
     sll_real64                              :: rmin
     sll_real64                              :: rmax
     sll_real64, dimension(:,:), allocatable :: rho
     sll_real64, dimension(:),   allocatable :: c  
     sll_real64, dimension(:),   allocatable :: Te
     sll_real64, dimension(:),   allocatable :: f
     sll_real64, dimension(:),   allocatable :: g    
     sll_real64                              :: Zi
     type(sll_fft_plan), pointer             :: fft_plan
     type(sll_fft_plan), pointer             :: inv_fft_plan
     type(layout_3D_t),  pointer             :: layout_fft
     type(layout_3D_t),  pointer             :: layout_lin_syst
  end type qns_2d_with_finite_diff_plan_par

contains


  function new_qns_2d_with_finite_diff_plan_par(start_layout, bc, rmin, rmax, rho, c, Te, f, g, Zi) result (plan)

    type(layout_3D_t), pointer                      :: start_layout
    character(len=100)                              :: bc ! Boundary_conditions
    sll_real64                                      :: rmin
    sll_real64                                      :: rmax
    sll_real64, dimension(:,:)                      :: rho
    sll_real64, dimension(:)                        :: c, Te, f, g    
    sll_real64                                      :: Zi
    sll_comp64, dimension(:),   allocatable         :: x
    sll_int32                                       :: nr, ntheta, ierr
    sll_int64                                       :: colsz
    type(qns_2d_with_finite_diff_plan_par), pointer :: plan

    nr = size(rho,1)
    ntheta = size(rho,2)

    SLL_ALLOCATE(plan, ierr)
    SLL_ALLOCATE(plan%rho(nr,ntheta), ierr)
    SLL_ALLOCATE(plan%c(nr), ierr)
    SLL_ALLOCATE(plan%Te(nr), ierr)
    SLL_ALLOCATE(plan%f(ntheta), ierr)
    SLL_ALLOCATE(plan%g(ntheta), ierr)
    SLL_ALLOCATE( x(ntheta), ierr )

    plan%bc     = bc
    plan%nr     = nr
    plan%ntheta = ntheta
    plan%rmin   = rmin
    plan%rmax   = rmax
    plan%rho    = rho
    plan%c      = c
    plan%Te     = Te
    plan%f      = f
    plan%g      = g
    plan%Zi     = Zi 

    ! For FFTs in theta-direction
    plan%fft_plan => new_plan_c2c_1d( ntheta, x, x, FFT_FORWARD )

    ! For inverse FFTs in theta-direction
    plan%inv_fft_plan => new_plan_c2c_1d( ntheta, x, x, FFT_INVERSE )

    ! Layout for FFTs-Inv_FFT in theta-direction
    plan%layout_fft => new_layout_3D( sll_world_collective )
    colsz  = sll_get_collective_size(sll_world_collective)
    call initialize_layout_with_distributed_3D_array( nr, ntheta, 1, int(colsz), 1, 1, plan%layout_fft )

    ! Layout for Linear systems in r-direction
    plan%layout_lin_syst => new_layout_3D( sll_world_collective )
    call initialize_layout_with_distributed_3D_array( nr, ntheta, 1, 1, int(colsz), 1, plan%layout_lin_syst )

    SLL_DEALLOCATE_ARRAY( x, ierr )

  end function new_qns_2d_with_finite_diff_plan_par


 subroutine solve_qn_2d_with_finite_diff_par(plan, phi)

    type(qns_2d_with_finite_diff_plan_par), pointer :: plan
    sll_real64, dimension(:,:)                      :: phi

  end subroutine solve_qn_2d_with_finite_diff_par


  subroutine delete_new_qns_2d_with_finite_diff_plan_par(plan)

       type (qns_2d_with_finite_diff_plan_par), pointer :: plan
       sll_int32                                        :: ierr

       ! Fixme: some error checking, whether the poisson pointer is associated
       ! for instance
       SLL_ASSERT( associated(plan) )

       call delete_fft_plan1d(plan%fft_plan)
       call delete_fft_plan1d(plan%inv_fft_plan)

       call delete_layout_3D( plan%layout_fft )
       call delete_layout_3D( plan%layout_lin_syst )

       SLL_DEALLOCATE_ARRAY(plan%rho, ierr)
       SLL_DEALLOCATE_ARRAY(plan%c, ierr)
       SLL_DEALLOCATE_ARRAY(plan%Te, ierr)
       SLL_DEALLOCATE_ARRAY(plan%f, ierr)
       SLL_DEALLOCATE_ARRAY(plan%g, ierr)
       SLL_DEALLOCATE(plan, ierr)

  end subroutine delete_new_qns_2d_with_finite_diff_plan_par


  subroutine dirichlet_matrix_resh_par(plan, j, a_resh)

    type(qns_2d_with_finite_diff_plan_par), pointer :: plan ! Matrix is sequential
    sll_real64                                      :: dr, dtheta, r, rmin, rmax, Zi
    sll_int32                                       :: i, j, nr, ierr
    sll_real64, dimension(:)                        :: a_resh
    sll_real64, dimension(:), allocatable           :: c, Te ! C & Te are the vector of the Cr & Te(i) respectively

    nr = plan%nr
    rmin = plan%rmin
    rmax = plan%rmax
    dr = (rmax-rmin)/(nr+1)
    dtheta = 2*sll_pi / plan%ntheta
    SLL_ALLOCATE( c(nr), ierr )
    c = plan%c
    SLL_ALLOCATE( Te(nr), ierr )
    Te = plan%Te
    Zi = plan%Zi        

    a_resh = 0.d0

    do i=1,nr
       r = rmin + i*dr
       if (i>1) a_resh(3*(i-1)+1) = -(1/dr**2 - c(i)/(2*dr))
       a_resh(3*(i-1)+2) = 2/dr**2 + 2/(r*dtheta)**2 + 1/(Zi*Te(i)) - 2*cos(j*dtheta)/(r*dtheta)**2
       if (i<nr) a_resh(3*(i-1)+3) = -( 1/dr**2 + c(i)/(2*dr) )
    enddo

    SLL_DEALLOCATE_ARRAY( c, ierr )
    SLL_DEALLOCATE_ARRAY( Te, ierr )

  end subroutine dirichlet_matrix_resh_par


  subroutine neumann_matrix_resh_par(plan, j, a_resh)
    type(qns_2d_with_finite_diff_plan_par), pointer :: plan ! Matrix is sequential
    sll_real64                                      :: dr, dtheta, Zi, rmin, rmax, r
    sll_int32                                       :: i, j, nr, ierr
    sll_real64, dimension(:)                        :: a_resh
    sll_real64, dimension(:), allocatable           :: c, Te 
    ! c & Te are the vector of the Cr & Te(i) respectively

    nr = plan%nr
    rmin = plan%rmin
    rmax = plan%rmax
    dr = (rmax-rmin)/(nr-1)
    dtheta = 2*sll_pi / plan%ntheta
    SLL_ALLOCATE( c(nr), ierr )
    c = plan%c
    SLL_ALLOCATE( Te(nr), ierr )
    Te = plan%Te
    Zi = plan%Zi

    a_resh = 0.d0

    a_resh(2) = 2/dr**2 + 2/(rmin*dtheta)**2 * (-cos(j*dtheta)+1) + 1/(Zi*Te(1))
    a_resh(3) = -2/dr**2

    do i=2,nr-1
       r = rmin + (i-1)*dr
       a_resh(3*(i-1)+1) = -(1/dr**2 - c(i)/(2*dr))
       a_resh(3*(i-1)+2) = 2/dr**2 + 2/(r*dtheta)**2 + 1/(Zi*Te(i)) - 2*cos(j*dtheta)/(r*dtheta)**2
       a_resh(3*(i-1)+3) = -( 1/dr**2 + c(i)/(2*dr) )
    enddo

    a_resh(3*(nr-1)+1) = -2/dr**2
    a_resh(3*(nr-1)+2) = 2/dr**2 + 2/(rmax*dtheta)**2 * (-cos(j*dtheta)+1) + 1/(Zi*Te(nr))

    SLL_DEALLOCATE_ARRAY( c, ierr )
    SLL_DEALLOCATE_ARRAY( Te, ierr )

  end subroutine neumann_matrix_resh_par


end module sll_qns_2d_with_finite_diff_par
