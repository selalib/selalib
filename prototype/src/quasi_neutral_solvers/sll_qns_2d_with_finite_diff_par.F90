!***************************************************************************
!
! Selalib 2012     
! Module: sll_qns_2d_with_finite_diff.F90
!
!> @brief 
!> Selalib 2D (r, theta) quasi-neutral solver with finite differences
!> Some arrays are here in 3D for remap utilities
!> Start date: March 13, 2012
!> Last modification: March 23, 2012
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
     character(len=100)                        :: bc!Boundary_conditions
     sll_int32                                 :: nr!Number of points in r-direction
     sll_int32                                 :: ntheta!Number of points in theta-direction
     sll_real64                                :: rmin
     sll_real64                                :: rmax
     sll_real64, dimension(:,:), allocatable   :: rho
     sll_real64, dimension(:),   allocatable   :: c  
     sll_real64, dimension(:),   allocatable   :: Te
     sll_real64, dimension(:),   allocatable   :: f
     sll_real64, dimension(:),   allocatable   :: g    
     sll_real64                                :: Zi
     type(sll_fft_plan), pointer               :: fft_plan
     type(sll_fft_plan), pointer               :: inv_fft_plan
     type(layout_3D_t),  pointer               :: layout_fft
     type(layout_3D_t),  pointer               :: layout_lin_sys
     sll_comp64, dimension(:,:,:), allocatable :: array_fft
     sll_comp64, dimension(:,:,:), allocatable :: array_lin_sys
     type(remap_plan_3D_t), pointer            :: rmp3_1
     type(remap_plan_3D_t), pointer            :: rmp3_2
  end type qns_2d_with_finite_diff_plan_par

contains


  function new_qns_2d_with_finite_diff_plan_par(start_layout, bc, rmin, rmax, rho, c, &
                                                            Te, f, g, Zi) result (plan)

    type(layout_3D_t), pointer                      :: start_layout
    character(len=100)                              :: bc ! Boundary_conditions
    sll_real64                                      :: rmin
    sll_real64                                      :: rmax
    sll_real64, dimension(:,:)                      :: rho
    sll_real64, dimension(:)                        :: c, Te, f, g    
    sll_real64                                      :: Zi
    sll_comp64, dimension(:),   allocatable         :: x
    sll_int32                                       :: nr, ntheta
    sll_int32                                       :: nr_loc, ntheta_loc
    sll_int32                                       :: ierr
    sll_int64                                       :: colsz
    type(qns_2d_with_finite_diff_plan_par), pointer :: plan

    colsz  = sll_get_collective_size(sll_world_collective)
    nr_loc = size(rho,1)
    nr     = nr_loc*int(colsz)
    ntheta = size(rho,2)
    ntheta_loc = ntheta/int(colsz)

    SLL_ALLOCATE(plan, ierr)
    SLL_ALLOCATE(plan%rho(nr_loc,ntheta_loc), ierr)
    SLL_ALLOCATE(plan%c(nr_loc), ierr)
    SLL_ALLOCATE(plan%Te(nr_loc), ierr)
    SLL_ALLOCATE(plan%f(ntheta_loc), ierr)
    SLL_ALLOCATE(plan%g(ntheta_loc), ierr)
    SLL_ALLOCATE( x(ntheta_loc), ierr )

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
    call initialize_layout_with_distributed_3D_array( nr, ntheta, 1, &
                                   int(colsz), 1, 1, plan%layout_fft )

    ! Layout for Linear systems in r-direction
    plan%layout_lin_sys => new_layout_3D( sll_world_collective )
    call initialize_layout_with_distributed_3D_array( nr, ntheta, 1, &
                               1, int(colsz), 1, plan%layout_lin_sys )

    SLL_ALLOCATE(plan%array_fft(nr_loc,ntheta,1), ierr)
    SLL_ALLOCATE(plan%array_lin_sys(nr,ntheta_loc,1), ierr)

    plan%rmp3_1 => NEW_REMAPPER_PLAN_3D(plan%layout_fft, plan%layout_lin_sys,plan%array_fft)
    plan%rmp3_2 => NEW_REMAPPER_PLAN_3D(plan%layout_lin_sys,plan%layout_fft,plan%array_lin_sys)

    SLL_DEALLOCATE_ARRAY( x, ierr )

  end function new_qns_2d_with_finite_diff_plan_par


 subroutine solve_qn_2d_with_finite_diff_par(plan, phi)

    type(qns_2d_with_finite_diff_plan_par), pointer :: plan
    sll_real64                                      :: dr, dtheta 
    sll_int32                                       :: nr, ntheta
    sll_int32                                       :: nr_loc, ntheta_loc
    sll_int32                                       :: i, j, ierr
    sll_real64, dimension(:,:)                      :: phi
    sll_comp64, dimension(plan%ntheta)              :: f, g
    sll_int32, dimension(plan%nr)                   :: ipiv
    sll_real64, dimension(3*plan%nr)                :: a_resh ! 3*n
    sll_real64, dimension(7*plan%nr)                :: cts ! 7*n allocation
    sll_int64                                       :: colsz ! collective size
    sll_int32, dimension(1:3)                       :: global
    sll_int32                                       :: ind

    colsz  = sll_get_collective_size(sll_world_collective)

    nr         = plan%nr
    ntheta     = plan%ntheta
    nr_loc     = nr/int(colsz)
    ntheta_loc = ntheta/int(colsz)

    if (plan%bc=='neumann') then
       dr = (plan%rmax-plan%rmin)/(nr-1)
    else ! 'dirichlet'
       dr = (plan%rmax-plan%rmin)/(nr+1)
    endif
    dtheta = 2*sll_pi/ntheta

    ! FFTs (in theta-direction)

    plan%array_fft(:,:,1) = cmplx(plan%rho, 0_f64, kind=f64)
    f = cmplx(plan%f, 0_f64, kind=f64)
    g = cmplx(plan%g, 0_f64, kind=f64)

    call apply_fft_c2c_1d( plan%fft_plan, f, f )
    call apply_fft_c2c_1d( plan%fft_plan, g, g )

    do i=1,nr_loc
       call apply_fft_c2c_1d( plan%fft_plan, plan%array_fft(i,:,1), plan%array_fft(i,:,1) )
       global = local_to_global_3D( plan%layout_fft, (/i, 1, 1/))
       ind = global(1)
       if (ind==1) then
          if (plan%bc=='neumann') then
             plan%array_fft(i,:,1) = plan%array_fft(i,:,1) + (plan%c(ind)-2/dr)*f 
          else ! 'dirichlet'
             plan%array_fft(i,:,1) = plan%array_fft(i,:,1) + (1/dr**2 - plan%c(ind)/(2*dr))*f
          endif
       elseif(ind==nr) then
          if (plan%bc=='neumann') then
             plan%array_fft(i,:,1) = plan%array_fft(i,:,1) + (plan%c(ind)+2/dr)*g
          else ! 'dirichlet'
             plan%array_fft(i,:,1) = plan%array_fft(i,:,1) + (1/dr**2 + plan%c(ind)/(2*dr))*g
          endif
       endif 
    enddo

    ! Remapping to solve linear systems
    call apply_remap_3D( plan%rmp3_1, plan%array_fft, plan%array_lin_sys ) 

    ! Solve linear systems (r-direction)
    do j=1,ntheta_loc
       global = local_to_global_3D( plan%layout_lin_sys, (/1, j, 1/))
       ind = global(2)
       if (plan%bc=='neumann') then
          call neumann_matrix_resh_par(plan, ind-1, a_resh)
       else ! 'dirichlet'
          call dirichlet_matrix_resh_par(plan, ind-1, a_resh)
       endif 
       call setup_cyclic_tridiag( a_resh, nr, cts, ipiv )
       call solve_cyclic_tridiag(cts,ipiv,plan%array_lin_sys(:,j,1),nr,plan%array_lin_sys(:,j,1))         
    enddo

    ! Remapping to do inverse FFTs
    call apply_remap_3D( plan%rmp3_2, plan%array_lin_sys, plan%array_fft ) 

    ! Inverse FFTs (in the theta-direction)
    do i=1,nr/int(colsz)
       call apply_fft_c2c_1d( plan%inv_fft_plan, plan%array_fft(i,:,1), plan%array_fft(i,:,1) ) 
    enddo

    phi = real(plan%array_fft(:,:,1), f64)/ntheta

  end subroutine solve_qn_2d_with_finite_diff_par


  subroutine delete_new_qns_2d_with_finite_diff_plan_par(plan)

       type (qns_2d_with_finite_diff_plan_par), pointer :: plan
       sll_int32                                        :: ierr

       ! Fixme: some error checking, whether the poisson pointer is 
       ! associated for instance
       SLL_ASSERT( associated(plan) )

       call delete_fft_plan1d(plan%fft_plan)
       call delete_fft_plan1d(plan%inv_fft_plan)

       call delete_layout_3D( plan%layout_fft )
       call delete_layout_3D( plan%layout_lin_sys )

       SLL_DEALLOCATE_ARRAY(plan%rho, ierr)
       SLL_DEALLOCATE_ARRAY(plan%c, ierr)
       SLL_DEALLOCATE_ARRAY(plan%Te, ierr)
       SLL_DEALLOCATE_ARRAY(plan%f, ierr)
       SLL_DEALLOCATE_ARRAY(plan%g, ierr)

       SLL_DEALLOCATE_ARRAY(plan%array_fft,ierr)
       SLL_DEALLOCATE_ARRAY(plan%array_lin_sys,ierr)

       SLL_DEALLOCATE(plan, ierr)

  end subroutine delete_new_qns_2d_with_finite_diff_plan_par


  subroutine dirichlet_matrix_resh_par(plan, j, a_resh)

    type(qns_2d_with_finite_diff_plan_par), pointer :: plan!Matrix is sequential
    sll_real64                                      :: dr, dtheta, Zi
    sll_real64                                      :: r, rmin, rmax
    sll_int32                                       :: i, j, nr, ierr
    sll_real64, dimension(:)                        :: a_resh
    sll_real64, dimension(:), allocatable           :: c, Te 
    ! C & Te are the vector of the Cr & Te(i) respectively

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
       a_resh(3*(i-1)+2) = 2/dr**2 + 2/(r*dtheta)**2 + 1/(Zi*Te(i)) &
                                      - 2*cos(j*dtheta)/(r*dtheta)**2
       if (i<nr) a_resh(3*(i-1)+3) = -( 1/dr**2 + c(i)/(2*dr) )
    enddo

    SLL_DEALLOCATE_ARRAY( c, ierr )
    SLL_DEALLOCATE_ARRAY( Te, ierr )

  end subroutine dirichlet_matrix_resh_par


  subroutine neumann_matrix_resh_par(plan, j, a_resh)
    type(qns_2d_with_finite_diff_plan_par), pointer :: plan!Matrix is sequential
    sll_real64                                      :: dr, dtheta, Zi
    sll_real64                                      :: rmin, rmax, r
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
       a_resh(3*(i-1)+2) = 2/dr**2 + 2/(r*dtheta)**2 + 1/(Zi*Te(i)) - &
                                          2*cos(j*dtheta)/(r*dtheta)**2
       a_resh(3*(i-1)+3) = -( 1/dr**2 + c(i)/(2*dr) )
    enddo

    a_resh(3*(nr-1)+1) = -2/dr**2
    a_resh(3*(nr-1)+2) = 2/dr**2 + 2/(rmax*dtheta)**2 * &
                       (-cos(j*dtheta)+1) + 1/(Zi*Te(nr))

    SLL_DEALLOCATE_ARRAY( c, ierr )
    SLL_DEALLOCATE_ARRAY( Te, ierr )

  end subroutine neumann_matrix_resh_par


end module sll_qns_2d_with_finite_diff_par
