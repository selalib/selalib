!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************


module sll_charge_to_density_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_accumulators.h"
  use sll_logical_meshes
  implicit none
  
contains

!!$  function new_charge_to_meshdens(   &
!!$                 all_charge, m2d ) result(ro)
!!$
!!$    type(sll_logical_mesh_2d), intent(in) :: m2d
!!$    type(charge_accumulator_cell),dimension(1:m2d%num_cells2*m2d%num_cells1), &
!!$    intent(in) :: all_charge
!!$    sll_real64, dimension(:,:), pointer  :: ro
!!$    sll_int32  :: ierr
!!$
!!$    SLL_ALLOCATE( ro(1:m2d%num_cells1+1, 1:m2d%num_cells2+1), ierr)
!!$
!!$    call sll_charge_to_meshdensity( all_charge, m2d, ro )
!!$
!!$  end function new_charge_to_meshdens

  subroutine sll_convert_charge_to_rho_2d_per_per( charge, rho )
    type(sll_charge_accumulator_2d), pointer :: charge
    sll_real64, dimension(:,:), intent(out)  :: rho
    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: ncx
    sll_int32  :: ncy
    sll_real64 :: deltax
    sll_real64 :: deltay
    sll_real64 :: factor

    ncx    = charge%mesh%num_cells1
    ncy    = charge%mesh%num_cells2
    deltax = charge%mesh%delta_eta1
    deltay = charge%mesh%delta_eta2
    factor = 1.0_f64/(deltax*deltay)

    ! loop over the node indices of rho. Edges might need special treatment,
    ! thus do the non-edge nodes first.
    do j=2,ncy
       do i=2,ncx
          rho(i,j) = SLL_GET_CHARGE_ACC_VALUE( charge, i,   j,  q_sw ) + &
                     SLL_GET_CHARGE_ACC_VALUE( charge, i-1, j,  q_se ) + &
                     SLL_GET_CHARGE_ACC_VALUE( charge, i-1, j-1,q_ne ) + &
                     SLL_GET_CHARGE_ACC_VALUE( charge, i,   j-1,q_nw )
          rho(i,j) = rho(i,j)*factor 
       end do
    end do

    do j=2,ncy
       rho(1,j) = SLL_GET_CHARGE_ACC_VALUE( charge,   1,   j,q_sw ) + &
                  SLL_GET_CHARGE_ACC_VALUE( charge,   1, j-1,q_nw ) + &  
                  SLL_GET_CHARGE_ACC_VALUE( charge, ncx,   j,q_se ) + &  
                  SLL_GET_CHARGE_ACC_VALUE( charge, ncx, j-1,q_ne )
       rho(1,j) = rho(1,j)*factor
       rho(ncx+1,j) = rho(1,j)
    end do

    do i=2,ncx
       rho(i,1) = SLL_GET_CHARGE_ACC_VALUE( charge, i-1,   1,q_se ) + &
                  SLL_GET_CHARGE_ACC_VALUE( charge,   i,   1,q_sw ) + &  
                  SLL_GET_CHARGE_ACC_VALUE( charge, i-1, ncy,q_ne ) + &  
                  SLL_GET_CHARGE_ACC_VALUE( charge,   i, ncy,q_nw )
       rho(i,1) = rho(i,1)*factor
       rho(i,ncy+1) = rho(i,1)
    end do

    rho(1,1) = SLL_GET_CHARGE_ACC_VALUE( charge,   1,  1,q_sw ) + &
               SLL_GET_CHARGE_ACC_VALUE( charge,   1,ncy,q_nw ) + &
               SLL_GET_CHARGE_ACC_VALUE( charge, ncx,  1,q_se ) + &
               SLL_GET_CHARGE_ACC_VALUE( charge, ncx,ncy,q_ne ) 
    rho(1,1) = rho(1,1)*factor
    rho(    1,ncy+1) = rho(1,1)
    rho(ncx+1,    1) = rho(1,1)
    rho(ncx+1,ncy+1) = rho(1,1)
!!$    density = 0.0_f64
!!$    do k = 1, ncx * ncy
!!$       i = mod( k-1, ncx ) + 1
!!$       j = int((k-1)/ncx ) + 1
!!$       density(i  ,j  ) = density(i,  j)   + all_charge(k)%q_sw
!!$       density(i+1,j  ) = density(i+1,j)   + all_charge(k)%q_se
!!$       density(i  ,j+1) = density(i,  j+1) + all_charge(k)%q_nw
!!$       density(i+1,j+1) = density(i+1,j+1) + all_charge(k)%q_ne
!!$    enddo
  end subroutine sll_convert_charge_to_rho_2d_per_per

  subroutine sll_convert_charge_to_rho_2d_per_per_CS( charge, rho )
    type(sll_charge_accumulator_2d_CS), pointer :: charge
    sll_real64, dimension(:,:), intent(out)  :: rho
    sll_int32  :: i, im1, ip2
    sll_int32  :: j, jm1, jp2
    sll_int32  :: k
    sll_int32  :: ncx
    sll_int32  :: ncy
    sll_real64 :: deltax
    sll_real64 :: deltay
    sll_real64 :: factor

    ncx    = charge%mesh%num_cells1
    ncy    = charge%mesh%num_cells2
    deltax = charge%mesh%delta_eta1
    deltay = charge%mesh%delta_eta2
    factor = 1.0_f64/(deltax*deltay)

    rho = 0.0_f64
    do j = 1, ncy
       jm1 = j - 1
       jp2 = j + 2
       if (j==1) then
          jm1 = ncy
       endif
       if (j==ncy) then
          jp2 = 2
       endif
       do i = 1, ncx
          k = i+(j-1)*ncx
          im1 = i - 1
          ip2 = i + 2
          if (i==1) then
             im1 = ncx
          endif
          if (i==ncx) then
             ip2 = 2
          endif
          rho(im1, j) = rho(im1, j) + charge%q_acc(k)%q_im1j * factor
          rho(i  , j) = rho(i,   j) + charge%q_acc(k)%q_ij * factor
          rho(i+1, j) = rho(i+1, j) + charge%q_acc(k)%q_ip1j * factor
          rho(ip2, j) = rho(ip2, j) + charge%q_acc(k)%q_ip2j * factor

          rho(im1, jm1) = rho(im1, jm1) + charge%q_acc(k)%q_im1jm1 * factor
          rho(i  , jm1) = rho(i  , jm1) + charge%q_acc(k)%q_ijm1 * factor
          rho(i+1, jm1) = rho(i+1, jm1) + charge%q_acc(k)%q_ip1jm1 * factor
          rho(ip2, jm1) = rho(ip2, jm1) + charge%q_acc(k)%q_ip2jm1 * factor

          rho(im1, j+1) = rho(im1, j+1) + charge%q_acc(k)%q_im1jp1 * factor
          rho(i  , j+1) = rho(i  , j+1) + charge%q_acc(k)%q_ijp1 * factor
          rho(i+1, j+1) = rho(i+1, j+1) + charge%q_acc(k)%q_ip1jp1 * factor
          rho(ip2, j+1) = rho(ip2, j+1) + charge%q_acc(k)%q_ip2jp1 * factor

          rho(im1, jp2) = rho(im1, jp2) + charge%q_acc(k)%q_im1jp2 * factor
          rho(i  , jp2) = rho(i  , jp2) + charge%q_acc(k)%q_ijp2 * factor
          rho(i+1, jp2) = rho(i+1, jp2) + charge%q_acc(k)%q_ip1jp2 * factor
          rho(ip2, jp2) = rho(ip2, jp2) + charge%q_acc(k)%q_ip2jp2 * factor
       enddo
    enddo

  end subroutine sll_convert_charge_to_rho_2d_per_per_CS


! ------------   COMMENTS  -------------
! The mesh elec field, e.g. E1(i,j), is defined for i=1,ncx+1 (the periodic BC
! is automatically considered : E1(1) = E1(ncx+1) and for j=1,ncy+1
! 
! Node (i,j) is the left_bottom corner of the cell i+(j-1)*ncx
! --------------------------------------
!

  subroutine sll_accumulate_field( E1, E2, E_accumulator)
    sll_real64, dimension(:,:), pointer, intent(in)        :: E1, E2
    type(electric_field_accumulator), pointer, intent(out) :: E_accumulator
    sll_int32  ::  i, j
    sll_int32  ::  ncx, ncy

    ncx  = E_accumulator%mesh%num_cells1
    ncy  = E_accumulator%mesh%num_cells2
    do j = 1, ncy
       do i = 1, ncx
          E_accumulator%e_acc(i+(j-1)*ncx)%Ex_sw = E1(i  ,j  )
          E_accumulator%e_acc(i+(j-1)*ncx)%Ex_se = E1(i+1,j  )
          E_accumulator%e_acc(i+(j-1)*ncx)%Ex_nw = E1(i  ,j+1)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ex_ne = E1(i+1,j+1)
          !
          E_accumulator%e_acc(i+(j-1)*ncx)%Ey_sw = E2(i  ,j  )
          E_accumulator%e_acc(i+(j-1)*ncx)%Ey_se = E2(i+1,j  )
          E_accumulator%e_acc(i+(j-1)*ncx)%Ey_nw = E2(i  ,j+1)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ey_ne = E2(i+1,j+1)
       enddo
    enddo

  end subroutine sll_accumulate_field


  subroutine sll_accumulate_field_CS( E1, E2, E_accumulator)
! ------   Recall : CS is for Cubic Splines  ------
    sll_real64, dimension(:,:), pointer, intent(in)        :: E1, E2
    type(electric_field_accumulator_CS), pointer, intent(out) :: E_accumulator
    sll_int32  ::  i, j, im1, ip2, jm1, jp2
    sll_int32  ::  ncx, ncy

    ncx  = E_accumulator%mesh%num_cells1
    ncy  = E_accumulator%mesh%num_cells2
    do j = 1, ncy
       jm1 = j - 1
       jp2 = j + 2
       if (j==1) then
          jm1 = ncy
       endif
       if (j==ncy) then
          jp2 = 2
       endif
       do i = 1, ncx
          im1 = i - 1
          ip2 = i + 2
          if (i==1) then
             im1 = ncx
          endif
          if (i==ncx) then
             ip2 = 2
          endif
          E_accumulator%e_acc(i+(j-1)*ncx)%Ex_im1j = E1(im1, j)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ex_ij   = E1(i,   j)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ex_ip1j = E1(i+1, j)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ex_ip2j = E1(ip2, j)

          E_accumulator%e_acc(i+(j-1)*ncx)%Ex_im1jm1 = E1(im1, jm1)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ex_ijm1   = E1(i,   jm1)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ex_ip1jm1 = E1(i+1, jm1)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ex_ip2jm1 = E1(ip2, jm1)

          E_accumulator%e_acc(i+(j-1)*ncx)%Ex_im1jp1 = E1(im1, j+1)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ex_ijp1   = E1(i,   j+1)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ex_ip1jp1 = E1(i+1, j+1)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ex_ip2jp1 = E1(ip2, j+1)

          E_accumulator%e_acc(i+(j-1)*ncx)%Ex_im1jp2 = E1(im1, jp2)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ex_ijp2   = E1(i,   jp2)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ex_ip1jp2 = E1(i+1, jp2)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ex_ip2jp2 = E1(ip2, jp2)
          ! ----- !
          ! ----- !
          E_accumulator%e_acc(i+(j-1)*ncx)%Ey_im1j = E2(im1, j)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ey_ij   = E2(i,   j)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ey_ip1j = E2(i+1, j)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ey_ip2j = E2(ip2, j)

          E_accumulator%e_acc(i+(j-1)*ncx)%Ey_im1jm1 = E2(im1, jm1)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ey_ijm1   = E2(i,   jm1)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ey_ip1jm1 = E2(i+1, jm1)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ey_ip2jm1 = E2(ip2, jm1)

          E_accumulator%e_acc(i+(j-1)*ncx)%Ey_im1jp1 = E2(im1, j+1)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ey_ijp1   = E2(i,   j+1)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ey_ip1jp1 = E2(i+1, j+1)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ey_ip2jp1 = E2(ip2, j+1)

          E_accumulator%e_acc(i+(j-1)*ncx)%Ey_im1jp2 = E2(im1, jp2)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ey_ijp2   = E2(i,   jp2)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ey_ip1jp2 = E2(i+1, jp2)
          E_accumulator%e_acc(i+(j-1)*ncx)%Ey_ip2jp2 = E2(ip2, jp2)
       enddo
    enddo

  end subroutine sll_accumulate_field_CS

end module sll_charge_to_density_module
