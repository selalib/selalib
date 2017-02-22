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

module sll_m_accumulators
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_cartesian_meshes, only: &
    sll_t_cartesian_mesh_2d

  implicit none

  public :: &
    sll_t_charge_accumulator_cell_2d, &
    sll_t_electric_field_accumulator, &
    sll_t_electric_field_accumulator_cs, &
    sll_t_field_accumulator_cell, &
    sll_t_field_accumulator_cs, &
    sll_f_new_charge_accumulator_2d, &
    sll_f_new_charge_accumulator_2d_cs, &
    sll_f_new_field_accumulator_2d, &
    sll_f_new_field_accumulator_cs_2d, &
    sll_s_reset_charge_accumulator_to_zero, &
    sll_s_reset_charge_accumulator_to_zero_cs, &
    sll_s_reset_field_accumulator_cs_to_zero, &
    sll_s_reset_field_accumulator_to_zero, &
    sll_t_charge_accumulator_2d, &
    sll_t_charge_accumulator_2d_cs, &
    sll_t_charge_accumulator_2d_cs_ptr, &
    sll_t_charge_accumulator_2d_ptr, &
    sll_o_delete, &
    sll_s_sum_accumulators, &
    sll_s_sum_accumulators_cs

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  ! The idea of having this data structure is to precompute the values of
  ! the electric field, store them redundantly on a cell-based structure and
  ! reduce the memory thrashing that would inevitably occur if we were to
  ! compute the values of the electric field from the potential repeatedly.
  ! These costs ought to be amortized if the number of points to be advected,
  ! or particles to be advanced (per cell) is relatively large compared with
  ! the cost of filling out this data structure.
  !
  ! The values of the electric field within the structure in 2D are named 
  ! like this:
  !                NW (north-west)       NE (north-east)
  !                  +--------------------+
  !                  |                    |
  !                  |                    |
  !                  |                    |
  !                  |                    |
  !                  |      cell(i)       |
  !                  |                    |
  !                  |                    |
  !                  |                    |
  !                  |                    |
  !                  +--------------------+
  !                SW (south-west)       SE (south-east)
  !
  ! Each corner represents a point with two electric field components defined:
  ! ex and ey (x-component and y-component of the field respectively).
  
  type sll_t_charge_accumulator_cell_2d
     sll_real64 :: q_sw
     sll_real64 :: q_se
     sll_real64 :: q_nw
     sll_real64 :: q_ne
  end type sll_t_charge_accumulator_cell_2d

  type sll_t_charge_accumulator_2d
     type(sll_t_cartesian_mesh_2d), pointer :: mesh
     type(sll_t_charge_accumulator_cell_2d), dimension(:), pointer :: q_acc
  end type sll_t_charge_accumulator_2d

  type sll_t_charge_accumulator_2d_ptr
     type(sll_t_charge_accumulator_2d), pointer :: q
  end type sll_t_charge_accumulator_2d_ptr


  type charge_accumulator_cell_2d_CS
     sll_real64 :: q_im1j
     sll_real64 :: q_ij
     sll_real64 :: q_ip1j
     sll_real64 :: q_ip2j
     
     sll_real64 :: q_im1jm1
     sll_real64 :: q_ijm1
     sll_real64 :: q_ip1jm1
     sll_real64 :: q_ip2jm1
     
     sll_real64 :: q_im1jp1
     sll_real64 :: q_ijp1
     sll_real64 :: q_ip1jp1
     sll_real64 :: q_ip2jp1
     
     sll_real64 :: q_im1jp2
     sll_real64 :: q_ijp2
     sll_real64 :: q_ip1jp2
     sll_real64 :: q_ip2jp2
  end type charge_accumulator_cell_2d_CS

  type sll_t_charge_accumulator_2d_cs
     type(sll_t_cartesian_mesh_2d), pointer :: mesh
     type(charge_accumulator_cell_2d_CS), dimension(:), pointer :: q_acc
  end type sll_t_charge_accumulator_2d_cs
  
  type sll_t_charge_accumulator_2d_cs_ptr
     type(sll_t_charge_accumulator_2d_cs), pointer :: q
  end type sll_t_charge_accumulator_2d_cs_ptr
  
  type sll_t_field_accumulator_cell
     sll_real64 :: Ex_sw
     sll_real64 :: Ex_se
     sll_real64 :: Ex_nw
     sll_real64 :: Ex_ne
     sll_real64 :: Ey_sw
     sll_real64 :: Ey_se
     sll_real64 :: Ey_nw
     sll_real64 :: Ey_ne
  end type sll_t_field_accumulator_cell

  type sll_t_electric_field_accumulator
     type(sll_t_cartesian_mesh_2d), pointer :: mesh
     type(sll_t_field_accumulator_cell), dimension(:), pointer :: e_acc
  end type sll_t_electric_field_accumulator

  type sll_t_field_accumulator_cs! needed for the Cubic Spline interpolation
     sll_real64 :: Ex_im1j,    Ey_im1j
     sll_real64 :: Ex_ij,      Ey_ij
     sll_real64 :: Ex_ip1j,    Ey_ip1j
     sll_real64 :: Ex_ip2j,    Ey_ip2j

     sll_real64 :: Ex_im1jm1,  Ey_im1jm1
     sll_real64 :: Ex_ijm1,    Ey_ijm1
     sll_real64 :: Ex_ip1jm1,  Ey_ip1jm1
     sll_real64 :: Ex_ip2jm1,  Ey_ip2jm1

     sll_real64 :: Ex_im1jp1,  Ey_im1jp1
     sll_real64 :: Ex_ijp1,    Ey_ijp1
     sll_real64 :: Ex_ip1jp1,  Ey_ip1jp1
     sll_real64 :: Ex_ip2jp1,  Ey_ip2jp1

     sll_real64 :: Ex_im1jp2,  Ey_im1jp2
     sll_real64 :: Ex_ijp2,    Ey_ijp2
     sll_real64 :: Ex_ip1jp2,  Ey_ip1jp2
     sll_real64 :: Ex_ip2jp2,  Ey_ip2jp2
  end type sll_t_field_accumulator_cs

  type sll_t_electric_field_accumulator_cs
     type(sll_t_cartesian_mesh_2d), pointer :: mesh
     type(sll_t_field_accumulator_cs), dimension(:), pointer :: e_acc
  end type sll_t_electric_field_accumulator_cs

  interface sll_o_delete
     module procedure delete_charge_accumulator_2d
  end interface sll_o_delete

  interface operator(+)
     module procedure q_acc_add, q_acc_add_CS
  end interface operator(+)

contains
  
  function sll_f_new_charge_accumulator_2d( mesh_2d ) result(acc)
    type(sll_t_cartesian_mesh_2d), pointer       :: mesh_2d
    type(sll_t_charge_accumulator_2d), pointer :: acc
    sll_int32  :: num_cells1
    sll_int32  :: num_cells2
    sll_int32  :: num_cells_total
    sll_int32  :: ierr

    if( .not. associated(mesh_2d) ) then
       print *, 'ERROR: sll_f_new_charge_accumulator_2d(), passed mesh is not ', &
            'associated. Exiting...'
       stop
    end if

    SLL_ALLOCATE( acc, ierr)
    acc%mesh        => mesh_2d
    num_cells1      = mesh_2d%num_cells1
    num_cells2      = mesh_2d%num_cells2
    num_cells_total = num_cells1*num_cells2
    SLL_ALLOCATE( acc%q_acc(num_cells_total), ierr)
    call sll_s_reset_charge_accumulator_to_zero( acc )

  end function sll_f_new_charge_accumulator_2d
  
  function q_acc_add( q1, q2 ) result(res)
    type(sll_t_charge_accumulator_cell_2d) :: res
    type(sll_t_charge_accumulator_cell_2d), intent(in) :: q1
    type(sll_t_charge_accumulator_cell_2d), intent(in) :: q2

    res%q_sw = q1%q_sw + q2%q_sw
    res%q_se = q1%q_se + q2%q_se
    res%q_nw = q1%q_nw + q2%q_nw
    res%q_ne = q1%q_ne + q2%q_ne
  end function q_acc_add

  function q_acc_add_CS( q1, q2 ) result(res)
    type(charge_accumulator_cell_2d_CS) :: res
    type(charge_accumulator_cell_2d_CS), intent(in) :: q1
    type(charge_accumulator_cell_2d_CS), intent(in) :: q2

    res%q_im1j = q1%q_im1j + q2%q_im1j
    res%q_ij   = q1%q_ij   + q2%q_ij
    res%q_ip1j = q1%q_ip1j + q2%q_ip1j
    res%q_ip2j = q1%q_ip2j + q2%q_ip2j

    res%q_im1jm1 = q1%q_im1jm1 + q2%q_im1jm1
    res%q_ijm1   = q1%q_ijm1   + q2%q_ijm1
    res%q_ip1jm1 = q1%q_ip1jm1 + q2%q_ip1jm1
    res%q_ip2jm1 = q1%q_ip2jm1 + q2%q_ip2jm1

    res%q_im1jp1 = q1%q_im1jp1 + q2%q_im1jp1
    res%q_ijp1   = q1%q_ijp1   + q2%q_ijp1
    res%q_ip1jp1 = q1%q_ip1jp1 + q2%q_ip1jp1
    res%q_ip2jp1 = q1%q_ip2jp1 + q2%q_ip2jp1

    res%q_im1jp2 = q1%q_im1jp2 + q2%q_im1jp2
    res%q_ijp2   = q1%q_ijp2   + q2%q_ijp2
    res%q_ip1jp2 = q1%q_ip1jp2 + q2%q_ip1jp2
    res%q_ip2jp2 = q1%q_ip2jp2 + q2%q_ip2jp2
  end function q_acc_add_CS

  subroutine sll_s_reset_charge_accumulator_to_zero( acc )
    type(sll_t_charge_accumulator_2d), pointer :: acc
    sll_int32 :: i
    sll_int32 :: num_cells

    num_cells = acc%mesh%num_cells1*acc%mesh%num_cells2

    do i=1,num_cells
       acc%q_acc(i)%q_sw = 0.0_f64
       acc%q_acc(i)%q_se = 0.0_f64
       acc%q_acc(i)%q_nw = 0.0_f64
       acc%q_acc(i)%q_ne = 0.0_f64
    end do
  end subroutine sll_s_reset_charge_accumulator_to_zero



  subroutine delete_charge_accumulator_2d( acc )
     type(sll_t_charge_accumulator_2d), pointer :: acc
     sll_int32  :: ierr
     
     if ( .not.associated(acc) ) then
        print*, 'delete_charge_accumulator_2d(): ',&
             'ERROR, pointer was not associated.'
        stop
     endif

     if( associated(acc%q_acc) ) then
        SLL_DEALLOCATE(acc%q_acc, ierr)
     end if
     SLL_DEALLOCATE(acc, ierr)
   end subroutine delete_charge_accumulator_2d


   function sll_f_new_charge_accumulator_2d_cs( mesh_2d ) result(acc)
     type(sll_t_cartesian_mesh_2d), pointer       :: mesh_2d
     type(sll_t_charge_accumulator_2d_cs), pointer :: acc
     sll_int32  :: num_cells1
     sll_int32  :: num_cells2
     sll_int32  :: num_cells_total
     sll_int32  :: ierr
     
     if( .not. associated(mesh_2d) ) then
        print *, 'ERROR: sll_f_new_charge_accumulator_2d(), passed mesh is not ', &
             'associated. Exiting...'
        stop
     end if
     
     SLL_ALLOCATE( acc, ierr)
     acc%mesh        => mesh_2d
     num_cells1      = mesh_2d%num_cells1
     num_cells2      = mesh_2d%num_cells2
     num_cells_total = num_cells1*num_cells2
     SLL_ALLOCATE( acc%q_acc(num_cells_total), ierr)
     call sll_s_reset_charge_accumulator_to_zero_cs( acc )
   end function sll_f_new_charge_accumulator_2d_cs

   subroutine sll_s_reset_charge_accumulator_to_zero_cs( acc )
     type(sll_t_charge_accumulator_2d_cs), pointer :: acc
     sll_int32 :: i
     sll_int32 :: num_cells
     
     num_cells = acc%mesh%num_cells1*acc%mesh%num_cells2
     
     do i=1,num_cells
        acc%q_acc(i)%q_im1j = 0.0_f64 ; acc%q_acc(i)%q_im1jm1 = 0.0_f64
        acc%q_acc(i)%q_ij   = 0.0_f64 ; acc%q_acc(i)%q_ijm1   = 0.0_f64
        acc%q_acc(i)%q_ip1j = 0.0_f64 ; acc%q_acc(i)%q_ip1jm1 = 0.0_f64
        acc%q_acc(i)%q_ip2j = 0.0_f64 ; acc%q_acc(i)%q_ip2jm1 = 0.0_f64

        acc%q_acc(i)%q_im1jp1 = 0.0_f64 ; acc%q_acc(i)%q_im1jp2 = 0.0_f64
        acc%q_acc(i)%q_ijp1   = 0.0_f64 ; acc%q_acc(i)%q_ijp2   = 0.0_f64
        acc%q_acc(i)%q_ip1jp1 = 0.0_f64 ; acc%q_acc(i)%q_ip1jp2 = 0.0_f64
        acc%q_acc(i)%q_ip2jp1 = 0.0_f64 ; acc%q_acc(i)%q_ip2jp2 = 0.0_f64 
     end do
   end subroutine sll_s_reset_charge_accumulator_to_zero_cs

   subroutine delete_charge_accumulator_2d_CS( acc )
     type(sll_t_charge_accumulator_2d_cs), pointer :: acc
     sll_int32  :: ierr
     
     if ( .not.associated(acc) ) then
        print*, 'delete_charge_accumulator_2d(): ',&
             'ERROR, pointer was not associated.'
        stop
     endif
     
     if( associated(acc%q_acc) ) then
        SLL_DEALLOCATE(acc%q_acc, ierr)
     end if
     SLL_DEALLOCATE(acc, ierr)
   end subroutine delete_charge_accumulator_2d_CS
   
   
   



   function sll_f_new_field_accumulator_2d( mesh_2d ) result(E_acc)
     type(sll_t_cartesian_mesh_2d), pointer        ::  mesh_2d
     type(sll_t_electric_field_accumulator), pointer ::  E_acc
     sll_int32  :: num_cells1
     sll_int32  :: num_cells2
     sll_int32  :: num_cells_total
     sll_int32  :: ierr
     
     if( .not. associated(mesh_2d) ) then
        print *, 'ERROR: sll_f_new_charge_accumulator_2d(), passed mesh is not ', &
             'associated. Exiting...'
        stop
     end if
     
     SLL_ALLOCATE( E_acc, ierr)
     E_acc%mesh        => mesh_2d
     num_cells1      = mesh_2d%num_cells1
     num_cells2      = mesh_2d%num_cells2
     num_cells_total = num_cells1*num_cells2
     SLL_ALLOCATE( E_acc%e_acc(num_cells_total), ierr)
     call sll_s_reset_field_accumulator_to_zero( E_acc )
   end function sll_f_new_field_accumulator_2d
   
  function sll_f_new_field_accumulator_cs_2d( mesh_2d ) result(E_acc)
    type(sll_t_cartesian_mesh_2d), pointer           ::  mesh_2d
    type(sll_t_electric_field_accumulator_cs), pointer ::  E_acc
    sll_int32  :: num_cells1
    sll_int32  :: num_cells2
    sll_int32  :: num_cells_total
    sll_int32  :: ierr

    if( .not. associated(mesh_2d) ) then
       print *, 'ERROR: sll_f_new_charge_accumulator_2d(), passed mesh is not ', &
            'associated. Exiting...'
       stop
    end if

    SLL_ALLOCATE( E_acc, ierr)
    E_acc%mesh        => mesh_2d
    num_cells1      = mesh_2d%num_cells1
    num_cells2      = mesh_2d%num_cells2
    num_cells_total = num_cells1*num_cells2
    SLL_ALLOCATE( E_acc%e_acc(num_cells_total), ierr)
    call sll_s_reset_field_accumulator_cs_to_zero( E_acc )

  end function sll_f_new_field_accumulator_cs_2d
  

  subroutine sll_s_reset_field_accumulator_to_zero( E_acc )
    type(sll_t_electric_field_accumulator), pointer :: E_acc
    sll_int32 :: i
    sll_int32 :: num_cells

    num_cells = E_acc%mesh%num_cells1*E_acc%mesh%num_cells2
    do i=1,num_cells
       E_acc%e_acc(i)%Ex_sw = 0.0_f64
       E_acc%e_acc(i)%Ex_se = 0.0_f64
       E_acc%e_acc(i)%Ex_nw = 0.0_f64
       E_acc%e_acc(i)%Ex_ne = 0.0_f64
       E_acc%e_acc(i)%Ey_sw = 0.0_f64
       E_acc%e_acc(i)%Ey_se = 0.0_f64
       E_acc%e_acc(i)%Ey_nw = 0.0_f64
       E_acc%e_acc(i)%Ey_ne = 0.0_f64
    end do
  end subroutine sll_s_reset_field_accumulator_to_zero

  subroutine sll_s_reset_field_accumulator_cs_to_zero( E_acc )
! -------------   CS is for Cubic Spline   -------------
    type(sll_t_electric_field_accumulator_cs), pointer :: E_acc
    sll_int32 :: i
    sll_int32 :: num_cells

    num_cells = E_acc%mesh%num_cells1*E_acc%mesh%num_cells2
    do i=1,num_cells
       E_acc%e_acc(i)%Ex_im1j = 0.0_f64 ;    E_acc%e_acc(i)%Ey_im1j = 0.0_f64
       E_acc%e_acc(i)%Ex_ij   = 0.0_f64 ;    E_acc%e_acc(i)%Ey_ij   = 0.0_f64
       E_acc%e_acc(i)%Ex_ip1j = 0.0_f64 ;    E_acc%e_acc(i)%Ey_ip1j = 0.0_f64
       E_acc%e_acc(i)%Ex_ip2j = 0.0_f64 ;    E_acc%e_acc(i)%Ey_ip2j = 0.0_f64
       E_acc%e_acc(i)%Ex_im1jm1 = 0.0_f64 ;  E_acc%e_acc(i)%Ey_im1jm1 = 0.0_f64
       E_acc%e_acc(i)%Ex_ijm1   = 0.0_f64 ;  E_acc%e_acc(i)%Ey_ijm1   = 0.0_f64
       E_acc%e_acc(i)%Ex_ip1jm1 = 0.0_f64 ;  E_acc%e_acc(i)%Ey_ip1jm1 = 0.0_f64
       E_acc%e_acc(i)%Ex_ip2jm1 = 0.0_f64 ;  E_acc%e_acc(i)%Ey_ip2jm1 = 0.0_f64
       E_acc%e_acc(i)%Ex_im1jp1 = 0.0_f64 ;  E_acc%e_acc(i)%Ey_im1jp1 = 0.0_f64
       E_acc%e_acc(i)%Ex_ijp1   = 0.0_f64 ;  E_acc%e_acc(i)%Ey_ijp1   = 0.0_f64
       E_acc%e_acc(i)%Ex_ip1jp1 = 0.0_f64 ;  E_acc%e_acc(i)%Ey_ip1jp1 = 0.0_f64
       E_acc%e_acc(i)%Ex_ip2jp1 = 0.0_f64 ;  E_acc%e_acc(i)%Ey_ip2jp1 = 0.0_f64
       E_acc%e_acc(i)%Ex_im1jp2 = 0.0_f64 ;  E_acc%e_acc(i)%Ey_im1jp2 = 0.0_f64
       E_acc%e_acc(i)%Ex_ijp2   = 0.0_f64 ;  E_acc%e_acc(i)%Ey_ijp2   = 0.0_f64
       E_acc%e_acc(i)%Ex_ip1jp2 = 0.0_f64 ;  E_acc%e_acc(i)%Ey_ip1jp2 = 0.0_f64
       E_acc%e_acc(i)%Ex_ip2jp2 = 0.0_f64 ;  E_acc%e_acc(i)%Ey_ip2jp2 = 0.0_f64
    end do
  end subroutine sll_s_reset_field_accumulator_cs_to_zero

  subroutine sll_s_sum_accumulators( tab, n_threads, n_cells )
    sll_int32, intent(in)  :: n_threads, n_cells
    type(sll_t_charge_accumulator_2d_ptr), dimension(:), pointer, intent(inout) :: tab
    sll_int32  :: i, j   
    
    do i = 1, n_cells  
       do j = 2, n_threads
          tab(1)%q%q_acc(i) = tab(1)%q%q_acc(i) + tab(j)%q%q_acc(i)
       enddo
    enddo
  end subroutine sll_s_sum_accumulators

  subroutine sll_s_sum_accumulators_cs( tab, n_threads, n_cells )
    sll_int32, intent(in)  :: n_threads, n_cells
    type(sll_t_charge_accumulator_2d_cs_ptr), dimension(:), pointer, intent(inout) :: tab
    sll_int32  :: i, j   
    
    do i = 1, n_cells  
       do j = 2, n_threads
          tab(1)%q%q_acc(i) = tab(1)%q%q_acc(i) + tab(j)%q%q_acc(i)
       enddo
    enddo
  end subroutine sll_s_sum_accumulators_cs
  
   ! The above are the only routines that should live here. Something like taking
   ! a list of particules and accumulating the charge, will change the charge but
   ! not the particles so it could legitimately live here. This would introduce
   ! a dependency on the particle group module. At the same time, a particle
   ! pusher would take the accumulators and change the particle list, thus it 
   ! would be tempting to put it in the particle group module, but this would 
   ! introduce the opposite dependency. Hence

!!$  subroutine sll_accumulate_charge_on_cell( &
!!$                particle, &
!!$                charge )
!!$
!!$    type(sll_particle_2d), intent(in) :: particle
!!$    type(charge_accumulator_cell), intent(inout) :: charge
!!$    sll_int64    ::  j
!!$    
!!$    charge%q_sw = charge%q_sw + &
!!$         particle%q * (1._f64 - particle%dx) * (1._f64 - particle%dy)
!!$
!!$    charge%q_se = charge%q_se + &
!!$         particle%q * particle%dx * (1._f64 - particle%dy)
!!$
!!$    charge%q_nw = charge%q_nw + &
!!$         particle%q * (1._f64 - particle%dx) * particle%dy
!!$
!!$    charge%q_ne = charge%q_ne + &
!!$         particle%q * particle%dx * particle%dy
!!$
!!$  end subroutine sll_accumulate_charge_on_cell


!!$  EDWIN CG
!!$  ! The following is a pseudo-accumulator. It is written for use while some
!!$  ! problems with the parallel loading of the previous accumulator are resolved.
!!$  ! This one is only useful to precompute the values of the electric field in
!!$  ! a structure that can have the same layout as the potential array.
!!$  type :: efield_2d_point
!!$     sll_real64 :: ex
!!$     sll_real64 :: ey
!!$  end type efield_2d_point



end module sll_m_accumulators
