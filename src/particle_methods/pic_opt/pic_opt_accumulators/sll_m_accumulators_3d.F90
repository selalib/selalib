module sll_m_accumulators_3d
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_cartesian_meshes, only: &
  sll_t_cartesian_mesh_3d

implicit none

public :: &
  sll_t_charge_accumulator_cell_3d, &
  sll_t_electric_field_accumulator_3d, &
  sll_t_field_accumulator_3d, &
  sll_f_new_charge_accumulator_3d, &
  sll_s_charge_accumulator_3d_init, &
  sll_f_new_field_accumulator_3d, &
  sll_s_field_accumulator_3d_init, &
  sll_s_reset_charge_accumulator_to_zero, &
  sll_s_reset_field_accumulator_to_zero, &
  sll_t_charge_accumulator_3d, &
  sll_s_sum_accumulators

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
  
  type sll_t_charge_accumulator_cell_3d
     real(kind=f64) :: q_sw
     real(kind=f64) :: q_se
     real(kind=f64) :: q_nw
     real(kind=f64) :: q_ne
  end type sll_t_charge_accumulator_cell_3d

  type sll_t_charge_accumulator_3d
     type(sll_t_cartesian_mesh_3d), pointer :: mesh
     type(sll_t_charge_accumulator_cell_3d), dimension(:), pointer :: q_acc
  end type sll_t_charge_accumulator_3d

  type sll_t_charge_accumulator_3d_ptr
     type(sll_t_charge_accumulator_3d), pointer :: q
  end type sll_t_charge_accumulator_3d_ptr


  type charge_accumulator_cell_3d_CS
     real(kind=f64) :: q_im1j
     real(kind=f64) :: q_ij
     real(kind=f64) :: q_ip1j
     real(kind=f64) :: q_ip2j
     
     real(kind=f64) :: q_im1jm1
     real(kind=f64) :: q_ijm1
     real(kind=f64) :: q_ip1jm1
     real(kind=f64) :: q_ip2jm1
     
     real(kind=f64) :: q_im1jp1
     real(kind=f64) :: q_ijp1
     real(kind=f64) :: q_ip1jp1
     real(kind=f64) :: q_ip2jp1
     
     real(kind=f64) :: q_im1jp2
     real(kind=f64) :: q_ijp2
     real(kind=f64) :: q_ip1jp2
     real(kind=f64) :: q_ip2jp2
  end type charge_accumulator_cell_3d_CS

  type sll_t_charge_accumulator_3d_cs
     type(sll_t_cartesian_mesh_3d), pointer :: mesh
     type(charge_accumulator_cell_3d_CS), dimension(:), pointer :: q_acc
  end type sll_t_charge_accumulator_3d_cs
  
  type sll_t_charge_accumulator_3d_cs_ptr
     type(sll_t_charge_accumulator_3d_cs), pointer :: q
  end type sll_t_charge_accumulator_3d_cs_ptr
  
  type sll_t_field_accumulator_cell
     real(kind=f64) :: Ex_sw
     real(kind=f64) :: Ex_se
     real(kind=f64) :: Ex_nw
     real(kind=f64) :: Ex_ne
     real(kind=f64) :: Ey_sw
     real(kind=f64) :: Ey_se
     real(kind=f64) :: Ey_nw
     real(kind=f64) :: Ey_ne
  end type sll_t_field_accumulator_cell

  type sll_t_electric_field_accumulator
     type(sll_t_cartesian_mesh_3d), pointer :: mesh
     type(sll_t_field_accumulator_cell), dimension(:), pointer :: e_acc
  end type sll_t_electric_field_accumulator

  type sll_t_field_accumulator_cs! needed for the Cubic Spline interpolation
     real(kind=f64) :: Ex_im1j,    Ey_im1j
     real(kind=f64) :: Ex_ij,      Ey_ij
     real(kind=f64) :: Ex_ip1j,    Ey_ip1j
     real(kind=f64) :: Ex_ip2j,    Ey_ip2j

     real(kind=f64) :: Ex_im1jm1,  Ey_im1jm1
     real(kind=f64) :: Ex_ijm1,    Ey_ijm1
     real(kind=f64) :: Ex_ip1jm1,  Ey_ip1jm1
     real(kind=f64) :: Ex_ip2jm1,  Ey_ip2jm1

     real(kind=f64) :: Ex_im1jp1,  Ey_im1jp1
     real(kind=f64) :: Ex_ijp1,    Ey_ijp1
     real(kind=f64) :: Ex_ip1jp1,  Ey_ip1jp1
     real(kind=f64) :: Ex_ip2jp1,  Ey_ip2jp1

     real(kind=f64) :: Ex_im1jp2,  Ey_im1jp2
     real(kind=f64) :: Ex_ijp2,    Ey_ijp2
     real(kind=f64) :: Ex_ip1jp2,  Ey_ip1jp2
     real(kind=f64) :: Ex_ip2jp2,  Ey_ip2jp2
  end type sll_t_field_accumulator_cs

  type sll_t_electric_field_accumulator_cs
     type(sll_t_cartesian_mesh_3d), pointer :: mesh
     type(sll_t_field_accumulator_cs), dimension(:), pointer :: e_acc
  end type sll_t_electric_field_accumulator_cs

  interface sll_o_delete
     module procedure delete_charge_accumulator_3d
  end interface sll_o_delete

  interface operator(+)
     module procedure q_acc_add, q_acc_add_CS
  end interface operator(+)

contains
  
  subroutine sll_s_charge_accumulator_3d_init( acc, mesh_3d ) 
    type(sll_t_charge_accumulator_3d)         :: acc
    type(sll_t_cartesian_mesh_3d),     target :: mesh_3d
    integer(kind=i32)  :: num_cells1
    integer(kind=i32)  :: num_cells2
    integer(kind=i32)  :: num_cells_total
    integer(kind=i32)  :: ierr

    acc%mesh        => mesh_3d
    num_cells1      = mesh_3d%num_cells1
    num_cells2      = mesh_3d%num_cells2
    num_cells_total = num_cells1*num_cells2
    allocate( acc%q_acc(num_cells_total), stat= ierr);        call sll_s_test_error_code( ierr, 'Memory allocation Failure.', "sll_m_accumulators.F90",   197);
    call sll_s_reset_charge_accumulator_to_zero( acc )

  end subroutine sll_s_charge_accumulator_3d_init
  
  function sll_f_new_charge_accumulator_3d( mesh_3d ) result(acc)
    type(sll_t_cartesian_mesh_3d), pointer     :: mesh_3d
    type(sll_t_charge_accumulator_3d), pointer :: acc
    integer(kind=i32)  :: ierr

    if( .not. associated(mesh_3d) ) then
       print *, 'ERROR: sll_f_new_charge_accumulator_3d(), passed mesh is not ', &
            'associated. Exiting...'
       stop
    end if

    allocate( acc, stat= ierr);        call sll_s_test_error_code( ierr, 'Memory allocation Failure.', "sll_m_accumulators.F90",   213);
    call sll_s_charge_accumulator_3d_init( acc, mesh_3d ) 

  end function sll_f_new_charge_accumulator_3d
  
  function q_acc_add( q1, q2 ) result(res)
    type(sll_t_charge_accumulator_cell_3d) :: res
    type(sll_t_charge_accumulator_cell_3d), intent(in) :: q1
    type(sll_t_charge_accumulator_cell_3d), intent(in) :: q2

    res%q_sw = q1%q_sw + q2%q_sw
    res%q_se = q1%q_se + q2%q_se
    res%q_nw = q1%q_nw + q2%q_nw
    res%q_ne = q1%q_ne + q2%q_ne
  end function q_acc_add

  function q_acc_add_CS( q1, q2 ) result(res)
    type(charge_accumulator_cell_3d_CS) :: res
    type(charge_accumulator_cell_3d_CS), intent(in) :: q1
    type(charge_accumulator_cell_3d_CS), intent(in) :: q2

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
    type(sll_t_charge_accumulator_3d) :: acc
    integer(kind=i32) :: i
    integer(kind=i32) :: num_cells

    num_cells = acc%mesh%num_cells1*acc%mesh%num_cells2

    do i=1,num_cells
       acc%q_acc(i)%q_sw = 0.0_f64
       acc%q_acc(i)%q_se = 0.0_f64
       acc%q_acc(i)%q_nw = 0.0_f64
       acc%q_acc(i)%q_ne = 0.0_f64
    end do
  end subroutine sll_s_reset_charge_accumulator_to_zero



  subroutine delete_charge_accumulator_3d( acc )
     type(sll_t_charge_accumulator_3d), pointer :: acc
     integer(kind=i32)  :: ierr
     
     if ( .not.associated(acc) ) then
        print*, 'delete_charge_accumulator_3d(): ',&
             'ERROR, pointer was not associated.'
        stop
     endif

     if( associated(acc%q_acc) ) then
        deallocate(acc%q_acc, stat= ierr);        call sll_s_test_error_code( ierr, 'Error in memory deallocation.', "sll_m_accumulators.F90",  283);   nullify(acc%q_acc);
     end if
     deallocate(acc, stat= ierr);        call sll_s_test_error_code( ierr, 'Error in memory deallocation.', "sll_m_accumulators.F90",  285);   nullify(acc);
   end subroutine delete_charge_accumulator_3d


   subroutine sll_s_charge_accumulator_3d_cs_init( acc, mesh_3d ) 
     type(sll_t_cartesian_mesh_3d),        target  :: mesh_3d
     type(sll_t_charge_accumulator_3d_cs)          :: acc
     integer(kind=i32)  :: num_cells1
     integer(kind=i32)  :: num_cells2
     integer(kind=i32)  :: num_cells_total
     integer(kind=i32)  :: ierr
     
     acc%mesh        => mesh_3d
     num_cells1      = mesh_3d%num_cells1
     num_cells2      = mesh_3d%num_cells2
     num_cells_total = num_cells1*num_cells2
     allocate( acc%q_acc(num_cells_total), stat= ierr);        call sll_s_test_error_code( ierr, 'Memory allocation Failure.', "sll_m_accumulators.F90",   301);
     call sll_s_reset_charge_accumulator_to_zero_cs( acc )
   end subroutine sll_s_charge_accumulator_3d_cs_init

   function sll_f_new_charge_accumulator_3d_cs( mesh_3d ) result(acc)
     type(sll_t_cartesian_mesh_3d), pointer       :: mesh_3d
     type(sll_t_charge_accumulator_3d_cs), pointer :: acc
     integer(kind=i32)  :: ierr
     
     if( .not. associated(mesh_3d) ) then
        print *, 'ERROR: sll_f_new_charge_accumulator_3d(), passed mesh is not ', &
             'associated. Exiting...'
        stop
     end if
     
     allocate( acc, stat= ierr);        call sll_s_test_error_code( ierr, 'Memory allocation Failure.', "sll_m_accumulators.F90",   316);
     call sll_s_charge_accumulator_3d_cs_init(acc, mesh_3d)

   end function sll_f_new_charge_accumulator_3d_cs

   subroutine sll_s_reset_charge_accumulator_to_zero_cs( acc )
     type(sll_t_charge_accumulator_3d_cs) :: acc
     integer(kind=i32) :: i
     integer(kind=i32) :: num_cells
     
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

   subroutine delete_charge_accumulator_3d_CS( acc )
     type(sll_t_charge_accumulator_3d_cs), pointer :: acc
     integer(kind=i32)  :: ierr
     
     if ( .not.associated(acc) ) then
        print*, 'delete_charge_accumulator_3d(): ',&
             'ERROR, pointer was not associated.'
        stop
     endif
     
     if( associated(acc%q_acc) ) then
        deallocate(acc%q_acc, stat= ierr);        call sll_s_test_error_code( ierr, 'Error in memory deallocation.', "sll_m_accumulators.F90",  352);   nullify(acc%q_acc);
     end if
     deallocate(acc, stat= ierr);        call sll_s_test_error_code( ierr, 'Error in memory deallocation.', "sll_m_accumulators.F90",  354);   nullify(acc);
   end subroutine delete_charge_accumulator_3d_CS
   
   
   



   subroutine sll_s_field_accumulator_3d_init(E_acc, mesh_3d ) 
     type(sll_t_electric_field_accumulator)         ::  E_acc
     type(sll_t_cartesian_mesh_3d),          target ::  mesh_3d
     integer(kind=i32)  :: num_cells1
     integer(kind=i32)  :: num_cells2
     integer(kind=i32)  :: num_cells_total
     integer(kind=i32)  :: ierr
     
     E_acc%mesh      => mesh_3d
     num_cells1      = mesh_3d%num_cells1
     num_cells2      = mesh_3d%num_cells2
     num_cells_total = num_cells1*num_cells2
     allocate( E_acc%e_acc(num_cells_total), stat= ierr);        call sll_s_test_error_code( ierr, 'Memory allocation Failure.', "sll_m_accumulators.F90",   374);
     call sll_s_reset_field_accumulator_to_zero( E_acc )

   end subroutine sll_s_field_accumulator_3d_init


   function sll_f_new_field_accumulator_3d( mesh_3d ) result(E_acc)
     type(sll_t_cartesian_mesh_3d), pointer        ::  mesh_3d
     type(sll_t_electric_field_accumulator), pointer ::  E_acc
     integer(kind=i32)  :: num_cells1
     integer(kind=i32)  :: num_cells2
     integer(kind=i32)  :: num_cells_total
     integer(kind=i32)  :: ierr
     
     if( .not. associated(mesh_3d) ) then
        print *, 'ERROR: sll_f_new_charge_accumulator_3d(), passed mesh is not ', &
             'associated. Exiting...'
        stop
     end if
     
     allocate( E_acc, stat= ierr);        call sll_s_test_error_code( ierr, 'Memory allocation Failure.', "sll_m_accumulators.F90",   394);
     call sll_s_field_accumulator_3d_init(E_acc, mesh_3d )

   end function sll_f_new_field_accumulator_3d
   
  subroutine sll_s_field_accumulator_cs_3d_init( E_acc, mesh_3d ) 
    type(sll_t_cartesian_mesh_3d), target     ::  mesh_3d
    type(sll_t_electric_field_accumulator_cs) ::  E_acc
    integer(kind=i32)  :: num_cells1
    integer(kind=i32)  :: num_cells2
    integer(kind=i32)  :: num_cells_total
    integer(kind=i32)  :: ierr

    E_acc%mesh      => mesh_3d
    num_cells1      = mesh_3d%num_cells1
    num_cells2      = mesh_3d%num_cells2
    num_cells_total = num_cells1*num_cells2
    allocate( E_acc%e_acc(num_cells_total), stat= ierr);        call sll_s_test_error_code( ierr, 'Memory allocation Failure.', "sll_m_accumulators.F90",   411);
    call sll_s_reset_field_accumulator_cs_to_zero( E_acc )

  end subroutine sll_s_field_accumulator_cs_3d_init
  
  function sll_f_new_field_accumulator_cs_3d( mesh_3d ) result(E_acc)
    type(sll_t_cartesian_mesh_3d), pointer             ::  mesh_3d
    type(sll_t_electric_field_accumulator_cs), pointer ::  E_acc
    integer(kind=i32)  :: num_cells1
    integer(kind=i32)  :: num_cells2
    integer(kind=i32)  :: num_cells_total
    integer(kind=i32)  :: ierr

    if( .not. associated(mesh_3d) ) then
       print *, 'ERROR: sll_f_new_charge_accumulator_3d(), passed mesh is not ', &
            'associated. Exiting...'
       stop
    end if

    allocate( E_acc, stat= ierr);        call sll_s_test_error_code( ierr, 'Memory allocation Failure.', "sll_m_accumulators.F90",   430);

    call sll_s_field_accumulator_cs_3d_init( E_acc, mesh_3d ) 

  end function sll_f_new_field_accumulator_cs_3d
  

  subroutine sll_s_reset_field_accumulator_to_zero( E_acc )
    type(sll_t_electric_field_accumulator) :: E_acc
    integer(kind=i32) :: i
    integer(kind=i32) :: num_cells

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
    type(sll_t_electric_field_accumulator_cs) :: E_acc
    integer(kind=i32) :: i
    integer(kind=i32) :: num_cells

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
    integer(kind=i32), intent(in)  :: n_threads, n_cells
    type(sll_t_charge_accumulator_3d_ptr), dimension(:), pointer, intent(inout) :: tab
    integer(kind=i32)  :: i, j   
    
    do i = 1, n_cells  
       do j = 2, n_threads
          tab(1)%q%q_acc(i) = tab(1)%q%q_acc(i) + tab(j)%q%q_acc(i)
       enddo
    enddo
  end subroutine sll_s_sum_accumulators

  subroutine sll_s_sum_accumulators_cs( tab, n_threads, n_cells )
    integer(kind=i32), intent(in)  :: n_threads, n_cells
    type(sll_t_charge_accumulator_3d_cs_ptr), dimension(:), pointer, intent(inout) :: tab
    integer(kind=i32)  :: i, j   
    
    do i = 1, n_cells  
       do j = 2, n_threads
          tab(1)%q%q_acc(i) = tab(1)%q%q_acc(i) + tab(j)%q%q_acc(i)
       enddo
    enddo
  end subroutine sll_s_sum_accumulators_cs
  


end module sll_m_accumulators_3d
