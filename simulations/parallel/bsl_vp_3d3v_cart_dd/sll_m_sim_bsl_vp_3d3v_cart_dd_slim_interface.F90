module sll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface
#include "sll_working_precision.h"
#include "sll_memory.h"

use, intrinsic :: ISO_C_Binding

use sll_m_sim_bsl_vp_3d3v_cart_dd_slim, only: &
    sll_t_sim_bsl_vp_3d3v_cart_dd_slim

use sll_m_sim_6d_utilities, only : &
    sll_s_compute_charge_density_6d_dd_slim, &
    sll_s_time_history_diagnostics, &
    sll_s_additional_time_history_diagnostics_dd,&
    sll_s_compute_momentum_energy_6d_dd, &
    sll_s_plot_diagnostic_time, &
    sll_s_plot_diagnostic_time_dd

use sll_m_poisson_3d_periodic_par, only : &
    sll_s_poisson_3d_periodic_par_compute_e_from_phi, &
    sll_s_poisson_3d_periodic_par_solve

use sll_m_remapper, only : &
    sll_o_apply_remap_3d

use sll_m_collective, only : &
    sll_s_collective_bcast_3d_real64

#ifdef _OPENMP
  use omp_lib
#endif

implicit none

private :: &
    sim_bsl_vp_3d3v_cart_dd_slim_init_f
public :: &
    sim_bsl_vp_3d3v_cart_dd_slim_allocate, &
    sim_bsl_vp_3d3v_cart_dd_slim_init, &
    sim_bsl_vp_3d3v_cart_dd_slim_run, &
!    sim_bsl_vp_3d3v_cart_dd_slim_write_diagnostics_init, &
!    sim_bsl_vp_3d3v_cart_dd_slim_write_diagnostics, &
    sim_bsl_vp_3d3v_cart_dd_slim_delete
contains

subroutine sim_bsl_vp_3d3v_cart_dd_slim_allocate(sim, sim_cptr)
    type(c_ptr), intent(out) :: sim_cptr
    type(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), pointer, intent(out) :: sim
    sll_int32 :: ierr
    SLL_ALLOCATE( sim, ierr )
    ! cf. https://stackoverflow.com/questions/57731079/pointer-casting-in-fortran
    sim_cptr = c_loc(sim)
end subroutine sim_bsl_vp_3d3v_cart_dd_slim_allocate

! inspired by https://stackoverflow.com/questions/5609502/fortran-accepting-string-from-c
subroutine sim_bsl_vp_3d3v_cart_dd_slim_init(sim_cptr, filename_in) bind(C,name='sim_bsl_vp_3d3v_cart_dd_slim_init')
    type(c_ptr), intent(out) :: sim_cptr
    type(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), pointer :: sim
    character(kind=c_char), dimension(*), intent(in) :: filename_in
    character(len=256) :: filename
    integer :: length, i
    call sim_bsl_vp_3d3v_cart_dd_slim_allocate(sim, sim_cptr)

    ! find string length
    length=0
    do
       if (filename_in(length+1) == C_NULL_CHAR) exit
       length = length + 1
    end do

    write (filename,*) (filename_in(i),i=1,length) ! inserts a leading whitespace
    ! print *, 'In Fortran, got string: ', (filename_in(i),i=1,length), '(',length,').'
    ! print *, 'In Fortran, got string: ',filename, '().'

    call sim_bsl_vp_3d3v_cart_dd_slim_init_f(sim, length, filename(2:))
end subroutine sim_bsl_vp_3d3v_cart_dd_slim_init

subroutine sim_bsl_vp_3d3v_cart_dd_slim_init_f(sim, length, filename)
    type(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), pointer, intent(in) :: sim
    integer, intent(in) :: length
    character(len=length), intent(in) :: filename

    ! print *, 'that is,', filename, '(',length,').'
    call sim%init_from_file(filename)
end subroutine sim_bsl_vp_3d3v_cart_dd_slim_init_f

subroutine sim_bsl_vp_3d3v_cart_dd_slim_run(sim_cptr) bind(C,name='sim_bsl_vp_3d3v_cart_dd_slim_run')
    type(c_ptr), intent(in) :: sim_cptr
    type(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), pointer :: sim

    call c_f_pointer(sim_cptr, sim)

    call sim%run()
    sim%ctest = .false.
end subroutine sim_bsl_vp_3d3v_cart_dd_slim_run

subroutine sim_bsl_vp_3d3v_cart_dd_slim_delete(sim_cptr) bind(C,name='sim_bsl_vp_3d3v_cart_dd_slim_delete')
    type(c_ptr), intent(in) :: sim_cptr
    type(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), pointer :: sim

    call c_f_pointer(sim_cptr, sim)

    call sim%delete()
end subroutine sim_bsl_vp_3d3v_cart_dd_slim_delete

subroutine sim_bsl_vp_3d3v_cart_dd_slim_get_distribution(sim_cptr, discotec_pointer) bind(C,name='sim_bsl_vp_3d3v_cart_dd_slim_get_distribution')
    type(c_ptr), intent(in) :: sim_cptr
    type(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), pointer :: sim

    type(C_ptr), intent( out ) :: discotec_pointer
    sll_real64, pointer:: distribution(:,:,:,:,:,:)

    call c_f_pointer(sim_cptr, sim)

    call sim%get_distribution(distribution)
    discotec_pointer = c_loc(distribution)
end subroutine sim_bsl_vp_3d3v_cart_dd_slim_get_distribution

subroutine sim_bsl_vp_3d3v_cart_dd_slim_set_distribution(sim_cptr, discotec_pointer) bind(C,name='sim_bsl_vp_3d3v_cart_dd_slim_set_distribution')
    type(c_ptr), intent(in) :: sim_cptr
    type(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), pointer :: sim
    type(C_ptr), value, intent( in ) :: discotec_pointer
    sll_real64, pointer :: distribution(:,:,:,:,:,:)
    sll_int32 :: local_size(6)

    call c_f_pointer(sim_cptr, sim)
    call sim%get_local_size(local_size)

    call c_f_pointer(discotec_pointer, distribution, local_size)
    call sim%set_distribution(distribution)

end subroutine sim_bsl_vp_3d3v_cart_dd_slim_set_distribution

subroutine sim_bsl_vp_3d3v_cart_dd_slim_get_local_size(sim_cptr, local_size) bind(C,name='sim_bsl_vp_3d3v_cart_dd_slim_get_local_size')
    type(c_ptr), intent(in) :: sim_cptr
    type(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), pointer :: sim
    integer(c_int32_t), dimension(6), intent( out ) :: local_size
    call c_f_pointer(sim_cptr, sim)
    call sim%get_local_size(local_size)
end subroutine sim_bsl_vp_3d3v_cart_dd_slim_get_local_size


subroutine sim_bsl_vp_3d3v_cart_dd_slim_print_etas(sim_cptr) bind(C,name='sim_bsl_vp_3d3v_cart_dd_slim_print_etas')
    type(c_ptr), intent(in) :: sim_cptr
    type(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), pointer :: sim
    sll_int32, dimension(6) :: local_size
    sll_int32 :: n

    call c_f_pointer(sim_cptr, sim)

    call sim%get_local_size(local_size)
    do n=1, local_size(1)
        print *, 'etas x1 ',sim%etas(1)%vals(n)
    end do
    do n=1, local_size(2)
        print *, 'etas x2 ',sim%etas(2)%vals(n)
    end do
    do n=1, local_size(3)
        print *, 'etas x3 ',sim%etas(3)%vals(n)
    end do
    do n=1, local_size(4)
        print *, 'etas v1 ',sim%etas(4)%vals(n)
    end do
    do n=1, local_size(5)
        print *, 'etas v2 ',sim%etas(5)%vals(n)
    end do
    do n=1, local_size(6)
        print *, 'etas v3 ',sim%etas(6)%vals(n)
    end do
end subroutine sim_bsl_vp_3d3v_cart_dd_slim_print_etas


subroutine sim_bsl_vp_3d3v_cart_dd_slim_write_diagnostics_init(sim_cptr) bind(C,name='sim_bsl_vp_3d3v_cart_dd_slim_write_diagnostics_init')
    type(c_ptr), intent(in) :: sim_cptr
    type(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), pointer :: sim
    call c_f_pointer(sim_cptr, sim)

    ! Solve Poisson
    call sll_s_compute_charge_density_6d_dd_slim(sim%f6d, &
                sim%decomposition, &
                sim%rho, &
                sim%collective_3d_velocity, &
                sim%mesh6d%volume_eta456)

    ! Compute electric fields
    if (all(sim%topology_3d_velocity%coords == 0)) then
        ! call the existing Poisson solver
        call sll_s_poisson_3d_periodic_par_solve( sim%poisson, sim%rho, sim%phi )
        call sll_s_poisson_3d_periodic_par_compute_e_from_phi( sim%poisson, sim%phi, &
            sim%ex, sim%ey, sim%ez)
    end if

    call sll_s_collective_bcast_3d_real64(sim%collective_3d_velocity, sim%ex, 0)
    call sll_s_collective_bcast_3d_real64(sim%collective_3d_velocity, sim%ey, 0)
    call sll_s_collective_bcast_3d_real64(sim%collective_3d_velocity, sim%ez, 0)
!         sim%ex = 0.
!         sim%ey = 0.01*(erf(itime/1000.-2.)+1.)
!         sim%ez = 0.

    call sll_s_time_history_diagnostics(&
            sim%decomposition%local%mn, &
            sim%decomposition%local%mx, &
            sim%decomposition%local%mn, &  ! previously: `sim%decomposition%local%lo`
            real(sim%first_time_step-1, f64)*sim%delta_t, &
            sim%mesh6d%volume, &
            sim%mesh6d%volume_eta123, &
            sim%etas,&
            sim%f6d, &
            sim%rho, &
            sim%phi, &
            sim%ex, &
            sim%ey, &
            sim%ez, &
            sim%thdiag_file_id, &
            sim%topology_3d_velocity%coords)
end subroutine sim_bsl_vp_3d3v_cart_dd_slim_write_diagnostics_init

subroutine sim_bsl_vp_3d3v_cart_dd_slim_write_diagnostics(sim_cptr, itime) bind(C,name='sim_bsl_vp_3d3v_cart_dd_slim_write_diagnostics')
    type(c_ptr), intent(in) :: sim_cptr
    integer(c_int32_t), intent(in) :: itime ! the number of the time step that we are doing the diagnostic for
    type(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), pointer :: sim
    call c_f_pointer(sim_cptr, sim)

    ! Solve Poisson
    call sll_s_compute_charge_density_6d_dd_slim(sim%f6d, &
                sim%decomposition, &
                sim%rho, &
                sim%collective_3d_velocity, &
                sim%mesh6d%volume_eta456)

    ! Compute electric fields
    if (all(sim%topology_3d_velocity%coords == 0)) then
        ! call the existing Poisson solver
        call sll_s_poisson_3d_periodic_par_solve( sim%poisson, sim%rho, sim%phi )
        call sll_s_poisson_3d_periodic_par_compute_e_from_phi( sim%poisson, sim%phi, &
            sim%ex, sim%ey, sim%ez)
    end if

    call sll_s_collective_bcast_3d_real64(sim%collective_3d_velocity, sim%ex, 0)
    call sll_s_collective_bcast_3d_real64(sim%collective_3d_velocity, sim%ey, 0)
    call sll_s_collective_bcast_3d_real64(sim%collective_3d_velocity, sim%ez, 0)

    call sll_s_time_history_diagnostics(&
            sim%decomposition%local%mn, &
            sim%decomposition%local%mx, &
            sim%decomposition%local%mn, &  ! previously: `sim%decomposition%local%lo`
            real(itime, f64)*sim%delta_t, &
            sim%mesh6d%volume, &
            sim%mesh6d%volume_eta123, &
            sim%etas,&
            sim%f6d, &
            sim%rho, &
            sim%phi, &
            sim%ex, &
            sim%ey, &
            sim%ez, &
            sim%thdiag_file_id, &
            sim%topology_3d_velocity%coords)

end subroutine sim_bsl_vp_3d3v_cart_dd_slim_write_diagnostics

subroutine sim_bsl_vp_3d3v_cart_dd_slim_advect_v(sim_cptr, delta_t) bind(C,name='sim_bsl_vp_3d3v_cart_dd_slim_advect_v')
    type(c_ptr), intent(in) :: sim_cptr
    type(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), pointer :: sim
    real(c_double), intent( in ) :: delta_t
    call c_f_pointer(sim_cptr, sim)
    call sim%advect_v(delta_t)
end subroutine sim_bsl_vp_3d3v_cart_dd_slim_advect_v

subroutine sim_bsl_vp_3d3v_cart_dd_slim_advect_x(sim_cptr) bind(C,name='sim_bsl_vp_3d3v_cart_dd_slim_advect_x')
    type(c_ptr), intent(in) :: sim_cptr
    type(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), pointer :: sim
!    real(c_double), intent( in ) :: delta_t
    call c_f_pointer(sim_cptr, sim)
    call sim%advect_x()
end subroutine sim_bsl_vp_3d3v_cart_dd_slim_advect_x

end module sll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface
