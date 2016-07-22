module sll_m_pic_diagnostics
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_arbitrary_degree_splines, only: &
    sll_s_uniform_b_splines_at_x

  use sll_m_particle_group_base, only: &
       sll_c_particle_group_base
  
  use sll_m_particle_mesh_coupling_base, only: &
    sll_c_particle_mesh_coupling

  use sll_m_maxwell_1d_base, only: &
    sll_c_maxwell_1d_base
  
  use sll_m_collective, only: &
       sll_o_collective_allreduce, &
    sll_s_collective_reduce_real64, &
    sll_v_world_collective

  use sll_mpi, only: &
       MPI_SUM

  implicit none

  public :: &
       sll_s_pic_diagnostics_Hpi, &
       sll_s_pic_diagnostics_eval_derivative_spline, &
       sll_s_pic_diagnostics_transfer, &
       sll_s_pic_diagnostics_vvb, &
       sll_s_pic_diagnostics_poynting
  
  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
contains

  !> compute v(index)-part of kinetic energy
  subroutine sll_s_pic_diagnostics_Hpi ( particle_group,  index, kinetic )
    class(sll_c_particle_group_base), intent(in)  :: particle_group !< particle group
    sll_int32,                        intent(in)  :: index !< velocity component
    sll_real64,                       intent(out) :: kinetic(1) !< value of \a index part of kinetic energy
    

    sll_real64 :: kinetic_local(1)
    sll_int32  :: i_part
    sll_real64 :: vi(3)
    sll_real64 :: wi(1)

    kinetic_local(1) = 0.0_f64
    do i_part = 1, particle_group%n_particles
       vi = particle_group%get_v(i_part)
       wi = particle_group%get_mass(i_part)
       ! Kinetic energy
       kinetic_local(1) = kinetic_local(1) + &
            (vi(index)**2)*wi(1)
    end do
    kinetic = 0.0_f64
    call sll_s_collective_reduce_real64( sll_v_world_collective, kinetic_local, 1, &
         MPI_SUM, 0, kinetic )
    
    
  end subroutine sll_s_pic_diagnostics_Hpi

  !> Compute the spline coefficient of the derivative of some given spline expansion
  subroutine sll_s_pic_diagnostics_eval_derivative_spline( position, xmin, delta_x, n_grid, field_dofs, degree, derivative )
    sll_real64, intent( in    ) :: position(:) !< particle position
    sll_real64, intent( in    ) :: xmin !< lower boundary of the domain
    sll_real64, intent( in    ) :: delta_x !< time step 
    sll_int32,  intent( in    ) :: n_grid !< number of grid points
    sll_real64, intent( in    ) :: field_dofs(:) !< coefficients of spline representation of the field
    sll_int32,  intent( in    ) :: degree !< degree of spline
    sll_real64, intent(   out ) :: derivative !< value of the derivative
    
    sll_int32 :: i1, der_degree, ind, index
    sll_real64 :: spline_val(degree)
    sll_real64 :: xi(3)
    
    der_degree = degree-1
    
    xi(1) = (position(1) - xmin)/delta_x
    index = ceiling(xi(1))
    xi(1) = xi(1) - real(index-1, f64)
    index = index - der_degree

    call sll_s_uniform_b_splines_at_x( der_degree, xi(1), spline_val )
    
    derivative = 0.0_f64

    do i1 = 1, degree
       ind = modulo(index+i1-2, n_grid)+1
       derivative = derivative + spline_val(i1)*&
            (field_dofs(ind)-field_dofs(modulo(ind-2, n_grid)+1))
    end do

    derivative = derivative/delta_x
    

  end subroutine sll_s_pic_diagnostics_eval_derivative_spline


  !> Compute \sum(particles)w_p( v_1,p e_1(x_p) + v_2,p e_2(x_p))
  subroutine sll_s_pic_diagnostics_transfer ( particle_group, kernel_smoother_0, kernel_smoother_1, efield_dofs, transfer)
    class(sll_c_particle_group_base), intent( in   )  :: particle_group   
    class(sll_c_particle_mesh_coupling) :: kernel_smoother_0  !< Kernel smoother (order p+1)
    class(sll_c_particle_mesh_coupling) :: kernel_smoother_1  !< Kernel smoother (order p)   
    sll_real64, intent( in    ) :: efield_dofs(:,:) !< coefficients of efield
    sll_real64, intent(   out ) :: transfer(1) !< result

    sll_int32 :: i_part
    sll_real64 :: xi(3), vi(3), wi, efield(2), transfer_local(1)

    transfer_local = 0.0_f64
    do i_part = 1, particle_group%n_particles
       xi = particle_group%get_x( i_part )
       wi = particle_group%get_charge( i_part )
       vi = particle_group%get_v( i_part )

       call kernel_smoother_1%evaluate &
            (xi(1), efield_dofs(:,1), efield(1))
       call kernel_smoother_0%evaluate &
            (xi(1), efield_dofs(:,2), efield(2))

       transfer_local(1) = transfer_local(1) + (vi(1) * efield(1) + vi(2) * efield(2))*wi
       
    end do

    call sll_o_collective_allreduce( sll_v_world_collective, transfer_local, 1, MPI_SUM, transfer )
    
  end subroutine sll_s_pic_diagnostics_transfer

  !> Compute \sum(particles) w_p v_1,p b(x_p) v_2,p
  subroutine sll_s_pic_diagnostics_vvb ( particle_group, kernel_smoother_1, bfield_dofs, vvb )
    class(sll_c_particle_group_base), intent( in   )  :: particle_group   !< particle group object
    class(sll_c_particle_mesh_coupling), intent( inout ) :: kernel_smoother_1  !< Kernel smoother (order p)  
    sll_real64,           intent( in    ) :: bfield_dofs(:) !< coefficients of bfield
    sll_real64,           intent(   out ) :: vvb(1) !< result

    sll_int32 :: i_part
    sll_real64 :: xi(3), vi(3), wi, bfield, vvb_local(1)

    vvb_local = 0.0_f64
    do i_part = 1, particle_group%n_particles
       xi = particle_group%get_x( i_part )
       wi = particle_group%get_charge( i_part )
       vi = particle_group%get_v( i_part )

       call kernel_smoother_1%evaluate &
            ( xi(1), bfield_dofs, bfield )

       vvb_local = vvb_local + wi * vi(1) * vi(2) * bfield
     
    end do

    call sll_o_collective_allreduce( sll_v_world_collective, vvb_local, 1, MPI_SUM, vvb )

  end subroutine sll_s_pic_diagnostics_vvb

  !> Compute e^T M_0^{-1}  R^T b
  subroutine sll_s_pic_diagnostics_poynting ( maxwell_solver, degree, efield_dofs, bfield_dofs, scratch, poynting )
    class(sll_c_maxwell_1d_base) :: maxwell_solver !< maxwell solver object
    sll_int32, intent( in    ) :: degree !< degree of finite element
    sll_real64, intent( in    ) :: efield_dofs(:) !< coefficients of efield
    sll_real64, intent( in    ) :: bfield_dofs(:) !< coefficients of bfield
    sll_real64, intent(   out ) :: scratch(:) !< scratch data 
    sll_real64, intent(   out ) :: poynting !< value of  e^T M_0^{-1}  R^T b

    scratch = 0.0_f64
    ! Multiply B by M_0^{-1}  R^T
    call maxwell_solver%compute_e_from_b ( 1.0_f64, bfield_dofs, scratch )

    poynting =  maxwell_solver%inner_product( efield_dofs, scratch, degree )

  end subroutine sll_s_pic_diagnostics_poynting

  


end module sll_m_pic_diagnostics
