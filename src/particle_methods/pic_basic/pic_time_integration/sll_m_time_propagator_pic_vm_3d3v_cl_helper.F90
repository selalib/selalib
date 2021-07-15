!> @ingroup pic_time_integration
!> @author Benedikt Perse, IPP
!> @brief Particle pusher based on antisymmetric splitting with AVF for 3d3v Vlasov-Maxwell with coordinate transformation.
!> @details MPI parallelization by domain cloning. General boundaries. Spline DoFs numerated by the point the spline starts.
module sll_m_time_propagator_pic_vm_3d3v_cl_helper
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_time_propagator_base, only: &
       sll_c_time_propagator_base

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  use sll_m_particle_mesh_coupling_base_3d, only: &
       sll_c_particle_mesh_coupling_3d


  implicit none

  public :: &
       sll_t_time_propagator_pic_vm_3d3v_cl_helper, &
       sll_p_boundary_particles_periodic, &
       sll_p_boundary_particles_singular, &
       sll_p_boundary_particles_reflection, &
       sll_p_boundary_particles_absorption, &
       sll_s_compute_particle_boundary_simple, &
       sll_s_compute_particle_boundary_trafo, &
       sll_s_compute_particle_boundary_trafo_current, &
       sll_s_compute_matrix_inverse

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: sll_p_boundary_particles_periodic = 0
  sll_int32, parameter :: sll_p_boundary_particles_singular = 1
  sll_int32, parameter :: sll_p_boundary_particles_reflection = 2
  sll_int32, parameter :: sll_p_boundary_particles_absorption = 3

  !> Hamiltonian splitting type for Vlasov-Maxwell 3d3v
  type  :: sll_t_time_propagator_pic_vm_3d3v_cl_helper
     class(sll_c_particle_mesh_coupling_3d), pointer :: particle_mesh_coupling !< Particle mesh coupling
     type( sll_t_mapping_3d ), pointer         :: map !< coordinate transformation

     sll_real64, allocatable :: j_dofs_local(:)!< MPI-processor local part of one component of \a j_dofs
     sll_real64, allocatable :: efield_dofs_work(:) !< scratch data

     sll_int32 :: spline_degree(3) !< Degree of the spline for j,B. Here 3.
     sll_real64 :: Lx(3) !< Size of the domain
     sll_real64 :: x_min(3) !< Lower bound for x domain
     sll_real64 :: x_max(3)
     sll_int32 :: boundary_particles = sll_p_boundary_particles_periodic


     sll_real64, pointer :: rhob(:) => null()
     sll_int32 :: counter_left = 0
     sll_int32 :: counter_right = 0

  end type sll_t_time_propagator_pic_vm_3d3v_cl_helper

contains

  !> invert matrix 
  subroutine sll_s_compute_matrix_inverse( x, y, bf, jm, sign)
    sll_real64, intent( in    ) :: x(:) !< Inputvariable
    sll_real64, intent(   out ) :: y(:) !< Outputvariable
    sll_real64, intent( in    ) :: bf(3) !< bfield
    sll_real64, intent( in    ) :: jm(3,3) !< Inverse of jacobian matrix
    sll_real64, intent( in    ) :: sign !< sign
    !local variables 
    sll_real64 :: a,b,c,d,e,f,g,h,i
    sll_real64 :: det

    a = 1-sign*( jm(1,1)*(jm(2,1)*bf(3)-jm(3,1)*bf(2)) + jm(2,1)*(jm(3,1)*bf(1)-jm(1,1)*bf(3)) + jm(3,1)*(jm(1,1)*bf(2)-jm(2,1)*bf(1)) )
    d =  -sign*( jm(1,2)*(jm(2,1)*bf(3)-jm(3,1)*bf(2)) + jm(2,2)*(jm(3,1)*bf(1)-jm(1,1)*bf(3)) + jm(3,2)*(jm(1,1)*bf(2)-jm(2,1)*bf(1)) )
    g =  -sign*( jm(1,3)*(jm(2,1)*bf(3)-jm(3,1)*bf(2)) + jm(2,3)*(jm(3,1)*bf(1)-jm(1,1)*bf(3)) + jm(3,3)*(jm(1,1)*bf(2)-jm(2,1)*bf(1)) )

    b =  -sign*( jm(1,1)*(jm(2,2)*bf(3)-jm(3,2)*bf(2)) + jm(2,1)*(jm(3,2)*bf(1)-jm(1,2)*bf(3)) + jm(3,1)*(jm(1,2)*bf(2)-jm(2,2)*bf(1)) )
    e = 1-sign*( jm(1,2)*(jm(2,2)*bf(3)-jm(3,2)*bf(2)) + jm(2,2)*(jm(3,2)*bf(1)-jm(1,2)*bf(3)) + jm(3,2)*(jm(1,2)*bf(2)-jm(2,2)*bf(1)) )
    h =  -sign*( jm(1,3)*(jm(2,2)*bf(3)-jm(3,2)*bf(2)) + jm(2,3)*(jm(3,2)*bf(1)-jm(1,2)*bf(3)) + jm(3,3)*(jm(1,2)*bf(2)-jm(2,2)*bf(1)) )

    c =  -sign*( jm(1,1)*(jm(2,3)*bf(3)-jm(3,3)*bf(2)) + jm(2,1)*(jm(3,3)*bf(1)-jm(1,3)*bf(3)) + jm(3,1)*(jm(1,3)*bf(2)-jm(2,3)*bf(1)) )
    f =  -sign*( jm(1,2)*(jm(2,3)*bf(3)-jm(3,3)*bf(2)) + jm(2,2)*(jm(3,3)*bf(1)-jm(1,3)*bf(3)) + jm(3,2)*(jm(1,3)*bf(2)-jm(2,3)*bf(1)) )
    i = 1-sign*( jm(1,3)*(jm(2,3)*bf(3)-jm(3,3)*bf(2)) + jm(2,3)*(jm(3,3)*bf(1)-jm(1,3)*bf(3)) + jm(3,3)*(jm(1,3)*bf(2)-jm(2,3)*bf(1)) )

    det = a*e*i + b*f*g + c*d*h - c*e*g - a*f*h - b*d*i

    y(1) = ( x(1)*(e*i - f*h) + x(2)*(c*h - b*i) + x(3)*(b*f - c*e) )/det
    y(2) = ( x(1)*(f*g - d*i) + x(2)*(a*i - c*g) + x(3)*(c*d - a*f) )/det
    y(3) = ( x(1)*(d*h - e*g) + x(2)*(b*g - a*h) + x(3)*(a*e - b*d) )/det

  end subroutine sll_s_compute_matrix_inverse


  !> Compute particle boundary
  subroutine compute_particle_boundary( self, xold, xnew, vi  )
    class(sll_t_time_propagator_pic_vm_3d3v_cl_helper), intent( inout ) :: self !< time propagator object 
    sll_real64,                                           intent( inout ) :: xold(3)
    sll_real64,                                           intent( inout ) :: xnew(3)
    sll_real64,                                           intent( inout ) :: vi(3)
    !local variables
    sll_real64 :: xbar

    if(xnew(1) < self%x_min(1) .or. xnew(1) > self%x_max(1) )then
       if(xnew(1) < self%x_min(1)  )then
          xbar = self%x_min(1)
          self%counter_left = self%counter_left+1
       else if(xnew(1) > self%x_max(1))then
          xbar = self%x_max(1)
          self%counter_right = self%counter_right+1
       end if


       select case(self%boundary_particles)
       case(sll_p_boundary_particles_reflection) 
          xnew(1) = 2._f64*xbar-xnew(1)
          vi(1) = -vi(1)
       case(sll_p_boundary_particles_absorption)
       case( sll_p_boundary_particles_periodic) 
          xnew(1) = self%x_min(1) + modulo(xnew(1)-self%x_min(1), self%Lx(1))
       case default
          print*,'error: boundary case missing', self%boundary_particles
       end select
    end if
    xnew(2:3) = self%x_min(2:3) + modulo(xnew(2:3)-self%x_min(2:3), self%Lx(2:3))

  end subroutine compute_particle_boundary


  !> compute particle boundary with coordinate transformation
  subroutine compute_particle_boundary_trafo( self, xold, xnew, vi  )
    class(sll_t_time_propagator_pic_vm_3d3v_cl_helper), intent( inout ) :: self !< time propagator object 
    sll_real64,                                           intent( inout ) :: xold(3)
    sll_real64,                                           intent( inout ) :: xnew(3)
    sll_real64,                                           intent( inout ) :: vi(3)
    !local variables
    sll_real64 :: xmid(3), xt(3), xbar, dx
    sll_real64 :: jmatrix(3,3)

    if(xnew(1) < 0._f64 .or. xnew(1) > 1._f64 )then
       if(xnew(1) < 0._f64  )then
          xbar = 0._f64
          self%counter_left = self%counter_left+1
       else if(xnew(1) > 1._f64)then
          xbar = 1._f64
          self%counter_right = self%counter_right+1
       end if
       dx = (xbar- xold(1))/(xnew(1)-xold(1))
       xmid = xold + dx * (xnew-xold)
       xmid(1) = xbar

       select case(self%boundary_particles)
       case(sll_p_boundary_particles_singular)
          if(xnew(1) < 0._f64 )then
             xnew(1) = -xnew(1)
             xnew(2) = xnew(2) + 0.5_f64
          else if(xnew(1) > 1._f64 )then
             jmatrix = self%map%jacobian_matrix_inverse_transposed(xmid)
             vi = vi-2._f64*(vi(1)*jmatrix(1,1)+vi(2)*jmatrix(2,1)+vi(3)*jmatrix(3,1))*jmatrix(:,1)/sum(jmatrix(:,1)**2)
             xnew(1) = 2._f64 - xnew(1)
             xnew(2) = 1._f64 - xnew(2)
          end if
       case(sll_p_boundary_particles_reflection)
          xt = xmid
          xt(2:3) = modulo(xt(2:3),1._f64)
          jmatrix = self%map%jacobian_matrix_inverse_transposed(xt)
          vi = vi-2._f64*(vi(1)*jmatrix(1,1)+vi(2)*jmatrix(2,1)+vi(3)*jmatrix(3,1))*jmatrix(:,1)/sum(jmatrix(:,1)**2)
          xnew(1) = 2._f64*xbar-xnew(1)
       case(sll_p_boundary_particles_absorption)
       case( sll_p_boundary_particles_periodic)
          xnew(1) = modulo(xnew(1), 1._f64)
       case default
          print*,'error: boundary case missing', self%boundary_particles
       end select
    else if(xnew(1) < -1._f64 .or.  xnew(1) > 2._f64 ) then
       print*, xnew
       SLL_ERROR("particle boundary", "particle out of bound")
    end if
    xnew(2:3) = modulo(xnew(2:3), 1._f64)

  end subroutine compute_particle_boundary_trafo


  !> compute particle boundary and current with coordinate transformation
  subroutine compute_particle_boundary_trafo_current( self, xold, xnew, vi, wi )
    class(sll_t_time_propagator_pic_vm_3d3v_cl_helper), intent( inout ) :: self !< time propagator object 
    sll_real64,                                           intent( in    ) :: xold(3)
    sll_real64,                                           intent( inout ) :: xnew(3)
    sll_real64,                                           intent( inout ) :: vi(3)
    sll_real64,                                           intent( in    ) :: wi(1)
    !local variables
    sll_real64 :: xmid(3), xt(3), vt(3), xbar, dx
    sll_real64 :: jmatrix(3,3)

    if(xnew(1) < 0._f64 .or. xnew(1) > 1._f64 )then
       if(xnew(1) < 0._f64  )then
          xbar = 0._f64
          self%counter_left = self%counter_left+1
       else if(xnew(1) > 1._f64)then
          xbar = 1._f64
          self%counter_right = self%counter_right+1
       end if

       dx = (xbar- xold(1))/(xnew(1)-xold(1))
       xmid = xold + dx * (xnew-xold)
       xmid(1) = xbar
       vt = (xmid - xold)*wi(1)
       call self%particle_mesh_coupling%add_current( xold, xmid, vt, self%j_dofs_local )
       select case(self%boundary_particles)
       case(sll_p_boundary_particles_singular)
          call self%particle_mesh_coupling%add_charge(xmid, wi(1), self%spline_degree, self%rhob)
          xmid(2) = xbar + (1._f64-2._f64*xbar)*xmid(2) + 0.5_f64-0.5_f64*xbar
          call self%particle_mesh_coupling%add_charge(xmid, -wi(1), self%spline_degree, self%rhob)
          xnew(1) = 2._f64*xbar-xnew(1)
          xnew(2) = xbar + (1._f64-2._f64*xbar)*xnew(2) + 0.5_f64-0.5_f64*xbar

          if(xnew(1) > 1._f64 )then
             xt = xmid
             xt(2:3) = modulo(xt(2:3),1._f64)
             jmatrix = self%map%jacobian_matrix_inverse_transposed(xt)
             vi = vi-2._f64*(vi(1)*jmatrix(1,1)+vi(2)*jmatrix(2,1)+vi(3)*jmatrix(3,1))*jmatrix(:,1)/sum(jmatrix(:,1)**2)
          end if
       case(sll_p_boundary_particles_reflection)
          xt = xmid
          xt(2:3) = modulo(xt(2:3),1._f64)
          jmatrix = self%map%jacobian_matrix_inverse_transposed(xt)
          vi = vi-2._f64*(vi(1)*jmatrix(1,1)+vi(2)*jmatrix(2,1)+vi(3)*jmatrix(3,1))*jmatrix(:,1)/sum(jmatrix(:,1)**2)
          xnew(1) = 2._f64*xbar-xnew(1)
       case(sll_p_boundary_particles_absorption)
          call self%particle_mesh_coupling%add_charge(xmid, wi(1), self%spline_degree, self%rhob)
       case( sll_p_boundary_particles_periodic)
          xnew(1) = modulo(xnew(1), 1._f64)
          xmid(1) = 1._f64-xbar
       case default
          print*,'error: boundary case missing', self%boundary_particles
       end select
       if(xnew(1) >= 0._f64 .and. xnew(1) <= 1._f64) then
          vt = (xnew - xmid)*wi(1)
          call self%particle_mesh_coupling%add_current( xmid, xnew, vt, self%j_dofs_local )
       end if
    else if(xnew(1) < -1._f64 .or. xnew(1) > 2._f64)then
       print*, xnew
       SLL_ERROR("particle boundary", "particle out of bound")
    else
       vt = (xnew - xold)*wi(1)
       call self%particle_mesh_coupling%add_current( xold, xnew, vt, self%j_dofs_local )
    end if
    xnew(2:3) = modulo(xnew(2:3), 1._f64)

  end subroutine compute_particle_boundary_trafo_current


  !> compute particle boundary and current with coordinate transformation and evaluate efield
  subroutine compute_particle_boundary_current_trafo_evaluate( self, xold, xnew, vi, wi, sign )
    class(sll_t_time_propagator_pic_vm_3d3v_cl_helper), intent( inout ) :: self !< time propagator object 
    sll_real64,                                           intent( in    ) :: xold(3)
    sll_real64,                                           intent( inout ) :: xnew(3)
    sll_real64,                                           intent( inout ) :: vi(3)
    sll_real64,                                           intent( in    ) :: wi(1)
    sll_real64,                                           intent( in    ) :: sign
    !local variables
    sll_real64 :: xmid(3), xt(3), vt(3), xbar, dx
    sll_real64 :: jmatrix(3,3), jmat(3,3), efield(3)
    sll_int32 :: j

    if(xnew(1) < 0._f64 .or. xnew(1) > 1._f64 )then
       if(xnew(1) < 0._f64  )then
          xbar = 0._f64
          self%counter_left = self%counter_left+1
       else if(xnew(1) > 1._f64)then
          xbar = 1._f64
          self%counter_right = self%counter_right+1
       end if
       dx = (xbar- xold(1))/(xnew(1)-xold(1))
       xmid = xold + dx * (xnew-xold)
       xmid(1) = xbar
       vt = (xmid - xold)*wi(1)
       call self%particle_mesh_coupling%add_current_evaluate( xold, xmid, vt, self%efield_dofs_work, &
            self%j_dofs_local, efield )
       jmat = self%map%jacobian_matrix_inverse_transposed(xold)
       xt = xmid
       xt(2:3) = modulo(xt(2:3), 1._f64)
       jmatrix = self%map%jacobian_matrix_inverse_transposed(xt)
       do j = 1, 3
          vi(j) = vi(j) + dx*sign *0.5_f64*((jmatrix(j,1)+jmat(j,1))*efield(1) + (jmatrix(j,2)+jmat(j,2))*efield(2) + (jmatrix(j,3)+ jmat(j,3))*efield(3))
       end do
       select case(self%boundary_particles)
       case(sll_p_boundary_particles_singular)
          call self%particle_mesh_coupling%add_charge(xmid, wi(1), self%spline_degree, self%rhob)
          xmid(2) = xbar + (1._f64-2._f64*xbar)*xmid(2) + 0.5_f64-0.5_f64*xbar
          xt(2:3) = modulo(xt(2:3), 1._f64)
          jmatrix = self%map%jacobian_matrix_inverse_transposed(xt)
          call self%particle_mesh_coupling%add_charge(xmid, -wi(1), self%spline_degree, self%rhob)
          xnew(1) = 2._f64*xbar-xnew(1)
          xnew(2) = xbar + (1._f64-2._f64*xbar)*xnew(2) + 0.5_f64-0.5_f64*xbar
          if(xnew(1) > 1._f64 )then
             vi = vi-2._f64*(vi(1)*jmatrix(1,1)+vi(2)*jmatrix(2,1)+vi(3)*jmatrix(3,1))*jmatrix(:,1)/sum(jmatrix(:,1)**2)
          end if
       case(sll_p_boundary_particles_reflection)
          vi = vi-2._f64*(vi(1)*jmatrix(1,1)+vi(2)*jmatrix(2,1)+vi(3)*jmatrix(3,1))*jmatrix(:,1)/sum(jmatrix(:,1)**2)
          xnew(1) = 2._f64*xbar-xnew(1)
       case(sll_p_boundary_particles_absorption)
          call self%particle_mesh_coupling%add_charge(xmid, wi(1), self%spline_degree, self%rhob)
       case( sll_p_boundary_particles_periodic)
          xnew(1) = modulo(xnew(1), 1._f64)
          xmid(1) = 1._f64-xbar
       case default
          print*,'error: boundary case missing', self%boundary_particles
       end select
       vt = (xnew - xmid)*wi(1)
       call self%particle_mesh_coupling%add_current_evaluate( xmid, xnew, vt, self%efield_dofs_work, &
            self%j_dofs_local, efield )
       xnew(2:3) = modulo(xnew(2:3), 1._f64)
       jmat = self%map%jacobian_matrix_inverse_transposed(xnew)
       do j = 1, 3
          vi(j) = vi(j) + (1._f64-dx)*sign *0.5_f64*((jmatrix(j,1)+jmat(j,1))*efield(1) + (jmatrix(j,2)+jmat(j,2))*efield(2) + (jmatrix(j,3)+ jmat(j,3))*efield(3))
       end do
    else if(xnew(1) < -1._f64 .or. xnew(1) > 2._f64)then
       print*, xnew
       SLL_ERROR("particle boundary", "particle out of bound")
    else   
       vt = (xnew - xold)*wi(1)
       call self%particle_mesh_coupling%add_current_evaluate( xold, xnew, vt, self%efield_dofs_work, &
            self%j_dofs_local, efield )
       jmat = self%map%jacobian_matrix_inverse_transposed(xold)
       xnew(2:3) = modulo(xnew(2:3), 1._f64)
       jmatrix = self%map%jacobian_matrix_inverse_transposed(xnew)
       do j = 1, 3
          vi(j) = vi(j) + sign *0.5_f64*((jmatrix(j,1)+jmat(j,1))*efield(1) + (jmatrix(j,2)+jmat(j,2))*efield(2) + (jmatrix(j,3)+ jmat(j,3))*efield(3))
       end do
    end if




  end subroutine compute_particle_boundary_current_trafo_evaluate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  !> compute new position
  subroutine sll_s_compute_particle_boundary_simple( boundary_particles, counter_left, counter_right, xold, xnew  )
    sll_int32,  intent( in    ) :: boundary_particles 
    sll_int32,  intent( inout ) :: counter_left 
    sll_int32,  intent( inout ) :: counter_right 
    sll_real64, intent( inout ) :: xold(3)
    sll_real64, intent( inout ) :: xnew(3)
    !local variables
    sll_real64 :: xbar

    if(xnew(1) < -1._f64 .or.  xnew(1) > 2._f64 ) then
       print*, xnew
       SLL_ERROR("particle boundary", "particle out of bound")
    else if(xnew(1) < 0._f64 .or. xnew(1) > 1._f64 )then
       if(xnew(1) < 0._f64  )then
          xbar = 0._f64
          counter_left = counter_left+1
       else if(xnew(1) > 1._f64)then
          xbar = 1._f64
          counter_right = counter_right+1
       end if
       select case(boundary_particles)
       case(sll_p_boundary_particles_singular)
          if(xnew(1) < 0._f64 )then
             xnew(1) = -xnew(1)
             xnew(2) = xnew(2) + 0.5_f64
          else if(xnew(1) > 1._f64 )then
             xnew(1) = 2._f64 - xnew(1)
             xnew(2) = 1._f64 - xnew(2)
          end if
       case(sll_p_boundary_particles_reflection)
          xnew(1) = 2._f64*xbar-xnew(1)
       case(sll_p_boundary_particles_absorption)
       case( sll_p_boundary_particles_periodic)
          xnew(1) = modulo(xnew(1), 1._f64)
       case default
          xnew(1) = modulo(xnew(1), 1._f64)
       end select
    end if
    xnew(2:3) = modulo(xnew(2:3), 1._f64)

  end subroutine sll_s_compute_particle_boundary_simple


  !> Compute particle boundary with coordinate transformation
  subroutine sll_s_compute_particle_boundary_trafo( boundary_particles, counter_left, counter_right, map, xold, xnew, vi  )
    sll_int32,  intent( in    ) :: boundary_particles 
    sll_int32,  intent( inout ) :: counter_left 
    sll_int32,  intent( inout ) :: counter_right
    type( sll_t_mapping_3d ), intent( inout ) :: map
    sll_real64,               intent( inout ) :: xold(3)
    sll_real64,               intent( inout ) :: xnew(3)
    sll_real64,               intent( inout ) :: vi(3)
    !local variables
    sll_real64 :: xmid(3), xt(3), xbar, dx
    sll_real64 :: jmatrix(3,3)

    if(xnew(1) < 0._f64 .or. xnew(1) > 1._f64 )then
       if(xnew(1) < 0._f64  )then
          xbar = 0._f64
          counter_left = counter_left+1
       else if(xnew(1) > 1._f64)then
          xbar = 1._f64
          counter_right = counter_right+1
       end if
       dx = (xbar- xold(1))/(xnew(1)-xold(1))
       xmid = xold + dx * (xnew-xold)
       xmid(1) = xbar

       select case(boundary_particles)
       case(sll_p_boundary_particles_singular)
          if(xnew(1) < 0._f64 )then
             xnew(1) = -xnew(1)
             xnew(2) = xnew(2) + 0.5_f64
          else if(xnew(1) > 1._f64 )then
             jmatrix = map%jacobian_matrix_inverse_transposed(xmid)
             vi = vi-2._f64*(vi(1)*jmatrix(1,1)+vi(2)*jmatrix(2,1)+vi(3)*jmatrix(3,1))*jmatrix(:,1)/sum(jmatrix(:,1)**2)
             xnew(1) = 2._f64 - xnew(1)
             xnew(2) = 1._f64 - xnew(2)
          end if
       case(sll_p_boundary_particles_reflection)
          xt = xmid
          xt(2:3) = modulo(xt(2:3),1._f64)
          jmatrix = map%jacobian_matrix_inverse_transposed(xt)
          vi = vi-2._f64*(vi(1)*jmatrix(1,1)+vi(2)*jmatrix(2,1)+vi(3)*jmatrix(3,1))*jmatrix(:,1)/sum(jmatrix(:,1)**2)
          xnew(1) = 2._f64*xbar-xnew(1)
       case(sll_p_boundary_particles_absorption)
       case( sll_p_boundary_particles_periodic)
          xnew(1) = modulo(xnew(1), 1._f64)
       case default
          xnew(1) = modulo(xnew(1), 1._f64)
       end select
    else if(xnew(1) < -1._f64 .or.  xnew(1) > 2._f64 ) then
       print*, xnew
       SLL_ERROR("particle boundary", "particle out of bound")
    end if
    xnew(2:3) = modulo(xnew(2:3), 1._f64)

  end subroutine sll_s_compute_particle_boundary_trafo


  !> Compute particle boundary and current with coordinate transformation
  subroutine sll_s_compute_particle_boundary_trafo_current(boundary_particles, counter_left, counter_right, map, particle_mesh_coupling, j_dofs_local, spline_degree, rhob, xold, xnew, vi, wi )
    sll_int32,  intent( in    ) :: boundary_particles 
    sll_int32,  intent( inout ) :: counter_left 
    sll_int32,  intent( inout ) :: counter_right
    type( sll_t_mapping_3d ), intent( inout ) :: map
    class(sll_c_particle_mesh_coupling_3d), intent( inout ) :: particle_mesh_coupling
    sll_real64, intent( inout ) :: j_dofs_local(:)
    sll_int32, intent( in    ) :: spline_degree(3)
    sll_real64, intent( inout ) :: rhob(:)
    sll_real64,                                           intent( in    ) :: xold(3)
    sll_real64,                                           intent( inout ) :: xnew(3)
    sll_real64,                                           intent( inout ) :: vi(3)
    sll_real64,                                           intent( in    ) :: wi(1)
    !local variables
    sll_real64 :: xmid(3), xt(3), vh(3), xbar, dx
    sll_real64 :: jmatrix(3,3)

    if(xnew(1) < 0._f64 .or. xnew(1) > 1._f64 )then
       if(xnew(1) < 0._f64  )then
          xbar = 0._f64
          counter_left = counter_left+1
       else if(xnew(1) > 1._f64)then
          xbar = 1._f64
          counter_right = counter_right+1
       end if

       dx = (xbar- xold(1))/(xnew(1)-xold(1))
       xmid = xold + dx * (xnew-xold)
       xmid(1) = xbar
       vh = (xmid - xold)*wi(1)
       call particle_mesh_coupling%add_current( xold, xmid, vh, j_dofs_local )
       select case(boundary_particles)
       case(sll_p_boundary_particles_singular)
          call particle_mesh_coupling%add_charge(xmid, wi(1), spline_degree, rhob)
          xmid(2) = xbar + (1._f64-2._f64*xbar)*xmid(2) + 0.5_f64-0.5_f64*xbar
          call particle_mesh_coupling%add_charge(xmid, -wi(1), spline_degree, rhob)
          xnew(1) = 2._f64*xbar-xnew(1)
          xnew(2) = xbar + (1._f64-2._f64*xbar)*xnew(2) + 0.5_f64-0.5_f64*xbar

          if(xnew(1) > 1._f64 )then
             xt = xmid
             xt(2:3) = modulo(xt(2:3),1._f64)
             jmatrix = map%jacobian_matrix_inverse_transposed(xt)
             vi = vi-2._f64*(vi(1)*jmatrix(1,1)+vi(2)*jmatrix(2,1)+vi(3)*jmatrix(3,1))*jmatrix(:,1)/sum(jmatrix(:,1)**2)
          end if
       case(sll_p_boundary_particles_reflection)
          xt = xmid
          xt(2:3) = modulo(xt(2:3),1._f64)
          jmatrix = map%jacobian_matrix_inverse_transposed(xt)
          vi = vi-2._f64*(vi(1)*jmatrix(1,1)+vi(2)*jmatrix(2,1)+vi(3)*jmatrix(3,1))*jmatrix(:,1)/sum(jmatrix(:,1)**2)
          xnew(1) = 2._f64*xbar-xnew(1)
       case(sll_p_boundary_particles_absorption)
          call particle_mesh_coupling%add_charge(xmid, wi(1), spline_degree, rhob)
       case( sll_p_boundary_particles_periodic)
          xnew(1) = modulo(xnew(1), 1._f64)
          xmid(1) = 1._f64-xbar
       case default
          xnew(1) = modulo(xnew(1), 1._f64)
          xmid(1) = 1._f64-xbar
       end select
    else if(xnew(1) < -1._f64 .or.  xnew(1) > 2._f64 ) then
       print*, xnew
       SLL_ERROR("particle boundary", "particle out of bound")
    else
       xmid = xold
    end if
    vh = (xnew - xmid)*wi(1)
    call particle_mesh_coupling%add_current( xmid, xnew, vh, j_dofs_local )
    xnew(2:3) = modulo(xnew(2:3), 1._f64)

  end subroutine sll_s_compute_particle_boundary_trafo_current


end module sll_m_time_propagator_pic_vm_3d3v_cl_helper
