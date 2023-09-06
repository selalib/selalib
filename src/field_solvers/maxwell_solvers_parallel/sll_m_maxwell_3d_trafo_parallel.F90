!> @ingroup maxwell_solvers
!> @brief
!> Module interface to solve Maxwell's equations with coordinate transformation in 3D
!> The linear systems are solved using PLAF
!> @details
!> 
!> @author
!> Benedikt Perse

module sll_m_maxwell_3d_trafo_parallel
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
       sll_o_collective_allreduce, &
       sll_s_collective_reduce_real64, &
       sll_f_get_collective_rank, &
       sll_f_get_collective_size, &
       sll_v_world_collective

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  use sll_m_matrix_csr, only: &
       sll_t_matrix_csr

  use sll_m_maxwell_3d_trafo, only : &
       sll_t_maxwell_3d_trafo

  use mpi, only: &
       mpi_sum

  use sll_m_profile_functions, only: &
       sll_t_profile_functions

  use sll_m_spline_fem_utilities_3d, only: &
       sll_s_spline_fem_mass3d, &
       sll_s_spline_fem_mixedmass3d

  use sll_m_spline_fem_utilities_3d_helper, only: &
       sll_s_spline_fem_sparsity_mass3d, &
       sll_s_spline_fem_sparsity_mixedmass3d


  implicit none

  public :: &
       sll_t_maxwell_3d_trafo_parallel

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, extends(sll_t_maxwell_3d_trafo) :: sll_t_maxwell_3d_trafo_parallel


   contains
     procedure :: &
          init => init_3d_trafo_parallel

  end type sll_t_maxwell_3d_trafo_parallel

contains

  subroutine init_3d_trafo_parallel( self, domain, n_dofs, s_deg_0, map, mass_tolerance, poisson_tolerance, solver_tolerance, adiabatic_electrons, profile )
    class(sll_t_maxwell_3d_trafo_parallel),   intent( inout ) :: self      !< Maxwell solver object
    sll_real64, intent(in) :: domain(3,2) 
    sll_int32,                       intent( in    ) :: n_dofs(3) !< number of degrees of freedom (here number of cells and grid points)
    sll_int32,                       intent( in    ) :: s_deg_0(3)   !< highest spline degree
    type(sll_t_mapping_3d), target, intent( inout    ) :: map       !< coordinate transformation
    sll_real64, intent(in), optional :: mass_tolerance !< tolerance for mass solver
    sll_real64, intent(in), optional :: poisson_tolerance !< tolerance for Poisson solver
    sll_real64, intent(in), optional :: solver_tolerance !< tolerance for Schur complement solver
    logical, intent(in), optional :: adiabatic_electrons  !< flag if adiabatic electrons are used
    type(sll_t_profile_functions), intent(in), optional :: profile !< temperature and density profiles
    ! local variables
    sll_int32  :: mpi_rank
    sll_int32  :: mpi_total
    sll_int32  :: begin(3)
    sll_int32  :: limit(3)
    sll_int32  :: j
    sll_int32  :: deg1(3), deg2(3)
    sll_real64, allocatable :: nullspace(:,:)
    type(sll_t_matrix_csr) :: mass0_local
    type(sll_t_matrix_csr) :: mass1_local(3,3)
    type(sll_t_matrix_csr) :: mass2_local(3,3) 

    self%Lx = map%Lx

    if (present( mass_tolerance) ) then
       self%mass_solver_tolerance = mass_tolerance
    else
       self%mass_solver_tolerance = 1d-12
    end if
    if (present( poisson_tolerance) ) then
       self%poisson_solver_tolerance = poisson_tolerance
    else
       self%poisson_solver_tolerance = 1d-12
    end if
    if (present( solver_tolerance) ) then
       self%solver_tolerance = solver_tolerance
    else
       self%solver_tolerance = 1d-12
    end if

    if( present( adiabatic_electrons ) ) then
       self%adiabatic_electrons = adiabatic_electrons
    end if

    if( present( profile ) ) then
       self%profile = profile
    end if

    self%map => map
    self%n_cells = n_dofs
    self%n_dofs = n_dofs
    self%n_total = product(n_dofs)
    self%n_total0 = self%n_total
    self%n_total1 = self%n_total

    self%delta_x = 1._f64 / real(n_dofs,f64)
    self%s_deg_0 = s_deg_0
    self%s_deg_1 = s_deg_0 - 1
    self%volume = product(self%delta_x)

    ! Allocate scratch data
    allocate( self%work0(1:self%n_total) )
    allocate( self%work(1:self%n_total*3) )
    allocate( self%work2(1:self%n_total*3) )

    mpi_total = sll_f_get_collective_size(sll_v_world_collective)
    mpi_rank = sll_f_get_collective_rank(sll_v_world_collective)


    call sll_s_compute_mpi_decomposition( mpi_rank, mpi_total, n_dofs, begin, limit )
    ! Assemble the sparse mass matrices
    !0-form
    deg1=self%s_deg_0
    if(self%adiabatic_electrons) then
       call sll_s_spline_fem_mass3d( n_dofs, deg1, -1, mass0_local, profile_m0, begin, limit )
    else
       call sll_s_spline_fem_mass3d( n_dofs, deg1, 0, mass0_local, profile_0, begin, limit )
    end if
    do j=1, 3
       deg1=self%s_deg_0
       deg1(j)=self%s_deg_1(j)
       call sll_s_spline_fem_mass3d( n_dofs, deg1, j, mass1_local(j,j), profile_1, begin, limit )
       deg1 =  self%s_deg_1
       deg1(j) = self%s_deg_0(j)
       call sll_s_spline_fem_mass3d( n_dofs, deg1, j, mass2_local(j,j), profile_2, begin, limit )
    end do
    if(self%map%flag2d)then
       deg1=self%s_deg_0
       deg1(1)=self%s_deg_1(1)
       deg2=self%s_deg_0
       deg2(2) = self%s_deg_1(2)
       call sll_s_spline_fem_mixedmass3d( n_dofs, deg1, deg2, [1,2], mass1_local(1,2), profile_m1, begin, limit )
       call sll_s_spline_fem_mixedmass3d( n_dofs, deg2, deg1, [2,1], mass1_local(2,1), profile_m1, begin, limit )

       deg1 =  self%s_deg_1
       deg1(1) = self%s_deg_0(1)
       deg2 = self%s_deg_1
       deg2(2) = self%s_deg_0(2)
       call sll_s_spline_fem_mixedmass3d( n_dofs, deg1, deg2, [1,2], mass2_local(1,2), profile_m2, begin, limit )
       call sll_s_spline_fem_mixedmass3d( n_dofs, deg2, deg1, [2,1], mass2_local(2,1), profile_m2, begin, limit )
       if(self%map%flag3d)then
          do j = 2, 3
             deg1=self%s_deg_0
             deg1(j)=self%s_deg_1(j)
             deg2=self%s_deg_0
             deg2(modulo(j,3)+1) = self%s_deg_1(modulo(j,3)+1)
             call sll_s_spline_fem_mixedmass3d( n_dofs, deg1, deg2, [j,modulo(j,3)+1], mass1_local(j,modulo(j,3)+1), profile_m1, begin, limit )
             call sll_s_spline_fem_mixedmass3d( n_dofs, deg2, deg1, [modulo(j,3)+1,j], mass1_local(modulo(j,3)+1,j), profile_m1, begin, limit )

             deg1 =  self%s_deg_1
             deg1(j) = self%s_deg_0(j)
             deg2 = self%s_deg_1
             deg2(modulo(j,3)+1) = self%s_deg_0(modulo(j,3)+1)
             call sll_s_spline_fem_mixedmass3d( n_dofs, deg1, deg2, [j,modulo(j,3)+1], mass2_local(j,modulo(j,3)+1), profile_m2, begin, limit )
             call sll_s_spline_fem_mixedmass3d( n_dofs, deg2, deg1, [modulo(j,3)+1,j], mass2_local(modulo(j,3)+1,j), profile_m2, begin, limit )
          end do
       end if
    end if


    !initialize global matrices
    deg1=self%s_deg_0
    call sll_s_spline_fem_sparsity_mass3d( deg1, n_dofs, self%mass0 )
    do j=1, 3
       deg1=self%s_deg_0
       deg1(j)=self%s_deg_1(j)
       call sll_s_spline_fem_sparsity_mass3d( deg1, n_dofs, self%mass1(j,j) )
       deg1 =  self%s_deg_1
       deg1(j) = self%s_deg_0(j)
       call sll_s_spline_fem_sparsity_mass3d( deg1, n_dofs, self%mass2(j,j) )
    end do
    if(self%map%flag2d)then
       deg1=self%s_deg_0
       deg1(1)=self%s_deg_1(1)
       deg2=self%s_deg_0
       deg2(2) = self%s_deg_1(2)
       call sll_s_spline_fem_sparsity_mixedmass3d( deg1, deg2, n_dofs, self%mass1(1,2) )
       call sll_s_spline_fem_sparsity_mixedmass3d( deg2, deg1, n_dofs, self%mass1(2,1) )
       deg1 =  self%s_deg_1
       deg1(1) = self%s_deg_0(1)
       deg2 = self%s_deg_1
       deg2(2) = self%s_deg_0(2)
       call sll_s_spline_fem_sparsity_mixedmass3d( deg1, deg2, n_dofs, self%mass2(1,2)  ) 
       call sll_s_spline_fem_sparsity_mixedmass3d( deg2, deg1, n_dofs, self%mass2(2,1) ) 
       if(self%map%flag3d)then
          do j = 2, 3
             deg1=self%s_deg_0
             deg1(j)=self%s_deg_1(j)
             deg2=self%s_deg_0
             deg2(modulo(j,3)+1) = self%s_deg_1(modulo(j,3)+1)
             call sll_s_spline_fem_sparsity_mixedmass3d( deg1, deg2, n_dofs, self%mass1(j,modulo(j,3)+1) )
             call sll_s_spline_fem_sparsity_mixedmass3d( deg2, deg1, n_dofs, self%mass1(modulo(j,3)+1,j) )
             deg1 =  self%s_deg_1
             deg1(j) = self%s_deg_0(j)
             deg2 = self%s_deg_1
             deg2(modulo(j,3)+1) = self%s_deg_0(modulo(j,3)+1)
             call sll_s_spline_fem_sparsity_mixedmass3d( deg1, deg2, n_dofs, self%mass2(j,modulo(j,3)+1)  ) 
             call sll_s_spline_fem_sparsity_mixedmass3d( deg2, deg1, n_dofs, self%mass2(modulo(j,3)+1,j) )
          end do
       end if
    end if

    ! MPI to sum up contributions from each processor
    call sll_o_collective_allreduce( sll_v_world_collective, mass0_local%arr_a(1,:), &
         self%mass0%n_nnz, MPI_SUM, self%mass0%arr_a(1,:))
    do j=1, 3
       call sll_o_collective_allreduce( sll_v_world_collective, mass1_local(j,j)%arr_a(1,:), &
            self%mass1(j,j)%n_nnz, MPI_SUM, self%mass1(j,j)%arr_a(1,:))
       call sll_o_collective_allreduce( sll_v_world_collective, mass2_local(j,j)%arr_a(1,:), &
            self%mass2(j,j)%n_nnz, MPI_SUM, self%mass2(j,j)%arr_a(1,:))
    end do
    if(self%map%flag2d)then
       call sll_o_collective_allreduce( sll_v_world_collective, mass1_local(1,2)%arr_a(1,:), &
            self%mass1(1,2)%n_nnz, MPI_SUM, self%mass1(1,2)%arr_a(1,:))
       call sll_o_collective_allreduce( sll_v_world_collective, mass1_local(2,1)%arr_a(1,:), &
            self%mass1(2,1)%n_nnz, MPI_SUM, self%mass1(2,1)%arr_a(1,:))
       call sll_o_collective_allreduce( sll_v_world_collective, mass2_local(1,2)%arr_a(1,:), &
            self%mass2(1,2)%n_nnz, MPI_SUM, self%mass2(1,2)%arr_a(1,:))
       call sll_o_collective_allreduce( sll_v_world_collective, mass2_local(2,1)%arr_a(1,:), &
            self%mass2(2,1)%n_nnz, MPI_SUM, self%mass2(2,1)%arr_a(1,:))
       if(self%map%flag3d)then
          do j = 2, 3
             call sll_o_collective_allreduce( sll_v_world_collective, mass1_local(j,modulo(j,3)+1)%arr_a(1,:), &
                  self%mass1(j,modulo(j,3)+1)%n_nnz, MPI_SUM, self%mass1(j,modulo(j,3)+1)%arr_a(1,:))
             call sll_o_collective_allreduce( sll_v_world_collective, mass1_local(modulo(j,3)+1,j)%arr_a(1,:), &
                  self%mass1(modulo(j,3)+1,j)%n_nnz, MPI_SUM, self%mass1(modulo(j,3)+1,j)%arr_a(1,:))
             call sll_o_collective_allreduce( sll_v_world_collective, mass2_local(j,modulo(j,3)+1)%arr_a(1,:), &
                  self%mass2(j,modulo(j,3)+1)%n_nnz, MPI_SUM, self%mass2(j,modulo(j,3)+1)%arr_a(1,:))
             call sll_o_collective_allreduce( sll_v_world_collective, mass2_local(modulo(j,3)+1,j)%arr_a(1,:), &
                  self%mass2(modulo(j,3)+1,j)%n_nnz, MPI_SUM, self%mass2(modulo(j,3)+1,j)%arr_a(1,:))
          end do
       end if
    end if

    call self%mass1_operator%create( 3, 3 )
    call self%mass2_operator%create( 3, 3 )

    do j=1,3
       call self%mass1_operator%set( j, j, self%mass1(j,j) )
       call self%mass2_operator%set( j, j, self%mass2(j,j) )
    end do
    if(self%map%flag2d)then
       call self%mass1_operator%set( 1, 2, self%mass1(1, 2) )
       call self%mass1_operator%set( 2, 1, self%mass1(2, 1) )
       call self%mass2_operator%set( 1, 2, self%mass2(1, 2) )
       call self%mass2_operator%set( 2, 1, self%mass2(2, 1) )
       if(self%map%flag3d)then
          do j = 2, 3
             call self%mass1_operator%set( j, modulo(j,3)+1, self%mass1(j, modulo(j,3)+1) )
             call self%mass1_operator%set( modulo(j,3)+1, j, self%mass1(modulo(j,3)+1, j) )
             call self%mass2_operator%set( j, modulo(j,3)+1, self%mass2(j, modulo(j,3)+1) )
             call self%mass2_operator%set( modulo(j,3)+1, j, self%mass2(modulo(j,3)+1, j) )
          end do
       end if
    end if

    call self%preconditioner_fft%init( self%Lx, n_dofs, s_deg_0 )
    call self%mass0_solver%create( self%mass0 )
    call self%mass1_solver%create( self%mass1_operator, self%preconditioner_fft%inverse_mass1_3d)
    call self%mass2_solver%create( self%mass2_operator, self%preconditioner_fft%inverse_mass2_3d)
    self%mass0_solver%atol = self%mass_solver_tolerance
    self%mass1_solver%atol = self%mass_solver_tolerance/maxval(self%Lx)
    self%mass2_solver%atol = self%mass_solver_tolerance/maxval(self%Lx)**2
    !self%mass1_solver%verbose = .true.
    !self%mass2_solver%verbose = .true.

    call self%poisson_matrix%create( self%mass1_operator, self%n_dofs, self%delta_x )
    ! Penalized Poisson operator
    allocate(nullspace(1,1:self%n_total))
    nullspace(1,:) = 1.0_f64
    call self%poisson_operator%create( self%poisson_matrix, nullspace, 1 )
    ! Poisson solver
    call self%poisson_solver%create( self%poisson_operator )
    self%poisson_solver%null_space = .true.
    self%poisson_solver%atol = self%poisson_solver_tolerance
    !self%poisson_solver%verbose = .true.
    self%poisson_solver%n_maxiter=40000


    ! Only for Schur complement eb solver
    call self%linear_op_schur_eb%create( self%mass1_operator, self%mass2_operator, self%n_total, self%n_dofs, self%delta_x   )
    call self%linear_solver_schur_eb%create( self%linear_op_schur_eb, self%preconditioner_fft%inverse_mass1_3d)
    self%linear_solver_schur_eb%atol = self%solver_tolerance/maxval(self%Lx)
    !self%linear_solver_schur_eb%verbose = .true.
    !self%linear_solver_schur_eb%n_maxiter = 2000

  contains
    function profile_m0( x, component)
      sll_real64 :: profile_m0
      sll_real64, intent(in) :: x(3)
      sll_int32, optional, intent(in)  :: component(:)

      profile_m0 = self%map%jacobian( x ) * self%profile%rho_0( x(1) )/ self%profile%T_e( x(1) )

    end function profile_m0

    function profile_0(x, component)
      sll_real64 :: profile_0
      sll_real64, intent(in) :: x(3)
      sll_int32, optional, intent(in)  :: component(:)

      profile_0 = self%map%jacobian( x ) 

    end function profile_0

    function profile_1(x, component)
      sll_real64 :: profile_1
      sll_real64, intent(in) :: x(3)
      sll_int32, optional, intent(in)  :: component(:)

      profile_1 = self%map%metric_inverse_single_jacobian( x, component(1), component(1) )

    end function profile_1

    function profile_2(x, component)
      sll_real64 :: profile_2
      sll_real64, intent(in) :: x(3)
      sll_int32, optional, intent(in)  :: component(:)

      profile_2 = self%map%metric_single_jacobian( x, component(1), component(1) )

    end function profile_2

    function profile_m1(x, component)
      sll_real64 :: profile_m1
      sll_real64, intent(in) :: x(3)
      sll_int32, optional, intent(in)  :: component(:)

      profile_m1 = self%map%metric_inverse_single_jacobian( x, component(1), component(2) )

    end function profile_m1

    function profile_m2(x, component)
      sll_real64 :: profile_m2
      sll_real64, intent(in) :: x(3)
      sll_int32, optional, intent(in)  :: component(:)

      profile_m2 = self%map%metric_single_jacobian( x, component(1), component(2) )

    end function profile_m2

  end subroutine init_3d_trafo_parallel



  subroutine sll_s_compute_mpi_decomposition( mpi_rank, mpi_total, n_cells, begin, limit )
    sll_int32, intent( in    ) :: mpi_rank
    sll_int32, intent( in    ) :: mpi_total
    sll_int32, intent( in    ) :: n_cells(3)
    sll_int32, intent(   out ) :: begin(3)
    sll_int32, intent(   out ) :: limit(3)

    if( mpi_total <= n_cells(3) )then
       begin(3)=1+mpi_rank*n_cells(3)/mpi_total
       limit(3)=(mpi_rank+1)*n_cells(3)/mpi_total
       if( mpi_rank == mpi_total-1 ) then
          limit(3)=n_cells(3)
       end if
       begin(2)=1
       limit(2)=n_cells(2)
    else if( mpi_total > n_cells(3) .and. mpi_total < n_cells(2)*n_cells(3) ) then
       if(modulo(real(n_cells(3)*n_cells(2),f64)/real(mpi_total,f64),1._f64)==0 )then
          begin(3)=1+mpi_rank*n_cells(3)/mpi_total
          limit(3)=begin(3)
          begin(2)=1+modulo((mpi_rank*n_cells(2)*n_cells(3))/mpi_total, n_cells(2))
          limit(2)=modulo(((mpi_rank+1)*n_cells(2)*n_cells(3))/mpi_total-1,n_cells(2))+1
       else
          SLL_WARNING('maxwell_3d_trafo_parallel:init', 'amount of mpi processes not optimal')
          if( mpi_rank < n_cells(3) ) then
             begin(3)=1+mpi_rank
             limit(3)=begin(3)
             begin(2)=1
             limit(2)=n_cells(2)
          else
             begin(3)=1
             limit(3)=0
             begin(2)=1
             limit(2)=0
          end if
       end if
    else
       if( mpi_rank < n_cells(2)*n_cells(3) )then
          begin(3)=1+mpi_rank/n_cells(2)
          limit(3)=begin(3)
          begin(2)=1+modulo(mpi_rank,n_cells(2))
          limit(2)=begin(2)
       else
          begin(3)=1
          limit(3)=0
          begin(2)=1
          limit(2)=0
       end if
    end if
    begin(1)=1
    limit(1)=n_cells(1)

  end subroutine sll_s_compute_mpi_decomposition


end module sll_m_maxwell_3d_trafo_parallel
