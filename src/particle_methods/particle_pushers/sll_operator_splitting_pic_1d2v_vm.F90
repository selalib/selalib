!> @ingroup particle_methods
!> @author Katharina Kormann
!> @brief Particle pusher based on operator splitting for 2x2v Vlasov-Poisson.
!> @details MPI parallelization by domain cloning. Periodic boundaries.
module sll_m_operator_splitting_pic_1d2v_vm
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_errors.h"


 
  use sll_module_pic_base
  use sll_m_kernel_smoother_base
  use sll_collective
  
  use sll_m_pic_maxwell_base

  implicit none

  !> Operator splitting type for Vlasov-Maxwell 1d2v
  type :: sll_operator_splitting_pic_1d2v_vm
     class(sll_pic_maxwell_base), pointer :: maxwell_pic
     !class(sll_maxwell_1d_base), pointer  :: maxwell _solver      !< Maxwell solver
     !class(sll_kernel_smoother_base), pointer :: kernel_smoother  !< Kernel smoother
     class(sll_particle_group_base), pointer  :: particle_group    !< Particle group

     sll_int32 :: n_fmodes


     sll_real64, allocatable :: efield_dofs(:,:) !< DoFs describing the two components of the electric field
     sll_real64, allocatable :: bfield_dofs(:) !< DoFs describing the magnetic field
     sll_real64, allocatable :: efield(:,:)    !< Efield at particle positions
     sll_real64, allocatable :: bfield(:)      !< Bfield at particle positions
     sll_real64, allocatable :: j_dofs(:,:)    !< DoFs for kernel representation of current density. 
     sll_real64, allocatable :: j_dofs_local(:)!< MPI-processor local part of one component of \a j_dofs

   contains
     procedure :: operatorHf => operatorHf_pic_1d2v_vm  !< Operator for H_f part
     procedure :: operatorHE => operatorHE_pic_1d2v_vm  !< Operator for H_E part
     procedure :: operatorHB => operatorHB_pic_1d2v_vm  !< Operator for H_B part

  end type sll_operator_splitting_pic_1d2v_vm

contains

  !---------------------------------------------------------------------------!
  !> Push Hf
  subroutine operatorHf_pic_1d2v_vm(this, dt)
    class(sll_operator_splitting_pic_1d2v_vm), intent(inout) :: this !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step

    !local variables
    sll_int32 :: i_part
    sll_real64 :: x_new(3)


    ! Update current density
    call this%maxwell_pic%compute_shape_factors(this%particle_group)
    this%j_dofs_local = 0.0_f64
    call this%maxwell_pic%accumulate_j_from_klimontovich(this%particle_group, &
         this%j_dofs_local, 1)
    ! MPI to sum up contributions from each processor
    this%j_dofs = 0.0_f64
    call sll_collective_allreduce( sll_world_collective, this%j_dofs_local, &
         this%maxwell_pic%n_dofs, MPI_SUM, this%j_dofs(:,1))
    this%j_dofs_local = 0.0_f64
    call this%maxwell_pic%accumulate_j_from_klimontovich(this%particle_group, &
         this%j_dofs_local, 2)
    ! MPI to sum up contributions from each processor
    call sll_collective_allreduce( sll_world_collective, this%j_dofs_local, &
         this%maxwell_pic%n_dofs, MPI_SUM, this%j_dofs(:,2))
    ! Update electric field (old part)
    this%efield_dofs(2:this%n_fmodes+1,:) = this%efield_dofs(2:this%n_fmodes+1,:) - &
         this%j_dofs(this%n_fmodes+2:2*this%n_fmodes+1,:)
    this%efield_dofs(this%n_fmodes+2:2*this%n_fmodes+1,:) = &
         this%efield_dofs(this%n_fmodes+2:2*this%n_fmodes+1,:) + &
         this%j_dofs(2:this%n_fmodes+1,:)


    ! X_new = X_old + dt * V
    do i_part=1,this%particle_group%n_particles
       x_new = this%particle_group%get_x(i_part) + dt * this%particle_group%get_v(i_part) 
       call this%particle_group%set_x(i_part, x_new)
    end do

    ! Update current density
    call this%maxwell_pic%compute_shape_factors(this%particle_group)
    this%j_dofs_local = 0.0_f64
    call this%maxwell_pic%accumulate_j_from_klimontovich(this%particle_group, &
         this%j_dofs_local, 1)
    ! MPI to sum up contributions from each processor
    this%j_dofs = 0.0_f64
    call sll_collective_allreduce( sll_world_collective, this%j_dofs_local, &
         this%maxwell_pic%n_dofs, MPI_SUM, this%j_dofs(:,1))
    this%j_dofs_local = 0.0_f64
    call this%maxwell_pic%accumulate_j_from_klimontovich(this%particle_group, &
         this%j_dofs_local, 2)
    ! MPI to sum up contributions from each processor
    call sll_collective_allreduce( sll_world_collective, this%j_dofs_local, &
         this%maxwell_pic%n_dofs, MPI_SUM, this%j_dofs(:,2))

    ! Update electric field (new part)
    this%efield_dofs(2:this%n_fmodes+1,:) = this%efield_dofs(2:this%n_fmodes+1,:) + &
         this%j_dofs(this%n_fmodes+2:2*this%n_fmodes+1,:)
    this%efield_dofs(this%n_fmodes+2:2*this%n_fmodes+1,:) = &
         this%efield_dofs(this%n_fmodes+2:2*this%n_fmodes+1,:) - &
         this%j_dofs(2:this%n_fmodes+1,:)

  end subroutine operatorHf_pic_1d2v_vm
  
  !---------------------------------------------------------------------------!
  ! Push H_E
  subroutine operatorHE_pic_1d2v_vm(this, dt)
    class(sll_operator_splitting_pic_1d2v_vm), intent(inout) :: this !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step

    !local variables
    sll_int32 :: i_part
    sll_real64 :: v_new(3)

    ! Evaluate efield at particle positions
    call this%maxwell_pic%evaluate_kernel_function(this%particle_group, &
         this%efield_dofs(:,1), this%efield(:,1))
    call this%maxwell_pic%evaluate_kernel_function(this%particle_group, &
         this%efield_dofs(:,2), this%efield(:,1))

    ! V_new = V_old + dt * E
    do i_part=1,this%particle_group%n_particles
       v_new = this%particle_group%get_v(i_part)
       v_new(1:2) = v_new(1:2) + dt * this%efield(i_part,:) 
       call this%particle_group%set_v(i_part, v_new)
    end do
    
    ! Update bfield
    call this%maxwell_pic%compute_B_from_E( &
         dt, this%efield_dofs(:,2), this%bfield_dofs)
    

  end subroutine operatorHE_pic_1d2v_vm
  

  !---------------------------------------------------------------------------!
  ! Push H_B
  subroutine operatorHB_pic_1d2v_vm(this, dt)
    class(sll_operator_splitting_pic_1d2v_vm), intent(inout) :: this !< time splitting object 
    sll_real64, intent(in) :: dt   !< time step

    !local variables
    sll_int32 :: i_part
    sll_real64 :: v_old(3), v_new(3)

    ! Evaluate bfield at particle positions
    call this%maxwell_pic%evaluate_kernel_function(this%particle_group, &
         this%bfield_dofs, this%bfield)

    ! V_new = V_old + dt * E
    do i_part=1,this%particle_group%n_particles
       v_old = this%particle_group%get_v(i_part)
       v_new(1) = v_old(1) * cos(this%bfield(i_part)*dt) + &
            v_old(2) * sin(this%bfield(i_part)*dt)
       v_new(2) = v_old(2) * cos(this%bfield(i_part)*dt) - &
            v_old(1) * sin(this%bfield(i_part)*dt) 
       call this%particle_group%set_v(i_part, v_new)
    end do
    
    ! Update efield2
    call this%maxwell_pic%compute_E_from_B(&
         dt, this%bfield_dofs, this%efield_dofs(:,2))
    

  end subroutine operatorHB_pic_1d2v_vm



  !---------------------------------------------------------------------------!
  !> Constructor.
  function sll_new_splitting_pic_1d2v_vm(maxwell_pic, particle_group) result(this)
    class(sll_operator_splitting_pic_1d2v_vm), pointer :: this !< time splitting object 
    class(sll_pic_maxwell_base),pointer, intent(in) :: maxwell_pic !< Maxwell solver
    class(sll_particle_group_base),pointer, intent(in) :: particle_group !< Particle group

    !local variables
    sll_int32 :: ierr

    SLL_ALLOCATE(this, ierr)

    this%maxwell_pic => maxwell_pic
    this%particle_group => particle_group

    SLL_ALLOCATE(this%efield_dofs(this%maxwell_pic%n_dofs,2), ierr)
    SLL_ALLOCATE(this%bfield_dofs(this%maxwell_pic%n_dofs), ierr)
    SLL_ALLOCATE(this%efield(this%particle_group%n_particles,2), ierr)
    SLL_ALLOCATE(this%bfield(this%particle_group%n_particles), ierr)
    SLL_ALLOCATE(this%j_dofs(this%maxwell_pic%n_dofs,2), ierr)
    SLL_ALLOCATE(this%j_dofs_local(this%maxwell_pic%n_dofs), ierr)

    ! Initialize the fields
    ! TODO: SHOULD NOT BE HARDCODED HERE
    this%bfield = 0.0_f64
    this%bfield(2) = 0.0001_f64
    this%efield_dofs(:,2) = 0.0_f64
    ! Compute E1 from Poisson
    ! Accumulate rho to j_dofs
    call this%maxwell_pic%compute_shape_factors(this%particle_group)
    this%j_dofs_local = 0.0_f64
    call this%maxwell_pic%accumulate_rho_from_klimontovich(this%particle_group, &
         this%j_dofs_local)
    ! MPI to sum up contributions from each processor
    this%j_dofs = 0.0_f64
    call sll_collective_allreduce( sll_world_collective, this%j_dofs_local, &
         this%maxwell_pic%n_dofs, MPI_SUM, this%j_dofs(:,1))
    !this%j_dofs(:,1) = this%j_dofs_local
   !print*, 'j0', this%j_dofs(:,1)
    ! Solve Poisson for rho
    call this%maxwell_pic%compute_E_from_rho(this%efield_dofs(:,1), this%j_dofs(:,1))
    !print*, 'E10', this%efield_dofs(:,1)

    this%n_fmodes = (this%maxwell_pic%n_dofs-1)/2
    
  end function sll_new_splitting_pic_1d2v_vm




end module sll_m_operator_splitting_pic_1d2v_vm
