! TODO: Use input from file to initialize and compare

program test_hamiltonian_splitting_cef_pic_vm_1d2v
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_m_particle_group_base
  use sll_m_particle_initializer
  use sll_m_particle_group_1d2v
  use sll_m_kernel_smoother_base
  use sll_m_kernel_smoother_spline_1d
  use sll_m_hamiltonian_splitting_cef_pic_vm_1d2v
  use sll_m_maxwell_1d_base
  use sll_m_maxwell_1d_fem
  use sll_m_constants, only : &
       sll_pi
  use sll_m_collective, only : &
       sll_boot_collective, sll_halt_collective

  ! Abstract particle group
  class(sll_particle_group_base), pointer :: particle_group

  ! Arrays for the fields
  sll_real64, pointer :: efield(:,:), efield_ref(:,:)
  sll_real64, pointer :: bfield(:), bfield_ref(:)

  ! Abstract kernel smoothers
  class(sll_c_kernel_smoother), pointer :: kernel_smoother_0     
  class(sll_c_kernel_smoother), pointer :: kernel_smoother_1
  
  ! Maxwell solver 
  ! Abstract 
  class(sll_maxwell_1d_base), pointer :: maxwell_solver
  
  ! Specific operator splitting
  class(sll_t_hamiltonian_splitting_cef_pic_vm_1d2v), pointer :: propagator
  
  ! Parameters
  sll_int32  :: n_particles
  sll_real64 :: eta_min, eta_max, domain(3)
  sll_int32  :: num_cells
  sll_real64 :: delta_t
  sll_int32  :: degree_smoother
  sll_int64  :: rnd_seed

  ! Helper 
  sll_int32  :: i_part
  sll_real64 :: xi(3)
  logical    :: passed
  sll_real64 :: error

  ! Reference
  sll_real64, allocatable :: particle_info_ref(:,:)
   
  call sll_boot_collective()

  ! Set parameters
  n_particles = 2!10
  eta_min = 0.0_f64
  eta_max = 4.0_f64*sll_pi
  num_cells = 10
  delta_t = 0.1_f64
  degree_smoother = 3
  passed = .TRUE.
  rnd_seed = 10

  domain = [eta_min, eta_max, eta_max - eta_min]
  
  ! Initialize
  particle_group => sll_new_particle_group_1d2v(n_particles, &
       n_particles ,1.0_f64, 1.0_f64, 1)

  ! Initial particle information   
  ! Data produce with following call
  !call sll_particle_initialize_sobol_landau_1d2v(particle_group, &
  !     [0.1_f64, 0.5_f64], domain(1),domain(3), &
  !     [1.0_f64, 1.0_f64] , rnd_seed)
  

  SLL_ALLOCATE(particle_info_ref(n_particles,4), i_part)
  particle_info_ref = 0.0_f64
  particle_info_ref = reshape([11.780972450961723D0,        5.4977871437821380D0,       -1.5341205443525459D0,       0.15731068461017067D0,       0.15731068461017067D0,       -1.5341205443525459D0,        6.8636759376074723D0,        5.7026946767517002D0   ], [n_particles, 4])

  ! Initialize particles from particle_info_ref
  xi = 0.0_f64
  do i_part =1, n_particles
     xi(1) = particle_info_ref(i_part, 1) 
     call particle_group%set_x(i_part, xi)
     xi(1:2) = particle_info_ref(i_part, 2:3)
     call particle_group%set_v(i_part, xi)
     xi(1) = particle_info_ref(i_part, 4)
     call particle_group%set_weights(i_part, xi(1))
  end do

  ! Initialize kernel smoother    
  kernel_smoother_1 => sll_new_smoother_spline_1d(&
       domain(1:2), [num_cells], &
       n_particles, degree_smoother-1, SLL_GALERKIN) 
  kernel_smoother_0 => &
       sll_new_smoother_spline_1d(domain(1:2), [num_cells], &
       n_particles, degree_smoother, SLL_GALERKIN) 
  
  ! Initialize Maxwell solver
  maxwell_solver => sll_new_maxwell_1d_fem([eta_min, eta_max], num_cells, &
       degree_smoother)

  SLL_ALLOCATE(efield(kernel_smoother_0%n_dofs,2),ierr)
  SLL_ALLOCATE(bfield(kernel_smoother_0%n_dofs),ierr)
  SLL_ALLOCATE(efield_ref(kernel_smoother_0%n_dofs,2),ierr)
  SLL_ALLOCATE(bfield_ref(kernel_smoother_0%n_dofs),ierr)

  efield = 1.0_f64
  bfield = 1.0_f64

  propagator => sll_new_hamiltonian_splitting_cef_pic_vm_1d2v(maxwell_solver, &
            kernel_smoother_0, kernel_smoother_1, particle_group, &
            efield, bfield, &
            eta_min, eta_max-eta_min)

  call propagator%operatorHp1(delta_t)


  ! Compare to reference
  ! Particle information after operatorV application 
  particle_info_ref = reshape([  11.627560396526469D0,        5.5135182122431550D0,       -1.5341205443525459D0,       0.15731068461017067D0,       0.15731068461017067D0,       -1.5341205443525459D0,        6.8636759376074723D0,        5.7026946767517002D0     ], [n_particles, 4])
  ! Compare computed values to reference values
  do i_part=1,n_particles
     xi = particle_group%get_x(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,1))> 1D-14) then
        passed = .FALSE.
     end if
     xi = particle_group%get_v(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,2))> 1D-14) then
        passed = .FALSE.
     elseif (abs(xi(2)-particle_info_ref(i_part,3))> 1D-14) then
        passed = .FALSE.
     end if
     xi(1:1) = particle_group%get_charge(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,4))> 1D-14) then
        passed = .FALSE.
     end if
  end do

  

  call propagator%operatorHE(delta_t)
  ! Compare to reference
  ! Particle information after operatorV application 
  particle_info_ref = reshape([   11.627560396526469D0,        5.5135182122431550D0,       -1.4007142828624286D0,       0.25441596552957291D0,       0.25438068275455678D0,       -1.4104377250793367D0,        6.8636759376074723D0,        5.7026946767517002D0    ], [n_particles, 4])
  ! Compare computed values to reference values
    do i_part=1,n_particles
     xi = particle_group%get_x(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,1))> 1D-14) then
        passed = .FALSE.
     end if
     xi = particle_group%get_v(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,2))> 1D-14) then
        passed = .FALSE.
     elseif (abs(xi(2)-particle_info_ref(i_part,3))> 1D-14) then
        passed = .FALSE.
     end if
     xi(1:1) = particle_group%get_charge(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,4))> 1D-14) then
        passed = .FALSE.
     end if
  end do


  call propagator%operatorHB(delta_t)

  ! Compare to reference
  ! Particle information after operatorV application 
  particle_info_ref = reshape([  11.627560396526469D0,        5.5135182122431550D0,       -1.3683192894303602D0,       0.11227480640273069D0,       0.39295337655538259D0,       -1.4287954464130088D0,        6.8636759376074723D0,        5.7026946767517002D0    ], [n_particles, 4])
  ! Compare computed values to reference values
    do i_part=1,n_particles
     xi = particle_group%get_x(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,1))> 1D-14) then
        passed = .FALSE.
     end if
     xi = particle_group%get_v(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,2))> 1D-14) then
        passed = .FALSE.
     elseif (abs(xi(2)-particle_info_ref(i_part,3))> 1D-14) then
        passed = .FALSE.
     end if
     xi(1:1) = particle_group%get_charge(i_part)
     if (abs(xi(1)-particle_info_ref(i_part,4))> 1D-14) then
        passed = .FALSE.
     end if
  end do

  bfield_ref = [ 0.99819008852299340D0,       0.99241234971657222D0,       0.98838969093169826D0,        1.0023321692551928D0,        1.0118169182044585D0,        1.0057521726101171D0,        1.0016358135685306D0,        1.0014153888891870D0,       0.99955401601159288D0,       0.99850139228965729D0 ] 

  efield_ref = reshape([1.0139333322572233D0,       0.99694980433046210D0,       0.98105828960960473D0,       0.96702117187072789D0,       0.98564985298325469D0,        1.0000855014667056D0,        1.0477798381501124D0,        1.2387410309599800D0,        1.3810186780785210D0,        1.1543013648190512D0,        1.0152805152929967D0,        1.1106903616633781D0,        1.2546884206799935D0,        1.2259034757779683D0,        1.0789193162085262D0,        1.0062750756095262D0,       0.98537165227992685D0,       0.96780452240085313D0,       0.97331093901524457D0,       0.99202671628832617D0   ],[num_cells,2])

  error = maxval(bfield-bfield_ref)

  if (error> 1E-14) then
     passed = .FALSE.
  end if

  error = maxval(efield-efield_ref)

  if (error> 1E-14) then
     passed = .FALSE.
  end if



  if (passed .EQV. .TRUE.) then
     print*, 'PASSED'
  else
     print*, 'FAILED'
     stop
  end if
  
  call sll_halt_collective()
end program test_hamiltonian_splitting_cef_pic_vm_1d2v
