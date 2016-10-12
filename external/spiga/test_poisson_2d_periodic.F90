!
! file: ex2.f90
!
!
! usage:
!   > mpirun -np 1 bin/test_poisson_2d_periodic parameters/parameters_2d.nml solvers/solver_mgmres_nullspace.nml
!
! authors:
!   ahmed ratnani  - ratnaniahmed@gmail.com
!

module testcase_poisson_2d
use spi_m_fem_bspline
use spi_m_space_bspline 
use dsc_m_fem_tensor
use dsc_m_space_tensor 
use spi_m_field_2d
use spi_m_poisson_2d
use spi_m_mapping_2d
use spi_m_mapping_gallery
use spi_m_parameters_discretization
implicit none
  type(spi_t_parameters_discretization), target :: parameters
  type(spi_t_fem_bspline), target :: fem_u 
  type(spi_t_fem_bspline), target :: fem_v 
  type(dsc_t_fem_tensor), target :: fem 
  type(spi_t_space_bspline), target :: trial_space_u 
  type(spi_t_space_bspline), target :: trial_space_v 
  type(dsc_t_space_tensor), target :: trial_space
  type(spi_t_poisson_2d), target :: pde 
  type(spi_t_mapping_2d), target :: mapping
  real(spi_rk), dimension(2) :: P_11
  real(spi_rk), dimension(2) :: P_12
  real(spi_rk), dimension(2) :: P_21
  real(spi_rk), dimension(2) :: P_22
  integer :: p_u
  integer :: p_v
  integer :: n_elements_u
  integer :: n_elements_v
  integer :: n_procs_u
  integer :: n_procs_v

contains

  ! ............................................
  subroutine test1 ()
  implicit none
    ! local
    integer :: i_err
    character(len=10) :: label 
    character(len=1024) :: filename
    integer :: i_global
    integer, parameter :: filestream = 111
    character(len=256)            :: filename_solver
    character(len=256)            :: filename_params
    logical                       :: file_exists, equal, empty
    real :: start, finish

    ! ... 
    call spi_initialize(i_err)
    ! ... 

    ! ............................................
    ! Check that input argument was given
    ! ............................................
    if (command_argument_count() /= 2 ) then
      write(*,*) "ERROR: exactly 2 input argument is required"
      stop
    end if
    ! ............................................

    ! ............................................
    ! Read name of reference file from input argument
    ! ............................................
    call get_command_argument( 1, filename_params )
    call get_command_argument( 2, filename_solver )
    ! ............................................

    ! ............................................
    ! Check that file exists    
    ! ............................................
    inquire( file=trim( filename_solver ), exist=file_exists )
    if (.not. file_exists) then
      write(*,*) &
        "ERROR: reference file '"//trim( filename_solver )//"' does not exist"
      stop
    end if
    inquire( file=trim( filename_params ), exist=file_exists )
    if (.not. file_exists) then
      write(*,*) &
        "ERROR: reference file '"//trim( filename_params )//"' does not exist"
      stop
    end if
    ! ............................................

    ! ... 
    call parameters % create(filename=filename_params)
    call parameters % print_info()

    p_u = parameters % p_u
    p_v = parameters % p_v
    n_elements_u = parameters % n_elements_u
    n_elements_v = parameters % n_elements_v
    n_procs_u = parameters % n_procs_u
    n_procs_v = parameters % n_procs_v
    ! ... 

    ! ... 
    P_11 = (/ 0.0_spi_rk, 0.0_spi_rk /)
    P_21 = (/ 1.0_spi_rk, 0.0_spi_rk /)
    P_12 = (/ 0.0_spi_rk, 1.0_spi_rk /)
    P_22 = (/ 1.0_spi_rk, 1.0_spi_rk /)

    call spi_mapping_bilinear(mapping, P_11, P_12, P_21, P_22)

!    call mapping % read_from_file('mapping.nml')
!    call mapping % export('output/test_poisson_2d_mapping.nml')
    ! ... 

    ! ... fem with periodic boundary conditions 
    call fem_u % create( p_u, &
                       & n_elements=n_elements_u, &
                       & type_bc_min=spi_bc_periodic, & 
                       & type_bc_max=spi_bc_periodic)

    call fem_v % create( p_v, &
                       & n_elements=n_elements_v, &
                       & type_bc_min=spi_bc_periodic, & 
                       & type_bc_max=spi_bc_periodic)

    call fem % create(fem_u, fem_v)
    ! ... 

    ! ... 
    call trial_space_u % create(fem_u)
    call trial_space_v % create(fem_v)

    call trial_space % create(fem, trial_space_u, trial_space_v) 
    ! ... 

    ! ... 
    call pde % create(trial_space=trial_space, &
                    & parameters=parameters          , &
                    & filename_solver=filename_solver, &
                    & is_distributed=.true.          )
    ! ... 

    ! ... 
    pde % ptr_field % ptr_function => function_f
    call pde % assemble(map=mapping)
    ! ... 

    ! ... save the matrix in the mm format
    write (label, "(A5,I1)") "proc-", pde % ptr_ddm % i_current_color
    filename = trim(label) // "_poisson_2d_periodic_matrix.mm"
    call pde % matrix % export(filename)

    filename = trim(label) // "_poisson_2d_periodic_rhs.txt"
    call pde % rhs % export (filename, i_format=plf_vector_format_txt)
    ! ...

    ! ...
    call pde % solve()
    write (label, "(A5,I1)") "proc-", pde % ptr_ddm % i_current_color

    filename = trim(label) // "_poisson_2d_periodic_field.nml"
    call pde % ptr_field % export (filename)
    ! ...

    ! ... 
    pde % ptr_field % ptr_function => function_u
    call pde % assemble(map=mapping)
    ! ... 

    ! ... 
    print *, "L2-error: "
    call pde % norms % print_info()
    ! ... 

    ! ...
    call pde % free()
    call trial_space % free()
    call trial_space_v % free()
    call trial_space_u % free()
    call fem % free()
    call fem_v % free()
    call fem_u % free()
    call parameters % free() 
    call mapping % free() 
    ! ...

    ! ...
    call spi_finalize(i_err)
    ! ...

  end subroutine test1
  ! ............................................

  ! ............................................................... 
  subroutine function_f(ptr_field, x, t, f)
  implicit none
    class(dsc_t_field_abstract), pointer :: ptr_field
    real(kind=spi_rk), dimension(:,:), intent(in) :: x 
    real(kind=spi_rk), intent(in) :: t
    real(kind=spi_rk), dimension(:,:), intent(inout) :: f 
    ! local
    real(kind=spi_rk) :: k1 
    real(kind=spi_rk) :: k2
    
    k1 = 2.0 * spi_pi 
    k2 = 2.0 * spi_pi 
    f(1,:) = (k1**2 + k2**2) * cos(k1*x(1,:)) * cos(k2*x(2,:))
  
  end subroutine function_f 
  ! ............................................................... 

  ! ............................................................... 
  subroutine function_u(ptr_field, x, t, f)
  implicit none
    class(dsc_t_field_abstract), pointer :: ptr_field
    real(kind=spi_rk), dimension(:,:), intent(in) :: x
    real(kind=spi_rk), intent(in) :: t
    real(kind=spi_rk), dimension(:,:), intent(inout) :: f 
    ! local
    real(kind=spi_rk) :: k1 
    real(kind=spi_rk) :: k2
    
    k1 = 2.0 * spi_pi 
    k2 = 2.0 * spi_pi 
    f(1,:) = cos(k1*x(1,:)) * cos(k2*x(2,:))
  
  end subroutine function_u 
  ! ............................................................... 

end module testcase_poisson_2d   
  
! ............................................
program main
use testcase_poisson_2d
implicit none

  call test1()

end program main
! ............................................
