!
! file: ex1.f90
!
!
! usage:
!   > ./ex1
!
! authors:
!   ahmed ratnani  - ratnaniahmed@gmail.com
!

! ............................................
program main
use spl_m_mapping_2d
use spl_m_mapping_gallery
use jrk_m_global
use jrk_m_assembler_1d1d
use jrk_m_gallery_1d1d
use dsc_m_space_abstract
use dsc_m_field_abstract
use plf_m_matrix_csr
use plf_m_vector
use plf_m_linear_solver_driver
use jrk_m_context
use jrk_m_context_utilities
use plf_m_ddm_parameters

implicit none
  integer :: i_err
  class(dsc_t_space_abstract), allocatable :: trial_space
  class(dsc_t_space_abstract), allocatable :: test_space
  type(jrk_t_assembler_1d1d), target :: assembler
  type(dsc_t_field)          , target      :: phi
  type(jrk_t_context)        , target      :: context
  type(jrk_t_context)        , target      :: test_context
  type(plf_t_vector), target :: rhs
  type(plf_t_matrix_csr), target :: matrix
  type(plf_t_vector), target :: norms
  type(spl_t_mapping_2d), target :: mapping
  type(plf_t_ddm_parameters), target :: ddm_params
  type(plf_t_linear_solver_driver)  , target :: linear_solver

  real(jrk_rk), dimension(2) :: P_11
  real(jrk_rk), dimension(2) :: P_12
  real(jrk_rk), dimension(2) :: P_21
  real(jrk_rk), dimension(2) :: P_22
  real(plf_rk), dimension(:), allocatable :: x
  real(plf_rk), dimension(:), allocatable :: y
  integer, parameter :: n_norms = 3
  integer, parameter :: filestream = 111
  character(len=256)            :: filename_solver
  character(len=256)            :: filename_params
  logical                       :: file_exists, equal, empty

  ! ...
  call jrk_initialize(i_err)
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
  P_11 = (/ 0.0_jrk_rk, 0.0_jrk_rk /)
  P_21 = (/ 1.0_jrk_rk, 0.0_jrk_rk /)
  P_12 = (/ 0.0_jrk_rk, 1.0_jrk_rk /)
  P_22 = (/ 1.0_jrk_rk, 1.0_jrk_rk /)

  call spl_mapping_bilinear(mapping, P_11, P_12, P_21, P_22)
!  call mapping % read_from_file('mapping.nml')
!  call mapping % export('output/test_projector_2d_mapping.nml')
  ! ... 

  ! ...
  call context % create(filename=filename_params)
  call context % print_info()
  ! ...

  call test_context % create(filename=filename_params)
  call test_context % print_info()
  ! ...
  
  ! ...
  call ddm_params % read_from_file(context % ptr_parameters % ddm_filename)
  call ddm_params % print_info()
  ! ...
  
  ! ...
  call jrk_allocate_space (trial_space, context)
  call jrk_allocate_space (test_space, test_context)
  ! ...

  ! ... TODO duplicate must be deferred
!  call trial_space % duplicate(test_space)
  ! ...

  ! ... allocation is postponed until we give the ddm
  call matrix % create(n_block_rows=1, n_block_cols=1)
  ! ...

  ! ... create fields with the trial space
  call phi % create(space=trial_space)
  ! ...

  ! ... allocation is postponed until we give the ddm
  call rhs % create(n_blocks=1)
  ! ...

  ! ...
  call norms % create(n_size=n_norms)
  ! ...

  ! ...
  call assembler % create( trial_space                 &
                       & , test_space                  &
                       & , n_fields   = 1              &
                       & , n_rhs      = 1              &
                       & , n_matrices = 1              &
                       & , norms = norms               &
                       & , mapping = mapping           &
                       & , ddm_parameters = ddm_params )
  ! ...

  ! ...
  call assembler % set_matrix(matrix, i_matrix=1)
  call assembler % set_rhs(rhs, i_rhs=1)
  call assembler % set_field(phi, i_field=1)
  ! ...

  ! ...
  assembler % ptr_matrix_contribution => build_matrix_stiffness
  assembler % ptr_rhs_contribution    => build_rhs

  call assembler % assemble(mapping=mapping)
  ! ...

  ! ...
  call matrix % export("matrix.mm")
  call rhs % export("rhs.txt", i_format=plf_vector_format_txt)
  ! ...

  ! ...
  call linear_solver % create(matrix, filename_solver)
  ! ...

  ! ...
  allocate(x(trial_space % ptr_numbering % n_vertex))
  allocate(y(test_space % ptr_numbering % n_vertex))

  x = 0.0_plf_rk
  y = 0.0_plf_rk
  call rhs % get(y)
  call linear_solver % solve(y, x)
  call phi % set(x)
  ! ...

  ! ...
  assembler % ptr_fields_evaluation   => evaluate_fields
  assembler % ptr_norms_evaluation    => evaluate_norms

  call assembler % assemble(mapping=mapping)
  ! ...

  ! ...
  call phi % export("phi.nml")

  print *, ">> norms : "
  call norms % print_info()
  ! ...

  ! ...
  deallocate(x)
  deallocate(y)

  call linear_solver % free()
  call assembler % free()
  call norms % free()
  call rhs % free()
  call phi % free()
  call test_space % free()
  call trial_space % free()
  call mapping % free()
  call test_context % free()
  call context % free()
  ! ...

  ! ...
  call jrk_finalize(i_err)
  ! ...

contains
  ! ..................................................
  subroutine evaluate_fields(self)
  implicit none
    class(jrk_t_assembler_abstract), intent(inout) :: self
    ! local
    integer :: i_trial_basis
    integer :: i_point
    integer :: n_points
    real(plf_rk) :: vj_0 
    real(plf_rk) :: vj_x 
    real(plf_rk) :: vj_y 

    ! ...
    do i_point = 1, self % n_points
    
      vj_0  = self % arr_vj(1, 1, i_point)
      vj_x  = self % arr_vj(1, 2, i_point)
      vj_y  = self % arr_vj(1, 3, i_point)
    
      self % field_values(1, 1, :, i_point) = self % field_values(1, 1, :, i_point) &
                                    & + vj_0 * self % field_coeffs(:)  

      self % field_values(1, 2, :, i_point) = self % field_values(1, 2, :, i_point) &
                                    & + vj_x * self % field_coeffs(:)  

      self % field_values(1, 3, :, i_point) = self % field_values(1, 3, :, i_point) &
                                    & + vj_y * self % field_coeffs(:)  
    end do
    ! ...

  end subroutine evaluate_fields
  ! ..................................................

  ! ..................................................
  subroutine evaluate_norms(self)
  implicit none
    class(jrk_t_assembler_abstract), intent(inout) :: self
    ! local
    integer :: i_diagnostic
    integer :: i_field
    integer :: i_point
    integer :: n_points
    real(plf_rk) :: wvol 
    real(plf_rk) :: norm_l2 
    real(plf_rk) :: norm_h1 
    real(plf_rk) :: f_0
    real(plf_rk) :: f_x
    real(plf_rk) :: f_y
    real(plf_rk) :: field_0
    real(plf_rk) :: field_x
    real(plf_rk) :: field_y

    ! ...
    n_points = self % n_points
    ! ...

    ! ...
    i_diagnostic = 0
    do i_field = 1, self % n_fields
      ! ...
      norm_l2 = 0.0_plf_rk
      norm_h1 = 0.0_plf_rk
      do i_point = 1, n_points

        f_0 = 0.0_plf_rk
        f_x = 0.0_plf_rk
        f_y = 0.0_plf_rk

        field_0 = self % field_values(1, 1, i_field, i_point) 
        field_x = self % field_values(1, 2, i_field, i_point) 
        field_y = self % field_values(1, 3, i_field, i_point) 

        wvol = self % weights(i_point) * self % jacobians(i_point) 

        norm_l2 = norm_l2 + (field_0 - f_0)**2 * wvol
        norm_h1 = norm_h1 + ((field_x - f_x)**2 + (field_y - f_y)**2) * wvol
      end do
      ! ...

      ! ...
      i_diagnostic = i_diagnostic + 1
      call self % ptr_norms % add_value(i_diagnostic, norm_l2)
      ! ...

      ! ...
      i_diagnostic = i_diagnostic + 1
      call self % ptr_norms % add_value(i_diagnostic, norm_h1)
      ! ...
    end do
    ! ...

  end subroutine evaluate_norms
  ! ..................................................

  ! ..................................................
  subroutine build_rhs(self)
  implicit none
    class(jrk_t_assembler_abstract), intent(inout) :: self
    ! local
    integer :: i_point
    integer :: n_points
    real(plf_rk) :: f_value
    real(plf_rk) :: vi_0 
    real(plf_rk) :: wvol 
    real(kind=plf_rk), parameter :: pi = 3.1415926535897931 
    real(kind=plf_rk), parameter :: k1 = 2.0 * pi 
    real(kind=plf_rk), parameter :: k2 = 2.0 * pi 
    real(kind=plf_rk), dimension(self % n_points) :: f_values

    ! ...
    n_points = self % n_points
    ! ...

    ! ...
    self % rhs_contribution = 0.0_plf_rk  
    ! ...

    ! ...
    f_values = sin(k1 * self % arr_x (1, :)) &
           & * sin(k2 * self % arr_x (2, :))
    ! ...

    ! ...
    do i_point = 1, n_points
      vi_0  = self % arr_vi(1, 1, i_point)

      f_value = f_values(i_point) 

      wvol = self % weights(i_point) * self % jacobians(i_point) 
   
      self % rhs_contribution = self % rhs_contribution + f_value * vi_0 * wvol
    end do
    ! ...

  end subroutine build_rhs
  ! ..................................................

end program main
! ............................................
