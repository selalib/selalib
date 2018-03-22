module sll_m_poisson_2d_fem_ssm_assembler
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision, only: f64

  use sll_m_elliptic_2d_fem_ssm_weak_form, only: sll_c_elliptic_2d_fem_ssm_weak_form

  implicit none

  public :: sll_t_poisson_2d_fem_ssm_assembler

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type :: sll_t_poisson_2d_fem_ssm_assembler

    integer :: n1
    integer :: n2

    class(sll_c_elliptic_2d_fem_ssm_weak_form), pointer :: weak_form

  contains

    procedure :: init        => s_poisson_2d_fem_ssm_assembler__init
    procedure :: add_element => s_poisson_2d_fem_ssm_assembler__add_element

  end type sll_t_poisson_2d_fem_ssm_assembler

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Initializer
  subroutine s_poisson_2d_fem_ssm_assembler__init( self, n1, n2, weak_form )
    class(sll_t_poisson_2d_fem_ssm_assembler)         , intent(inout) :: self
    integer                                           , intent(in   ) :: n1
    integer                                           , intent(in   ) :: n2
    class(sll_c_elliptic_2d_fem_ssm_weak_form), target, intent(in   ) :: weak_form

    self % n1 = n1
    self % n2 = n2

    self % weak_form => weak_form

  end subroutine

  ! Add element in stiffness and mass matrices
  subroutine s_poisson_2d_fem_ssm_assembler__add_element( &
    self        , &
    k1          , &
    k2          , &
    data_1d_eta1, &
    data_1d_eta2, &
    int_volume  , &
    inv_metric  , &
    A           , &
    M )
    class(sll_t_poisson_2d_fem_ssm_assembler), intent(in   ) :: self
    integer                                  , intent(in   ) :: k1
    integer                                  , intent(in   ) :: k2
    real(wp)                                 , intent(in   ) :: data_1d_eta1(:,:,:,:)
    real(wp)                                 , intent(in   ) :: data_1d_eta2(:,:,:,:)
    real(wp)                                 , intent(in   ) :: int_volume(:,:,:,:)
    real(wp)                                 , intent(in   ) :: inv_metric(:,:,:,:,:,:)
    real(wp)                                 , intent(inout) :: A(:,:)
    real(wp)                                 , intent(inout) :: M(:,:)

    integer  :: i, j, i1, i2, j1, j2, p1, p2
    real(wp) :: Aij, Mij
    integer  :: il1, il2, jl1, jl2

    ! Extract degree of B-splines
    p1 = size( data_1d_eta1, 3 ) - 1
    p2 = size( data_1d_eta2, 3 ) - 1

    associate( n1 => self % n1, n2 => self % n2 )

      ! Cycle over basis functions in 2D trial space (same as test space)
      do jl2 = 1, p2+1

        j2 = modulo( k2 + jl2 - 2, n2 )

        do jl1 = 1, p1+1

          j1 = k1 + jl1 - 2

          ! Second matrix global index
          j = j1 * n2 + j2 + 1

          ! Cycle over basis functions in 2D test space
          do il2 = 1, p2+1

            i2 = modulo( k2 + il2 - 2, n2 )

            do il1 = 1, p1+1

              i1 = k1 + il1 - 2

              ! First matrix global index
              i = i1 * n2 + i2 + 1

              call self % weak_form % element( &
                test_values_and_derivs_eta1  = data_1d_eta1(:,:,il1,k1)  , &
                test_values_and_derivs_eta2  = data_1d_eta2(:,:,il2,k2)  , &
                trial_values_and_derivs_eta1 = data_1d_eta1(:,:,jl1,k1)  , &
                trial_values_and_derivs_eta2 = data_1d_eta2(:,:,jl2,k2)  , &
                int_volume                   = int_volume(:,:,k1,k2)    , &
                inv_metric                   = inv_metric(:,:,k1,k2,:,:), &
                Aij = Aij, &
                Mij = Mij )

              A(i,j) = A(i,j) + Aij
              M(i,j) = M(i,j) + Mij

            end do
          end do ! End cycle over basis functions in 2D test space

        end do
      end do ! End cycle over basis functions in 2D trial space

    end associate

  end subroutine s_poisson_2d_fem_ssm_assembler__add_element

end module sll_m_poisson_2d_fem_ssm_assembler
