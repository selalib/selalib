!>  @brief
!> module for evaluation of a polynomial
!> in pp form using Horner's algorithm
!
!> @authors:
!>  Celine Caldini-Queiros  - celine.caldini@gmail.com
!
module sll_m_horner
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  implicit none

  public :: sll_f_horner_1d_eval, &
       sll_f_horner_2d_eval, &
       sll_f_horner_3d_eval

  private
  
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains

  !.................................................
  !> @brief     horner algorithm in 1D
  !> @description evaluation +  derivative of the polynomials
  !>
  !> @param[in] coeffs: pp coeffs 
  !> @param[in] x : point of evaluation
  !> @param[out] res : results
  function sll_f_horner_1d_eval(coeffs, x, jderiv) result (res)
  implicit none
    real(8), dimension (:), intent (in) :: coeffs
    real(8), intent (in)  :: x
    integer, intent (in)  :: jderiv 
    real(8)                             :: res
    integer :: i
    integer :: k
    real(8) :: fmmjdr 

    k = ubound(coeffs,1)
    fmmjdr = real(k - 1 - jderiv, 8)
    res = coeffs(k)

    do i = 1,k-1-jderiv
       res = (res / fmmjdr) * x + coeffs (k-i) 
       fmmjdr = fmmjdr - 1._8
    end do
  end function sll_f_horner_1d_eval
  !.................................................
  
  !.................................................
  !> @brief     horner algorithm in 2D
  !>
  !> @param[in] coeffs: pp coeffs 
  !> @param[in] x : point of evaluation (2D)
  !> @param[out] res : results
  function sll_f_horner_2d_eval (coeffs, x,jderiv) result (res)
  implicit none
    real(8), dimension (:,:), intent (in) :: coeffs
    real(8), dimension(:), intent (in) :: x
    integer, dimension(:), intent (in)  :: jderiv 
    real(8)                            :: res

    !local
    real(8), dimension(ubound(coeffs,2)) :: coeffi
    integer :: j
    integer :: k_v

    k_v = ubound(coeffs,2)
    coeffi(:) = 0.0_8

    do j = 1 , k_v
      coeffi(j) = sll_f_horner_1d_eval(coeffs(:,j),x(1),jderiv(1))
    end do
    
    res = sll_f_horner_1d_eval(coeffi,x(2),jderiv(2)) 

  end function sll_f_horner_2d_eval
  !.................................................
 
  !...............................................
  !> @brief     horner algorithm in 3D
  !>
  !> @param[in] coeffs: pp coeffs 
  !> @param[in] x : point of evaluation (3D)
  !> @param[out] res : results
  function sll_f_horner_3d_eval (coeffs, x,jderiv) result (res)
  implicit none
    real(8), dimension (:,:,:), intent (in) :: coeffs
    real(8), dimension(:), intent (in) :: x
    integer, dimension(:), intent (in)  :: jderiv 
    real(8):: res
    
    !local
    real(8), dimension(ubound(coeffs,1),ubound(coeffs,1)) :: coeffij
    real(8), dimension(ubound(coeffs,1)) :: coeffi
    integer :: i
    integer :: j 
    integer :: k_u
    integer :: k_v
    integer :: k_z

    k_u = ubound(coeffs,1)
    k_v = ubound(coeffs,2)
    k_z = ubound(coeffs,3)

    do i = 1 , k_z   
       coeffi(i) = sll_f_horner_2d_eval(coeffs(:,:,i),x(1:2),jderiv(1:2))
    end do
    res = sll_f_horner_1d_eval(coeffi,x(3),jderiv(3)) 


  end function sll_f_horner_3d_eval
  !...............................................

end module sll_m_horner
