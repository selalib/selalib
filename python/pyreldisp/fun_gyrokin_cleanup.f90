!-------------------------------------------------------------------------------
!**BEGIN PROLOGUE Function_Input_Module 
!**DATE WRITTEN   1999/06/01  (year/month/day)
!**AUTHORS
!
!  Peter Kravanja, Marc Van Barel
!
!    Department of Computer Science
!    Katholieke Universiteit Leuven
!    Celestijnenlaan 200 A
!    B-3001 Heverlee
!    BELGIUM
!    E-mail: Peter.Kravanja@na-net.ornl.gov
!            Marc.VanBarel@cs.kuleuven.ac.be
!
!  Omiros Ragos, Michael N. Vrahatis, and Filareti A. Zafiropoulos
!
!    Department of Mathematics
!    University of Patras
!    GR-261.10 Patras
!    GREECE
!    E-mail: {ragos,vrahatis,phikapa}@math.upatras.gr
!
!**DESCRIPTION
!
!  The user can specify the function (as well as its first derivative)
!  that is to be considered, the region where this function is analytic
!  and the values of the parameters NEWTONZ and NEWTONF, which are used
!  in the iterative refinement procedure that improves the computed 
!  approximations for the zeros. 
!
!**END PROLOGUE Function_Input_Module
!-------------------------------------------------------------------------------
!
!  Modified to add plasma dispersion functions for 
!    4D drift-kinetic cylindrical case
!   based on the complex error function. 
!
!  - Eric Sonnendrucker: 2008/12/11
!  - Modified by Virginie Grandgirard: 2013/12/11
!-------------------------------------------------------------------------------

module function_input_module

  use precision_module

  implicit none

  integer, parameter :: max=61
  double precision   :: ordre_grandeur = 1
  complex (kind=dp), allocatable :: vector(:,:)

  public

  !-----------------------------------------------------------------------------
  ! If ICON = 3 or 4 (as specified in Zeal_Input_Module), then specify 
  ! the values of NEWTONZ and NEWTONF.
  ! These variables are used as follows. The modified Newton's method,
  ! which takes into account the multiplicity of a zero and converges
  ! quadratically, is used to refine the calculated approximations for
  ! the zeros.
  !
  ! The iteration stops if:
  ! the relative distance between two successive approximations is
  ! at most NEWTONZ
  ! or
  ! the absolute value of the function at the last approximation is
  ! at most NEWTONF
  ! or
  ! if a maximum number of iterations is exceeded.
  !-----------------------------------------------------------------------------
  real (kind=dp) :: NEWTONZ = 5.0E-08_DP
  real (kind=dp) :: NEWTONF = 1.0E-14_DP
  real (kind=dp) :: PI = 3.1415926535897931_DP

  integer        :: Zi
  real (kind=dp) :: kzeta, ktheta
  integer        :: NNr
  real (kind=dp) :: dr
  real (kind=dp) :: B0

  real (kind=dp), allocatable :: Ti     (:)
  real (kind=dp), allocatable :: dlogTi (:)
  real (kind=dp), allocatable :: Te     (:)
  real (kind=dp), allocatable :: dlogn0 (:)
  real (kind=dp), allocatable :: ddlogn0(:)
  real (kind=dp), allocatable :: rmesh  (:)

  real (kind=dp), allocatable :: btheta(:)
  real (kind=dp), allocatable :: bz    (:)

  ! matrix of operator
  complex (kind=dp), allocatable :: A(:,:) 
  ! derivative of A (diagonal)
  complex (kind=dp), allocatable :: Ap(:)   

  contains


  !-----------------------------------------------------------------------------
  ! Define the function F and its derivative DF.
  !-----------------------------------------------------------------------------
  subroutine FDF(omega,F,DF)

    complex (kind=dp), intent(in   ) :: omega
    complex (kind=dp), intent(  out) :: F, DF

    complex (kind=dp) :: D, Dprime

    call matrix(omega)
    call matrix_derivative(omega)
    call det2(D, Dprime)
    F  = D     /ordre_grandeur
    DF = Dprime/ordre_grandeur

  end subroutine FDF


  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine integral(omega,r,int)
    
    complex (kind=dp), intent(in   ) :: omega
    integer          , intent(in   ) :: r
    complex (kind=dp), intent(  out) :: int

    complex (kind=dp) :: zeta, dzeta, rho
    real    (kind=dp) :: kstar

    kstar = (kzeta*bz(r)+(ktheta*btheta(r)/rmesh(r)))*sqrt(2.0_dp*Ti(r))
    rho   = omega/kstar

    call FriedConte(rho,zeta,dzeta)

    int = (1+rho*zeta)/Ti(r) & 
          -ktheta/(rmesh(r)*B0*kstar)*((dlogn0(r)-0.5_DP*dlogTi(r))*zeta &
                                       +dlogTi(r)*rho*(1+rho*zeta))
  end subroutine integral


  !-----------------------------------------------------------------------------
  ! A is the tridiagonal matrix of equation (35) of
  ! D. Coulette, N. Besse / Journal of Computational Physics 248 (2013) 1-32
  !
  ! A is filled for a given value of \omega as follows:
  ! A(1,i) is the upper diagonal term in column i
  ! A(2,i) is the       diagonal term in column i
  ! A(3,i) is the lower diagonal term in column i 
  !-----------------------------------------------------------------------------
  subroutine matrix(omega)

    complex (kind=dp), intent(in) :: omega

    complex (kind=dp) :: zeta, dzeta, int
    integer           :: i
    
    do i = 1, NNr 

      call integral(omega,i,int)

      A(2,i) = 2._dp+dr**2 &
               *((ktheta**2-0.25_dp)/rmesh(i)**2+1._dp/(Zi*Te(i)) &
                 +0.5_dp*(0.5_dp*dlogn0(i)**2 &
                          + dlogn0(i)/rmesh(i)+ddlogn0(i))+int)
      A(1,i) = -1.0_dp
      A(3,i) = -1.0_dp

    enddo

  end subroutine matrix


  !-----------------------------------------------------------------------------
  ! Ap is the element-wise derivative of A with respect to \omega.
  ! The only non-vanishing derivatives are on the diagonal elements.
  !-----------------------------------------------------------------------------
  subroutine matrix_derivative(omega)

    complex (kind=dp), intent(in) :: omega

    complex (kind=dp) :: zeta, dzeta, rho
    real    (kind=dp) :: kstar
    integer           :: i

    do i = 1, NNr  

      kstar = (kzeta*bz(i)+(ktheta*btheta(i)/rmesh(i)))*sqrt(2.0_dp*Ti(i))
      rho   = omega/kstar

      call FriedConte(rho,zeta,dzeta)

      Ap(i) = dr*dr &
              *((zeta+rho*dzeta)/Ti(i)-ktheta/(B0*rmesh(i)*kstar) &
                *(dlogn0(i)*dzeta+dlogTi(i) &
                  *(1._dp+2._dp*rho*zeta+(rho**2-0.5_dp)*dzeta)))/kstar
    enddo

  end subroutine matrix_derivative


  !-----------------------------------------------------------------------------
  ! Implement equation (36) of
  ! D. Coulette, N. Besse / Journal of Computational Physics 248 (2013) 1-32
  !
  ! At step i:
  ! D     is D_{i  }(\omega) of equation (36);
  ! D_im1 is D_{i-1}(\omega) of equation (36);
  ! D_im2 is D_{i-2}(\omega) of equation (36).
  !-----------------------------------------------------------------------------
  subroutine det(D)

    complex (kind=dp), intent(out) :: D

    complex (kind=dp) :: D_im1, D_im2
    integer           :: i

    D_im1 = A(2,1)
    D_im2 = 1

    do i = 2, NNr-1

      D = A(2,i)*D_im1-D_im2

      D_im2 = D_im1
      D_im1 = D

    enddo

  end subroutine det


  !-----------------------------------------------------------------------------
  ! Implement equations (36) and (37) of
  ! D. Coulette, N. Besse / Journal of Computational Physics 248 (2013) 1-32
  !
  ! At step i:
  ! D     is  D_{i  }(\omega) of equation (36);
  ! D_im1 is  D_{i-1}(\omega) of equation (36);
  ! D_im2 is  D_{i-2}(\omega) of equation (36);
  ! Dprime     is D'_{i  }(\omega) of equation (37);
  ! Dprime_im1 is D'_{i-1}(\omega) of equation (37);
  ! Dprime_im2 is D'_{i-2}(\omega) of equation (37);
  !-----------------------------------------------------------------------------
  subroutine det2(D, Dprime)

    complex (kind=dp), intent(out) :: D, Dprime
    
    complex (kind=dp) :: D_im1, D_im2
    complex (kind=dp) :: Dprime_im1, Dprime_im2
    integer           :: i

    D_im1 = A(2,1)
    D_im2 = 1
    Dprime_im1 = Ap(1)
    Dprime_im2 = 0

    do i = 2, NNr-1

      D      = A(2,i)*D_im1-D_im2
      Dprime = Ap(i)*D_im1+A(2,i)*Dprime_im1-Dprime_im2

      D_im2 = D_im1
      D_im1 = D
      Dprime_im2 = Dprime_im1
      Dprime_im1 = Dprime

    enddo

  end subroutine det2


  !-----------------------------------------------------------------------------
  subroutine init

    complex (kind=dp) :: D, omega
    real    (kind=dp) :: tab(20)
    integer           :: i

    if (.not.allocated(A) ) allocate(A (3,NNr))
    if (.not.allocated(Ap)) allocate(Ap(NNr  ))

    do i = 1, 20
      omega = cmplx(-1+i*0.1,0.0)
      call matrix(omega)
      call det(D)
      tab(i) = real(D)
    enddo

    ordre_grandeur = maxval(tab)

  end subroutine init


  !-----------------------------------------------------------------------------
  ! Compute the eigenvector of A corresponding to the eigenvalue \omega
  ! by singular value decomposition (SVD) of A.
  !
  ! Use zgbbrd and zbdsqr subroutines from LAPACK.
  !-----------------------------------------------------------------------------
  subroutine kernel(omega)

    complex (kind=dp), intent(in) :: omega

    complex (kind=dp) :: PT    (NNr,NNr)
    complex (kind=dp) :: Q     (NNr,NNr)
    complex (kind=dp) :: C     (NNr,NNr)
    complex (kind=dp) :: WORK  (NNr    )
    real    (kind=dp) :: D     (NNr    )
    real    (kind=dp) :: E     (NNr-1  )
    real    (kind=dp) :: RWORK (NNr    )
    real    (kind=dp) :: RRWORK(4*NNr  )
    integer           :: i, j, ii, INFO
    
    call matrix(omega)

    call zgbbrd('B',NNr,NNr,0,1,1,A,3,D,E,Q,NNr,PT,NNr,C,NNr,WORK,RWORK,INFO)
    if ( INFO .ne. 0 ) then
      write(*,'(a)') 'zeal failed because there is an error in the compute &
                     &of singular value decomposition'
      stop
    end if

    call zbdsqr('U',NNr,NNr,NNr,0,D,E,PT,NNr,Q,NNr,C,1,RRWORK,INFO)
    if ( INFO .ne. 0 ) then
      write(*,'(a)') 'zeal failed because there is an error in the compute &
                     &of singular value decomposition'
      stop
    end if

    ! Singular values are given in D in decreasing order, so
    ! the kernel corresponds to the last singular values.
    ii = 0
    do while ( D(NNr-ii) .lt. 1.e-10 )
      ii = ii+1
    enddo

    write (*,'(a,i0)') 'Dimension of kernel = ', ii
    write (*,'(a,f16.14,a,f16.14,a)') &
      'omega = (', real(omega), '+', imag(omega), 'j)'

    if (ii .eq. 0) then
      write(*,'(a,f16.14,f16.14,a)') 'zealpy failed because A(', omega, ') &
                                     &has 0-dimensional kernel'
      stop
    else
      if (.not.allocated(vector)) allocate(vector(NNr, ii))

      do i = 1, NNr
        do j = 1, ii
          vector(i,j) = PT(NNr-ii+j, i)
          ! print*,'vect propre ',i,vector(i,j)
        enddo
      enddo
    endif

  end subroutine kernel


  !-----------------------------------------------------------------------------
  ! Given a rectangular region specified by LV and H (see the module
  ! Zeal_Input_Module), decide whether the function is analytic
  ! inside this region or not.
  !-----------------------------------------------------------------------------
  function valreg(lv,h)

    logical valreg
    real (kind=dp), intent(in) :: lv(2), h(2)

    valreg = .true.
    !  The following statement can be used for functions that have a
    !  branch cut along the non-positive real axis.
    !
    !   valreg = .not. ( lv(2)*(lv(2)+h(2)) <= zero .and. lv(1) <= zero )

  end function valreg

end module function_input_module
