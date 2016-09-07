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
!  Modified to add plasma dispersion functions for 4D drift-kinetic cylindrical
!  case based on the complex error function. 
!
!  - Eric Sonnendrucker: 2008/12/11
!  - Modified by Virginie Grandgirard: 2013/12/11
!  - Modified by Yaman Gü¢lü and Edoardo Zoni: 2016/08
!
!  The equation solved is the quasi-neutrality equation
!
!  ((-1/(r*n0))*d/dr(r*n0*d/dr)+m**2/r**2+1/Ti+I)phi=0,
!
!  where phi=phi(omega,r) is the Laplace transform of the Fourier coefficient
!  of the field for a given mode (m,n). The integral I=I(omega,r) can be
!  expressed analytically in terms of the plasma dispersion function Z:
!
!  I=(1+rho*Z(rho))/Ti-m/(r*B0*kstar)
!    *((dlogn0(r)-0.5_DP*dlogTi(r))*Z(rho)+dlogTi(r)*rho*(1+rho*Z(rho))),
!
!  with rho:=omega/kstar. The above equation is solved via the replacement
!  psi:=sqrt(r*n0)*phi. The resulting "Schroedinger-type" equation is
!
!  (-d**2/dr**2+m**2/r**2+1/Ti+1/(sqrt(r*n0))*d**2/dr**2(sqrt(r*n0))+I)psi=0.
!
!  The discrete problem E(omega)psi=0, resulting from the discretization
!  of the above radial differential operator on a uniform grid, is solved
!  through two steps:
!
!  1) solve the dispersion equation D(omega)=det(E(omega))=0;
!  2) find the matching eigenvector by singular value decomposition of E(omega).
!
!  For references, see:
!
!  1) sections 2.3.1 and 3.3 of
!     D. Coulette, N. Besse / Journal of Computational Physics 248 (2013) 1-32;
!
!  2) appendix B of
!     Latu, Mehrenberger, Gü¢lü, Ottaviani, Sonnendrücker /
!     Field-aligned interpolation for semi-Lagrangian gyrokinetic simulations
!     (submitted on May 9, 2016).  
!-------------------------------------------------------------------------------

module function_input_module

  use precision_module

  implicit none

  integer, parameter :: max=61
  double precision   :: order_magnitude = 1
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
  complex (kind=dp), allocatable :: Aprime(:)   

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
    F  = D     /order_magnitude
    DF = Dprime/order_magnitude

  end subroutine FDF


  !-----------------------------------------------------------------------------
  ! Compute the integral I calling the plasma dispersion function
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
  ! Aprime is the element-wise derivative of A with respect to \omega.
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

      Aprime(i) = dr*dr &
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
    Dprime_im1 = Aprime(1)
    Dprime_im2 = 0

    do i = 2, NNr-1

      D      = A(2,i)*D_im1-D_im2
      Dprime = Aprime(i)*D_im1+A(2,i)*Dprime_im1-Dprime_im2

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
    if (.not.allocated(Aprime)) allocate(Aprime(NNr  ))

    do i = 1, 20
      omega = cmplx(-1+i*0.1,0.0)
      call matrix(omega)
      call det(D)
      tab(i) = real(D)
    enddo

    order_magnitude = maxval(tab)

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
    integer           :: i, j
    integer           :: dim_kernel, INFO
    
    call matrix(omega)

    call zgbbrd('B',NNr,NNr,0,1,1,A,3,D,E,Q,NNr,PT,NNr,C,NNr,WORK,RWORK,INFO)
    if ( INFO .ne. 0 ) then
      write(*,'(a)') 'zeal failed because there is an error in the compute &
                     &of singular value decomposition'
      stop
    end if

    ! Previous call was
    ! call zbdsqr('U',NNr,NNr,NNr,0,D,EE,PT,NNr,Q,NNr,C,1,RRWORK,INFO)
    ! with EE of dimension NNr: wrong dimension for zbdsqr subroutine.
    call zbdsqr('U',NNr,NNr,NNr,0,D,E,PT,NNr,Q,NNr,C,1,RRWORK,INFO)
    if ( INFO .ne. 0 ) then
      write(*,'(a)') 'Zeal failed because there is an error in the compute &
                     &of singular value decomposition'
      stop
    end if

    ! Singular values are given in D in decreasing order, so
    ! the kernel corresponds to the last singular values.
    dim_kernel = 0
    do while ( D(NNr-dim_kernel) .lt. 1.e-10 )
      dim_kernel = dim_kernel+1
    enddo

    write (*,'(a,i0)') 'Dimension of kernel = ', dim_kernel
    write (*,'(a,f20.18,a,f20.18,a)') &
          'omega = (', real(omega), '+', imag(omega), 'j)'

    if (dim_kernel .eq. 0) then
      write(*,'(/,a,f20.18,a,f20.18,a,/,a,/)') 'Zealpy failed because A(', &
           real(omega), '+', imag(omega), 'j)', 'has 0-dimensional kernel.'
      stop
    else
      if (.not.allocated(vector)) allocate(vector(NNr, dim_kernel))

      do i = 1, NNr
        do j = 1, dim_kernel
          vector(i,j) = PT(NNr-dim_kernel+j, i)
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
  function valreg(LV,H)

    logical valreg
    real (kind=dp), intent(in) :: LV(2), H(2)

    valreg = .true.
    !  The following statement can be used for functions that have a
    !  branch cut along the non-positive real axis.
    !
    !   valreg = .not. ( LV(2)*(LV(2)+H(2)) <= zero .and. LV(1) <= zero )

  end function valreg

end module function_input_module
