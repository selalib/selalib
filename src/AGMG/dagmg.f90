! COPYRIGHT (c) 2012 Universite Libre de Bruxelles (ULB)
!
! This file is part of AGMG software package
! Release 3.2.1-aca built on "Mar 20 2014" by Yvan Notay
!
! ALL USAGE OF AGMG IS SUBJECT TO LICENSE. PLEASE REFER TO THE FILE "LICENSE".
! IF YOU OBTAINED A DCOPY OF THIS SOFTWARE WITHOUT THIS FILE,
! PLEASE CONTACT ynotay@ulb.ac.be
!
! In particular, if you have a free academic license:
!
! (1) You must be a member of an educational, academic or research institution.
!     The license agreement automatically terminates once you no longer fulfill
!     this requirement.
!
! (2) You are obliged to cite AGMG in any publication or report as:
!     "Yvan Notay, AGMG software and documentation;
!      see http://homepages.ulb.ac.be/~ynotay/AGMG".
!
! (3) You may not make available to others the software in any form, either
!     as source or as a precompiled object.
!
! (4) You may not use AGMG for the benefit of any third party or for any
!     commercial purposes. Note that this excludes the use within the
!     framework of a contract with an industrial partner.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! See the accompanying userguide for more details on how to use the software,
! and the README file for installation instructions.
!
! See the Web pages <http://homepages.ulb.ac.be/~ynotay/AGMG> for
! release information and possible upgrade.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DICLAIMER:
!    AGMG is provided on an "AS IS" basis, without any explicit or implied
!    WARRANTY; see the file "LICENSE" for more details.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   If you use AGMG for research, please observe that your work benefits
!   our past research efforts that allowed the development of AGMG.
!   Hence, even if you do not see it directly, the results obtained thanks
!   to the use of AGMG depend on the results in publications [1-3] below,
!   where the main algorithms used in AGMG are presented and justified.
!   It is then a normal duty to cite these publications (besides citing
!   AGMG itself) in any scientific work depending on the usage of AGMG,
!   as you would do with any former research result you are using.
!
! [1] Y. Notay, An aggregation-based algebraic multigrid method,
!    Electronic Transactions on Numerical Analysis, vol. 37, pp. 123-146, 2010
!
! [2] A. Napov and Y. Notay, An algebraic multigrid method with guaranteed
!    convergence rate, SIAM J. Sci. Comput., vol. 34, pp. A1079-A1109, 2012.
!
! [3] Y. Notay, Aggregation-based algebraic multigrid for convection-diffusion
!    equations, SIAM J. Sci. Comput., vol. 34, pp. A2288-A2316, 2012.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This file provides dagmg and dependencies; dagmg is a sequential
! implementation for real matrices in double precision
! of the method presented in [1], where the used algorithms are described
! in detail. From realease 3.x, the coarsening has been modified according
! to the results in [2,3], see the release notes in the README file
! for more details.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! SUBROUTINE dagmg (MAIN DRIVER): see bottom of this file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! PARAMETERS DEFINITON -  may be tuned by expert users
!-----------------------------------------------------------------------
  MODULE dagmg_mem
    SAVE
!
! INTEGER
!
!  maxlev  is the maximal number of levels
!          (should be large enough - much larger than needed is armless).
!  real_len is the length of 1 REAL(kind(0.0d0)) in byte
!        (used only to display information on memory usage).
!  nsmooth  is the number of pre- and post-smoothing steps
!       (see smoothtype for details) .
!
!  smoothtype indicates which smoother is used:
!
!    if smoothtype==1, the smoother is based on Gauss-Seidel sweeps;
!    if smoothtype==0, the smoother is based on SOR sweeps with automatic
!       estimation of the relaxation parameter (often reduces to Gauss-Seidel);
!    if smoothtype==-1, the smoother is based on SOR sweeps with relaxation
!       parameter specified in constant omega defined below;
!
!     Scheme used in all three cases: nsmooth sweeps for both pre- and
!     post-smoothing, with:
!        pre-smoothing: Forward sweep, then Backward sweep, then Forward, etc
!       post-smoothing: Backward sweep, then Forward sweep, then Backward, etc;
!
!    if smoothtype==2, the smoother is ILU(0), with
!      (nsmooth/2) pre- and (nsmooth/2+mod(nsmooth,1)) post-smoothing steps.
!
!  nstep   is the maximum number of coarsening steps;
!          nstep<0 means that coarsening is stopped only according to
!          the matrix size, see parameter maxcoarsesize.
!  nlvcyc  is the number of coarse levels (from bottom) on which V-cycle
!          formulation is enforced (Rmk: K-cycle always allowed at first
!          coarse level).
!  npass   is the maximal number of pairwise aggregation passes for each
!          coarsening step, according to the algorithms in [2,3].
!  maxcoarsesize: when the size of the coarse grid matrix is less than or
!                 equal to maxcoarsesize*N^(1/3),  it is factorized exactly
!                 and coarsening is stopped;
!         maxcoarsesizeslow: in case of slow coarsening,
!                 exact factorization can be performed when the size of
!                 the coarse grid matrix is less than or equal to
!                 maxcoarsesizeslow*N^(1/3).
!         (N is the number of rows of the input matrix)
    INTEGER, PARAMETER :: maxlev=40, real_len=8
    INTEGER, PARAMETER :: nsmooth=1, smoothtype=0, nstep=-1, nlvcyc=0
    INTEGER, PARAMETER :: npass=2,maxcoarsesize=40,maxcoarsesizeslow=400
!
! REAL
!
!  omega is the relaxation parameter for SOR sweeps (smoothtype=-1)
!  resi is the threshold t for the relative residual error in inner FCG & GCR
!       iterations, see Algorithm 3.2 in [1]
!  trspos is a threshold: if a row has a positive offdiagonal entry larger
!         than trspos times the diagonal entry, the corresponding node is
!         transferred unaggregated to the coarse grid
!  kaptg_ is the threshold used to accept or not a tentative aggregate
!         when applying the coarsening algorithms from [2,3];
!         kaptg_blocdia is used for control based on bloc diagonal smoother [2];
!         kaptg_dampJac is used for control based on Jacobi smoother [3].
!  checkdd is the threshold to keep outside aggregation nodes where
!         the matrix is strongly diagonally dominant (based on mean of row
!         and column);
!         In fact, AGMG uses the maximum of |checkdd| and of the default value
!            according to kaptg_ as indicated in [2,3]
!            (hence |checkdd| < 1 ensures that one uses this default value)
!         checkdd <0 : consider |checkdd|, but base the test on minus the
!                sum of offdiagonal elements, without taking the absolute value.
!  targetcoarsefac is the target coarsening factor (parameter tau in the main
!         coarsening algorithm in [2,3]): further pairwise aggregation passes
!         are omitted once the number of nonzero entries has been reduced by a
!         factor of at least targetcoarsefac.
!  fracnegrcsum: if, at some level, more than fracnegrcsum*nl nodes,
!         where nl is the total number of nodes at that level, have
!         negative mean row and column sum, then the aggregation algorithm
!         of [2,3] is modified, exchanging all diagonal entries for the mean
!         row and column sum (that is, the algorithm is applied to a
!         modified matrix with mean row and colum sum enforced to be zero).
    REAL(kind(0.0d0)), PARAMETER :: omega=0.8, resi=0.2, trspos=0.45
    REAL(kind(0.0d0)), PARAMETER :: kaptg_blocdia=8, kaptg_dampJac=10
    REAL(kind(0.0d0)), PARAMETER :: checkdd=0.5
    REAL(kind(0.0d0)), PARAMETER :: targetcoarsefac=2.0**npass
    REAL(kind(0.0d0)), PARAMETER :: fracnegrcsum=0.25
!!!!!!!!!!!!!!!!!!!! END of PARAMETERS DEFINITION -----------------
!!!!!!!!!!!!!!!!!!! Internal variables declaration
!
    TYPE InnerData
       REAL(kind(0.0d0)), DIMENSION(:), POINTER :: a
       INTEGER, DIMENSION(:), POINTER :: ja
       INTEGER, DIMENSION(:), POINTER :: ia
       INTEGER, DIMENSION(:), POINTER :: il
       INTEGER, DIMENSION(:), POINTER :: iu
       REAL(kind(0.0d0)), DIMENSION(:), POINTER :: p
       INTEGER, DIMENSION(:), POINTER :: idiag
       INTEGER, DIMENSION(:), POINTER :: ind
       INTEGER, DIMENSION(:), POINTER :: iext
       INTEGER, DIMENSION(:), POINTER :: ilstout
       INTEGER, DIMENSION(:), POINTER :: lstout
       INTEGER, DIMENSION(:), POINTER :: ilstin
       INTEGER, DIMENSION(:), POINTER :: lstin
       INTEGER, DIMENSION(:), POINTER :: iblockl
    END TYPE InnerData
!
    TYPE(InnerData) :: dt(maxlev)
    REAL(kind(0.0d0)), ALLOCATABLE :: scald(:)
    INTEGER :: nn(maxlev),kstat(2,maxlev)=0,innermax(maxlev)
    INTEGER :: nlev,nwrkcum,iout,nrst,nbblock
    INTEGER :: maxcoarset,maxcoarseslowt
    REAL(kind(0.0d0)) :: memi=0.0,memax=0.0,memr=0.0,mritr,rlenilen,smoothtp,omeg
    REAL(kind(0.0d0)) :: wcplex(4),fracnz(maxlev),flop=0.0d0
    LOGICAL :: spd,wfo,wff,woo,allzero,trans,transint,zerors,gcrin,coasing
    REAL(kind(0.0d0)), PARAMETER :: cplxmax=3.0, xsi=0.6d0
    REAL(kind(0.0d0)), PARAMETER :: repsmach=SQRT(EPSILON(1.0d0)),epsmach=EPSILON(1.0d0)
    INTEGER :: nlc(2),nlcp(2),nlc1(2),icum,imult
    REAL(kind(0.0d0)) :: ngl(2),nglp(2),nlctot(2),ngl1(2),ngltot(2),RELRESL1
    INTEGER, DIMENSION(:), POINTER :: iextL1,ilstoutL1,lstoutL1,ilstinL1,lstinL1
    INTEGER, PARAMETER :: IRANK=-9999
    REAL(kind(0.0d0)), PARAMETER ::    &
         checkddJ=MAX(ABS(checkdd),kaptg_dampJac/(kaptg_dampJac-2))
    REAL(kind(0.0d0)), PARAMETER ::    &
         checkddB=MAX(ABS(checkdd),(kaptg_blocdia+1)/(kaptg_blocdia-1))
  END MODULE dagmg_mem
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! END of Internal variables declaration
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! TIMING
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_mestime(id,cputm,eltm)
    IMPLICIT NONE
    INTEGER, SAVE :: cpt_init(10)=-1,cpt_fin,cpt_max,freq,cpt
    REAL, SAVE :: t1(10), t2
    REAL(kind(0.0d0)) :: cputm,eltm
    INTEGER :: id
    IF (id>0) THEN
       !Next line may be uncommented if FORTRAN 95 function
       !CPU_TIME is implemented
       !   CALL CPU_TIME(t2)
       CALL SYSTEM_CLOCK(cpt_fin,freq,cpt_max)
       !
       cpt = cpt_fin - cpt_init(id)
       IF (cpt_fin < cpt_init(id)) cpt = cpt + cpt_max
       eltm = dble(cpt) / freq
       cputm = dble(t2 - t1(id))
       !
    ELSE
       !
       CALL SYSTEM_CLOCK(cpt_init(-id),freq,cpt_max)
       !Next line may be uncommented if FORTRAN 95 function
       !CPU_TIME is implemented
       !   CALL CPU_TIME(t1(-id))
       !
    END IF
    RETURN
  END SUBROUTINE dagmg_mestime
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! END of TIMING
!-----------------------------------------------------------------------
  MODULE dagmg_PARDISSO
    TYPE PARDISO_HANDLE
       INTEGER(KIND=8) DUMMY
    END TYPE PARDISO_HANDLE
  END MODULE dagmg_PARDISSO
  MODULE dagmg_PARDISO
    USE dagmg_PARDISSO
    INTERFACE
       SUBROUTINE PARDISO( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA  &
            , PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
         USE dagmg_PARDISSO
         TYPE(PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
         INTEGER,          INTENT(IN)    :: MAXFCT
         INTEGER,          INTENT(IN)    :: MNUM
         INTEGER,          INTENT(IN)    :: MTYPE
         INTEGER,          INTENT(IN)    :: PHASE
         INTEGER,          INTENT(IN)    :: N
         INTEGER,          INTENT(IN)    :: IA(*)
         INTEGER,          INTENT(IN)    :: JA(*)
         INTEGER,          INTENT(IN)    :: PERM(*)
         INTEGER,          INTENT(IN)    :: NRHS
         INTEGER,          INTENT(INOUT) :: IPARM(*)
         INTEGER,          INTENT(IN)    :: MSGLVL
         INTEGER,          INTENT(OUT)   :: ERROR
         REAL(kind(0.0d0)), INTENT(IN)    :: A(*)
         REAL(kind(0.0d0)), INTENT(INOUT) :: B(*)
         REAL(kind(0.0d0)), INTENT(OUT)   :: X(*)
       END SUBROUTINE PARDISO
    END INTERFACE
  END MODULE dagmg_PARDISO
!-----------------------------------------------------------------------
    SUBROUTINE PARDISO( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA    &
         , PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
      USE dagmg_PARDISSO
      TYPE(PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
      INTEGER,          INTENT(IN)    :: MAXFCT
      INTEGER,          INTENT(IN)    :: MNUM
      INTEGER,          INTENT(IN)    :: MTYPE
      INTEGER,          INTENT(IN)    :: PHASE
      INTEGER,          INTENT(IN)    :: N
      INTEGER,          INTENT(IN)    :: IA(*)
      INTEGER,          INTENT(IN)    :: JA(*)
      INTEGER,          INTENT(IN)    :: PERM(*)
      INTEGER,          INTENT(IN)    :: NRHS
      INTEGER,          INTENT(INOUT) :: IPARM(*)
      INTEGER,          INTENT(IN)    :: MSGLVL
      INTEGER,          INTENT(OUT)   :: ERROR
      REAL(kind(0.0d0)), INTENT(IN)    :: A(*)
      REAL(kind(0.0d0)), INTENT(INOUT) :: B(*)
      REAL(kind(0.0d0)), INTENT(OUT)   :: X(*)
      ERROR=1
      X(1:N)=0.0d0
      STOP
    END SUBROUTINE PARDISO
!-----------------------------------------------------------------------
MODULE dagmg_ALLROUTINES
CONTAINS
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_relmem
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l
    !
    !
    DO l=1,nlev-1
       IF(nn(l) > 0) THEN
            DEALLOCATE(dt(l)%a,dt(l)%ja,dt(l)%il,dt(l)%iu)
            IF (smoothtp .NE. 1) DEALLOCATE(dt(l)%p)
            IF (nn(l+1) > 0) DEALLOCATE(dt(l)%ind)
       END IF
    END DO
    IF (nn(nlev) > 0 .AND. nlev > 1)         &
              DEALLOCATE(dt(nlev)%a,dt(nlev)%ja,dt(nlev)%ia)
    memi=0
    memr=0
    memax=0
    RETURN
  END SUBROUTINE dagmg_relmem
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_applyprec( N,f,X,a,ja,ia,ijb)
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: N
    REAL(kind(0.0d0)), OPTIONAL :: f(N), X(N)
    INTEGER, OPTIONAL :: ja(*), ia(N+1), ijb
    REAL(kind(0.0d0)), OPTIONAL :: a(*)
    REAL(kind(0.0d0)), ALLOCATABLE, SAVE :: S(:)
    REAL(kind(0.0d0)) :: tau=1.35
    INTEGER :: iter, nvec, i
    !
       iter=ijb-2
       nvec=min(2,iter)
       mritr=nwrkcum+nvec*N
       ALLOCATE (S(nwrkcum+nvec*N))
    IF (iter > 1) THEN
       CALL dagmg_prec_matv(1,f,X,S(N+1),S(2*N+1),.TRUE.,a,ja,ia)
       X(1:N)=tau*X(1:N)
       DO i=2,iter-1
          CALL dagmg_prec_matv(1,f,S,S(N+1),S(2*N+1),.TRUE.,a,ja,ia)
          X(1:N)=X(1:N)+tau*S(1:N)
          f(1:N)=f(1:N)-tau*S(N+1:2*N)
       END DO
       CALL dagmg_prec_matv(1,f,S,S(N+1),S(2*N+1),.FALSE.,a,ja,ia)
       X(1:N)=X(1:N)+tau*S(1:N)
    ELSE
       CALL dagmg_prec_matv(1,f,X,S,S(N+1),.FALSE.,a,ja,ia)
    END IF
    kstat(2,1)=kstat(2,1)+iter
    !
     DEALLOCATE(S)
    RETURN
  END SUBROUTINE dagmg_applyprec
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_partroword(n, a, ja, ia, idiag, w, iw, iext)
!
!
!
    IMPLICIT NONE
    INTEGER :: n, ja(*), ia(n+1), idiag(n)
    INTEGER, OPTIONAL :: iext(*)
    INTEGER, TARGET :: iw(*)
    REAL(kind(0.0d0)) :: a(*), p
    REAL(kind(0.0d0)), TARGET :: w(*)
    INTEGER :: i, j, k, ipos, nzu, nze
    DO i = 1, n
       ipos = ia(i)
       nzu = 0
       nze = 0
       DO k = ia(i), ia(i+1) - 1
          j = ja(k)
          IF (j.GT.i) THEN
             nzu = nzu + 1
             w(nzu) = a(k)
             iw(nzu) = j
          ELSEIF (j.LT.i) THEN
             a(ipos) = a(k)
             ja(ipos) = j
             ipos = ipos + 1
          ELSE
             p = a(k)
          ENDIF
       ENDDO
       !
       idiag(i) = ipos
       a(ipos) = p
       ja(ipos) = i
       !
       DO k = 1, nzu
          ipos = ipos + 1
          a(ipos) = w(k)
          ja(ipos) = iw(k)
       ENDDO
       !
    ENDDO
    RETURN
  END SUBROUTINE dagmg_partroword
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_csrdlu(n,a,ja,ia,idiag,ao,jao,iup,iext,iextop)
!
!
    IMPLICIT NONE
    INTEGER :: n, ja(*), jao(n+1:*)
    INTEGER, TARGET :: ia(n+1), idiag(n+1)
    INTEGER, OPTIONAL, TARGET :: iup(n+1), iext(n+1),iextop(n+1)
    REAL(kind(0.0d0)) :: a(*),ao(*)
    INTEGER :: i,j,ili,iui,iei,ipos,ili0,iui0,iei0,nl,nu,k,next
    INTEGER, POINTER :: il(:), iu(:), iexto(:)
    nl=0
    nu=0
    il => idiag(1:n+1)
    IF (PRESENT(iup)) THEN
       iu => iup(1:n+1)
       DO i=1,n
          ili=idiag(i)-ia(i)
          iui=ia(i+1)-idiag(i)-1
          iup(i)=iui
          idiag(i)=ili
          nl=nl+ili
          nu=nu+iui
       END DO
    ELSE
       iu => ia(1:n+1)
       DO i=1,n
          ili=idiag(i)-ia(i)
          iui=ia(i+1)-idiag(i)-1
          ia(i)=iui
          idiag(i)=ili
          nl=nl+ili
          nu=nu+iui
       END DO
    END IF
    ipos=0
    ili=n+1
    iui=ili+nl
    iei=iui+nu
       DO i=1,n
          ili0=ili
          iui0=iui
          DO j=1,il(i)
             ipos=ipos+1
             ao(ili)=a(ipos)
             jao(ili)=ja(ipos)
             ili=ili+1
          END DO
          ipos=ipos+1
          ao(i)=a(ipos)
          DO j=1,iu(i)
             ipos=ipos+1
             ao(iui)=a(ipos)
             jao(iui)=ja(ipos)
             iui=iui+1
          END DO
          iu(i)=iui0
          il(i)=ili0
       END DO
       iu(n+1)=iui
       il(n+1)=ili
    RETURN
  END SUBROUTINE dagmg_csrdlu
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_csrdluT(n,a,ja,ia,idiag,ao,jao,il,iu)
!
    IMPLICIT NONE
    INTEGER :: n, ja(*), ia(n+1), idiag(n), jao(n+1:*), il(n+1), iu(n+1)
    REAL(kind(0.0d0)) :: a(*),ao(*)
    INTEGER :: i,j,ili,iui,iei,ipos,ili0,iui0,iei0,nl,nu,k,next
    nl=0
    nu=0
       il(2:n+1)=0
       iu(2:n+1)=0
       DO i=1,n
          iui=idiag(i)-ia(i)
          ili=ia(i+1)-idiag(i)-1
          nl=nl+ili
          nu=nu+iui
          DO k=ia(i),idiag(i)-1
             iu(ja(k)+1)=iu(ja(k)+1)+1
          END DO
          DO k=idiag(i)+1,ia(i+1)-1
             il(ja(k)+1)=il(ja(k)+1)+1
          END DO
       END DO
    ipos=0
    ili=n+1
    iui=ili+nl
    iei=iui+nu
       il(1)=ili
       iu(1)=iui
       DO i=1,n
          il(i+1) = il(i) + il(i+1)
          iu(i+1) = iu(i) + iu(i+1)
       END DO
       DO i=1,n
          DO k=ia(i),idiag(i)-1
             j = ja(k)
             next = iu(j)
             ao(next) = a(k)
             jao(next) = i
             iu(j) = next+1
          END DO
          ao(i)=a(idiag(i))
          DO k=idiag(i)+1,ia(i+1)-1
             j = ja(k)
             next = il(j)
             ao(next) = a(k)
             jao(next) = i
             il(j) = next+1
          END DO
       END DO
       DO i=n,1,-1
          il(i+1) = il(i)
          iu(i+1) = iu(i)
       END DO
       il(1)=ili
       iu(1)=iui
    RETURN
  END SUBROUTINE dagmg_csrdluT
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_csrmv(n, a, ja, ifja, ia, x, ifx, y, iad, flop, lstin)
!
!
!
!
!
    IMPLICIT NONE
    INTEGER :: n, ifja, ja(ifja:*), ia(n+1), iad, ifx, i, k, j
    INTEGER, OPTIONAL :: lstin(0:*)
    REAL(kind(0.0d0)) :: a(*), x(ifx:*), y(n), t
    REAL(kind(0.0d0)) :: flop
    !
       IF (iad .LT. -1) THEN
          DO i=1,n
             t=y(i)
             DO k=ia(i),ia(i+1)-1
                t = t - a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       ELSE IF (iad .GT. 1) THEN
          DO i=1,n
             t=y(i)
             DO k=ia(i),ia(i+1)-1
                t = t + a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       ELSE IF (iad .EQ. -1) THEN
          DO i=1,n
             t=0.0d0
             DO k=ia(i),ia(i+1)-1
                t = t - a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       ELSE
          DO i=1,n
             t=0.0d0
             DO k=ia(i),ia(i+1)-1
                t = t + a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       END IF
       !
    flop=flop+dble(2*(ia(n+1)-ia(1)))
    RETURN
  END SUBROUTINE dagmg_csrmv
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_mcsrmv(n, a, ja, ifja, ia, x, ifx, y, iad, flop)
!
!
!
!
    IMPLICIT NONE
    INTEGER :: n, ifja, ja(ifja:*), ia(n+1), iad, ifx, i, k, j
    REAL(kind(0.0d0)) :: a(*), x(ifx:*), y(n), t
    REAL(kind(0.0d0)) :: flop
    !
       IF (iad .LT. -1) THEN
          DO i=1,n
             t=y(i)-a(i)*x(i)
             DO k=ia(i),ia(i+1)-1
                t = t - a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
          flop=flop+dble(2*(ia(n+1)-ia(1)+n))
       ELSE IF (iad .GT. 1) THEN
          DO i=1,n
             t=y(i)+a(i)*x(i)
             DO k=ia(i),ia(i+1)-1
                t = t + a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
          flop=flop+dble(2*(ia(n+1)-ia(1)+n))
       ELSE IF (iad .EQ. -1) THEN
          DO i=1,n
             t=-a(i)*x(i)
             DO k=ia(i),ia(i+1)-1
                t = t - a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
          flop=flop+dble(2*(ia(n+1)-ia(1))+n)
       ELSE
          DO i=1,n
             t=a(i)*x(i)
             DO k=ia(i),ia(i+1)-1
                t = t + a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
          flop=flop+dble(2*(ia(n+1)-ia(1))+n)
       END IF
       !
    RETURN
  END SUBROUTINE dagmg_mcsrmv
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_csrlsolve(n, a, ja, ifja, ia, p, x, y, iunit, flop)
!
!
!
    IMPLICIT NONE
    INTEGER :: n, ifja, ja(ifja:*), ia(n+1), iunit, i, k
    REAL(kind(0.0d0)) :: a(*), p(n), x(n), y(n), t
    REAL(kind(0.0d0)) :: flop
!
    IF (iunit .LT. 0) THEN
       x(1) = y(1)
       DO i=2,n
          t = 0.0d0
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = y(i) + p(i)*t
       END DO
       flop=flop+dble(n+2*(ia(n+1)-ia(1)))
    ELSE IF (iunit .GT. 0) THEN
       x(1) = p(1)*y(1)
       DO i=2,n
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = p(i)*t
       END DO
       flop=flop+dble(n+2*(ia(n+1)-ia(1)))
    ELSE
       x(1) = y(1)
       DO i=2,n
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = t
       END DO
       flop=flop+dble(2*(ia(n+1)-ia(1)))
    END IF
    RETURN
  END SUBROUTINE dagmg_csrlsolve
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_csrusolve(n, a, ja, ifja, ia, p, x, y, iunit, flop)
!
!
!
    IMPLICIT NONE
    INTEGER :: n, ifja, ja(ifja:*), ia(n+1), iunit, i, k
    REAL(kind(0.0d0)) :: a(*), p(n), x(n), y(n), t
    REAL(kind(0.0d0)) :: flop
!
    IF (iunit .LT. 0) THEN
       x(n) = y(n)
       DO i=n-1,1,-1
          t = 0.0d0
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = y(i) + p(i)*t
       END DO
       flop=flop+dble(n+2*(ia(n+1)-ia(1)))
    ELSE IF (iunit .GT. 0) THEN
       x(n) = p(n)*y(n)
       DO i=n-1,1,-1
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = p(i)*t
       END DO
       flop=flop+dble(n+2*(ia(n+1)-ia(1)))
    ELSE
       x(n) = y(n)
       DO i=n-1,1,-1
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = t
       END DO
       flop=flop+dble(2*(ia(n+1)-ia(1)))
    END IF
    RETURN
  END SUBROUTINE dagmg_csrusolve
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_FlexCG(N,f,X,ITER,RESID,a,ja,ia,init)
    !
    !
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER       :: N, ITER, init
    REAL(kind(0.0d0)) :: f(N), X(N)
    INTEGER       :: ja(*), ia(N+1)
    REAL(kind(0.0d0)) :: a(*)
    INTEGER       :: MAXIT, ierr, kk
    REAL(kind(0.0d0)) :: TOL, BNORM, RESID, dum0, RESID0
    REAL(kind(0.0d0)) :: ALPHA, BET0, RHO, dum(6), td
    REAL(kind(0.0d0)) :: BET1
    REAL(kind(0.0d0)), ALLOCATABLE :: SD(:)
    REAL(kind(0.0d0)), ALLOCATABLE :: S(:),fsc(:)
    REAL(kind(0.0d0)), external :: DDOT
    REAL(kind(0.0d0)), external :: DNRM2
    INTEGER , parameter :: IONE=1
    !
    mritr=nwrkcum+4*N
    ALLOCATE (S(2*N+1:nwrkcum+5*N),SD(N))
    !
    flop=0.0d0
    kstat=0
    IF (wfo) THEN
       WRITE(iout,940) IRANK
    END IF
    IF (wff) THEN
       IF (  trans) THEN
          WRITE(iout,941)
       END IF
       IF (nlev > 1) THEN
        IF (smoothtp == -1) THEN
          WRITE(iout,943) nsmooth
          WRITE(iout,944) omeg
        ELSEIF (smoothtp == 1) THEN
          WRITE(iout,945) nsmooth
          WRITE(iout,947)
        ELSE
          WRITE(iout,946) nsmooth
          WRITE(iout,947)
        END IF
       END IF
    END IF
    !
    TOL = RESID
    MAXIT = ITER
    RESID = 1.0d0
    ITER = 0
    dum(3) = DNRM2(N, f, IONE)**2
    flop=flop+dble(2*N)
      IF (init.EQ.1) THEN
          CALL dagmg_matv(N,x,SD,a,ja,ia,trans )
          f(1:n)=f(1:n)-SD(1:N)
       dum(2) = DNRM2(N, f, IONE)**2
       BNORM=SQRT(dum(3))
       RESID=SQRT(dum(2))
       RESID0=RESID
       IF (BNORM.EQ.0.0d0) THEN
          !
          !
          IF(wff) THEN
             WRITE(iout,998)
             WRITE(iout,999)
          END IF
          X(1:N)=0.0d0
          RETURN
       END IF
       RESID=RESID/BNORM
       IF(wff.AND. (MAXIT <= 0 .OR. RESID <= TOL)) THEN
          WRITE(iout,900) 0, resid*bnorm, resid
       END IF
    END IF
    DO WHILE ( ITER < MAXIT .AND. RESID > TOL )
       ITER = ITER + 1
       !
       !
       CALL dagmg_prec_matv(1,f,S(1+3*N),S(1+4*N),S(1+5*N),.FALSE.,a,ja,ia)
       IF ( ITER > 1 ) THEN
          dum(1) = - DDOT(N, SD, IONE, S(1+3*N), IONE)
          BET0=dum(1)
          BET1=BET0/RHO
          CALL DAXPY(N, BET1, S(1+2*N), IONE, S(1+3*N), IONE)
          flop=flop+dble(4*N)
       ENDIF
       CALL DCOPY(N, S(1+3*N), IONE, S(1+2*N), IONE)
       CALL dagmg_matv(N,S(1+2*N),SD,a,ja,ia,trans )
       !
       dum(1) =  DDOT(N, S(1+2*N), IONE, SD, IONE)
       dum(2) =  DDOT(N,S(1+2*N),IONE,f,IONE)
       flop=flop+dble(4*N)
       !
       IF (ITER==1) THEN
          IF (init == 0) THEN
             BNORM=SQRT(dum(3))
             RESID0=BNORM
             IF (BNORM.EQ.0.0d0) THEN
                IF(wff) THEN
                   WRITE(iout,998)
                   WRITE(iout,999)
                END IF
                X(1:N)=0.0d0
                ITER=0
                RETURN
             END IF
          END IF
          IF(wff) THEN
             WRITE(iout,900) 0, resid*bnorm, resid
          END IF
       ELSE
       END IF
       !
       RHO=dum(1)
       ALPHA=dum(2)/RHO
       !
       IF (ITER == 1 .AND. init == 0) THEN
          CALL DCOPY(N,S(1+2*N),IONE,X,IONE)
          CALL DSCAL(N,ALPHA,X,IONE)
          flop=flop+dble(N)
       ELSE
          CALL DAXPY(N, ALPHA, S(1+2*N), IONE, X, IONE)
          flop=flop+dble(2*N)
       END IF
       CALL DAXPY(N, -ALPHA, SD, IONE, f, IONE)
       !
       dum0 = DNRM2(N,f,IONE)**2
       RESID=dum0
       RESID=SQRT(RESID)
       RESID=RESID/BNORM
       IF (wff) THEN
          WRITE(iout,900) iter, resid*bnorm, resid
       END IF
       flop=flop+dble(4*N)
    END DO
    !
    IF( resid > tol ) THEN
       IF (woo) THEN
          WRITE(iout,'()')
          WRITE(iout,950) iter
          WRITE(iout,951)
          WRITE(iout,'()')
       END IF
       iter=-iter
    ELSE
       IF (wff) THEN
          WRITE(iout,952) iter
          WRITE(iout,'()')
       END IF
    END IF
    RESID=RESID*BNORM/RESID0
    !
    kstat(2,1)=ABS(iter)
    !
    DEALLOCATE(S,SD)
    IF(ALLOCATED(fsc)) DEALLOCATE(fsc)
    RETURN
900 FORMAT('****  Iter=',i5,'        Resid=',e9.2,                 &
         '        Relat. res.=',e9.2)
940 FORMAT(i3,                                                     &
         '*SOLUTION: flexible conjugate gradient iterations (FCG(1))')
941 FORMAT('****     Rmk: solve system with the transpose of the input matrix')
943 FORMAT('****     AMG preconditioner with',i2,             &
         ' SOR pre- and post-smoothing sweeps')
944 FORMAT('****         at each level  (relaxation factor: omega=',f5.2,')')
945 FORMAT('****     AMG preconditioner with',i2,             &
         ' Gauss-Seidel pre- and post-smoothing sweeps')
946 FORMAT('****     AMG preconditioner with',i2,             &
         ' ILU(0) pre- plus post-smoothing step(s)')
947 FORMAT('****         at each level')
950 FORMAT('**** !!!   WARNING!!!',I5,' ITERATIONS WERE')
951 FORMAT('**** INSUFFICIENT TO ACHIEVE THE DESIRED ACCURACY')
952 FORMAT('****  - Convergence reached in',I5,' iterations -')
998 FORMAT('**** The norm of the right hand side is zero:')
999 FORMAT('****     set x equal to the zero vector and exit')
    !
  END SUBROUTINE dagmg_FlexCG
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_GCR(N,f,X,ITER,RESID,a,ja,ia,init)
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER       :: N, ITER, init
    REAL(kind(0.0d0)) ::  f(N), X(N)
    INTEGER       :: ja(*), ia(N+1)
    REAL(kind(0.0d0)) :: a(*)
    INTEGER   :: MAXIT,i,itm,irst,ierr,j,k
    REAL(kind(0.0d0)) ::  ALPHA, BET0
    REAL(kind(0.0d0)) ::  RESID, RESID2, BNORM2,  RHO, TRS, RESID0
    REAL(kind(0.0d0)) ::  TOL, TOL2BNORM2, TOLT, dum0
    REAL(kind(0.0d0)) :: dum(6),xd,t
    REAL(kind(0.0d0)), ALLOCATABLE :: Sc(:),Siv(:),SiR(:)
    REAL(kind(0.0d0)), ALLOCATABLE :: Su(:),R(:),fsc(:),Sw(:)
    REAL(kind(0.0d0)), external :: DDOT
    REAL(kind(0.0d0)), external :: DNRM2
    INTEGER , parameter :: IONE=1
    INTEGER  :: itm1, m, info
    !
    mritr=nwrkcum+N*2*nrst+((nrst+1)*nrst)/2+nrst
    ALLOCATE (Su(N*nrst),Sc(N*nrst)                  &
             ,SiR(((nrst+1)*nrst)/2),Siv(nrst)       &
             ,R(nwrkcum))
    !
    flop=0.0d0
    kstat=0
    IF (wfo) THEN
       WRITE(iout,938) IRANK,nrst
    END IF
    IF (wff) THEN
       IF (  trans) THEN
          WRITE(iout,941)
       END IF
       IF (nlev > 1) THEN
        IF (smoothtp == -1) THEN
          WRITE(iout,943) nsmooth
          WRITE(iout,944) omeg
        ELSEIF (smoothtp == 1) THEN
          WRITE(iout,945) nsmooth
          WRITE(iout,947)
        ELSE
          WRITE(iout,946) nsmooth
          WRITE(iout,947)
       END IF
      END IF
    END IF
    !
    m=MIN(nrst,ITER)
    TRS=EPSILON(1.0d0)
    TRS=SQRT(SQRT(TRS))
    TOL = RESID
    TOL2BNORM2 = TOL
    MAXIT = ITER
    RESID2= 1.0d0
    ITER = 0
    dum(3) = DNRM2(N, f, IONE)**2
    flop=flop+dble(2*N)
    !
      IF (init==1) THEN
          CALL dagmg_matv(N,x,Sc,a,ja,ia,trans )
          f(1:n)=f(1:n)-Sc(1:N)
       dum(2) = DNRM2(N, f, IONE)**2
       BNORM2=dum(3)
       RESID2=dum(2)
       IF (BNORM2.EQ.0.0d0) THEN
          !
          !
          IF(wff) THEN
             WRITE(iout,998)
             WRITE(iout,999)
          END IF
          X(1:N)=0.0d0
          RETURN
       END IF
       TOL2BNORM2 = TOL*TOL*BNORM2
       IF (wff.AND. (MAXIT <= 0 .OR. RESID2 <= TOL2BNORM2)) THEN
          WRITE(iout,900) 0, 0,SQRT(resid2),SQRT(resid2/bnorm2)
       END IF
       RESID0=SQRT(RESID2)
    END IF
    itm  = -1
    irst = 0
    DO WHILE ( ITER < MAXIT .AND. RESID2 > TOL2BNORM2 )
       itm  = itm  + 1
       ITER = ITER + 1
       !
       IF (itm == m) THEN
          CALL DTPTRS('U','N','U',m,IONE,SiR,Siv,m,info)
          IF (irst == 0 .AND. init == 0) THEN
             CALL DGEMV('N',N,m,1.0d0,Su,       &
                  N,Siv,IONE,0.0d0,X,IONE)
             flop=flop+dble(2*m*N+m*(m+1))
          ELSE
             CALL DGEMV('N',N,m,1.0d0,Su,        &
                  N,Siv,IONE,1.0d0,X,IONE)
             flop=flop+dble((2*m+1)*N+m*(m+1))
          END IF
          itm=0
          irst=irst+1
       END IF
       !
       !
       CALL dagmg_prec_matv(1,f,Su(1+itm*N),Sc(1+itm*N),R,.FALSE.,a,ja,ia)
       CALL dagmg_matv(N, Su(1+itm*N), Sc(1+itm*N)             &
            , a, ja, ia, trans )
       !
       !
       IF (itm > 0) THEN
          DO i=0,itm-1
             dum(1)=DDOT(N,Sc(1+i*N),IONE,Sc(1+itm*N),IONE)
             bet0=dum(1)
             bet0=bet0/SiR(1+i+(i*(i+1))/2)
             SiR(1+i+(itm*(itm+1))/2)=bet0
             CALL DAXPY(N,-bet0,Sc(1+i*N),IONE,Sc(1+itm*N),IONE)
             flop=flop+dble(4*N)
          END DO
       END IF
       !
       !
       dum(1)=DNRM2(N,Sc(1+itm*N),IONE)**2
       dum(2)=DDOT(N, Sc(1+itm*N), IONE, f, IONE)
       IF (ITER == 1) THEN
          IF (init == 0) THEN
             BNORM2=dum(3)
             RESID2=BNORM2
             RESID0=SQRT(BNORM2)
             IF (BNORM2.EQ.0.0d0) THEN
                !
                !
                IF(wff) THEN
                   WRITE(iout,998)
                   WRITE(iout,999)
                END IF
                X(1:N)=0.0d0
                RETURN
             END IF
             TOL2BNORM2=TOL*TOL*BNORM2
          END IF
          IF (wff) THEN
             WRITE(iout,900) 0, 0,SQRT(resid2),SQRT(resid2/bnorm2)
          END IF
          TOLT = MAX(TOL2BNORM2,TRS*RESID2)
       ELSE
       END IF
       rho=dum(1)
       alpha=dum(2)
       SiR(1+itm+(itm*(itm+1))/2)=rho
       bet0=alpha/rho
       Siv(1+itm)=bet0
       !
       CALL DAXPY(N, -bet0, Sc(1+itm*N), IONE, f, IONE)
       flop=flop+dble(6*N)
       !
       RESID2 = RESID2 - alpha*alpha/rho
       IF (RESID2 <= TOLT) THEN
          RESID2 = DNRM2(N,f,IONE)**2
          flop=flop+dble(2*N)
          TOLT = MAX(TOL2BNORM2,TRS*RESID2)
       END IF
       IF (wff)THEN
          WRITE(iout,900) iter,irst,SQRT(ABS(resid2)),SQRT(ABS(resid2/bnorm2))
       END IF
       !
    END DO
    !
    IF (itm >= 0) THEN
       itm1=itm+1
       CALL DTPTRS('U','N','U',itm1, IONE,SiR,Siv,m,info)
       IF (irst == 0 .AND. init == 0) THEN
          CALL DGEMV('N',N,itm1,1.0d0,Su,        &
               N,Siv,IONE,0.0d0,X,IONE)
          flop=flop+dble(2*(itm+1)*N+(itm+1)*(itm+2))
       ELSE
          CALL DGEMV('N',N,itm1,1.0d0,Su,        &
               N,Siv,IONE,1.0d0,X,IONE)
          flop=flop+dble((2*(itm+1)+1)*N+(itm+1)*(itm+2))
       END IF
    END IF
    !
    RESID=SQRT(RESID2/BNORM2)
    IF( resid > tol ) THEN
       IF (woo) THEN
          WRITE(iout,'()')
          WRITE(iout,950) iter
          WRITE(iout,951)
          WRITE(iout,'()')
       END IF
       iter=-iter
    ELSE
       IF (wff) THEN
          WRITE(iout,952) iter
          WRITE(iout,'()')
       END IF
    END IF
    RESID=RESID*SQRT(BNORM2)/RESID0
    !
    DEALLOCATE (Su,Sc,R,Siv,SiR)
    IF(ALLOCATED(fsc)) DEALLOCATE(fsc)
    IF(ALLOCATED(Sw)) DEALLOCATE(Sw)
    !
    kstat(2,1)=ABS(iter)
    !
    RETURN
900 FORMAT('****  Iter=',i5,' (',i2,' rest.)        Resid=',e9.2,    &
         '        Relat. res.=',e9.2)
938 FORMAT(i3,'*SOLUTION: GCR iterations (GCR(',i2,'))')
941 FORMAT('****     Rmk: solve system with the transpose of the input matrix')
943 FORMAT('****     AMG preconditioner with',i2,             &
         ' SOR pre- and post-smoothing sweeps')
944 FORMAT('****         at each level  (relaxation factor: omega=',f5.2,')')
945 FORMAT('****     AMG preconditioner with',i2,             &
         ' Gauss-Seidel pre- and post-smoothing sweeps')
946 FORMAT('****     AMG preconditioner with',i2,             &
         ' ILU(0) pre- plus post-smoothing step(s)')
947 FORMAT(  '****         at each level')
950 FORMAT('**** !!!   WARNING!!!',I5,' ITERATIONS WERE')
951 FORMAT('**** INSUFFICIENT TO ACHIEVE THE DESIRED ACCURACY')
952 FORMAT('****  - Convergence reached in',I5,' iterations -')
998 FORMAT('**** The norm of the right hand side is zero:')
999 FORMAT('****     set x equal to the zero vector and exit')
    !
  END SUBROUTINE dagmg_GCR
!--------------------------------------------------------------------
  SUBROUTINE dagmg_matv(n, x, y, a, ja, ia, transpose,         &
       iext, lstout, ilstout, lstin, ilstin)
!
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: n, ja(*), ia(n+1), i, kk, k1, k2, ier, idum(1)
    INTEGER, OPTIONAL :: iext(n), lstout(*), ilstout(*), lstin(0:*), ilstin(*)
    REAL(kind(0.0d0)) :: y(n), a(*), t, xx
    REAL(kind(0.0d0)) :: x(n)
    LOGICAL :: transpose
    !
    !
    IF (.NOT.transpose) THEN
       DO i = 1, n
          k1 = ia(i)
          xx = x(ja(k1))
          t = a(k1) * xx
          k2 = ia(i+1)-1
          DO kk = k1+1, k2
             xx = x(ja(kk))
             t = t + a(kk)*xx
          ENDDO
          y(i) = t
       ENDDO
    ELSE
       y(1:n)=0.0d0
       DO i = 1, n
          xx = x(i)
          DO kk = ia(i), ia(i+1)-1
             y(ja(kk)) = y(ja(kk)) + a(kk)*xx
          ENDDO
       ENDDO
    END IF
    !
    flop = flop + dble(2 *(ia(n+1)-1)-n)
    !
    RETURN
  END SUBROUTINE dagmg_matv
!-----------------------------------------------------------------------
  RECURSIVE SUBROUTINE dagmg_prec_matv(l,B,X,AX,R,matv,a,ja,ia)
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l
    REAL(kind(0.0d0)), OPTIONAL ::  B(*), X(*), AX(*), R(*), A(*)
    LOGICAL, OPTIONAL :: matv
    INTEGER, OPTIONAL :: ja(*), ia(*)
    REAL(kind(0.0d0)) ::  dum(1)
    LOGICAL :: update
    INTEGER :: is,n,nnext,XN,BN,AXN,RN,idum(1)
    INTEGER, PARAMETER :: smoothilu=max(smoothtype,1)-1
    INTEGER, PARAMETER :: npostsmooth=nsmooth+smoothilu*mod(nsmooth,2)
    INTEGER, PARAMETER ::  npresmooth=nsmooth-smoothilu*mod(nsmooth,2)
    n=nn(l)
    nnext=nn(l+1)
    IF (l == nlev) THEN
       !
       !
       X(1:N)=B(1:N)
       CALL dagmg_DIRseq(n,X,2,a,ja,ia)
       !
       RETURN
       !
    END IF
    !
    !
    !
    IF (npresmooth == 0) THEN
       !
       X(1:n)=0.0d0
       !
    ELSEIF (npresmooth == 1) THEN
       !
       CALL dagmg_fwGS(l,B,X,AX,.FALSE.,R)
       !
    ELSE
       !
       update=.FALSE.
       DO is=2,npresmooth,2
          CALL dagmg_fwbwsolve1(l,B,X,AX,update,R)
          update=.TRUE.
       END DO
       IF (mod(npresmooth,2) .EQ. 1) THEN
          !
          CALL dagmg_fwGS(l,B,X,AX,.TRUE.,R)
          !
       END IF
    END IF
    !
    !
    IF (nnext > 0) THEN
       !
       !
       XN=1
       BN=XN+nnext
       IF (l+1 == nlev) BN=XN
       IF (npresmooth == 0) THEN
         CALL dagmg_restag(N,nnext,B,R(BN),dt(l)%ind,flop)
       ELSE
         CALL dagmg_restag(N,nnext,AX,R(BN),dt(l)%ind,flop)
       END IF
       !
       IF (l+1 == nlev) THEN
          !
          CALL dagmg_DIRseq(nnext,R,2)
          !
       ELSE IF ( innermax(l+1) <= 1 ) THEN
          !
          AXN=BN+nnext
          RN=AXN+nnext
          CALL dagmg_prec_matv(l+1,R(BN),R(XN),R(AXN),R(RN),.FALSE.)
          !
          kstat(1,l+1)=MAX(kstat(1,l+1),1)
          kstat(2,l+1)=kstat(2,l+1)+1
          !
       ELSE
          !
          CALL dagmg_inner_iter(nnext,R(XN),R(BN),l+1)
          !
       END IF
       !
       CALL dagmg_prolag(N,nnext,X,R(XN),dt(l)%ind,flop)
       !
    END IF
    !
    !
    IF (npostsmooth == 1) THEN
       !
       CALL dagmg_bwGS(l,B,X,AX,matv)
       !
    ELSE
       IF (mod(npostsmooth,2) .EQ. 1) THEN
          !
          CALL dagmg_bwGS(l,B,X,AX,.FALSE.)
          !
       END IF
       !
       update=.FALSE.
       DO is=2,npostsmooth,2
          IF (is .GE. npostsmooth-1) update=matv
          CALL dagmg_fwbwsolve2(l,B,X,AX,update,R)
       END DO
    END IF
    RETURN
  END SUBROUTINE dagmg_prec_matv
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_fwGS(l,B,X,AX,update,R)
!
!
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l, n, idum(1), ier
    REAL(kind(0.0d0)) ::  B(*), X(*), AX(*), R(*)
    LOGICAL :: update
    n=nn(l)
    !
    IF (update) THEN
       !
       !
       IF (smoothtp==1) THEN
         CALL dagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%a   &
                              ,R,AX,1,flop)
       ELSE
         CALL dagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%p   &
                              ,R,AX,1,flop)
       END IF
       X(1:n)=X(1:n)+R(1:n)
       flop=flop+dble(n)
       !
       !
       !
       IF (smoothtp == 1) THEN
         CALL  dagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,R,1,AX,-1,flop)
       ELSE
         CALL  dagmg_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,R,1,AX,-1,flop)
       ENDIF
       !
    ELSE
       !
       !
       IF (smoothtp==1) THEN
         CALL dagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%a,  &
                               X,B,1,flop)
       ELSE
         CALL dagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%p,  &
                               X,B,1,flop)
       END IF
       !
       IF (smoothtp == 1) THEN
         CALL  dagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,X,1,AX,-1,flop)
       ELSE
         CALL  dagmg_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,X,1,AX,-1,flop)
       ENDIF
       !
    END IF
    RETURN
  END SUBROUTINE dagmg_fwGS
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_bwGS(l,B,X,AX,matv)
!
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l, n, idum(1), ier
    REAL(kind(0.0d0)) ::  B(*), X(*), AX(*)
    LOGICAL :: matv
    n=nn(l)
    !
    !
    AX(1:n)=B(1:n)
    IF (smoothtp==1) THEN
      CALL  dagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,-2,flop)
    ELSE
      CALL  dagmg_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,-2,flop)
    END IF
    !
    !
    IF (smoothtp==1) THEN
      CALL dagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%a,X,AX,1,flop)
    ELSE
      CALL dagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%p,X,AX,1,flop)
    END IF
    !
    IF (.NOT.matv) RETURN
    !
    !
    IF (smoothtp==1) THEN
      CALL  dagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,2,flop)
    ELSE
      CALL  dagmg_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,2,flop)
    END IF
    RETURN
  END SUBROUTINE dagmg_bwGS
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_fwbwsolve1(l,B,X,AX,update,R)
!
!
!
!
!
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l, n, idum(1), ier
    REAL(kind(0.0d0)) ::  B(*), X(*), AX(*), R(*)
    LOGICAL :: update
    !
    n=nn(l)
    IF (.NOT.update) THEN
       !
       !
       IF (smoothtp == 1) THEN
         CALL dagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%a    &
                              ,R,B,-1,flop)
       ELSE
         CALL dagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%p    &
                              ,R,B,-1,flop)
       ENDIF
       !
       !
       IF (smoothtp == 1) THEN
         CALL dagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%a    &
                              ,X,R,1,flop)
       ELSE
         CALL dagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%p    &
                              ,X,R,1,flop)
       ENDIF
       !
       !
       AX(1:n)=B(1:n)-R(1:n)
       flop=flop+dble(n)
       !
       !
       IF (smoothtp == 1) THEN
         CALL  dagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,-2,flop)
       ELSE
         CALL  dagmg_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,-2,flop)
       ENDIF
    ELSE
       !
       !
       IF (smoothtp == 1) THEN
         CALL dagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%a    &
                              ,R,AX,-1,flop)
       ELSE
         CALL dagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%p    &
                              ,R,AX,-1,flop)
       ENDIF
       !
       !
       IF (smoothtp == 1) THEN
         CALL dagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%a    &
                              ,R(n+1),R,1,flop)
       ELSE
         CALL dagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%p    &
                              ,R(n+1),R,1,flop)
       ENDIF
       !
       !
       X(1:n)=X(1:n)+R(n+1:2*n)
       !
       !
       AX(1:n)=AX(1:n)-R(1:n)
       flop=flop+dble(2*n)
       !
       !
       IF (smoothtp == 1) THEN
         CALL  dagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,R(n+1),1,AX  &
                           ,-2,flop)
       ELSE
         CALL  dagmg_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,R(n+1),1,AX  &
                           ,-2,flop)
       ENDIF
       !
    END IF
    !
    RETURN
  END SUBROUTINE dagmg_fwbwsolve1
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_fwbwsolve2(l,B,X,AX,matv,R)
!
!
!
!
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l, n, idum(1), ier
    REAL(kind(0.0d0)) ::  B(*), X(*), AX(*), R(*)
    LOGICAL :: matv
    !
    n=nn(l)
    !
    IF (smoothtp == 1) THEN
      CALL  dagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,X,1,AX,0,flop)
    ELSE
      CALL  dagmg_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,X,1,AX,0,flop)
    ENDIF
    !
    R(1:n)=B(1:n)-AX(1:n)
    !
    !
    IF (smoothtp == 1) THEN
      CALL dagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%a    &
                           ,R(n+1),R,1,flop)
    ELSE
      CALL dagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%p    &
                           ,R(n+1),R,1,flop)
    ENDIF
    IF (matv) THEN
      IF (smoothtp == 1) THEN
         AX(1:n)=AX(1:n)+R(n+1:2*n)/dt(l)%a(1:n)
      ELSE
         AX(1:n)=AX(1:n)+R(n+1:2*n)/dt(l)%p(1:n)
      ENDIF
      flop=flop+dble(2*n)
    END IF
    R(n+1:2*n) = R(n+1:2*n) - X(1:n)
    IF (smoothtp == 1) THEN
      CALL dagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%a    &
                           ,R,R(n+1),-1,flop)
    ELSE
      CALL dagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%p    &
                           ,R,R(n+1),-1,flop)
    ENDIF
    X(1:n)=X(1:n) + R(1:n)
    flop=flop+dble(3*n)
    !
    IF (.NOT.matv) RETURN
    IF (smoothtp == 1) THEN
      CALL  dagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,2,flop)
    ELSE
      CALL  dagmg_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,2,flop)
    ENDIF
    RETURN
  END SUBROUTINE dagmg_fwbwsolve2
!-----------------------------------------------------------------------
  RECURSIVE SUBROUTINE dagmg_inner_iter(n,X,R,l)
    !
    !
    !
    !
    !
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER   :: N, ITER, l, ierr, k
    REAL(kind(0.0d0)) ::  RESID, BNORM, det
    REAL(kind(0.0d0)) :: X(N), R(N,*)
    REAL(kind(0.0d0)) :: alpha1,alpha2,bet0,bet1,bet2,rho1,rho2,gamm0,gamm1,zeta
    REAL(kind(0.0d0)) :: dum(10),y1,y2
    REAL(kind(0.0d0)), external :: DDOT
    REAL(kind(0.0d0)), external :: DNRM2
    INTEGER , parameter :: IONE=1
    !
    !
    dum(3)=DNRM2(N, R, IONE)**2
    !
    !
    IF (dum(3) .EQ. 0.0d0) THEN
       X(1:N)=0.0d0
       flop=flop+dble(2*N)
       RETURN
    END IF
    ITER = 1
    !
    !
    CALL dagmg_prec_matv(l,R,X,R(1,2),R(1,3),.TRUE.)
    !
    !
    dum(1) = DDOT(N,X,IONE,R(1,2),IONE)
    dum(2) = DDOT(N,X,IONE,R,IONE)
    IF (ABS(dum(1)) .LE. repsmach*ABS(dum(2))) THEN
       flop=flop+dble(2*N)
       GOTO 100
    END IF
    BNORM = dum(3)
    rho1=dum(1)
    alpha1=dum(2)
    !
    bet0=alpha1/rho1
    !
    !
    CALL DAXPY(N, -bet0, R(1,2), IONE, R, IONE)
    RESID = DNRM2(N,R,IONE)**2
    IF (RESID <= resi*resi*BNORM) THEN
       CALL DSCAL( N, bet0, X, IONE )
       flop=flop+dble(11*N)
       kstat(1,l)=MAX(kstat(1,l),iter)
       kstat(2,l)=kstat(2,l)+iter
       RETURN
    END IF
    !
    ITER = 2
    !
    CALL dagmg_prec_matv(l,R,R(1,3),R(1,4),R(1,5),.TRUE.)
    !
    !
    !
    dum(1) = DDOT(N,R(1,3),IONE,R(1,2),IONE)
    dum(2) = DDOT(N,R(1,3),IONE,R,IONE)
    dum(3) = DDOT(N,R(1,3),IONE,R(1,4),IONE)
    IF (.NOT.spd) dum(4) = DDOT(N,X,IONE,R(1,4),IONE)
    gamm0 =dum(1)
    alpha2=dum(2)
    rho2  =dum(3)
    IF (spd) THEN
       gamm1 = gamm0
    ELSE
       gamm1 = dum(4)
    END IF
    !
    bet1=rho2-gamm0*gamm1/rho1
    bet2=(alpha1-alpha2*gamm1/bet1)/rho1
    zeta=alpha2/bet1
    !
    IF (ABS(bet1).LE.repsmach*ABS(alpha2) .OR.   &
        ABS(bet2)*repsmach.GE.1.0d0)            THEN
       flop=flop+dble(6*N)
       IF (.NOT.spd) flop=flop+dble(2*N)
      GOTO 200
    END IF
    CALL DSCAL(N, bet2, X, IONE)
    CALL DAXPY(N, zeta, R(1,3), IONE, X, IONE)
    !
    !
    kstat(1,l)=MAX(kstat(1,l),iter)
    kstat(2,l)=kstat(2,l)+iter
    !
    flop=flop+dble(19*N)
    IF (.NOT.spd) flop=flop+dble(2*N)
    !
    RETURN
    !
100 CONTINUE
    dum(1)=DNRM2(N,R(1,2),IONE)**2
    dum(2)=DDOT(N, R(1,2), IONE, R, IONE )
    BNORM = dum(3)
110 CONTINUE
    rho1=dum(1)
    alpha1=dum(2)
    !
    bet0=alpha1/rho1
    !
    RESID=BNORM-alpha1*bet0
    IF (RESID <= resi*resi*BNORM) THEN
       CALL DSCAL( N, bet0, X, IONE )
       flop=flop+dble(7*N)
       kstat(1,l)=MAX(kstat(1,l),iter)
       kstat(2,l)=kstat(2,l)+iter
       RETURN
    END IF
    !
    CALL DAXPY(N, -bet0, R(1,2), IONE, R, IONE)
    !
    ITER = 2
    !
    !
    CALL dagmg_prec_matv(l,R,R(1,3),R(1,4),R(1,5),.TRUE.)
    !
    !
    dum(1) = DDOT(N,R(1,4),IONE,R(1,2),IONE)
    dum(2) = DDOT(N,R(1,4),IONE,R,IONE)
    dum(3) = DNRM2(N,R(1,4),IONE)**2
    gamm0 =dum(1)
    alpha2=dum(2)
    rho2  =dum(3)
    gamm1 = gamm0
    !
    bet1=rho2-gamm0*gamm1/rho1
    bet2=(alpha1-alpha2*gamm1/bet1)/rho1
    zeta=alpha2/bet1
    !
    CALL DSCAL(N, bet2, X, IONE)
    CALL DAXPY(N, zeta, R(1,3), IONE, X, IONE)
    !
    !
    kstat(1,l)=MAX(kstat(1,l),iter)
    kstat(2,l)=kstat(2,l)+iter
    !
    flop=flop+dble(19*N)
    !
    RETURN
    !
200 CONTINUE
    !
    !
    !
    !
    !
    !
    !
    !
    dum(1) = DNRM2(N,R(1,2),IONE)**2
    dum(2) = DDOT(N,R(1,4),IONE,R(1,2),IONE)
    dum(3) = DNRM2(N,R(1,4),IONE)**2
    dum(4) = DDOT(N,R(1,2),IONE,R,IONE)
    dum(5) = DDOT(N,R(1,4),IONE,R,IONE)
     !
     dum(4) = dum(4)+bet0*dum(1)
     dum(5) = dum(5)+bet0*dum(2)
     det = dum(1)*dum(3)-dum(2)*dum(2)
     y1  = (dum(3)*dum(4)-dum(2)*dum(5))/det
     y2  = (-dum(2)*dum(4)+dum(1)*dum(5))/det
     CALL DSCAL( N, y1, X, IONE )
     CALL DAXPY( N, y2, R(1,3), IONE, X, IONE )
    kstat(1,l)=MAX(kstat(1,l),iter)
    kstat(2,l)=kstat(2,l)+iter
    !
    flop=flop+dble(23*N)
    !
    RETURN
    !
  END SUBROUTINE dagmg_inner_iter
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_prolag(n, nc, V, B, ind, flop)
!
!
    IMPLICIT NONE
    INTEGER :: n, nc, ind(n), k, i
    REAL(kind(0.0d0)) :: V(n), B(nc)
    REAL(kind(0.0d0)) :: flop
    !
    DO i = 1, n
       k = ind(i)
       IF (k.GT.0) THEN
          V(i) = V(i)+B(k)
       ENDIF
    ENDDO
    flop = flop + dble(n)
    RETURN
  END SUBROUTINE dagmg_prolag
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_restag(n, nc, V, B, ind, flop)
!
!
    IMPLICIT NONE
    INTEGER :: n, nc, ind(n), k, i
    REAL(kind(0.0d0)) :: V(n), B(nc)
    REAL(kind(0.0d0)) :: flop
    !
    B(1:nc)=0.0d0
    !
    DO i = 1, n
       k = ind(i)
       IF (k.GT.0) B(k) = B(k) + V(i)
    ENDDO
    flop = flop + dble(n)
    RETURN
  END SUBROUTINE dagmg_restag
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
  RECURSIVE SUBROUTINE dagmg_setupL1(n,a,ja,ia,listrank,ifl)
!
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: n
    INTEGER :: ja(*),ia(n+1)
    REAL(kind(0.0d0)) :: a(*)
    INTEGER :: ifl,listrank(ifl:*)
    INTEGER :: nc,ierr,i,j,k,nz
    LOGICAL :: slcoarse
    INTEGER, POINTER, DIMENSION(:) :: jap
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: ap
    REAL(kind(0.0d0)) :: eta,dum(3)
    CHARACTER(len=13) :: prtint
    REAL (kind(0.0d0)) :: fff(1)
    nn(1)=n
    nlc(1)=n
    nlc(2)=ia(n+1)-ia(1)
    ngl=nlc
    maxcoarset=FLOOR(maxcoarsesize*(n**(1.0d0/3)))
    maxcoarseslowt=FLOOR(maxcoarsesizeslow*(n**(1.0d0/3)))
    IF ( 0 == nstep  .OR. 1 == maxlev .OR. n <= maxcoarset ) nlev=1
       !
       wcplex=1.0d0
       nlctot=nlc
       ngltot=ngl
       ngl1=ngl
       nlc1=nlc
       icum=1
       fracnz(1)=1.0d0
       allzero=.FALSE.
       nglp=0.0d0
       !
       IF (wfo) THEN
          IF (n > 9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') dble(n)
          ELSE
             WRITE(prtint(1:12),'(i12)') n
          END IF
          WRITE(iout,'()')
          WRITE(iout,918) prtint(1:12)
          IF (nlc(2) > 9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') dble(nlc(2))
          ELSE
             WRITE(prtint(1:12),'(i12)') nlc(2)
          END IF
          WRITE(iout,919) prtint(1:12),dble(nlc(2))/dble(n)
          WRITE(iout,'()')
       END IF
    nlcp=nlc
    nglp=ngl
    !
    IF (1 /= nlev) THEN
       !
       !
       CALL dagmg_aggregation(1,n,a,ja,ia,nc)
       !
       !
       CALL dagmg_setup(2,nc)
       !
    ELSE
       !
       !
       DEALLOCATE(dt(1)%idiag)
       memi=memi-(n+1)
       CALL dagmg_DIRseq(n,fff,1,a,ja,ia)
       IF (wfo) THEN
          WRITE(iout,911) IRANK,flop/dble(2*nlc1(2))
       END IF
       wcplex(4)=0.0d0
       nwrkcum=1
       IF (ngl(1) > maxcoarset)   wcplex(4)=-1.0d0
       IF (ngl(1) > maxcoarseslowt) wcplex(4)=ngl(1)/(1000*ngl1(1)**(1.0d0/3))
       !
       !
       RETURN
    END IF
    !
    !
       !
       nz=ia(n+1)-ia(1)
       IF (transint) THEN
          ALLOCATE(ap(nz),jap(nz-n),dt(1)%il(n+1),dt(1)%iu(n+1))
          memi=memi+nz+n+2
          memr=memr+nz
          memax=MAX(memax,memr+memi*rlenilen)
          CALL dagmg_csrdluT(n,a,ja,ia,dt(1)%idiag              &
               ,ap,jap,dt(1)%il,dt(1)%iu)
          DEALLOCATE(dt(1)%idiag)
          memi=memi-(n+1)
       ELSE
          ALLOCATE(ap(nz),jap(nz-n),dt(1)%iu(n+1))
          dt(1)%il => dt(1)%idiag
          memi=memi+nz+1
          memr=memr+nz
          memax=MAX(memax,memr+memi*rlenilen)
       !
          CALL dagmg_csrdlu(n,a,ja,ia,dt(1)%idiag              &
               ,ap,jap,dt(1)%iu )
          NULLIFY(dt(1)%idiag)
       END IF
       dt(1)%a  => ap
       dt(1)%ja => jap
       NULLIFY(ap,jap)
       !
       !
       innermax(nlev)=0
       innermax(1)=1
       eta=xsi/((1-xsi)*(cplxmax-1))
       icum=1
       DO i=2,nlev-1
          innermax(i)=min(2,floor(xsi**(i-1)/(eta*fracnz(i)*icum)))
          IF (nlev-i.LE.nlvcyc .AND. i.GT.2) innermax(i)=1
          icum=icum*innermax(i)
          wcplex(2)=wcplex(2)+icum*fracnz(i)
          wcplex(3)=wcplex(3)+(2**(i-1))*fracnz(i)
       END DO
       wcplex(2)=wcplex(2)+(2**(nlev-1))*fracnz(nlev)
       wcplex(3)=wcplex(3)+(2**(nlev-1))*fracnz(nlev)
       IF (nsmooth.GT.1 .OR. smoothtype.GT.1)  THEN
          nwrkcum=2*nn(nlev-1)
       ELSE
          nwrkcum=nn(nlev)
       END IF
       nwrkcum=max(nwrkcum,1)
       DO i=nlev-2,1,-1
          nwrkcum=nwrkcum+3*nn(i+1)
          IF (innermax(i+1) > 1) nwrkcum=nwrkcum+2*nn(i+1)
          IF (nsmooth.GT.1 .OR. smoothtype.GT.1)  nwrkcum=max(2*nn(i),nwrkcum)
       END DO
       IF (wfo) THEN
          WRITE(iout,'()')
          WRITE(iout,954) nlctot(1)/dble(nlc1(1))
          WRITE(iout,955) nlctot(2)/dble(nlc1(2))
       END IF
       IF (wff) THEN
          WRITE(iout,956) wcplex(3)
          WRITE(iout,957) wcplex(2)
          WRITE(iout,'()')
       END IF
    RETURN
911 FORMAT(i3,'*','        Exact factorization:',f12.3,' work units (*)')
918 FORMAT('****','        Number of unknowns:', A12)
919 FORMAT('****','                 Nonzeros :', A12,                &
         ' (per row:',f7.2,')')
920 FORMAT('****','        Number of variables:',A12,              &
         '          (reduction ratio:',f5.2,')')
921 FORMAT('****','                   Nonzeros:',A12,              &
         ' (per row:',f4.1,    &
         '; red. ratio:',f5.2,')')
954 FORMAT('****','                  Grid complexity:',f9.2)
955 FORMAT('****','              Operator complexity:',f9.2)
956 FORMAT('****','  Theoretical Weighted complexity:',f9.2, &
                ' (K-cycle at each level)' )
957 FORMAT('****','    Effective Weighted complexity:',f9.2, &
                ' (V-cycle enforced where needed)' )
  END SUBROUTINE dagmg_setupL1
!-----------------------------------------------------------------------
  RECURSIVE SUBROUTINE dagmg_setup(l,n,listrank)
!
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l,n
    INTEGER, OPTIONAL :: listrank(n+1:*)
    INTEGER :: nc,ierr,i,j,k,nz
    LOGICAL :: slcoarse
    INTEGER, POINTER, DIMENSION(:) :: jap
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: ap
    LOGICAL, SAVE :: slowcoarse
    REAL(kind(0.0d0)) :: fw,eta,dum(3)
    CHARACTER(len=13) :: prtint
    REAL (kind(0.0d0)) :: fff(1)
    !
    nn(l)=n
    nlc(1)=n
    IF (n > 0) THEN
       nlc(2)=dt(l)%ia(n+1)-dt(l)%ia(1)
    ELSE
       nlc(2)=0
    END IF
    ngl=nlc
    IF (l==2) slowcoarse=.FALSE.
    slcoarse = 2*nlcp(1) < 3*nlc(1) .AND. 2*nlcp(2) < 3*nlc(2)
    IF( l == nstep+1  .OR. l == maxlev                        &
         .OR. ( nlc(1) <= maxcoarset)                         &
         .OR. ( nlcp(1) < 2*nlc(1) .AND. nlcp(2) < 2*nlc(2)   &
                            .AND. nlc(1) <= maxcoarseslowt )  &
         .OR. ( slowcoarse .AND. slcoarse )                   &
         .OR. nlcp(1) == nlc(1) )                       THEN
         nlev=l
    END IF
    slowcoarse=slcoarse
    fracnz(l)=ngl(2)/ngl1(2)
    IF (l==2) wcplex(1)=ngl1(2)/ngl(2)
       !
       nlctot=nlctot+nlc
       ngltot=ngltot+ngl
       IF (wfo) THEN
          WRITE(iout,'()')
          WRITE(iout,914) IRANK,l
          IF (n > 0) THEN
           IF (n > 9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') dble(n)
           ELSE
             WRITE(prtint(1:12),'(i12)') n
           END IF
           WRITE(iout,920) prtint(1:12),dble(nlcp(1))/dble(n)
           IF (nlc(2) > 9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') dble(nlc(2))
           ELSE
             WRITE(prtint(1:12),'(i12)') nlc(2)
           END IF
           WRITE(iout,921) prtint(1:12),dble(nlc(2))/dble(n),        &
                dble(nlcp(2))/dble(nlc(2))
          ELSE
           WRITE(iout,924) n
          END IF
       END IF
    nlcp=nlc
    nglp=ngl
    IF (n == 0) THEN
       nc=0
       allzero=.TRUE.
       !
       !
       RETURN
    END IF
    IF (l /= nlev) THEN
       !
       !
       CALL dagmg_aggregation(l,n,dt(l)%a,dt(l)%ja,dt(l)%ia,nc)
       !
       !
       CALL dagmg_setup(l+1,nc)
       !
    ELSE
       !
       !
       DEALLOCATE(dt(l)%idiag)
       memi=memi-(n+1)
       CALL dagmg_DIRseq(n,fff,1,dt(l)%a,dt(l)%ja,dt(l)%ia)
       IF (wfo) THEN
          WRITE(iout,911) IRANK,flop/dble(2*nlc1(2))
          IF (coasing) THEN
             WRITE(iout,912) IRANK
          END IF
       END IF
       wcplex(4)=0.0d0
       IF (ngl(1) > maxcoarset)   wcplex(4)=-1.0d0
       IF (ngl(1) > maxcoarseslowt) wcplex(4)=ngl(1)/(1000*ngl1(1)**(1.0d0/3))
       !
       !
       RETURN
    END IF
    !
    !
       !
       nz=dt(l)%ia(n+1)-dt(l)%ia(1)
       IF (transint) THEN
          ALLOCATE(ap(nz),jap(nz-n),dt(l)%il(n+1),dt(l)%iu(n+1))
          memi=memi+nz+n+2
          memr=memr+nz
          memax=MAX(memax,memr+memi*rlenilen)
          CALL dagmg_csrdluT(n,dt(l)%a,dt(l)%ja,dt(l)%ia,dt(l)%idiag &
               ,ap,jap,dt(l)%il,dt(l)%iu)
          memi=memi-size(dt(l)%ja)
          memr=memr-size(dt(l)%a)
          DEALLOCATE(dt(l)%a,dt(l)%ja)
          DEALLOCATE(dt(l)%idiag,dt(l)%ia)
          memi=memi-2*(n+1)
       ELSE
          ALLOCATE(ap(nz),jap(nz-n))
          dt(l)%iu => dt(l)%ia
          dt(l)%il => dt(l)%idiag
          memi=memi+nz-n
          memr=memr+nz
          memax=MAX(memax,memr+memi*rlenilen)
          CALL dagmg_csrdlu(n,dt(l)%a,dt(l)%ja,dt(l)%ia,dt(l)%idiag &
               ,ap,jap )
          memi=memi-size(dt(l)%ja)
          memr=memr-size(dt(l)%a)
          DEALLOCATE(dt(l)%a,dt(l)%ja)
          NULLIFY(dt(l)%idiag,dt(l)%ia)
       END IF
       dt(l)%a  => ap
       dt(l)%ja => jap
       NULLIFY(ap,jap)
    !
    RETURN
911 FORMAT(i3,'*','        Exact factorization:',f12.3,' work units (*)')
912 FORMAT(i3,'*','        Warning: coarsest grid matrix treated as singular')
914 FORMAT(i3,'*','                      Level:',I12)
918 FORMAT('****','        Number of unknowns:', A12)
919 FORMAT('****','                 Nonzeros :', A12,                &
         ' (per row:',f7.2,')')
920 FORMAT('****','        Number of variables:',A12,              &
         '          (reduction ratio:',f5.2,')')
921 FORMAT('****','                   Nonzeros:',A12,              &
         ' (per row:',f4.1,    &
         '; red. ratio:',f5.2,')')
924 FORMAT('****','        Number of variables:',i12)
  END SUBROUTINE dagmg_setup
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_smoothsetup
    USE dagmg_mem
    IMPLICIT NONE
    REAL(kind(0.0d0)) :: unmominv
    INTEGER :: l,i
    !
    !
    smoothtp=smoothtype
    IF (smoothtype == 0) THEN
       IF (spd) THEN
          smoothtp=1
       ELSE IF (omeg <= 0.97) THEN
          smoothtp=-1
       ELSE
          smoothtp=1
       END IF
    END IF
    IF (smoothtp == -1) THEN
       IF (smoothtype == -1) omeg=omega
       unmominv=1.0d0-1.0d0/omeg
       DO l=1,nlev-1
          IF (nn(l) .GT. 0) THEN
            ALLOCATE(dt(l)%p(nn(l)))
            memr=memr+nn(l)
            DO i=1,nn(l)
               dt(l)%p(i)=omeg/dt(l)%a(i)
               dt(l)%a(i)=unmominv*dt(l)%a(i)
            END DO
          END IF
       END DO
       memax=MAX(memax,memr+memi*rlenilen)
    ELSEIF (smoothtp == 1) THEN
       DO l=1,nlev-1
          IF (nn(l) .GT. 0) THEN
            DO i=1,nn(l)
               dt(l)%a(i)=1.0d0/dt(l)%a(i)
            END DO
          END IF
       END DO
    ELSEIF (smoothtp == 2) THEN
       DO l=1,nlev-1
          IF (nn(l) .GT. 0) THEN
            ALLOCATE(dt(l)%p(nn(l)))
            memr=memr+nn(l)
            CALL dagmg_ilu0(nn(l),dt(l)%a,dt(l)%ja,dt(l)%iu,dt(l)%il,dt(l)%p)
          END IF
       END DO
       memax=MAX(memax,memr+memi*rlenilen)
    END IF
    RETURN
  END SUBROUTINE dagmg_smoothsetup
!------------------------------------------------------------------
  SUBROUTINE dagmg_ilu0(n,a,ja,iu,il,p)
    IMPLICIT NONE
    INTEGER :: n, ja(n+1:*), il(n+1), iu(n+1)
    REAL(kind(0.0d0)) :: a(*), p(n), t
    INTEGER :: i,j,k,l
!
    DO i=1,n
      t=a(i)
      DO j=il(i),il(i+1)-1
         k=ja(j)
         DO l=iu(k),iu(k+1)-1
            IF (ja(l) .EQ. i) THEN
               t=t-a(j)*a(l)*p(k)
               EXIT
            ENDIF
         ENDDO
      ENDDO
      a(i)=a(i)-t
      p(i)=1.0d0/t
    ENDDO
!
    RETURN
  END SUBROUTINE dagmg_ilu0
!------------------------------------------------------------------
  SUBROUTINE dagmg_aggregation(l,n,a,ja,ia,nc,listrank)
!
!
!
!
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l,n,nc
    INTEGER :: ja(*),ia(n+1)
    INTEGER, OPTIONAL :: listrank(*)
    REAL (kind(0.0d0)) :: a(*), dum
    INTEGER :: ier,i,j,k,maxdg,np,kpass,nzc,m1,ndd,nzp,isize,nddp,npass1,nz,i0
    LOGICAL :: skipass
    !
    INTEGER, POINTER, DIMENSION(:) :: jan,ian,idiagn,iextn,ind2,lcg,lcgn
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: an
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: sinn
    !
    INTEGER, POINTER, DIMENSION(:) :: jap,iap,idiagp,iextp,lcg1
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: ap
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: sip
    !
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ldd,iw,iperm,riperm
    REAL(kind(0.0d0)), ALLOCATABLE, DIMENSION(:) :: si1,w
    REAL(kind(0.0d0)), ALLOCATABLE, DIMENSION(:) :: wc
    !
    !
    IF (l .EQ. 1) THEN
       IF (wfo) THEN
          WRITE (iout,901) IRANK
       END IF
       IF (wff) THEN
          IF (spd) THEN
             WRITE (iout,906)
          ELSE IF (  transint) THEN
             WRITE (iout,908)
          END IF
          IF (.not.spd) then
             WRITE (iout,902) 'Jacobi',kaptg_dampJac,checkddJ
          ELSE
             WRITE (iout,902) 'BlockD',kaptg_blocdia,checkddB
          END IF
          IF (checkdd < 0) THEN
             WRITE (iout,904)
          END IF
          WRITE(iout,903) npass,targetcoarsefac
          WRITE (iout,905) trspos
       END IF
    END IF
    !
    IF (l .EQ. 1) THEN
       ALLOCATE(si1(n),ind2(n),iperm(n),riperm(n))
       memi=memi+4*n+1
       memr=memr+n
       CALL dagmg_setCMK(n,ja,ia,dt(l)%idiag,riperm,iperm)
    ELSE
       ALLOCATE(si1(n),ind2(n),iperm(n))
       memi=memi+2*n
       memr=memr+n
       iperm(1:n)=1
    END IF
    memax=MAX(memax,memr+memi*rlenilen)
    !
    call dagmg_prepareagg(n,a,ja,ia,dt(l)%idiag,ind2,iperm,si1,ndd,l )
    !
    IF (ndd .EQ. n) THEN
       nc=0
       nzc=0
       DEALLOCATE(si1,iperm,ind2)
       memi=memi-2*n
       memr=memr-n
       IF (l .EQ. 1) THEN
          DEALLOCATE(riperm)
          memi=memi-n
          IF (smoothtype .EQ. 0) omeg=1.0d0
       END IF
       GOTO 999
    END IF
    IF (smoothtype.EQ.0 .AND. l.EQ.1 .AND. (.NOT.spd)) THEN
        ALLOCATE(w(2*n))
        CALL dagmg_omegaest(n,a,ja,ia,dt(l)%idiag,si1,w )
        DEALLOCATE(w)
    END IF
    !
    ALLOCATE(ldd(ndd),lcg(2*(n-ndd)))
    memi=memi+2*n-ndd
    memax=MAX(memax,memr+memi*rlenilen)
    !
    IF (dble(n) .GT. targetcoarsefac*(n-ndd)) THEN
       skipass=.TRUE.
       npass1=npass+1
    ELSE
       skipass=.FALSE.
       npass1=npass
    END IF
    !
    !
    IF (l > 1) THEN
       IF (spd) THEN
          CALL dagmg_findpairs_SI(l,n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,iperm )
       ELSE
          CALL dagmg_findpairs_GI(l,n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,iperm )
       END IF
       DEALLOCATE(iperm)
       memi=memi-n
    ELSE
       IF (spd) THEN
          CALL dagmg_findpairs_SI1(l,n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,riperm,iperm )
       ELSE
          CALL dagmg_findpairs_GI1(l,n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,riperm,iperm )
       END IF
       DEALLOCATE(iperm,riperm)
       memi=memi-2*n
    END IF
10  CONTINUE
    nz=ia(n+1)-1
    !
    IF (npass1.GT.1) THEN
       isize=nc
    ELSE
       isize=1
    END IF
    ALLOCATE(an(nz-2*(n-nc)+ndd),jan(nz-2*(n-nc)+ndd)             &
         ,ian(nc+1),idiagn(nc+1),sinn(isize),wc(nc),iw(2*nc))
    memi=memi+nz-2*(n-nc)+ndd+2*(nc+1)+2*nc+isize
    memr=memr+nz-2*(n-nc)+ndd+nc
    memax=MAX(memax,memr+memi*rlenilen)
    CALL dagmg_setcg(n,a,ja,ia,dt(l)%idiag,si1,ind2,lcg,nc,an,jan,ian      &
         ,idiagn,sinn,npass1.GT.1,maxdg,iw,wc )
    DEALLOCATE(wc,iw)
    memi=memi-2*nc
    memr=memr-nc
    nzc=ian(nc+1)-1
    IF (dble(nz).GT.targetcoarsefac*nzc .OR. npass1.LE.1) THEN
       DEALLOCATE(si1,sinn,lcg)
       IF(ALLOCATED(ldd)) DEALLOCATE(ldd)
       memi=memi-2*(n-ndd)
       memr=memr-n-isize
       dt(l)%ind  => ind2
       dt(l+1)%a    => an
       dt(l+1)%ja   => jan
       dt(l+1)%ia   => ian
       dt(l+1)%idiag=> idiagn
       NULLIFY(ind2,an,jan,ian,idiagn)
       GOTO 999
    END IF
    !
    DEALLOCATE(ind2)
    memi=memi-n
    !
    lcg1 => lcg
    NULLIFY(lcg)
    m1=1
    !
    !
    !
    DO kpass=2,npass1
       m1=2*m1
       np  = nc
       nzp = nzc
       ap     => an
       jap    => jan
       iap    => ian
       idiagp => idiagn
       sip    => sinn
       NULLIFY(an,jan,ian,idiagn,sinn)
       ALLOCATE(lcg(2*np),ind2(np),w(maxdg),iw(maxdg))
       memi=memi+maxdg+3*np
       memr=memr+maxdg
       memax=MAX(memax,memr+memi*rlenilen)
       !
       !
       ind2(1:np)=-1
       IF (spd) THEN
          CALL dagmg_findpairs_SF(l,np,ap,jap,iap,idiagp,sip  &
               ,ind2,lcg,nc                                  &
               ,m1,lcg1,a,ja,ia,dt(l)%idiag,si1,w,iw )
       ELSE
          CALL dagmg_findpairs_GF(l,np,ap,jap,iap,idiagp,sip     &
               ,ind2,lcg,nc                                     &
               ,m1,lcg1,a,ja,ia,dt(l)%idiag,si1,w,iw )
       END IF
       DEALLOCATE(w,iw)
       memi=memi-maxdg
       memr=memr-maxdg
       IF (kpass.NE.npass1) THEN
          isize=nc
       ELSE
          isize=1
          DEALLOCATE(si1)
          memr=memr-n
       END IF
       !
       !
       ALLOCATE(an(nzp-2*(np-nc)),jan(nzp-2*(np-nc))                 &
            ,ian(nc+1),idiagn(nc+1),sinn(isize),wc(nc),iw(2*nc))
       memi=memi+nzp-2*(np-nc)+2*(nc+1)+2*nc+isize
       memr=memr+nzp-2*(np-nc)+nc
       memax=MAX(memax,memr+memi*rlenilen)
       !
       !
       CALL dagmg_setcg(np,ap,jap,iap,idiagp,sip,ind2,lcg,nc,an     &
            ,jan,ian,idiagn,sinn,kpass.NE.npass1,maxdg,iw,wc )
       memi=memi-SIZE(jap)-2*(np+1)-np-2*nc
       memr=memr-SIZE(ap)-SIZE(sip)-nc
       DEALLOCATE(ap,jap,iap,idiagp,sip,ind2,wc,iw)
       !
       !
       ALLOCATE(lcgn(2*m1*nc))
       memi=memi+2*m1*nc
       memax=MAX(memax,memr+memi*rlenilen)
       CALL dagmg_lcgmix(nc,m1,lcg1,lcg,lcgn)
       memi=memi-SIZE(lcg)-SIZE(lcg1)
       DEALLOCATE(lcg,lcg1)
       lcg1 => lcgn
       NULLIFY(lcgn)
       nzc=ian(nc+1)-1
       IF ( kpass.NE.npass1 .AND. dble(nz).GT.targetcoarsefac*nzc ) THEN
          DEALLOCATE(si1)
          memr=memr-n
          EXIT
       END IF
    END DO
    !
    memr=memr-SIZE(sinn)
    DEALLOCATE(sinn)
    !
    ALLOCATE(dt(l)%ind(n))
    memi=memi+n
    memax=MAX(memax,memr+memi*rlenilen)
    CALL dagmg_setind(nc,ndd,ldd,lcg1,2*m1,dt(l)%ind)
    memi=memi-ndd-SIZE(lcg1)
    DEALLOCATE(lcg1,ldd)
    !
       dt(l+1)%a    => an
       dt(l+1)%ja   => jan
       dt(l+1)%ia   => ian
       dt(l+1)%idiag=> idiagn
       NULLIFY(an,jan,ian,idiagn)
999 CONTINUE
    !
    RETURN
901 FORMAT(i3,'*SETUP: Coarsening by multiple pairwise aggregation')
902 FORMAT('****       Quality threshold (',A6,'):',f6.2, &
         ' ;  Strong diag. dom. trs:',f5.2)
903 FORMAT('****         Maximal number of passes:',i3,     &
         '  ; Target coarsening factor:',f5.2)
904 FORMAT('****           Diag. dom. checked w.r.t. sum of offdiag', &
          ' (no absolute vaues)')
905 FORMAT('****',22x,'Threshold for rows with large pos. offdiag.:',f5.2)
906 FORMAT('****  Rmk: Setup performed assuming the matrix symmetric')
908 FORMAT('****  Rmk: Setup performed for the transpose of the input matrix')
  END SUBROUTINE dagmg_aggregation
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_setCMK(n,ja,ia,idiag,riperm,iperm,iext)
    !
    !
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),riperm(*),iperm(n)
    INTEGER, OPTIONAL :: iext(*)
    LOGICAL :: exc
    INTEGER :: i,j,jj,jk,jd,k,kk,j1,j2,i1,i2,ijs,ijs1,ijs2,dg,mindg,kdim
    INTEGER :: ifirst,ilast,kb,i0
    !
    ifirst=1
    ilast=n
    i2=ifirst
    mindg=n+1
    DO i = ifirst,ilast
       dg=ia(i+1)-ia(i)
       IF (dg .GT. 1) THEN
          iperm(i)=-dg
          IF (dg.LT.mindg) THEN
             mindg=dg
             jj=i
          END IF
       ELSE
          riperm(i2)=i
          iperm(i)=i2
          i2=i2+1
       END IF
    ENDDO
    !
    ijs=ifirst-1
    i1=i2
15  CONTINUE
    !
    IF (i2 .LE. ilast) THEN
      riperm(i2)=jj
      iperm(jj)=i2
    END IF
    !
    DO WHILE (i1.LE.i2 .AND. i2.LT.ilast)
       !
       !
       i=riperm(i1)
       ijs1=i2+1
       !
       j1 =ia(i)
       jd=idiag(i)
       j2 = ia (i+1)-1
       DO kk = j1,jd-1
          j=ja(kk)
          IF (iperm(j) .LT. 0) THEN
             i2=i2+1
             riperm(i2)=j
          END IF
       ENDDO
       DO kk = jd+1,j2
          j=ja(kk)
          IF (iperm(j) .LT. 0) THEN
             i2=i2+1
             riperm(i2)=j
          END IF
       ENDDO
       !
       ijs2=i2
       exc=.TRUE. .AND. ijs2.GT.ijs1
       DO WHILE(exc)
          exc=.FALSE.
          DO kk=ijs1+1,ijs2
             IF( iperm(riperm(kk)) .GT. iperm(riperm(kk-1)) )THEN
                j=riperm(kk)
                riperm(kk)=riperm(kk-1)
                riperm(kk-1)=j
                exc=.TRUE.
             END IF
          END DO
       END DO
       DO kk=ijs1,ijs2
          iperm(riperm(kk))=kk
       END DO
       !
       i1=i1+1
    END DO
    IF (i2 .LT. ilast) THEN
       !
       jj=0
       DO WHILE (jj .EQ. 0)
          ijs=ijs+1
          IF (ijs .GT. ilast) THEN
             mindg=mindg+1
             ijs=ifirst
          END IF
          ijs1=ijs
          IF (iperm(ijs1).LT.0 .AND. ia(ijs1+1)-ia(ijs1).EQ.mindg) &
               jj=ijs1
       END DO
       i2=i2+1
       GOTO 15
    END IF
    !
 !
 !
    RETURN
  END SUBROUTINE dagmg_setCMK
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_prepareagg(n,a,ja,ia,idiag,ind2,iperm,si,ndd,l,iext)
    !
    !
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind2(n),iperm(n)
    INTEGER, OPTIONAL :: iext(*)
    INTEGER :: ndd, l
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) , TARGET :: si(n)
    REAL(kind(0.0d0)) :: checkddl,oda,odm,ods,vald
    INTEGER :: i,j,jj,jk,jd,k,kk,j1,j2,i1,i2,ijs,ijs1,ijs2,dg,kdim,nnegrcs
    INTEGER :: ifirst,ilast,i0,kb
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: odmax,odabs,osi
    !
    IF (.NOT.spd) THEN
       checkddl=checkddJ
    ELSE
       checkddl=checkddB
    END IF
    !
    IF (.NOT.spd) THEN
       !
       !
       IF (checkdd > 0) THEN
          ALLOCATE(odmax(n),odabs(n))
          memi=memi+2*n
          odabs(1:n)=0.0d0
       ELSE
          ALLOCATE(odmax(n))
          memi=memi+n
       END IF
       osi => si(1:n)
       si(1:n)=0.0d0
       odmax(1:n)=0.0d0
       memax=MAX(memax,memr+memi*rlenilen)
       ifirst=1
       ilast=n
        DO i=ilast,ifirst,-1
          j =ia(i)
          jd=idiag(i)
          jj=ia(i+1)-1
          DO k = j,jd-1
            jk=ja(k)
             osi(jk)=osi(jk)+a(k)
             odmax(jk)=max(odmax(jk),a(k))
             IF (checkdd > 0) odabs(jk)=odabs(jk)+abs(a(k))
          ENDDO
          DO k = jd+1,jj
            jk=ja(k)
             osi(jk)=osi(jk)+a(k)
             odmax(jk)=max(odmax(jk),a(k))
             IF (checkdd > 0) odabs(jk)=odabs(jk)+abs(a(k))
          ENDDO
        ENDDO
    END IF
    !
    ndd=0
    nnegrcs=0
    !
    ifirst=1
    ilast=n
     DO i=ifirst,ilast
       j1 =ia(i)
       jd=idiag(i)
       j2 = ia (i+1)-1
       vald = a(jd)
       odm=0.0d0
       oda=0.0d0
       ods=0.0d0
       DO kk = j1,jd-1
          ods=ods+a(kk)
          odm=max(odm,a(kk))
          IF (checkdd > 0) oda=oda+abs(a(kk))
       ENDDO
       DO kk = jd+1,j2
          ods=ods+a(kk)
          odm=max(odm,a(kk))
          IF (checkdd > 0) oda=oda+abs(a(kk))
       ENDDO
       !
       IF (.NOT.spd) THEN
          ods=(osi(i)+ods)/2
          odm=max(odm,odmax(i))
          IF (checkdd > 0) oda=(oda+odabs(i))/2
       END IF
       !
       IF ((vald+ods) .LT. -repsmach*ABS(vald)) nnegrcs=nnegrcs+1
       !
       !
       si(i)=-ods
       IF ( (checkdd.GT.0 .AND. vald.GT.checkddl*oda)     &
            .OR. (checkdd.LT.0 .AND. vald.GT.checkddl*abs(ods)) ) THEN
          !
          ind2(i)=0
          ndd=ndd+1
       ELSE
          ind2(i)=-1
          IF (odm .GT. trspos*vald) iperm(i)=0
       ENDIF
     END DO
    !
    !
    zerors=.FALSE.
    IF (nnegrcs.GT.fracnegrcsum*n) THEN
       zerors=.TRUE.
       ndd=0
       ind2(1:n)=-1
    END IF
    !
    IF (.NOT.spd) THEN
       IF (checkdd > 0) THEN
          DEALLOCATE(odmax,odabs)
          memi=memi-2*n
       ELSE
          DEALLOCATE(odmax)
          memi=memi-n
       END IF
    END IF
    RETURN
  END SUBROUTINE dagmg_prepareagg
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_omegaest(n,a,ja,ia,idiag,si,w,iext)
    !
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n)
    INTEGER, OPTIONAL :: iext(*)
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n),w(2*n),t,tm
    INTEGER :: i,j,j1,j2,jd,ll,l2,kl,kk
    INTEGER :: ifirst,ilast,kb
    w(1:2*n)=0.0d0
    ifirst=1
    ilast=n
     DO i=ilast,ifirst,-1
       j1 =ia(i)
       jd=idiag(i)
       j2=ia(i+1)-1
       DO kk=j1,jd-1
        j=ja(kk)
          l2=ia(j+1)-1
          kl=0
          DO ll=idiag(j)+1,l2
            IF (ja(ll).EQ.i) THEN
              kl=ll
              EXIT
            END IF
          END DO
          IF (kl.EQ.0) THEN
             w(j)=w(j)+ABS(a(kk))
          ELSE
             w(j)=w(j)+ABS(a(kk)-a(kl))-ABS(a(kl))
          END IF
       ENDDO
       t=w(i)
       DO kk=jd+1,j2
          t=t+ABS(a(kk))
       ENDDO
       w(i)=t/(max(a(idiag(i)),si(i)))
     ENDDO
     tm=0.0d0
    ifirst=1
    ilast=n
     DO i=ifirst,ilast
       j1 =ia(i)
       jd=idiag(i)
       j2=ia(i+1)-1
       DO kk=jd+1,j2
        j=ja(kk)
          kl=0
          DO ll=ia(j),idiag(j)-1
            IF (ja(ll).EQ.i) THEN
               kl=ll
               EXIT
            END IF
          END DO
          IF (kl.EQ.0) THEN
             w(n+j)=w(n+j)+ABS(a(kk))*w(i)
          ELSE
             w(n+j)=w(n+j)+(ABS(a(kk)-a(kl))-ABS(a(kl)))*w(i)
          END IF
       ENDDO
       t=w(n+i)
       DO kk=j1,jd-1
          t=t+ABS(a(kk))*w(ja(kk))
       ENDDO
       t=t/(max(a(idiag(i)),si(i)))
       tm=max(t,tm)
     ENDDO
    IF (tm .LE. 1.0d0) THEN
      omeg=1.0d0
    ELSE
      omeg=2/(1.0d0+SQRT(tm))
    ENDIF
    RETURN
  END SUBROUTINE dagmg_omegaest
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_findpairs_GF(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,m1,lcg1,a1,ja1,ia1,idiag1,si1,rtent,jtent )
!
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: m1,ja1(*),ia1(*),idiag1(*),jtent(*),lcg1(m1,n)
    REAL(kind(0.0d0)) :: a1(*)
    REAL(kind(0.0d0)) :: si1(*),rtent(*)
!
!
!
!
!
!
!
!
!
!
!
!
!----------------
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_dampJac
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
!
      DO WHILE (nmark.LT.nt)
       isel=ijs
       ijs=ijs+1
       !
       !
       IF (ind(isel) .GE. 0) CYCLE
       !
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       IF (lcg1(2,isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       ntentleft=0
       !
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j) .GE. 0) CYCLE
          IF(lcg1(2,j).EQ.0) CYCLE
          kk=0
          IF (i .LT. idiag(isel)) THEN
             j2=ia(j+1)-1
             DO jk=idiag(j)+1,j2
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ELSE
             DO jk=ia(j),idiag(j)-1
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ENDIF
          vals=-a(i)/2
          IF(kk .NE. 0) vals=vals-a(kk)/2
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
             eta1=2*si(isel)
             eta2=2*si(j)
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
             eta1=2*a(idiag(isel))
             eta2=2*a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !
          IF (sig1 > 0.0d0) THEN
             del1=rsi
          ELSE
             del1=rsi+2*sig1
          END IF
          IF (sig2 > 0.0d0) THEN
             del2=rsj
          ELSE
             del2=rsj+2*sig2
          END IF
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=((eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=(eta1*eta2)/(eta1+eta2)
             valp=vals/valp
          END IF
          IF (valp > kaptg) CYCLE
          !
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          ntentleft=ntentleft+1
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          rtent(ntentleft)=tent
          jtent(ntentleft)=j
          CYCLE
9         CONTINUE
          rtent(ntentleft)=val
          jtent(ntentleft)=ipair
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       IF (ipair .EQ. 0) GOTO 25
20     CONTINUE
       CALL dagmg_checktentagg_GF
       IF (.NOT.acc) THEN
          ipair = 0
          IF (ntentleft .GT.0) THEN
             i=1
             j=1
             DO WHILE (i .LE. ntentleft)
                IF (jtent(j).GT.0) THEN
                   tent=rtent(j)
                   IF (ipair.EQ.0) GOTO 22
                   IF (16*(tent-val).LT.-1) GOTO 22
                   IF (16*(tent-val).LT.1 .AND. j.LT.ipair) GOTO 22
                   GOTO 23
22                 CONTINUE
                   val=tent
                   ipair=jtent(j)
                   ijtent=j
23                 CONTINUE
                   i=i+1
                END IF
                j=j+1
             END DO
             ntentleft=ntentleft-1
             jtent(ijtent)=0
             GOTO 20
          END IF
       END IF
       !
25     CONTINUE
       !
       IF (ipair .EQ. 0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
    RETURN
  CONTAINS
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_checktentagg_GF
!
!
!
!
!
    INTEGER, PARAMETER :: mm=max(2**(npass+1),8)
    REAL(kind(0.0d0)) :: W(mm,mm), sig(mm), AGe(mm), v(mm)
    REAL(kind(0.0d0)) :: alpha, alp, tmp, beta, f1, f2
    INTEGER :: j,jj,k,l,m,info, setdim1, setdim, l2, k2
    INTEGER :: set(mm), l1, wdthT
    REAL(kind(0.0d0)) :: T
    LOGICAL :: exc
!
    IF (m1.eq.2) THEN
       IF (lcg1(2,isel) .LT. 0) THEN
          IF (lcg1(2,ipair) .LT. 0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             setdim=2
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             set(3)=lcg1(2,ipair)
             setdim=3
          END IF
          l1=1
       ELSE
          IF (lcg1(2,ipair) .LT. 0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             setdim=3
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             set(4)=lcg1(2,ipair)
             setdim=4
          END IF
          l1=2
       END IF
    ELSE
       l1=m1
       IF (lcg1(m1,isel).LT.0) l1=-lcg1(m1,isel)
       set(1:l1)=lcg1(1:l1,isel)
       l2=m1
       IF (lcg1(m1,ipair).LT.0) l2=-lcg1(m1,ipair)
       set(l1+1:l1+l2)=lcg1(1:l2,ipair)
       setdim=l1+l2
    END IF
!
    exc=.TRUE.
    DO WHILE(exc)
       exc=.FALSE.
       DO l=2,SetDim
          IF( set(l)<set(l-1) )THEN
             jj=set(l)
             set(l)=set(l-1)
             set(l-1)=jj
             exc=.TRUE.
          END IF
       END DO
    END DO
!
    DO j=1,SetDim
       jj=Set(j)
       sig(j)=si1(jj)
       IF (zerors) THEN
          W(j,j)=sig(j)
          AGe(j)=0.0d0
       ELSE
          W(j,j)=a1(idiag1(jj))
          AGe(j)=W(j,j)-sig(j)
       END IF
       l2=j+1
       DO l=l2,SetDim
          W(j,l)=0.0d0
          W(l,j)=0.0d0
       END DO
       k2=ia1(jj+1)-1
       DO k=idiag1(jj)+1,k2
          DO l=l2,SetDim
             m=Set(l)
             IF(ja1(k)==m)THEN
                alpha=a1(k)/2
                W(j,l)=alpha
                W(l,j)=alpha
                EXIT
             END IF
          END DO
       END DO
       DO k=ia1(jj),idiag1(jj)-1
          DO l=1,j-1
             m=Set(l)
             IF(ja1(k)==m)THEN
                alpha=a1(k)/2
                W(j,l)=W(j,l)+alpha
                W(l,j)=W(j,l)
                EXIT
             END IF
          END DO
       END DO
    END DO
!
    DO j=1,SetDim
       DO k=1,SetDim
          IF (j.ne.k) THEN
             sig(j)=sig(j)+W(j,k)
          END IF
       ENDDO
       IF (sig(j) < 0.0d0)  AGe(j)=AGe(j)+2*sig(j)
   !
       v(j)=W(j,j)
   !
   !
   !
       W(j,j)=umdbndmum1*W(j,j)-abs(sig(j))
       IF (j .eq. 1) THEN
          beta=v(j)
          alp=abs(AGe(j))
       ELSE
          beta=beta+v(j)
          alp=max(alp,abs(AGe(j)))
       END IF
    END DO
!
!
    beta=dbndmum1/beta
    DO j=1,SetDim
       DO k=1,SetDim
          W(j,k)=W(j,k)+beta*v(j)*v(k)
       END DO
    END DO
!
!
    IF (alp.LT.repsmach*beta) THEN
       SetDim1=SetDim-1
    ELSE
       SetDim1=SetDim
    END IF
!
!
    acc=.FALSE.
!
    SELECT CASE (SetDim1)
    CASE (1)
       GOTO 11
    CASE (2)
       GOTO 12
    CASE (3)
       GOTO 13
    CASE (4)
       GOTO 14
    CASE (5)
       GOTO 15
    CASE (6)
       GOTO 16
    CASE (7)
       GOTO 17
    CASE (8)
       GOTO 18
    CASE DEFAULT
       CALL DPOTRF('U',SetDim1,W,mm,info)
       IF (info .NE. 0) RETURN
       GOTO 10
    END SELECT
18  CONTINUE
    IF (W(8,8) .LE. 0.0d0) RETURN
    W(7,7) = W(7,7) - (W(7,8)/W(8,8)) * W(7,8)
    T = W(6,8)/W(8,8)
    W(6,7) = W(6,7) - T * W(7,8)
    W(6,6) = W(6,6) - T * W(6,8)
    T = W(5,8)/W(8,8)
    W(5,7) = W(5,7) - T * W(7,8)
    W(5,6) = W(5,6) - T * W(6,8)
    W(5,5) = W(5,5) - T * W(5,8)
    T = W(4,8)/W(8,8)
    W(4,7) = W(4,7) - T * W(7,8)
    W(4,6) = W(4,6) - T * W(6,8)
    W(4,5) = W(4,5) - T * W(5,8)
    W(4,4) = W(4,4) - T * W(4,8)
    T = W(3,8)/W(8,8)
    W(3,7) = W(3,7) - T * W(7,8)
    W(3,6) = W(3,6) - T * W(6,8)
    W(3,5) = W(3,5) - T * W(5,8)
    W(3,4) = W(3,4) - T * W(4,8)
    W(3,3) = W(3,3) - T * W(3,8)
    T = W(2,8)/W(8,8)
    W(2,7) = W(2,7) - T * W(7,8)
    W(2,6) = W(2,6) - T * W(6,8)
    W(2,5) = W(2,5) - T * W(5,8)
    W(2,4) = W(2,4) - T * W(4,8)
    W(2,3) = W(2,3) - T * W(3,8)
    W(2,2) = W(2,2) - T * W(2,8)
    T = W(1,8)/W(8,8)
    W(1,7) = W(1,7) - T * W(7,8)
    W(1,6) = W(1,6) - T * W(6,8)
    W(1,5) = W(1,5) - T * W(5,8)
    W(1,4) = W(1,4) - T * W(4,8)
    W(1,3) = W(1,3) - T * W(3,8)
    W(1,2) = W(1,2) - T * W(2,8)
    W(1,1) = W(1,1) - T * W(1,8)
17  CONTINUE
    IF (W(7,7) .LE. 0.0d0) RETURN
    W(6,6) = W(6,6) - (W(6,7)/W(7,7)) * W(6,7)
    T = W(5,7)/W(7,7)
    W(5,6) = W(5,6) - T * W(6,7)
    W(5,5) = W(5,5) - T * W(5,7)
    T = W(4,7)/W(7,7)
    W(4,6) = W(4,6) - T * W(6,7)
    W(4,5) = W(4,5) - T * W(5,7)
    W(4,4) = W(4,4) - T * W(4,7)
    T = W(3,7)/W(7,7)
    W(3,6) = W(3,6) - T * W(6,7)
    W(3,5) = W(3,5) - T * W(5,7)
    W(3,4) = W(3,4) - T * W(4,7)
    W(3,3) = W(3,3) - T * W(3,7)
    T = W(2,7)/W(7,7)
    W(2,6) = W(2,6) - T * W(6,7)
    W(2,5) = W(2,5) - T * W(5,7)
    W(2,4) = W(2,4) - T * W(4,7)
    W(2,3) = W(2,3) - T * W(3,7)
    W(2,2) = W(2,2) - T * W(2,7)
    T = W(1,7)/W(7,7)
    W(1,6) = W(1,6) - T * W(6,7)
    W(1,5) = W(1,5) - T * W(5,7)
    W(1,4) = W(1,4) - T * W(4,7)
    W(1,3) = W(1,3) - T * W(3,7)
    W(1,2) = W(1,2) - T * W(2,7)
    W(1,1) = W(1,1) - T * W(1,7)
16  CONTINUE
    IF (W(6,6) .LE. 0.0d0) RETURN
    W(5,5) = W(5,5) - (W(5,6)/W(6,6)) * W(5,6)
    T = W(4,6)/W(6,6)
    W(4,5) = W(4,5) - T * W(5,6)
    W(4,4) = W(4,4) - T * W(4,6)
    T = W(3,6)/W(6,6)
    W(3,5) = W(3,5) - T * W(5,6)
    W(3,4) = W(3,4) - T * W(4,6)
    W(3,3) = W(3,3) - T * W(3,6)
    T = W(2,6)/W(6,6)
    W(2,5) = W(2,5) - T * W(5,6)
    W(2,4) = W(2,4) - T * W(4,6)
    W(2,3) = W(2,3) - T * W(3,6)
    W(2,2) = W(2,2) - T * W(2,6)
    T = W(1,6)/W(6,6)
    W(1,5) = W(1,5) - T * W(5,6)
    W(1,4) = W(1,4) - T * W(4,6)
    W(1,3) = W(1,3) - T * W(3,6)
    W(1,2) = W(1,2) - T * W(2,6)
    W(1,1) = W(1,1) - T * W(1,6)
15  CONTINUE
    IF (W(5,5) .LE. 0.0d0) RETURN
    W(4,4) = W(4,4) - (W(4,5)/W(5,5)) * W(4,5)
    T = W(3,5)/W(5,5)
    W(3,4) = W(3,4) - T * W(4,5)
    W(3,3) = W(3,3) - T * W(3,5)
    T = W(2,5)/W(5,5)
    W(2,4) = W(2,4) - T * W(4,5)
    W(2,3) = W(2,3) - T * W(3,5)
    W(2,2) = W(2,2) - T * W(2,5)
    T = W(1,5)/W(5,5)
    W(1,4) = W(1,4) - T * W(4,5)
    W(1,3) = W(1,3) - T * W(3,5)
    W(1,2) = W(1,2) - T * W(2,5)
    W(1,1) = W(1,1) - T * W(1,5)
14  CONTINUE
    IF (W(4,4) .LE. 0.0d0) RETURN
    W(3,3) = W(3,3) - (W(3,4)/W(4,4)) * W(3,4)
    T = W(2,4)/W(4,4)
    W(2,3) = W(2,3) - T * W(3,4)
    W(2,2) = W(2,2) - T * W(2,4)
    T = W(1,4)/W(4,4)
    W(1,3) = W(1,3) - T * W(3,4)
    W(1,2) = W(1,2) - T * W(2,4)
    W(1,1) = W(1,1) - T * W(1,4)
13  CONTINUE
    IF (W(3,3) .LE. 0.0d0) RETURN
    W(2,2) = W(2,2) - (W(2,3)/W(3,3)) * W(2,3)
    T = W(1,3)/W(3,3)
    W(1,2) = W(1,2) - T * W(2,3)
    W(1,1) = W(1,1) - T * W(1,3)
12  CONTINUE
    IF (W(2,2) .LE. 0.0d0) RETURN
    W(1,1) = W(1,1) - (W(1,2)/W(2,2)) * W(1,2)
11  CONTINUE
    IF (W(1,1) .LE. 0.0d0) RETURN
10  CONTINUE
!
    acc=.TRUE.
!
    RETURN
  END SUBROUTINE dagmg_checktentagg_GF
  END SUBROUTINE dagmg_findpairs_GF
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_findpairs_SF(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,m1,lcg1,a1,ja1,ia1,idiag1,si1,rtent,jtent )
!
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: m1,ja1(*),ia1(*),idiag1(*),jtent(*),lcg1(m1,n)
    REAL(kind(0.0d0)) :: a1(*)
    REAL(kind(0.0d0)) :: si1(*),rtent(*)
!
!
!
!
!
!
!
!
!
!
!
!
!----------------
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_blocdia
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
!
      DO WHILE (nmark.LT.nt)
       isel=ijs
       ijs=ijs+1
       !
       !
       IF (ind(isel) .GE. 0) CYCLE
       !
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       IF (lcg1(2,isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       ntentleft=0
       !
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j) .GE. 0) CYCLE
          IF(lcg1(2,j).EQ.0) CYCLE
          vals=-a(i)
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !
          IF (sig1 > 0.0d0) THEN
             del1=rsi
             eta1=rsi+2*sig1
          ELSE
             del1=rsi+2*sig1
             eta1=rsi
          END IF
          IF (eta1 < 0.0d0) CYCLE
          IF (sig2 > 0.0d0) THEN
             del2=rsj
             eta2=rsj+2*sig2
          ELSE
             del2=rsj+2*sig2
             eta2=rsj
          END IF
          IF (eta2 < 0.0d0) CYCLE
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
               valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=(vals+(eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=vals+(eta1*eta2)/(eta1+eta2)
             IF (vals < 0.0d0) CYCLE
             valp=vals/valp
          END IF
          IF (valp > kaptg) CYCLE
          !
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          ntentleft=ntentleft+1
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          rtent(ntentleft)=tent
          jtent(ntentleft)=j
          CYCLE
9         CONTINUE
          rtent(ntentleft)=val
          jtent(ntentleft)=ipair
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       IF (ipair .EQ. 0) GOTO 25
20     CONTINUE
       CALL dagmg_checktentagg_SF
       IF (.NOT.acc) THEN
          ipair = 0
          IF (ntentleft .GT.0) THEN
             i=1
             j=1
             DO WHILE (i .LE. ntentleft)
                IF (jtent(j).GT.0) THEN
                   tent=rtent(j)
                   IF (ipair.EQ.0) GOTO 22
                   IF (16*(tent-val).LT.-1) GOTO 22
                   IF (16*(tent-val).LT.1 .AND. j.LT.ipair) GOTO 22
                   GOTO 23
22                 CONTINUE
                   val=tent
                   ipair=jtent(j)
                   ijtent=j
23                 CONTINUE
                   i=i+1
                END IF
                j=j+1
             END DO
             ntentleft=ntentleft-1
             jtent(ijtent)=0
             GOTO 20
          END IF
       END IF
       !
25     CONTINUE
       !
       IF (ipair .EQ. 0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
    RETURN
  CONTAINS
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_checktentagg_SF
!
!
!
!
!
    INTEGER, PARAMETER :: mm=max(2**(npass+1),8)
    REAL(kind(0.0d0)) :: W(mm,mm), sig(mm), AGe(mm), v(mm)
    REAL(kind(0.0d0)) :: alpha, alp, tmp, beta, f1, f2
    INTEGER :: j,jj,k,l,m,info, setdim1, setdim, l2, k2
    INTEGER :: set(mm), l1, wdthT
    REAL(kind(0.0d0)) :: T
    LOGICAL :: exc
!
    IF (m1.eq.2) THEN
       IF (lcg1(2,isel) .LT. 0) THEN
          IF (lcg1(2,ipair) .LT. 0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             setdim=2
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             set(3)=lcg1(2,ipair)
             setdim=3
          END IF
          l1=1
       ELSE
          IF (lcg1(2,ipair) .LT. 0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             setdim=3
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             set(4)=lcg1(2,ipair)
             setdim=4
          END IF
          l1=2
       END IF
    ELSE
       l1=m1
       IF (lcg1(m1,isel).LT.0) l1=-lcg1(m1,isel)
       set(1:l1)=lcg1(1:l1,isel)
       l2=m1
       IF (lcg1(m1,ipair).LT.0) l2=-lcg1(m1,ipair)
       set(l1+1:l1+l2)=lcg1(1:l2,ipair)
       setdim=l1+l2
    END IF
!
    exc=.TRUE.
    DO WHILE(exc)
       exc=.FALSE.
       DO l=2,SetDim
          IF( set(l)<set(l-1) )THEN
             jj=set(l)
             set(l)=set(l-1)
             set(l-1)=jj
             exc=.TRUE.
          END IF
       END DO
    END DO
!
    DO j=1,SetDim
       jj=Set(j)
       sig(j)=si1(jj)
       IF (zerors) THEN
          W(j,j)=sig(j)
          AGe(j)=0.0d0
       ELSE
          W(j,j)=a1(idiag1(jj))
          AGe(j)=W(j,j)-sig(j)
       END IF
       l2=j+1
       DO l=l2,SetDim
          W(j,l)=0.0d0
          W(l,j)=0.0d0
       END DO
       k2=ia1(jj+1)-1
       DO k=idiag1(jj)+1,k2
          DO l=l2,SetDim
             m=Set(l)
             IF(ja1(k)==m)THEN
                alpha=a1(k)
                W(j,l)=alpha
                W(l,j)=alpha
                EXIT
             END IF
          END DO
       END DO
    END DO
!
    DO j=1,SetDim
       DO k=1,SetDim
          IF (j.ne.k) THEN
             sig(j)=sig(j)+W(j,k)
          END IF
       ENDDO
       IF (sig(j) < 0.0d0)  AGe(j)=AGe(j)+2*sig(j)
   !
   !
       W(j,j)=W(j,j)-abs(sig(j))
   !
   !
   !
   !
       tmp=2*abs(sig(j))
       W(j,j)=W(j,j)-bndmum1m1*tmp
   !
       v(j)=tmp+AGe(j)
       IF (j .eq. 1) THEN
          beta=v(j)
          alp=abs(AGe(j))
       ELSE
          beta=beta+v(j)
          alp=max(alp,abs(AGe(j)))
       END IF
    END DO
!
!
    beta=bndmum1m1/beta
    DO j=1,SetDim
       DO k=1,SetDim
          W(j,k)=W(j,k)+beta*v(j)*v(k)
       END DO
    END DO
!
!
    IF (alp.LT.repsmach*beta) THEN
       SetDim1=SetDim-1
    ELSE
       SetDim1=SetDim
    END IF
!
!
    acc=.FALSE.
!
    SELECT CASE (SetDim1)
    CASE (1)
       GOTO 11
    CASE (2)
       GOTO 12
    CASE (3)
       GOTO 13
    CASE (4)
       GOTO 14
    CASE (5)
       GOTO 15
    CASE (6)
       GOTO 16
    CASE (7)
       GOTO 17
    CASE (8)
       GOTO 18
    CASE DEFAULT
       CALL DPOTRF('U',SetDim1,W,mm,info)
       IF (info .NE. 0) RETURN
       GOTO 10
    END SELECT
18  CONTINUE
    IF (W(8,8) .LE. 0.0d0) RETURN
    W(7,7) = W(7,7) - (W(7,8)/W(8,8)) * W(7,8)
    T = W(6,8)/W(8,8)
    W(6,7) = W(6,7) - T * W(7,8)
    W(6,6) = W(6,6) - T * W(6,8)
    T = W(5,8)/W(8,8)
    W(5,7) = W(5,7) - T * W(7,8)
    W(5,6) = W(5,6) - T * W(6,8)
    W(5,5) = W(5,5) - T * W(5,8)
    T = W(4,8)/W(8,8)
    W(4,7) = W(4,7) - T * W(7,8)
    W(4,6) = W(4,6) - T * W(6,8)
    W(4,5) = W(4,5) - T * W(5,8)
    W(4,4) = W(4,4) - T * W(4,8)
    T = W(3,8)/W(8,8)
    W(3,7) = W(3,7) - T * W(7,8)
    W(3,6) = W(3,6) - T * W(6,8)
    W(3,5) = W(3,5) - T * W(5,8)
    W(3,4) = W(3,4) - T * W(4,8)
    W(3,3) = W(3,3) - T * W(3,8)
    T = W(2,8)/W(8,8)
    W(2,7) = W(2,7) - T * W(7,8)
    W(2,6) = W(2,6) - T * W(6,8)
    W(2,5) = W(2,5) - T * W(5,8)
    W(2,4) = W(2,4) - T * W(4,8)
    W(2,3) = W(2,3) - T * W(3,8)
    W(2,2) = W(2,2) - T * W(2,8)
    T = W(1,8)/W(8,8)
    W(1,7) = W(1,7) - T * W(7,8)
    W(1,6) = W(1,6) - T * W(6,8)
    W(1,5) = W(1,5) - T * W(5,8)
    W(1,4) = W(1,4) - T * W(4,8)
    W(1,3) = W(1,3) - T * W(3,8)
    W(1,2) = W(1,2) - T * W(2,8)
    W(1,1) = W(1,1) - T * W(1,8)
17  CONTINUE
    IF (W(7,7) .LE. 0.0d0) RETURN
    W(6,6) = W(6,6) - (W(6,7)/W(7,7)) * W(6,7)
    T = W(5,7)/W(7,7)
    W(5,6) = W(5,6) - T * W(6,7)
    W(5,5) = W(5,5) - T * W(5,7)
    T = W(4,7)/W(7,7)
    W(4,6) = W(4,6) - T * W(6,7)
    W(4,5) = W(4,5) - T * W(5,7)
    W(4,4) = W(4,4) - T * W(4,7)
    T = W(3,7)/W(7,7)
    W(3,6) = W(3,6) - T * W(6,7)
    W(3,5) = W(3,5) - T * W(5,7)
    W(3,4) = W(3,4) - T * W(4,7)
    W(3,3) = W(3,3) - T * W(3,7)
    T = W(2,7)/W(7,7)
    W(2,6) = W(2,6) - T * W(6,7)
    W(2,5) = W(2,5) - T * W(5,7)
    W(2,4) = W(2,4) - T * W(4,7)
    W(2,3) = W(2,3) - T * W(3,7)
    W(2,2) = W(2,2) - T * W(2,7)
    T = W(1,7)/W(7,7)
    W(1,6) = W(1,6) - T * W(6,7)
    W(1,5) = W(1,5) - T * W(5,7)
    W(1,4) = W(1,4) - T * W(4,7)
    W(1,3) = W(1,3) - T * W(3,7)
    W(1,2) = W(1,2) - T * W(2,7)
    W(1,1) = W(1,1) - T * W(1,7)
16  CONTINUE
    IF (W(6,6) .LE. 0.0d0) RETURN
    W(5,5) = W(5,5) - (W(5,6)/W(6,6)) * W(5,6)
    T = W(4,6)/W(6,6)
    W(4,5) = W(4,5) - T * W(5,6)
    W(4,4) = W(4,4) - T * W(4,6)
    T = W(3,6)/W(6,6)
    W(3,5) = W(3,5) - T * W(5,6)
    W(3,4) = W(3,4) - T * W(4,6)
    W(3,3) = W(3,3) - T * W(3,6)
    T = W(2,6)/W(6,6)
    W(2,5) = W(2,5) - T * W(5,6)
    W(2,4) = W(2,4) - T * W(4,6)
    W(2,3) = W(2,3) - T * W(3,6)
    W(2,2) = W(2,2) - T * W(2,6)
    T = W(1,6)/W(6,6)
    W(1,5) = W(1,5) - T * W(5,6)
    W(1,4) = W(1,4) - T * W(4,6)
    W(1,3) = W(1,3) - T * W(3,6)
    W(1,2) = W(1,2) - T * W(2,6)
    W(1,1) = W(1,1) - T * W(1,6)
15  CONTINUE
    IF (W(5,5) .LE. 0.0d0) RETURN
    W(4,4) = W(4,4) - (W(4,5)/W(5,5)) * W(4,5)
    T = W(3,5)/W(5,5)
    W(3,4) = W(3,4) - T * W(4,5)
    W(3,3) = W(3,3) - T * W(3,5)
    T = W(2,5)/W(5,5)
    W(2,4) = W(2,4) - T * W(4,5)
    W(2,3) = W(2,3) - T * W(3,5)
    W(2,2) = W(2,2) - T * W(2,5)
    T = W(1,5)/W(5,5)
    W(1,4) = W(1,4) - T * W(4,5)
    W(1,3) = W(1,3) - T * W(3,5)
    W(1,2) = W(1,2) - T * W(2,5)
    W(1,1) = W(1,1) - T * W(1,5)
14  CONTINUE
    IF (W(4,4) .LE. 0.0d0) RETURN
    W(3,3) = W(3,3) - (W(3,4)/W(4,4)) * W(3,4)
    T = W(2,4)/W(4,4)
    W(2,3) = W(2,3) - T * W(3,4)
    W(2,2) = W(2,2) - T * W(2,4)
    T = W(1,4)/W(4,4)
    W(1,3) = W(1,3) - T * W(3,4)
    W(1,2) = W(1,2) - T * W(2,4)
    W(1,1) = W(1,1) - T * W(1,4)
13  CONTINUE
    IF (W(3,3) .LE. 0.0d0) RETURN
    W(2,2) = W(2,2) - (W(2,3)/W(3,3)) * W(2,3)
    T = W(1,3)/W(3,3)
    W(1,2) = W(1,2) - T * W(2,3)
    W(1,1) = W(1,1) - T * W(1,3)
12  CONTINUE
    IF (W(2,2) .LE. 0.0d0) RETURN
    W(1,1) = W(1,1) - (W(1,2)/W(2,2)) * W(1,2)
11  CONTINUE
    IF (W(1,1) .LE. 0.0d0) RETURN
10  CONTINUE
!
    acc=.TRUE.
!
    RETURN
  END SUBROUTINE dagmg_checktentagg_SF
  END SUBROUTINE dagmg_findpairs_SF
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_findpairs_GI(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,ipc )
!
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: nddl
    INTEGER :: ldd(nddl)
    LOGICAL :: skipass
    INTEGER :: ipc(n)
!
!
!
!
!
!
!
!
!
!
!
!----------------
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_dampJac
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
!
      DO WHILE (nmark.LT.nt)
       isel=ijs
       ijs=ijs+1
       !
       !
       IF (ind(isel) .EQ. 0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       !
       IF (ind(isel) .GE. 0) CYCLE
       !
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       IF (ipc(isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       !
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j) .GE. 0) CYCLE
          IF(ipc(j).EQ.0) CYCLE
          kk=0
          IF (i .LT. idiag(isel)) THEN
             j2=ia(j+1)-1
             DO jk=idiag(j)+1,j2
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ELSE
             DO jk=ia(j),idiag(j)-1
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ENDIF
          vals=-a(i)/2
          IF(kk .NE. 0) vals=vals-a(kk)/2
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
             eta1=2*si(isel)
             eta2=2*si(j)
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
             eta1=2*a(idiag(isel))
             eta2=2*a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !
          IF (sig1 > 0.0d0) THEN
             del1=rsi
          ELSE
             del1=rsi+2*sig1
          END IF
          IF (sig2 > 0.0d0) THEN
             del2=rsj
          ELSE
             del2=rsj+2*sig2
          END IF
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=((eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=(eta1*eta2)/(eta1+eta2)
             valp=vals/valp
          END IF
          IF (valp > kaptg) CYCLE
          !
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       !
       IF (ipair .EQ. 0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
!
!
    IF (ifirst .AND. 3*nc.GT.2*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       imult=2*imult
       ifirst=.FALSE.
       IF (wfo) THEN
          WRITE (iout,901) IRANK,imult*kaptg_dampJac
       END IF
901    FORMAT(i3, &
       '*   Coarsening too slow: quality threshold increased to',f6.2)
       ind(1:n)=-1
       DO i=1,nddl
          ind(ldd(i))=0
       END DO
       GOTO 1
    ENDIF
    RETURN
  END SUBROUTINE dagmg_findpairs_GI
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_findpairs_SI(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,ipc )
!
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: nddl
    INTEGER :: ldd(nddl)
    LOGICAL :: skipass
    INTEGER :: ipc(n)
!
!
!
!
!
!
!
!
!
!
!
!----------------
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_blocdia
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
!
      DO WHILE (nmark.LT.nt)
       isel=ijs
       ijs=ijs+1
       !
       !
       IF (ind(isel) .EQ. 0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       !
       IF (ind(isel) .GE. 0) CYCLE
       !
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       IF (ipc(isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       !
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j) .GE. 0) CYCLE
          IF(ipc(j).EQ.0) CYCLE
          vals=-a(i)
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !
          IF (sig1 > 0.0d0) THEN
             del1=rsi
             eta1=rsi+2*sig1
          ELSE
             del1=rsi+2*sig1
             eta1=rsi
          END IF
          IF (eta1 < 0.0d0) CYCLE
          IF (sig2 > 0.0d0) THEN
             del2=rsj
             eta2=rsj+2*sig2
          ELSE
             del2=rsj+2*sig2
             eta2=rsj
          END IF
          IF (eta2 < 0.0d0) CYCLE
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
               valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=(vals+(eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=vals+(eta1*eta2)/(eta1+eta2)
             IF (vals < 0.0d0) CYCLE
             valp=vals/valp
          END IF
          IF (valp > kaptg) CYCLE
          !
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       !
       IF (ipair .EQ. 0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
!
!
    IF (ifirst .AND. 3*nc.GT.2*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       imult=2*imult
       ifirst=.FALSE.
       IF (wfo) THEN
          WRITE (iout,901) IRANK,imult*kaptg_blocdia
       END IF
901    FORMAT(i3, &
       '*   Coarsening too slow: quality threshold increased to',f6.2)
       ind(1:n)=-1
       DO i=1,nddl
          ind(ldd(i))=0
       END DO
       GOTO 1
    ENDIF
    RETURN
  END SUBROUTINE dagmg_findpairs_SI
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_findpairs_GI1(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,riperm,iperm )
!
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: nddl
    INTEGER :: ldd(nddl)
    LOGICAL :: skipass
    INTEGER :: iperm(n),riperm(n)
!
!
!
!
!
!
!
!
!
!
!
!----------------
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_dampJac
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
!
      DO WHILE (nmark.LT.nt)
       isel=ijs
       isel=riperm(ijs)
       ijs=ijs+1
       !
       !
       IF (ind(isel) .EQ. 0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       !
       IF (ind(isel) .GE. 0) CYCLE
       !
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       IF (iperm(isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       !
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j) .GE. 0) CYCLE
          IF(iperm(j).EQ.0) CYCLE
          kk=0
          IF (i .LT. idiag(isel)) THEN
             j2=ia(j+1)-1
             DO jk=idiag(j)+1,j2
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ELSE
             DO jk=ia(j),idiag(j)-1
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ENDIF
          vals=-a(i)/2
          IF(kk .NE. 0) vals=vals-a(kk)/2
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
             eta1=2*si(isel)
             eta2=2*si(j)
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
             eta1=2*a(idiag(isel))
             eta2=2*a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !
          IF (sig1 > 0.0d0) THEN
             del1=rsi
          ELSE
             del1=rsi+2*sig1
          END IF
          IF (sig2 > 0.0d0) THEN
             del2=rsj
          ELSE
             del2=rsj+2*sig2
          END IF
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=((eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=(eta1*eta2)/(eta1+eta2)
             valp=vals/valp
          END IF
          IF (valp > kaptg) CYCLE
          !
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. iperm(j).LT.iperm(ipair))  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       !
       IF (ipair .EQ. 0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
!
!
    IF (ifirst .AND. 3*nc.GT.2*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       imult=2*imult
       ifirst=.FALSE.
       IF (wfo) THEN
          WRITE (iout,901) IRANK,imult*kaptg_dampJac
       END IF
901    FORMAT(i3, &
       '*   Coarsening too slow: quality threshold increased to',f6.2)
       ind(1:n)=-1
       DO i=1,nddl
          ind(ldd(i))=0
       END DO
       GOTO 1
    ENDIF
    RETURN
  END SUBROUTINE dagmg_findpairs_GI1
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_findpairs_SI1(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,riperm,iperm )
!
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: nddl
    INTEGER :: ldd(nddl)
    LOGICAL :: skipass
    INTEGER :: iperm(n),riperm(n)
!
!
!
!
!
!
!
!
!
!
!
!----------------
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_blocdia
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
!
      DO WHILE (nmark.LT.nt)
       isel=ijs
       isel=riperm(ijs)
       ijs=ijs+1
       !
       !
       IF (ind(isel) .EQ. 0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       !
       IF (ind(isel) .GE. 0) CYCLE
       !
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       IF (iperm(isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       !
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j) .GE. 0) CYCLE
          IF(iperm(j).EQ.0) CYCLE
          vals=-a(i)
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !
          IF (sig1 > 0.0d0) THEN
             del1=rsi
             eta1=rsi+2*sig1
          ELSE
             del1=rsi+2*sig1
             eta1=rsi
          END IF
          IF (eta1 < 0.0d0) CYCLE
          IF (sig2 > 0.0d0) THEN
             del2=rsj
             eta2=rsj+2*sig2
          ELSE
             del2=rsj+2*sig2
             eta2=rsj
          END IF
          IF (eta2 < 0.0d0) CYCLE
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
               valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=(vals+(eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=vals+(eta1*eta2)/(eta1+eta2)
             IF (vals < 0.0d0) CYCLE
             valp=vals/valp
          END IF
          IF (valp > kaptg) CYCLE
          !
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. iperm(j).LT.iperm(ipair))  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       !
       IF (ipair .EQ. 0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
!
!
    IF (ifirst .AND. 3*nc.GT.2*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       imult=2*imult
       ifirst=.FALSE.
       IF (wfo) THEN
          WRITE (iout,901) IRANK,imult*kaptg_blocdia
       END IF
901    FORMAT(i3, &
       '*   Coarsening too slow: quality threshold increased to',f6.2)
       ind(1:n)=-1
       DO i=1,nddl
          ind(ldd(i))=0
       END DO
       GOTO 1
    ENDIF
    RETURN
  END SUBROUTINE dagmg_findpairs_SI1
!------------------------------------------------------------------
  SUBROUTINE dagmg_lcgmix(nc,m,lcg1,lcg,lcgn)
    INTEGER :: nc,m,lcg1(m,*),lcg(2,*),lcgn(2*m,*),i,l,l1,l2
    IF (m.eq.2) THEN
       DO i=1,nc
          IF(lcg(2,i) .EQ. 0) THEN
             lcgn(1,i)=lcg1(1,lcg(1,i))
             lcgn(2,i)=0
             lcgn(4,i)=-1
          ELSE IF(lcg(2,i) .LT. 0) THEN
             IF (lcg1(2,lcg(1,i)) .LT. 0) THEN
                lcgn(1,i)=lcg1(1,lcg(1,i))
                lcgn(2,i)=-1
                lcgn(4,i)=-1
             ELSE
                lcgn(1,i)=lcg1(1,lcg(1,i))
                lcgn(2,i)=lcg1(2,lcg(1,i))
                lcgn(4,i)=-2
             END IF
          ELSE
             IF (lcg1(2,lcg(1,i)) .LT. 0) THEN
                IF (lcg1(2,lcg(2,i)) .LT. 0) THEN
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(1,lcg(2,i))
                   lcgn(4,i)=-2
                ELSE
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(1,lcg(2,i))
                   lcgn(3,i)=lcg1(2,lcg(2,i))
                   lcgn(4,i)=-3
                END IF
             ELSE
                IF (lcg1(2,lcg(2,i)) .LT. 0) THEN
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(2,lcg(1,i))
                   lcgn(3,i)=lcg1(1,lcg(2,i))
                   lcgn(4,i)=-3
                ELSE
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(2,lcg(1,i))
                   lcgn(3,i)=lcg1(1,lcg(2,i))
                   lcgn(4,i)=lcg1(2,lcg(2,i))
                END IF
             END IF
          END IF
       END DO
    ELSE
       DO i=1,nc
          IF(lcg(2,i) .EQ. 0) THEN
             lcgn(1,i)=lcg1(1,lcg(1,i))
             lcgn(2,i)=0
             lcgn(2*m,i)=-1
          ELSE
             lcgn(2,i)=-1
             l1=m
             IF (lcg1(m,lcg(1,i)).LT.0) l1=-lcg1(m,lcg(1,i))
             lcgn(1:l1,i)=lcg1(1:l1,lcg(1,i))
             IF(lcg(2,i) .LT. 0) THEN
                l=l1
             ELSE
                l2=m
                IF (lcg1(m,lcg(2,i)).LT.0) l2=-lcg1(m,lcg(2,i))
                lcgn(l1+1:l1+l2,i)=lcg1(1:l2,lcg(2,i))
                l=l1+l2
             END IF
             IF(l .LT. 2*m) lcgn(2*m,i)=-l
          END IF
       END DO
    END IF
    RETURN
  END SUBROUTINE dagmg_lcgmix
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_setind(nc,ndd,ldd,lcg,m,ind)
    INTEGER :: nc,m,lcg(m,*),nll,ldd(ndd),ind(*),i,k,l
    DO i=1,ndd
       ind(ldd(i))=0
    END DO
    DO i=1,nc
       l=m
       IF (lcg(m,i) .LT. 0) l=-lcg(m,i)
       DO k=1,l
          ind(lcg(k,i))=i
       END DO
    END DO
    RETURN
  END SUBROUTINE dagmg_setind
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_setcg(n,a,ja,ia,idiag,si,ind,lcg           &
       ,nc,a2,ja2,ia2,idiag2,si2,ysi,maxdg,iw,w,iext,iext2)
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind(n),nc,lcg(2,nc)
    INTEGER :: ja2(*),ia2(nc+1),idiag2(nc),maxdg
    INTEGER, TARGET :: iw(2*nc)
    INTEGER, OPTIONAL :: iext(*),iext2(*)
    REAL(kind(0.0d0)) :: a(*),a2(*),w(nc),vald
    REAL(kind(0.0d0)) :: si(n),si2(*),sii
    LOGICAL :: ysi
    INTEGER :: nz,nzu,i,jj,jc,jcol,ki,kb,kf,jpos
    INTEGER, POINTER, DIMENSION(:) :: iw1, iw2
    !
    iw1 => iw(1:nc)
    iw2 => iw(nc+1:2*nc)
    !
    nz = 0
    iw1(1:nc)=0
    maxdg=0
    ia2(1)=1
    DO i = 1,nc
       sii=0.0d0
       vald=0.0d0
       nzu=0
       DO ki= 1,2
          jj = lcg(ki,i)
          IF (ki.EQ.1 .OR. jj.GT.0) THEN
             IF (ysi) sii=sii+si(jj)
             kf = ia(jj+1) - 1
             DO kb = ia(jj),kf
                jc = ja(kb)
                jcol = ind(jc)
                IF (jcol .GT. 0) THEN
                   IF (jcol .LT. i) THEN
                      jpos = iw1(jcol)
                      IF (jpos.EQ.0) THEN
                         nz = nz+1
                         ja2(nz) = jcol
                         iw1(jcol) = nz
                         a2(nz) = a(kb)
                      ELSE
                         a2(jpos) = a2(jpos) + a(kb)
                      ENDIF
                   ELSE IF (jcol .GT. i) THEN
                      jpos = iw1(jcol)
                      IF (jpos.EQ.0) THEN
                         nzu = nzu+1
                         iw2(nzu) = jcol
                         iw1(jcol) = nzu
                         w(nzu) = a(kb)
                      ELSE
                         w(jpos) = w(jpos) + a(kb)
                      ENDIF
                   ELSE
                      vald=vald+a(kb)
                      IF (ysi .AND. jc.NE.jj) sii=sii+a(kb)
                   ENDIF
                ENDIF
             ENDDO
          ENDIF
       ENDDO
       nz=nz+1
       a2(nz)=vald
       idiag2(i)=nz
       ja2(nz)=i
       a2(nz+1:nz+nzu)=w(1:nzu)
       ja2(nz+1:nz+nzu)=iw2(1:nzu)
       nz=nz+nzu
       maxdg=max(maxdg,nz-ia2(i))
       DO kb = ia2(i), nz
          iw1(ja2(kb))=0
       ENDDO
       IF (ysi) si2(i)=sii
       ia2(i+1)=nz+1
    ENDDO
    RETURN
  END SUBROUTINE dagmg_setcg
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_LAPACK(n,f,ijb,a,ja,ia)
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ijb
    REAL(kind(0.0d0)), TARGET :: f(n)
    REAL(kind(0.0d0)), OPTIONAL, TARGET :: a(*)
    INTEGER, OPTIONAL :: ia(n+1)
    INTEGER, OPTIONAL, TARGET :: ja(*)
    !
    REAL(kind(0.0d0)), ALLOCATABLE, SAVE :: ac(:,:)
    INTEGER, ALLOCATABLE, SAVE :: ipiv(:)
    REAL(kind(0.0d0)), SAVE :: iflop
    INTEGER :: i,kk,iext
    INTEGER  :: ierr
    INTEGER , parameter :: IONE=1
    !
    ierr=0
    IF (ijb == -2) THEN
       !
       DEALLOCATE (ac,ipiv)
       !
    ELSE IF (ijb == 1) THEN
       !
       ALLOCATE (ac(n,n),ipiv(n))
       memi=memi+n
       memr=memr+n*n
       memax=MAX(memax,memr+memi*rlenilen)
       ac=0.0d0
       DO i=1,n
          DO kk=ia(i),ia(i+1)-1
             ac(i,ja(kk))=a(kk)
          END DO
       END DO
       CALL DGETRF(N,N,ac,N,ipiv,ierr)
       IF (ierr /= 0) THEN
          WRITE(iout, *) ' FATAL ERROR in GETRF: ierror=',ierr
          STOP
       END IF
       iflop=2*dble(n)**2-n
       flop=(2*1.0d0)/(3*1.0d0)*(dble(n)**3)
       !
    ELSE IF (ijb == 2) THEN
       !
       CALL DGETRS('N',N,IONE,ac,N,ipiv,f,N,ierr)
       IF (ierr /= 0) THEN
          WRITE(iout, *) ' FATAL ERROR in GETRS: ierror=',ierr
          STOP
       END IF
       flop=flop+iflop
       !
    END IF
    !
    RETURN
  END SUBROUTINE dagmg_LAPACK
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_DIRseq(n,f,ijob,a,ja,ia)
    USE dagmg_mem
    USE dagmg_PARDISO
    IMPLICIT NONE
    INTEGER :: n,ijob
    REAL(kind(0.0d0)), TARGET :: f(n)
    REAL(kind(0.0d0)), OPTIONAL, TARGET :: a(*)
    INTEGER, OPTIONAL :: ia(n+1)
    INTEGER, OPTIONAL, TARGET :: ja(*)
    REAL(kind(0.0d0)) , SAVE :: iflop
    TYPE(PARDISO_HANDLE), ALLOCATABLE, SAVE :: PT(:)
    INTEGER, SAVE :: IPARM(64)
    REAL(kind(0.0d0)), ALLOCATABLE, SAVE :: X(:)
    INTEGER, SAVE :: MAXFCT, MNUM, MTYPE, PHASE, NRHS, MSGLVL, ERROR
!!!!!!!
!
!
      TYPE DMUMPS_ROOT_STRUC
        SEQUENCE
        INTEGER MBLOCK, NBLOCK, NPROW, NPCOL
        INTEGER MYROW, MYCOL
        INTEGER ROOT_SIZE, TOT_ROOT_SIZE
        INTEGER :: CNTXT_BLACS, truc
        INTEGER, DIMENSION(:), POINTER :: RG2L_ROW
        INTEGER, DIMENSION(:), POINTER :: RG2L_COL
        INTEGER , DIMENSION(:), POINTER :: IPIV
        INTEGER, DIMENSION( 9 ) :: DESCRIPTOR, DESCB
        LOGICAL yes, gridinit_done
        INTEGER LPIV, brol
!       Used to access Schur easily from root structure
        DOUBLE PRECISION, DIMENSION(:), POINTER :: SCHUR_POINTER
        INTEGER SCHUR_MLOC, SCHUR_NLOC, SCHUR_LLD, machin
!
!      Data for nullspace/QR
!
        DOUBLE PRECISION, DIMENSION(:), POINTER :: QR_TAU
        DOUBLE PRECISION     QR_RCOND
!
!      Givens rotations
!
        INTEGER MAXG, GIND
        DOUBLE PRECISION, DIMENSION(:),POINTER::GROW, GCOS, GSIN
!
!      RRRLU data
!
        INTEGER ELG_MAX,NULL_MAX
        INTEGER ELIND,EUIND,NLUPDATE,NUUPDATE
        INTEGER,DIMENSION(:),POINTER::PERM_ROW,PERM_COL
        INTEGER,DIMENSION(:),POINTER::ELROW, EUROW, PTREL, PTREU
        DOUBLE PRECISION, DIMENSION(:), POINTER :: ELELG, EUELG, DL
!
      END TYPE DMUMPS_ROOT_STRUC
      TYPE DMUMPS_STRUC
        SEQUENCE
!
! This structure contains all parameters
! for the interface to the user, plus internal
! information
!
! *****************
! INPUT PARAMETERS
! *****************
!    -----------------
!    MPI Communicator
!    -----------------
        INTEGER COMM
!    ------------------
!    Problem definition
!    ------------------
!    Solver (SYM=0 unsymmetric,SYM=1 symmetric Positive Definite,
!        SYM=2 general symmetric)
!    Type of parallelism (PAR=1 host working, PAR=0 host not working)
        INTEGER SYM, PAR
        INTEGER JOB
!    --------------------
!    Order of Input matrix
!    --------------------
        INTEGER N
!
!    ----------------------------------------
!    Assembled input matrix : User interface
!    ----------------------------------------
        INTEGER NZ
        DOUBLE PRECISION, DIMENSION(:), POINTER :: A
        INTEGER, DIMENSION(:), POINTER :: IRN, JCN
        DOUBLE PRECISION, DIMENSION(:), POINTER :: COLSCA, ROWSCA, pad0
!
!       ------------------------------------
!       Case of distributed assembled matrix
!       matrix on entry:
!       ------------------------------------
        INTEGER NZ_loc, pad1
        INTEGER, DIMENSION(:), POINTER :: IRN_loc, JCN_loc
        DOUBLE PRECISION, DIMENSION(:), POINTER :: A_loc, pad2
!
!    ----------------------------------------
!    Unassembled input matrix: User interface
!    ----------------------------------------
        INTEGER NELT, pad3
        INTEGER, DIMENSION(:), POINTER :: ELTPTR
        INTEGER, DIMENSION(:), POINTER :: ELTVAR
        DOUBLE PRECISION, DIMENSION(:), POINTER :: A_ELT, pad4
!
!    ---------------------------------------------
!    Symmetric permutation :
!               PERM_IN if given by user (optional)
!    ---------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: PERM_IN
!
!
! ******************
! INPUT/OUTPUT data
! ******************
!    --------------------------------------------------------
!    RHS / SOL_LOC
!    -------------
!       right-hand side and solution
!    -------------------------------------------------------
        DOUBLE PRECISION, DIMENSION(:), POINTER :: RHS, REDRHS
        DOUBLE PRECISION, DIMENSION(:), POINTER :: RHS_SPARSE
        DOUBLE PRECISION, DIMENSION(:), POINTER :: SOL_LOC
        INTEGER, DIMENSION(:), POINTER :: IRHS_SPARSE
        INTEGER, DIMENSION(:), POINTER :: IRHS_PTR
        INTEGER, DIMENSION(:), POINTER :: ISOL_LOC
        INTEGER LRHS, NRHS, NZ_RHS, LSOL_LOC, LREDRHS
        INTEGER pad5
!    ----------------------------
!    Control parameters,
!    statistics and output data
!    ---------------------------
        INTEGER ICNTL(40)
        INTEGER INFO(40)
        INTEGER INFOG(40)
        DOUBLE PRECISION COST_SUBTREES
        DOUBLE PRECISION CNTL(15)
        DOUBLE PRECISION RINFO(20)
        DOUBLE PRECISION RINFOG(20)
!    ---------------------------------------------------------
!    Permutations computed during analysis:
!       SYM_PERM: Symmetric permutation
!       UNS_PERM: Column permutations (optionnal)
!    ---------------------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: SYM_PERM, UNS_PERM
!
!    -------------------------------------
!    Case of distributed matrix on entry:
!    DMUMPS potentially provides mapping
!    -------------------------------------
        INTEGER, DIMENSION(:), POINTER :: MAPPING
!
!    -------------------------------
!    Deficiency and null space basis
!    -------------------------------
        DOUBLE PRECISION, DIMENSION(:,:), POINTER :: NULL_SPACE
        INTEGER Deficiency, pad6
!    -----
!    Schur
!    -----
        INTEGER NPROW, NPCOL, MBLOCK, NBLOCK
        INTEGER SCHUR_MLOC, SCHUR_NLOC, SCHUR_LLD
        INTEGER SIZE_SCHUR
        INTEGER, DIMENSION(:), POINTER :: LISTVAR_SCHUR
        DOUBLE PRECISION, DIMENSION(:), POINTER :: SCHUR
        DOUBLE PRECISION, DIMENSION(:), POINTER :: SCHUR_CINTERFACE
!    --------------
!    Version number
!    --------------
        CHARACTER(LEN=16) VERSION_NUMBER
!    -----------
!    Out-of-core
!    -----------
        CHARACTER(LEN=256) :: OOC_TMPDIR
        CHARACTER(LEN=64) :: OOC_PREFIX
!    ------------------------------------------
!    To save the matrix in matrix market format
!    ------------------------------------------
        CHARACTER(LEN=256) WRITE_PROBLEM
!
!
! **********************
! INTERNAL Working data
! *********************
        INTEGER INST_Number
!       For MPI
        INTEGER COMM_NODES, MYID_NODES, COMM_LOAD
        INTEGER  MYID, NPROCS, NSLAVES
        INTEGER ASS_IRECV
        INTEGER, DIMENSION(:), POINTER :: POIDS
        INTEGER LBUFR
        INTEGER LBUFR_BYTES
        INTEGER, DIMENSION(:), POINTER ::  BUFR
!       For analysis/facto/solve phases
        INTEGER MAXIS1, pad7
        INTEGER KEEP(500)
        INTEGER(8) KEEP8(150)
!       IS is used for the factors + workspace for contrib. blocks
        INTEGER, DIMENSION(:), POINTER :: IS
!       is1 (maxis1) contains working arrays computed
!       and used only during analysis
        INTEGER, DIMENSION(:), POINTER :: IS1
!       The following data/arrays are computed during the analysis
!       phase and used during the factorization and solve phases.
        INTEGER LNA
        INTEGER NBSA
        INTEGER,POINTER,DIMENSION(:)::STEP, NE_STEPS, ND_STEPS
!  Info for pruning tree
        INTEGER,POINTER,DIMENSION(:)::Step2node
!  ---------------------
        INTEGER,POINTER,DIMENSION(:)::FRERE_STEPS, DAD_STEPS
        INTEGER,POINTER,DIMENSION(:)::FILS, PTRAR, FRTPTR, FRTELT
        INTEGER,POINTER,DIMENSION(:)::NA, PROCNODE_STEPS
!       The two pointer arrays computed in facto and used by the solve
!          (except the factors) are PTLUST_S and PTRFAC.
        INTEGER, DIMENSION(:), POINTER :: PTLUST_S
        INTEGER(8), DIMENSION(:), POINTER :: PTRFAC
!       main real working arrays for factorization/solve phases
        DOUBLE PRECISION, DIMENSION(:), POINTER :: S
!       Information on mapping
        INTEGER, DIMENSION(:), POINTER :: PROCNODE
!       Input matrix ready for numerical assembly
!           -arrowhead format in case of assembled matrix
!           -element format otherwise
        INTEGER, DIMENSION(:), POINTER :: INTARR
        DOUBLE PRECISION, DIMENSION(:), POINTER :: DBLARR
!       Element entry: internal data
        INTEGER NELT_LOC, LELTVAR, NA_ELT, pad8
        INTEGER, DIMENSION(:), POINTER :: ELTPROC
!       Candidates and node partitionning
        INTEGER, DIMENSION(:,:), POINTER :: CANDIDATES
        INTEGER, DIMENSION(:),   POINTER :: ISTEP_TO_INIV2
        INTEGER, DIMENSION(:),   POINTER :: FUTURE_NIV2
        INTEGER, DIMENSION(:,:), POINTER :: TAB_POS_IN_PERE
        LOGICAL, DIMENSION(:), POINTER :: I_AM_CAND
!       For heterogeneous architecture
        INTEGER, DIMENSION(:), POINTER :: MEM_DIST
!       Compressed RHS
        INTEGER, DIMENSION(:),   POINTER :: POSINRHSCOMP
        DOUBLE PRECISION, DIMENSION(:), POINTER :: RHSCOMP
!       For C interface
!   Info on the subtrees to be used during factorization
        DOUBLE PRECISION, DIMENSION(:),   POINTER :: MEM_SUBTREE
        INTEGER, DIMENSION(:),   POINTER :: MY_ROOT_SBTR
        INTEGER, DIMENSION(:),   POINTER :: MY_FIRST_LEAF
        INTEGER, DIMENSION(:),   POINTER :: MY_NB_LEAF
        INTEGER, DIMENSION(:),   POINTER :: DEPTH_FIRST
        DOUBLE PRECISION, DIMENSION(:),   POINTER :: COST_TRAV
        INTEGER NBSA_LOCAL, zwave
        INTEGER(8) :: MAX_SURF_MASTER
        INTEGER :: LWK_USER, zozo
        DOUBLE PRECISION, DIMENSION(:), POINTER :: WK_USER
!    For simulating parallel out-of-core stack.
        DOUBLE PRECISION, DIMENSION(:),POINTER ::CB_SON_SIZE
!   Instance number used/managed by the C/F77 interface
        INTEGER INSTANCE_NUMBER
!    OOC management data that must persist from factorization to solve.
        INTEGER OOC_MAX_NB_NODES_FOR_ZONE
        INTEGER, DIMENSION(:,:),   POINTER :: OOC_INODE_SEQUENCE
        INTEGER(8),DIMENSION(:,:), POINTER :: OOC_SIZE_OF_BLOCK
        INTEGER(8), DIMENSION(:,:),   POINTER :: OOC_VADDR
        INTEGER,DIMENSION(:), POINTER :: OOC_TOTAL_NB_NODES
        INTEGER,DIMENSION(:), POINTER :: OOC_NB_FILES
        CHARACTER,DIMENSION(:,:), POINTER :: OOC_FILE_NAMES
        INTEGER,DIMENSION(:), POINTER :: OOC_FILE_NAME_LENGTH
!    Indices of nul pivots
        INTEGER,DIMENSION(:), POINTER :: PIVNUL_LIST
!    Internal control array
        DOUBLE PRECISION DKEEP(30)
!    Array needed to manage additionnal candidate processor
        INTEGER, DIMENSION(:,:), POINTER :: SUP_PROC
!   ------------------------
!   Root structure(internal)
!   ------------------------
        TYPE (DMUMPS_ROOT_STRUC) :: root
      END TYPE DMUMPS_STRUC
!!!!!!!
    TYPE(DMUMPS_STRUC), SAVE :: mumps_par
    LOGICAL, PARAMETER :: mklpardiso=.FALSE.
    REAL(kind(0.0d0)), ALLOCATABLE, TARGET, SAVE :: coak(:)
    INTEGER :: i, j, k, nnzi, idum(1), kzi
    REAL(kind(0.0d0)) :: td, dum(1)
    REAL(kind(0.0d0)) :: absmax, tt
    REAL(kind(0.0d0)), external :: DDOT
    REAL(kind(0.0d0)), external :: DNRM2
    INTEGER , parameter :: IONE=1
    MAXFCT=1
    MNUM=1
    MTYPE=11
    NRHS=1
    MSGLVL=0
    !
    IF (ijob == -2) THEN
       !
       IF (mklpardiso) THEN
          PHASE=-1
          CALL PARDISO( PT, MAXFCT, MNUM, MTYPE, PHASE, N      &
               , dum, idum, idum                               &
               , idum, NRHS, IPARM, MSGLVL, dum, dum, ERROR )
          IF (ERROR .NE. 0) THEN
             PRINT 999, ERROR, PHASE
             STOP
          END IF
          DEALLOCATE(PT,X)
       ELSE
          mumps_par%JOB = -2
          CALL DAGMG_MUMPS(mumps_par)
          IF (transint) THEN
             DEALLOCATE(mumps_par%JCN_loc)
             IF (nlev .EQ. 1) DEALLOCATE(mumps_par%IRN_loc,mumps_par%A_loc)
          ELSE
             DEALLOCATE(mumps_par%IRN_loc)
             IF (nlev .EQ. 1) DEALLOCATE(mumps_par%JCN_loc,mumps_par%A_loc)
          END IF
       END IF
       IF (coasing) DEALLOCATE(coak)
       !
    ELSE IF (ijob == 1) THEN
       !
       absmax=MAXVAL(ABS(a(1:ia(n+1)-1)))
       IF (mklpardiso) THEN
          ALLOCATE(PT(64))
          do i = 1, 64
            PT(i)%DUMMY = 0
          end do
          IPARM(1:64)=0
          IPARM(1)=1
          IPARM(2)=2
          IPARM(6)=1
          IPARM(10)=13
          IPARM(11)=1
          IF (transint) THEN
             IPARM(12)=2
          ELSE
             IPARM(12)=0
          END IF
          IPARM(18)=-1
          IPARM(19)=-1
          IPARM(20)=0
          ERROR  = 0
          PHASE=12
          CALL dagmg_sortrows(n,a,ja,ia)
          CALL PARDISO( PT, MAXFCT, MNUM, MTYPE, PHASE, N      &
               , a, ia, ja                                     &
               , idum, NRHS, IPARM, MSGLVL, dum, dum, ERROR )
          IF (ERROR .NE. 0) THEN
             PRINT 999, ERROR, PHASE
             STOP
          END IF
          flop=1e6*dble(IPARM(19))
          iflop=2*dble(IPARM(18))-n
          ALLOCATE(X(n))
       ELSE
          !
          mumps_par%COMM = 0
          mumps_par%JOB = -1
          mumps_par%SYM = 0
          mumps_par%PAR = 1
          CALL DAGMG_MUMPS(mumps_par)
          mumps_par%ICNTL(2)=-1
          mumps_par%ICNTL(3)=-1
          mumps_par%ICNTL(4)=0
          mumps_par%ICNTL(14)=80
          mumps_par%ICNTL(18)=3
          mumps_par%ICNTL(24)=1
          nnzi=ia(n+1)-1
          mumps_par%NZ_loc=nnzi
          IF (transint) THEN
             ALLOCATE( mumps_par%JCN_loc(nnzi) )
             do i=1,n
                do j=ia(i),ia(i+1)-1
                   mumps_par%JCN_loc(j)=i
                end do
             end do
             IF (nlev .GT. 1) THEN
               mumps_par%IRN_loc => dt(nlev)%ja(1:nnzi)
               mumps_par%A_loc => dt(nlev)%a(1:nnzi)
             ELSE
               ALLOCATE( mumps_par%IRN_loc(nnzi) , mumps_par%A_loc(nnzi) )
               mumps_par%IRN_loc = ja(1:nnzi)
               mumps_par%A_loc = a(1:nnzi)
               memi=memi+nnzi
               memr=memr+nnzi
             END IF
          ELSE
             ALLOCATE( mumps_par%IRN_loc(nnzi) )
             do i=1,n
                do j=ia(i),ia(i+1)-1
                   mumps_par%IRN_loc(j)=i
                 end do
             end do
             IF (nlev .GT. 1) THEN
               mumps_par%JCN_loc => dt(nlev)%ja(1:mumps_par%NZ_loc)
               mumps_par%A_loc => dt(nlev)%a(1:mumps_par%NZ_loc)
             ELSE
               ALLOCATE( mumps_par%JCN_loc(nnzi) , mumps_par%A_loc(nnzi) )
               mumps_par%JCN_loc = ja(1:nnzi)
               mumps_par%A_loc = a(1:nnzi)
               memi=memi+nnzi
               memr=memr+nnzi
             END IF
          END IF
          memi=memi+nnzi
          memax=MAX(memax,memr+memi*rlenilen)
          !
          mumps_par%N=n
          mumps_par%JOB = 4
          CALL DAGMG_MUMPS(mumps_par)
          flop=mumps_par%RINFO(3)
          i=mumps_par%INFO(27)
          IF (i .GT. 0) THEN
             iflop=2*dble(i)-n
          ELSE
             iflop=-2.0d6*dble(i)-n
          END IF
       END IF
       ALLOCATE(coak(n))
       coak(1:n)=1.0d0
       IF (mklpardiso) THEN
          PHASE=33
          CALL PARDISO( PT, MAXFCT, MNUM, MTYPE, PHASE, N             &
               , a, ia, ja                                            &
               , idum, NRHS, IPARM, MSGLVL, coak, X, ERROR )
          IF (ERROR .NE. 0) THEN
             PRINT 999, ERROR, PHASE
             STOP
          END IF
       ELSE
          mumps_par%JOB = 3
          mumps_par%RHS => coak(1:n)
          CALL DAGMG_MUMPS(mumps_par)
       END IF
       tt=DNRM2(N,coak,IONE)
       IF (  SQRT(dble(N))/(ABS(tt)*absmax) .LE. 10*epsmach  ) THEN
          coasing=.TRUE.
          coak(1:n)=(1.0d0/tt)*coak(1:n)
       ELSE
          coasing=.FALSE.
          DEALLOCATE(coak)
       ENDIF
    ELSE IF (ijob == 2) THEN
       IF (mklpardiso) THEN
          PHASE=33
          IF (PRESENT(A)) THEN
             CALL PARDISO( PT, MAXFCT, MNUM, MTYPE, PHASE, N      &
                  , a, ia, ja                                     &
                  , idum, NRHS, IPARM, MSGLVL, f, X, ERROR )
          ELSE
             CALL PARDISO( PT, MAXFCT, MNUM, MTYPE, PHASE, N      &
                  , dt(nlev)%a, dt(nlev)%ia, dt(nlev)%ja          &
                  , idum, NRHS, IPARM, MSGLVL, f, X, ERROR )
          END IF
          IF (ERROR .NE. 0) THEN
             PRINT 999, ERROR, PHASE
             STOP
          END IF
       ELSE
          mumps_par%JOB = 3
          mumps_par%RHS => f(1:n)
          CALL DAGMG_MUMPS(mumps_par)
       END IF
       flop=flop+iflop
       IF (coasing) THEN
          td=DDOT(N,coak,IONE,f,IONE)
          CALL DAXPY(N,-td,coak,IONE,f,IONE)
       END IF
    END IF
    !
    RETURN
999 FORMAT('Error detected in PARDISO; error code:', i2,' ; Phase:',i2)
  END SUBROUTINE dagmg_DIRseq
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_sortrows(n,a,ja,ia)
    IMPLICIT NONE
    INTEGER :: n, ia(n+1), ja(*), i, k, jj
    REAL(kind(0.0d0)) :: a(*), aa
    LOGICAL exc
    DO i=1,n
       exc=.TRUE.
       DO WHILE(exc)
          exc=.FALSE.
          DO k=ia(i)+1,ia(i+1)-1
             IF(ja(k).LT.ja(k-1)) THEN
                jj=ja(k-1)
                aa=a(k-1)
                ja(k-1)=ja(k)
                a(k-1)=a(k)
                ja(k)=jj
                a(k)=aa
                exc=.TRUE.
             END IF
          END DO
       END DO
    END DO
    RETURN
  END SUBROUTINE dagmg_sortrows
!!!!!!!!!!!!!!!!!!
END MODULE dagmg_ALLROUTINES
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! MAIN DRIVER
!-----------------------------------------------------------------------
  SUBROUTINE dagmg( n,a,ja,ia,f,x,ijob,iprint,nrest,iter,tol )
    USE dagmg_mem
    USE dagmg_ALLROUTINES
    IMPLICIT NONE
    INTEGER    :: n,ia(n+1),ja(*),ijob,iprint,nrest,iter
    REAL (kind(0.0d0)) :: a(*),f(n),x(n)
    REAL (kind(0.0d0)) :: tol
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Arguments
!  =========
!
!  N       (input) INTEGER.
!          The dimension of the matrix.
!
!  A       (input/output) REAL (kind(0.0d0)). Numerical values of the matrix.
!  IA      (input/output) INTEGER. Pointers for every row.
!  JA      (input/output) INTEGER. Column indices.
!
!              AGMG ASSUMES THAT ALL DIAGONAL ENTRIES ARE POSITIVE
!
!          Detailed description of the matrix format
!
!              On input, IA(I), I=1,...,N, refers to the physical start
!              of row I. That is, the entries of row I are located
!              in A(K), where K=IA(I),...,IA(I+1)-1. JA(K) carries the
!              associated column indices. IA(N+1) must be defined in such
!              a way that the above rule also works for I=N (that is,
!              the last valid entry in arrays A,JA should correspond to
!              index K=IA(N+1)-1). According what is written
!              above, AGMG assumes that some of these JA(K) (for
!              IA(I)<= K < IA(I+1) ) is equal to I with corresponding
!              A(K) carrying the value of the diagonal element,
!              which must be positive.
!
!              A,IA,JA are "output" parameters because on exit the
!              entries of each row may occur in a different order (The
!              matrix is mathematically the same, but stored in
!              different way).
!
!  F       (input/output) REAL (kind(0.0d0)).
!          On input, the right hand side vector f.
!          Overwritten on output.
!          Significant only if IJOB==0, 2, 3, 10, 12, 100, 102, 110, 112
!
!  X       (input/output) REAL (kind(0.0d0)).
!          On input and if IJOB== 10, 12, 110, 112: initial guess
!             (for other values of IJOB, the default is used: the zero vector).
!          On output, the computed solution.
!
! IJOB     (input) INTEGER. Tells AGMG what has to be done.
!          0: performs setup + solve + memory release, no initial guess
!         10: performs setup + solve + memory release, initial guess in x(1:n)
!          1: performs setup only
!             (preprocessing: prepares all parameters for subsequent solves)
!          2: solves only (based on previous setup), no initial guess
!         12: solves only (based on previous setup), initial guess in x(1:n)
!          3: the vector returned in x(1:n) is not the solution of the linear
!                 system, but the result of the action of the multigrid
!                 preconditioner on the right hand side in f(1:n)
!         -1: erases the setup and releases internal memory
!
!   IJOB == 100,110,101,102,112: same as, respectively, IJOB==0,10,1,2,12
!       but, use the TRANSPOSE of the input matrix in A, JA, IA.
!
!   !!! IJOB==2,3,12,102,112 require that one has previously called AGMG
!       with IJOB==1 or IJOB==101
!
!   !!! (change with respect to versions 2.x) !!!
!       The preconditioner defined when calling AGMG
!         with IJOB==1 or IJOB==101 is entirely kept in internal memory.
!       Hence the arrays A, JA and IA are not accessed upon subsequent calls
!         with IJOB==3.
!       Upon subsequent calls with IJOB==2,12,102,112, a matrix needs to
!            be supplied in arrays A, JA, IA, but it will be used to
!            perform matrix vector product within the main iterative
!            solution process (and only for this).
!            Hence the system is solved with this matrix which
!            may differ from the matrix in A, JA, IA that was supplied
!            upon the previous call with IJOB==1 or IJOB==101;
!            then AGMG attempts to solve a linear system with the "new"
!            matrix (supplied when IJOB==2,12,102 or 112) using the
!            preconditioner set up for the "old" one (supplied when
!            IJOB==1 or 101).
!         The same remarks apply to IJOB >= 100 or not: the value IJOB==1
!            or 101 determines whether the preconditioner set up and stored
!            in internal memory is based on the matrix or its transpose;
!            the value IJOB==2,12 or 102,112 is used to determine whether
!            the linear system to be solved is with the matrix or its
!            transpose, independently of the set up.
!            Hence one may set up a preconditioner for a matrix and use it
!            for its transpose.
!       These functionalities (set up a preconditioner and use it for another
!            matrix) are provided for the sake of generality but should be
!            used with care; in general, set up is fast with AGMG and hence
!            it is recommended to rerun it even if the matrix changes only
!            slightly.
!
! IPRINT   (input) INTEGER.
!              Indicates the unit number where information is to be printed
!              (N.B.: 5 is converted to 6). If nonpositive, only error
!              messages are printed on standard output.
!
! NREST    (input) INTEGER.
!             Restart parameter for GCR (an implementation of GMRES)
!             Nonpositive values are converted to NREST=10 (default)
!
! !!  If NREST==1, Flexible CG is used instead of GCR (when IJOB==0,10,2,
!             12,100,110,102,112) and also (IJOB==0,1,100,101) performs some
!             simplification during the set up based on the assumption
!             that the matrix supplied in A, JA, IA is symmetric (there is
!             then no more difference between IJOB==1 and IJOB==101).
!
! !!!  NREST=1 Should be used if and only if the matrix is really SYMMETRIC
! !!!         (and positive definite).
!
!  ITER    (input/output) INTEGER.
!          On input, the maximum number of iterations. Should be positive.
!          On output, actual number of iterations.
!            If this number of iteration was insufficient to meet convergence
!            criterion, ITER will be returned negative and equal to the
!            opposite of the number of iterations performed.
!          Significant only if IJOB==0, 2, 10, 12, 100, 102, 110, 112
!
!  TOL     (input) REAL (kind(0.0d0)).
!          The tolerance on residual norm. Iterations are stopped whenever
!               || A*x-f || <= TOL* || f ||
!          Should be positive and less than 1.0
!          Significant only if IJOB==0, 2, 10, 12, 100, 102, 110, 112
!
!!!!! Remark !!!! Except insufficient number of iterations to achieve
!                 convergence (characterized by a negative value returned
!                 in ITER), all other detected errors are fatal and lead
!                 to a STOP statement.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
    INTEGER, SAVE :: nza
    REAL(kind(0.0d0)), SAVE :: cputm=0.0d0,eltm=0.0d0,memh
    LOGICAL, SAVE :: preprocessed=.FALSE.,solve=.FALSE.
    INTEGER :: i,init,ijb,k
    REAL(kind(0.0d0)) :: cputmp,eltmp=0.0d0
    REAL(kind(0.0d0)) :: resid,wupd
    INTEGER, POINTER, DIMENSION(:) :: iw
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: ascal
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: w
    REAL(kind(0.0d0)) :: adum(1),fdum(1)
    INTEGER :: listrank(1),ifirstlistrank
!
    wfo=.TRUE.
    woo=IRANK<=0
    iout=iprint
    IF (iprint <= 0) THEN
       iout=6
       wfo=.FALSE.
       IF (iprint < 0) woo=.FALSE.
    ELSE IF (iprint == 5) THEN
       iout=6
    END IF
    ijb=mod(ijob,100)
!
    wff=wfo.AND.(IRANK<=0)
!
    IF (MOD(ijb,10) >= 2 .AND. .NOT.preprocessed) THEN
       WRITE (iout,1001) IRANK,ijob
       STOP
    END IF
!
    IF (ijb < 0 .AND. solve) GOTO 450
    IF (ijb < 0) GOTO 500
    trans=ijob.GE.100
    IF (MOD(ijb,10) >= 2) GOTO 300
    IF (preprocessed) THEN
       CALL dagmg_relmem
       IF (.NOT.allzero)                        &
               CALL dagmg_DIRseq(nn(nlev),fdum,-2)
       preprocessed=.FALSE.
       solve=.FALSE.
       eltm=0.0d0
       cputm=0.0d0
       flop=0.0d0
       kstat=0
    END IF
    CALL dagmg_mestime(-1,0.0d0,0.0d0)
    IF (HUGE(n) > 1.0e10) THEN
       rlenilen=dble(8)/dble(real_len)
    ELSE
       rlenilen=dble(4)/dble(real_len)
    END IF
    spd=nrest.EQ.1
    transint=trans.AND.(.NOT.spd)
    nza=ia(n+1)-ia(1)
    nlev=0
    imult=1
    IF (wfo) THEN
       WRITE(iout,900) IRANK
    END IF
    ALLOCATE(dt(1)%idiag(n+1),w(n),iw(n))
    CALL dagmg_partroword(n,a,ja,ia,dt(1)%idiag,w,iw)
    DEALLOCATE(w,iw)
    CALL dagmg_setupL1(n,a,ja,ia,listrank,ifirstlistrank)
    CALL dagmg_smoothsetup
    preprocessed=.TRUE.
    memh=memr+memi*rlenilen
!
    CALL dagmg_mestime(1,cputmp,eltmp)
    IF(wfo)THEN
       IF (nlev > 1) THEN
          WRITE(iout,960) IRANK, memax/nza, real_len    &
                        , memax*real_len/(2**20)
          IF(MOD(ijb,10) == 1) THEN
             WRITE(iout,961) IRANK, memh/nza, real_len  &
                           , memh*real_len/(2**20)
          END IF
       END IF
   !CPU_TIME: next line may be uncommented if implemented
       WRITE(iout,997) IRANK,eltmp
       WRITE(iout,'()')
    END IF
!
    IF (MOD(ijb,10) == 1) RETURN
    GOTO 310
!
!
300 CONTINUE
310 CONTINUE
    CALL dagmg_mestime(-2,0.0d0,0.0d0)
    resid=tol
    nrst=nrest
    init=MAX(IJB/10,0)
    IF (nrst <= 0) nrst=10
   !
    IF (MOD(ijb,10) >= 3) THEN
       IF (wfo) THEN
          WRITE(iout,901) IRANK
       END IF
       CALL dagmg_applyprec(n,f,x,a,ja,ia,MOD(ijb,10))
       CALL dagmg_mestime(2,cputmp,eltmp)
       cputm=cputm+cputmp
       eltm=eltm+eltmp
       solve=.TRUE.
       RETURN
   !
   !
    ELSE IF (nrst > 1) THEN
   !
       CALL dagmg_GCR(n,f,x,iter,resid,a,ja,ia,init)
   !
   !
    ELSE
   !
       CALL dagmg_FlexCG(n,f,x,iter,resid,a,ja,ia,init)
   !
   !
    END IF
!
    CALL dagmg_mestime(2,cputm,eltm)
    solve=.FALSE.
    IF (wfo) THEN
       IF (wff .AND. iter.NE.0) THEN
          DO i=2,nlev-1
             WRITE(iout,955) min(i,nlev),kstat(2,i-1),kstat(2,i),       &
                  dble(kstat(2,i))/dble(kstat(2,i-1)),kstat(1,i)
          END DO
       END IF
       WRITE(iout,'()')
       IF (ITER .NE. 0) THEN
         flop=flop*log(1.0d0/10)/log(RESID)
         WRITE(iout,952) IRANK,flop/dble(2*nza)
       END IF
       WRITE(iout,962) IRANK,(memh+mritr)/nza, real_len   &
                     , (memh+mritr)*real_len/(2**20)
   !CPU_TIME: next line may be uncommented if implemented
       WRITE(iout,999) IRANK,eltm
       WRITE(iout,'()')
    END IF
   IF (MOD(ijb,10) > 0) THEN
       solve=.FALSE.
       eltm=0.0d0
       cputm=0.0d0
       flop=0.0d0
       kstat=0
       RETURN
    ELSE
       GOTO 500
    END IF
450 CONTINUE
    IF (wfo) THEN
       WRITE(iout,'()')
       WRITE(iout,990) IRANK
       WRITE(iout,'()')
    END IF
    IF (wfo) THEN
       IF (wff .AND. iter.NE.0) THEN
          DO i=2,nlev-1
             WRITE(iout,955) i,kstat(2,i-1),kstat(2,i),         &
                  dble(kstat(2,i))/dble(kstat(2,i-1)),kstat(1,i)
          END DO
       END IF
       WRITE(iout,'()')
       WRITE(iout,953) IRANK,flop/dble(2*nza)
       WRITE(iout,963) IRANK,(memh+mritr)/nza, real_len   &
                     , (memh+mritr)*real_len/(2**20)
   !CPU_TIME: next line may be uncommented if implemented
       WRITE(iout,999) IRANK,eltm
       WRITE(iout,'()')
    END IF
!
!
!
500 CONTINUE
    CALL dagmg_relmem
    IF (.NOT.allzero)                              &
            CALL dagmg_DIRseq(nn(nlev),fdum,-2)
    preprocessed=.FALSE.
    solve=.FALSE.
    eltm=0.0d0
    cputm=0.0d0
    flop=0.0d0
    kstat=0
    IF (wfo) THEN
       WRITE (iout,902) IRANK
       WRITE (iout,903) IRANK
    END IF
!
    RETURN
900 FORMAT(i3,'*ENTERING AGMG **********************************',&
         '***************************')
901 FORMAT(i3,'*ONE APPLICATION OF AGMG PRECONDITIONER')
902 FORMAT(i3,            &
  ' (*) 1 work unit represents the cost of 1 (fine grid) residual evaluation ')
903 FORMAT(i3,'*LEAVING AGMG * (MEMORY RELEASED) ***************',&
         '***************************')
952 FORMAT(i3,'*','       Number of work units:',f9.2,             &
              ' per digit of accuracy (*)')
953 FORMAT(i3,'*',' Total number of work units:',f9.2, '(*)')
955 FORMAT('****     level',i2,'   #call=',i6,'   #cycle=',i6,    &
         '   mean=',f7.2,'    max=',i3)
960 FORMAT(  i3,'*','         memory used (peak):',f9.2,        &
         ' real(',i1,') words per nnz (',f8.2,' Mb)')
961 FORMAT(  i3,'*','     memory still allocated:',f9.2,        &
         ' real(',i1,') words per nnz (',f8.2,' Mb)')
962 FORMAT(  i3,'*','   memory used for solution:',f9.2,        &
         ' real(',i1,') words per nnz (',f8.2,' Mb)')
963 FORMAT(  i3,'*','                memory used:',f9.2,        &
         ' real(',i1,') words per nnz (',f8.2,' Mb)')
990 FORMAT(i3,'*GLOBAL STATISTICS for preconditioner application:')
996 FORMAT(i3,'*','           Setup time (CPU):   ',ES10.2,     &
         ' seconds')
997 FORMAT(i3,'*','       Setup time (Elapsed):   ',ES10.2,     &
         ' seconds')
998 FORMAT(i3,'*','        Solution time (CPU):   ',ES10.2,     &
         ' seconds')
999 FORMAT(i3,'*','    Solution time (Elapsed):   ',ES10.2,     &
         ' seconds')
1001 FORMAT(i3,'*',' FATAL ERROR: setup not done: ijob=',i3, &
         ' is not allowed')
1002 FORMAT(i3,'*',' FATAL ERROR: ijob=',i3, &
         ' (i.e. >= 100: work with transpose) is not allowed in the || case')
  END SUBROUTINE dagmg
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! END of MAIN DRIVER
!------------------ END of source file ---------------------------------
