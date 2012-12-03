!**************************************************************
!  Copyright Euratom-CEA
!  Authors : 
!     Virginie Grandgirard (virginie.grandgirard@cea.fr)
!     Chantal Passeron (chantal.passeron@cea.fr)
!     Guillaume Latu (guillaume.latu@cea.fr)
!     Xavier Garbet (xavier.garbet@cea.fr)
!     Philippe Ghendrih (philippe.ghendrih@cea.fr)
!     Yanick Sarazin (yanick.sarazin@cea.fr)
!  
!  This code GYSELA (for GYrokinetic SEmi-LAgrangian) 
!  is a 5D gyrokinetic global full-f code for simulating 
!  the plasma turbulence in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************
!-----------------------------------------------------------
! file : filter.f90
! date : 01/04/2010
!   Low-pass filter to avoid spurious gradients that 
!   develop during turbulent regime.
!   The transpose4D data structures are needed because
!   a transposition is needed to perform filtering.
!-----------------------------------------------------------
module filter_module
  use prec_const
  use globals, only : hffilterfref, hffiltertype, &
    Nr, Ntheta
  use MPIutils_module, only : f4D_transp
      
  implicit none
      
  ! for the equilibrium distribution function
  real(RKIND), dimension(:,:,:), pointer, private :: fmu_eq_priv
      
  ! Wavelet parameters
  integer, private, parameter :: Np = 4
  integer, private, parameter :: Nu = 2
  integer, private, parameter :: itcmax = 2
  real(8), private, parameter :: &
    filterPred(0:3) = (/ -0.0625, 0.5625, 0.5625, -0.0625 /)
  real(8), private, parameter :: &
    filterUpdate(0:1) = (/ .25, .25 /)
  integer, private, parameter :: INVALID = -1e8
      
  !******************************
  contains
  !******************************
  
  !--------------------------------------------------------
  ! Forward wavelet transform applied on fn
  !--------------------------------------------------------  
  subroutine FWT2D(fn, beg, dom)
    ! Distribution function
    real(RKIND), dimension (:,:),pointer  :: fn
    ! domain beginning indices
    integer, dimension (:), intent(in)   :: beg
    ! domain size
    integer, dimension (:), intent(in)   :: dom
      
    ! local variables
    integer     :: itc, itf
    integer     :: yb, ya, xb, xa, n
    integer     :: decpred, decup, bp, bu, ep, eu
    real(RKIND) :: tmp, val
      
    if (iand(itcmax, itcmax-1) .ne. 0) &
      STOP "FWT2D: False precondition itcmax is not a power of two"
    
    decpred = np/2-1
    decup   = nu/2-1
    bp      = -np/2+1
    ep      =  np/2
    bu      = -nu/2+1
    eu      =  nu/2
    itc     = 2
    itf     = 1
    do while (itc .le. itcmax)
      do xb = beg(1)+itf, beg(1)+dom(1)-1, itc
        do yb = beg(2), beg(2)+dom(2)-1, itf
          ya = yb
          tmp = fn(xb,yb)
          do n = bp, ep
            xa = beg(1) + &
              mod(xb + itc*n - itf + dom(1) - beg(1), dom(1)) 
            tmp = tmp - filterpred(n+decpred) * fn(xa,ya)
          end do
          fn(xb,yb) = tmp
        end do
      end do
      if (nu .ne. 0) then
        do xb = beg(1), beg(1)+dom(1)-1, itc
          do yb = beg(2), beg(2)+dom(2)-1, itf
            ya = yb
            tmp = fn(xb,yb)
            do n = bu, eu
              xa = beg(1) + &
                mod(xb + itc*n - itf + dom(1) - beg(1), dom(1)) 
              tmp = tmp + filterupdate(n+decup) * fn(xa,ya)
            end do
            fn(xb,yb) = tmp
          end do
        end do
      end if
      do xb = beg(1), beg(1)+dom(1)-1, itf
        xa = xb
        do yb = beg(2)+itf, beg(2)+dom(2)-1, itc
          tmp = fn(xb,yb)
          do n = bp, ep
            ya = beg(2) + &
              mod(yb + itc*n - itf + 2*dom(2) - beg(2), dom(2))
            tmp = tmp - filterpred(n+decpred) * fn(xa,ya)
          end do
          fn(xb,yb) = tmp
        end do
      end do
      if (nu .ne. 0) then
        do xb = beg(1), beg(1)+dom(1)-1, itf
          xa = xb
          do yb = beg(2), beg(2)+dom(2)-1, itc
            tmp = fn(xb,yb)
            do n = bu, eu
              ya = beg(2) + &
                mod(yb + itc*n - itf + 2*dom(2) - beg(2), dom(2))
              tmp = tmp + filterUpdate(n+decup) * fn(xa,ya)
            end do
            fn(xb,yb) = tmp
          end do
        end do
      end if
      itc = itc*2
      itf = itf*2
    enddo
  end subroutine FWT2D
      
  !--------------------------------------------------------
  ! Performs the thresholding of the wavelet transform 
  !--------------------------------------------------------  
  subroutine filter2D(fn, beg, dom, flooreps)
    ! distribution function (depending on variables [x, y])
    real(RKIND), dimension (:,:),pointer :: fn
    ! domain beginning indices
    integer, dimension (:), intent(in)   :: beg
    ! domain size
    integer, dimension (:), intent(in)   :: dom
    ! treshold for the adaptive representation
    real(RKIND), intent(in )             :: flooreps
    ! local variables
    real(RKIND) :: aux
    integer     :: itc, itf
    integer     :: yb, ya, xb, xa, n
    integer     :: decpred, decup, bp, bu, ep, eu
    real(RKIND) :: tmp, val
      
    if (iand(itcmax, itcmax-1) .ne. 0) &
      STOP "FWT2D: False precondition itcmax is not a power of two"
    itc = 0
    itf = 0
    do xb = beg(1), beg(1)+dom(1)-1, 1
      do yb = beg(2), beg(2)+dom(2)-1, 1
        if  ((mod(xb-beg(1),itcmax).eq.0) .and. &
          (mod(yb-beg(2),itcmax).eq.0)) then
          itc = itc+1
        else 
          if (abs(fn(xb,yb)) .lt. flooreps) then
            fn(xb,yb) = 0._RKIND
            itf       = itf+1
          end if
        end if
      end do
    end do
  end subroutine filter2D
      
  !--------------------------------------------------------
  ! Inverse wavelet transform applied on fn
  !--------------------------------------------------------  
  subroutine IWT2D(fn, beg, dom)
    ! distribution function
    real(RKIND), dimension (:,:), pointer :: fn
    ! domain begin indices
    integer, dimension (:), intent(in)    :: beg
    ! domain size
    integer, dimension (:), intent(in)    :: dom
    ! local variables
    integer    :: itc, itf
    integer    :: yb, ya, xb, xa, n
    integer    :: decpred, decup, bp, bu, ep, eu
    real(RKIND) :: tmp, val
      
    if (iand(itcmax, itcmax-1) .ne. 0) &
      STOP "compressBH: False precond itcmax is not a power of two"
      
    decpred = np/2-1
    decup   = nu/2-1
    bp      = -np/2+1
    ep      =  np/2
    bu      = -nu/2+1
    eu      =  nu/2
    itc     = itcmax
    itf     = itc/2
    do while (itc .gt. 1)
      if (nu .ne. 0) then
        do xb = beg(1), beg(1)+dom(1)-1, itf
          xa = xb
          do yb = beg(2), beg(2)+dom(2)-1, itc
            tmp = fn(xb,yb)
            do n = bu, eu
              ya = beg(2) + &
                mod(yb + itc*n - itf + dom(2) - beg(2), dom(2))
              tmp = tmp - filterUpdate(n+decup) * fn(xa,ya)
            end do
            fn(xb,yb) = tmp
          end do
        end do
      end if
      do xb = beg(1), beg(1)+dom(1)-1, itf
        xa = xb
        do yb = beg(2)+itf, beg(2)+dom(2)-1, itc
          tmp = fn(xb,yb)
          do n = bp, ep
            ya = beg(2) + &
              mod(yb + itc*n - itf + dom(2) - beg(2), dom(2))
            tmp = tmp + filterpred(n+decpred) * fn(xa,ya)
          end do
          fn(xb,yb) = tmp
        end do
      end do
      if (nu .ne. 0) then
        do xb = beg(1), beg(1)+dom(1)-1, itc
          do yb = beg(2), beg(2)+dom(2)-1, itf
            ya = yb
            tmp = fn(xb,yb)
            do n = bu, eu
              xa = beg(1) + &
                mod(xb + itc*n - itf + dom(1) - beg(1), dom(1)) 
              tmp = tmp - filterupdate(n+decup) * fn(xa,ya)
            end do
            fn(xb,yb) = tmp
          end do
        end do
      end if
      do xb = beg(1)+itf, beg(1)+dom(1)-1, itc
        do yb = beg(2), beg(2)+dom(2)-1, itf
          ya = yb
          tmp = fn(xb,yb)
          do n = bp, ep
            xa = beg(1) + &
              mod(xb + itc*n - itf + dom(1) - beg(1), dom(1)) 
            tmp = tmp + filterpred(n+decpred) * fn(xa,ya)
          end do
          fn(xb,yb) = tmp
        end do
      end do
      itc = itc/2
      itf = itf/2
    enddo
  end subroutine IWT2D
      
  !--------------------------------------------------------
  !  Definition of the size of the structures which will
  !  be used. This size depends on the reference function
  !  that has been chosen, among:
  !   1: f-feq
  !   2: f-<f>_theta
  !   3: f is duplicated in r (to be periodic)
  !   4: f-<f>_theta with zeroing of the function at rmin
  !      and rmax
  !--------------------------------------------------------  
  subroutine setDomainSize(beg, nbeg, dom, dec)
    integer, intent(out) :: dom(1:2)
    integer, intent(out) :: dec(1:2)
    integer, intent(out) :: beg(1:2)
    integer, intent(out) :: nbeg(1:2)
      
    select case (hffiltertype) 
    case (1)
      ! FFT filters
      select case (hffilterfref) 
      case(1,2,4)
        beg(1:2)  = (/ 0, 0 /)
        nbeg(1:2) = (/ 1, 1 /)
        dom(1:2)  = (/ (Nr+1), Ntheta /)
        dec(1:2)  = (/ 0, 0 /)
      case(3)
        beg(1:2)  = (/ 0, 0 /)
        nbeg(1:2) = (/ 1, 1 /)
        dom(1:2)  = (/ 2*(Nr+1), Ntheta /)
        dec(1:2)  = (/ 0, 0 /)
      end select
    case(2)
      ! Wavelet filters
      select case (hffilterfref) 
      case(1,2,4)
        beg(1:2)  = (/ 0, 0 /)
        nbeg(1:2) = (/ 1, 1 /)
        dom(1:2)  = (/ Nr+1, Ntheta /)
        dec(1:2)  = (/ 0, 0 /)
      case(3)
        beg(1:2)  = (/ 0, 0 /)
        nbeg(1:2) = (/ 1, 1 /)
        dom(1:2)  = (/ 2*(Nr+1), Ntheta /)
        dec(1:2)  = (/ 0, 0 /)
      end select
    end select
  end subroutine setDomainSize
      
  !--------------------------------------------------------
  !  Copy f4Dtransp into a 2D buffer in order to apply
  !  the forward transform, according to the choice of
  !  the reference function:
  !   1: f-feq
  !   2: f-<f>_theta
  !   3: f is duplicated in r (to be periodic)
  !   4: f-<f>_theta with zeroing of the function at rmin
  !      and rmax
  !--------------------------------------------------------  
  subroutine copyforw2D(nfn, sum1, iphi, ivpar)
    ! Distribution functions (depending on variables [x, y])
    real(RKIND), dimension (:,:), pointer :: nfn
    real(RKIND), dimension (:),   pointer :: sum1
    integer, intent(in)                   :: iphi, ivpar
      
    real(RKIND) :: tmp
    integer     :: i, j
    integer     :: xa, xb, ya, yb, xc, frac8
    integer     :: dom(1:2), dec(1:2), beg(1:2), nbeg(1:2)
      
    call setDomainSize(beg,nbeg,dom,dec)
    frac8 = dom(1)/16
    if ((hffilterfref .eq. 2) .or. (hffilterfref .eq. 4)) then
      do i = 0, dom(1)-1
        xb = nbeg(1)+dec(1)
        tmp = 0._RKIND
        do j = 0, dom(2)-1
          tmp = tmp + f4D_transp(beg(1)+i,beg(2)+j,iphi,ivpar)
        end do
        sum1(xb+i) = tmp/real(dom(2))
      end do
    end if
    do j = 0, dom(2)-1
      ya = j+beg(2)
      yb = j+nbeg(2)+dec(2)
      xa = beg(1)
      xb = nbeg(1)+dec(1)
      xc = nbeg(1)+dec(1)+dom(1)-1
      select case (hffilterfref) 
      case (1)
        do i = 0, dom(1)-1 
          nfn(xb+i,yb) = f4D_transp(xa+i,ya,iphi,ivpar) - &
            fmu_eq_priv(xa+i,ya,ivpar)
        end do
      case (2)
        do i = 0, dom(1)-1 
          nfn(xb+i,yb) = f4D_transp(xa+i,ya,iphi,ivpar) - &
            sum1(xb+i)
        end do
      case (3)
        do i = 0, dom(1)/2-1 
          tmp = f4D_transp(xa+i,ya,iphi,ivpar)
          nfn(xb+i,yb) = tmp
          nfn(xc-i,yb) = tmp
        end do
      case (4)
        do i = 0, frac8-1 
          nfn(xb+i,yb) = 0._RKIND
        end do
        do i = frac8, dom(1)-frac8-1
          nfn(xb+i,yb) = f4D_transp(xa+i,ya,iphi,ivpar) - &
            sum1(xb+i)
        end do
        do i = dom(1)-frac8, dom(1)-1
          nfn(xb+i,yb) = 0._RKIND
        end do
      end select
    end do
  end subroutine copyforw2D
      
  !--------------------------------------------------------
  !  Copy the 2D buffer into f4Dtransp in order to apply
  !  the inverse transform, according to the choice of
  !  the reference function:
  !   1: f-feq
  !   2: f-<f>_theta
  !   3: f is duplicated in r (to be periodic)
  !   4: f-<f>_theta with zeroing of the function at rmin
  !      and rmax
  !--------------------------------------------------------  
  subroutine copyback2D(nfn, sum1, iphi, ivpar)
    ! Distribution functions (depending on variables [x, y])
    real(RKIND), dimension (:,:), pointer :: nfn
    real(RKIND), dimension (:),   pointer :: sum1
    integer, intent(in)                   :: iphi, ivpar
      
    real(RKIND) :: tmp
    integer     :: i, j, frac8
    integer     :: xa, xb, ya, yb, xc
    integer     :: dom(1:2), dec(1:2), beg(1:2), nbeg(1:2)
      
    call setDomainSize(beg,nbeg,dom,dec)
    do j = 0, dom(2)-1
      ya = j+beg(2)
      yb = j+nbeg(2)+dec(2)
      xa = beg(1)
      xb = nbeg(1)+dec(1)
      xc = nbeg(1)+dec(1)+dom(1)-1
      select case (hffilterfref) 
      case (1)
        do i = 0, dom(1)-1 
          f4D_transp(xa+i,ya,iphi,ivpar) = nfn(xb+i,yb) + &
            fmu_eq_priv(xa+i,ya,ivpar)
        end do
      case (2,4)
        do i = 0, dom(1)-1 
          f4D_transp(xa+i,ya,iphi,ivpar) = nfn(xb+i,yb) + &
            sum1(xb+i)
        end do
      case (3)
        do i = 0, dom(1)/2-1 
          f4D_transp(xa+i,ya,iphi,ivpar) = 0.5_RKIND * &
            (nfn(xb+i,yb)+nfn(xc-i,yb))
        end do
      end select
    end do
  end subroutine Copyback2D
      
  !--------------------------------------------------------
  ! Applying of the wavelet filter according to the
  !  frequency 'hffilterfreq' and the reference function 
  !  'hffilterfref'
  !--------------------------------------------------------  
  subroutine waveletfilter()
    use globals
    use MPIutils_module
    use OMPutils_module
    use fft_module
      
    integer     :: ir, itheta, iphi, ivpar, sign, tid
    real(RKIND) :: fact
    real(RKIND), dimension(:,:), pointer :: pwav
    real(RKIND), dimension(:)  , pointer :: psum1
    integer :: beg(1:2), nbeg(1:2)
    integer :: dom(1:2), dec(1:2)
      
    call setDomainSize(beg,nbeg,dom,dec)
      
#ifdef _OPENMP
!$OMP PARALLEL private(tid,ir,itheta,ivpar,iphi,&
!$OMP      pwav,psum1) default(shared) 
!$OMP BARRIER
    tid = 1+omp_get_thread_num()
#else
    tid = 1
#endif
    pwav  => Romp1_1Nrpb_1Nthetapb(tid)%val
    psum1 => Romp1_1Nthetap1(tid)%val
!$OMP BARRIER
    do ivpar = lstart,lend
!$OMP DO SCHEDULE(static) 
      do iphi = kstart,kend
        call copyforw2D(pwav, psum1, iphi, ivpar)
        call FWT2D(pwav, nbeg, dom)
        call filter2D(pwav, nbeg, dom, 1.e30_RKIND)
        call IWT2D(pwav, nbeg, dom)
        call copyback2D(pwav, psum1, iphi, ivpar)
      end do
!$OMP END DO
!$OMP BARRIER
    end do
!$OMP END PARALLEL
  end subroutine waveletfilter
      
  !--------------------------------------------------------
  ! Applying of the FFT filter according to the
  !  frequency 'hffilterfreq' and the reference function 
  !  'hffilterfref'
  !--------------------------------------------------------  
  subroutine fftfilter()
    use globals
    use MPIutils_module
    use OMPutils_module
    use fft_module
    
    integer     :: ir, itheta, iphi, ivpar, isign
    real(RKIND) :: fact
    complex(CKIND), dimension(:,:), pointer :: pthc, prc
    real(RKIND)   , dimension(:,:), pointer :: prth
    real(RKIND)   , dimension(:)  , pointer :: psum1
    integer :: beg(1:2), nbeg(1:2)
    integer :: dom(1:2), dec(1:2)
    integer :: tid
      
    call setDomainSize(beg,nbeg,dom,dec)
      
#ifdef _OPENMP
!$OMP PARALLEL private(tid,ir,itheta,isign,ivpar,iphi,&
!$OMP      fact,pthc,prc,prth,psum1) default(shared) 
!$OMP BARRIER
    tid = 1+omp_get_thread_num()
#else
    tid = 1
#endif
    psum1 => Romp1_1Nthetap1(tid)%val
    prth  => Romp1_1Nrpb_1Nthetapb(tid)%val
    if (hffilterfref .eq. 3) then
      prc  => Comp1_1Nrpb_1Nthetapb(tid)%val
      pthc => Comp1_1Nthetapb_1Nrpb(tid)%val
    else
      prc  => Comp1_1Nrp1_1Ntheta(tid)%val
      pthc => Comp1_1Ntheta_1Nrp1(tid)%val
    end if
    do ivpar = lstart,lend
!$OMP DO SCHEDULE(static) 
      do iphi = kstart,kend
        call copyforw2D(prth, psum1, iphi, ivpar)
        fact = 1._RKIND/real(dom(1)*dom(2))
        do ir = 0, dom(1)-1
          do itheta = 0, dom(2)-1
            pthc(itheta+1,ir+1) = &
              cmplx(prth(1+ir,1+itheta),0.,CKIND)
          end do
        end do
        isign = 1
        call fourrow(pthc,isign)
        do itheta = 0, dom(2)-1
          do ir = 0, dom(1)-1
            prc(ir+1,itheta+1) = pthc(itheta+1,ir+1)
          end do
        end do
        call fourrow(prc,isign)
        do ir = 0, dom(1)-1
          do itheta = 0, dom(2)-1
            if (((itheta .gt. nint(dom(2)/4.)) .and. &
              (itheta .le. nint((3.*dom(2))/4.))) .or. &
              ((ir .gt. nint((dom(1)-1)/4.))) .and. &
              (ir .le. nint((3.*(dom(1)-1))/4.))) then
              prc(ir+1,itheta+1) = 0._CKIND
            else
              prc(ir+1,itheta+1) = fact * prc(ir+1,itheta+1)
            end if
          end do
        end do
        isign = -1
        call fourrow(prc,isign)
        do itheta = 0, dom(2)-1
          do ir = 0, dom(1)-1
            pthc(itheta+1,ir+1) = prc(ir+1,itheta+1)
          end do
        end do
        call fourrow(pthc,isign)
        do ir = 0, dom(1)-1
          do itheta = 0, dom(2)-1
            prth(1+ir,1+itheta) = real(pthc(itheta+1,ir+1),RKIND)
          end do
        end do
        call copyback2D(prth, psum1, iphi, ivpar)
      end do
!$OMP END DO
      isign =1
!$OMP BARRIER
    end do
!$OMP END PARALLEL
  end subroutine fftfilter
      
  !------------------------------------------------------------
  ! Define low-pass filters (FFT and wavelet filters).
  ! The cut-off frequency is nearby two times the cell length.
  !------------------------------------------------------------
  subroutine hffilter_loc(Sfmu_eq)
    use globals
    use fdistribu5d_class
    use MPIutils_module
    real(RKIND), dimension(:,:,:), pointer :: Sfmu_eq
      
    logical :: activatefilter
      
    fmu_eq_priv => Sfmu_eq
      
    if ((hffilterfreq .ne. 0) .and. &
      (hffiltertype .le. 2) .and. (hffiltertype .ge. 1)) then
      activatefilter = (mod(iter_glob,hffilterfreq) .eq. 0)
    else
      activatefilter = .false.
    end if
      
    if (activatefilter) then
      select case (hffiltertype) 
      case(1)
        call fftfilter()
      case(2)
        call waveletfilter()
      end select
    end if
  end subroutine hffilter_loc
      
  !-------------------------------------------------------------
  ! Define low-pass filters (FFT and wavelet filters).
  ! The cut-off frequency is nearby two times the cell length.
  !-------------------------------------------------------------
  subroutine hffilter(f,Sfmu_eq)
    use globals
    use fdistribu5d_class
    use MPIutils_module
    type(fdistribu5d)            , intent(inout) :: f
    real(RKIND), dimension(:,:,:), pointer       :: Sfmu_eq
      
    logical :: activatefilter
      
    fmu_eq_priv => Sfmu_eq
      
    if ((hffilterfreq .ne. 0) .and. &
      (hffiltertype .le. 2) .and. (hffiltertype .ge. 1)) then
      activatefilter = (mod(iter_glob,hffilterfreq) .eq. 0)
    else
      activatefilter = .false.
    end if
      
    if (activatefilter) then
      call ppbarrier()
      call pptranspose_forward(f%values)   
      select case (hffiltertype) 
      case(1)
        call fftfilter()
        if (pglobal_id .eq. 0) print *,"----> FFT Filter", &
          iter_glob
      case(2)
        call waveletfilter()
        if (pglobal_id .eq. 0) print *,"----> WAV Filter", &
          iter_glob
      end select
      call ppbarrier()
      call pptranspose_backward(f%values)   
    end if
 end subroutine hffilter
      
end module filter_module
