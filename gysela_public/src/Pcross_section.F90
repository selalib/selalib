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
      
!-------------------------------------------------
! file : Pcross_section.f90
! date : 18/04/2005
!  used for saving cross_sections of the 
!  electric potential and the
!  distribution function in all directions
!------------------------------------------------
module Pcross_section_module
  use prec_const
  use mem_alloc_module
  use fdistribu5d_class
      
  implicit none
  
  !-------------------------------------
  !  PRIVATE VARIABLES
  !-------------------------------------
  real(RKIND), dimension(:,:), pointer :: frtheta_loc_tmp
  real(RKIND), dimension(:,:), pointer :: frtheta_glob_tmp
  real(RKIND), dimension(:,:), pointer :: frvpar_loc_tmp
  real(RKIND), dimension(:,:), pointer :: frvpar_glob_tmp
  real(RKIND), dimension(:,:), pointer :: fthetavpar_loc_tmp
  real(RKIND), dimension(:,:), pointer :: fthetavpar_glob_tmp
      
  !-------------------------------------
  !  PUBLIC VARIABLES
  !-------------------------------------
  public
  ! -> Electric Potential in 2D for output saving 
  real(RKIND), dimension(:,:), pointer :: Phirth_diag
  real(RKIND), dimension(:,:), pointer :: Phirphi_diag
  real(RKIND), dimension(:,:), pointer :: Phithphi_diag
  ! -> distribution function in 2D for output saving
  real(RKIND), dimension(:,:), pointer :: frtheta_v0_diag
  real(RKIND), dimension(:,:), pointer :: frtheta_vmax_diag
  real(RKIND), dimension(:,:), pointer :: fphivpar_mu0_diag
  real(RKIND), dimension(:,:), pointer :: fphivpar_mumax_diag
  real(RKIND), dimension(:,:), pointer :: frvpar_mu0_diag
  real(RKIND), dimension(:,:), pointer :: frvpar_mumax_diag
  real(RKIND), dimension(:,:), pointer :: fthetavpar_mu0_diag
  real(RKIND), dimension(:,:), pointer :: fthetavpar_mumax_diag
  real(RKIND), dimension(:,:), pointer :: fvparmu_diag
      
  !******************************
  contains
  !******************************
      
  !-------------------------------------------------
  ! Constructor 
  !-------------------------------------------------
  subroutine new_Pcross_section
    use globals, only : istart, iend, jstart, jend, &
      Nbproc_r, Nbproc_theta, dom_r, dom_theta, &
      Nr, Ntheta, Nphi, Nvpar, Nmu, &
      plocal_id, istart_id, jstart_id, &
      cible, cible_mu, cible_istart, cible_jstart
      
!bmodif
!VG!    !*** for f(r,theta) cross-sections ***
!VG!    call glob_allocate(frtheta_loc_tmp,istart,iend,jstart,jend)
!VG!    if (plocal_id.eq.cible_mu) then
!VG!      call glob_allocate(frtheta_glob_tmp,0,dom_r-1,0,&
!VG!        Nbproc_r*(Ntheta+1),'frtheta_glob_tmp')
!VG!    else
!VG!      call glob_allocate(frtheta_glob_tmp,0,1,0,1, &
!VG!        'frtheta_glob_tmp')
!VG!    end if
!VG!    !*** for f(r,vpar) cross-sections ***
!VG!    call glob_allocate(frvpar_loc_tmp,istart,iend,0,Nvpar)
!VG!    if (jstart_id.eq.cible_jstart) then
!VG!      call glob_allocate(frvpar_glob_tmp,0,dom_r-1,0,&
!VG!        Nbproc_r*(Nvpar+1),'frvpar_glob_tmp')
!VG!    else
!VG!      call glob_allocate(frvpar_glob_tmp,0,1,0,&
!VG!        1,'frvpar_glob_tmp')
!VG!    end if
    ! -> solution to avoid temporary the problem 
    !    of MPI_GATHER instruction
    !*** for f(r,theta) cross-sections ***
    call glob_allocate(frtheta_loc_tmp,0,dom_r-1, &
      0,dom_theta-1,'frtheta_loc_tmp')
    if (plocal_id.eq.cible_mu) then
      call glob_allocate(frtheta_glob_tmp,0,dom_r-1,0,&
        Nbproc_r*(Ntheta+1),'frtheta_glob_tmp')
    else
      call glob_allocate(frtheta_glob_tmp,0,1,0,1, &
        'frtheta_glob_tmp')
    end if
    !*** for f(r,vpar) cross-sections ***
    call glob_allocate(frvpar_loc_tmp,0,dom_r-1,0,Nvpar, &
      'frvpar_loc_tmp')
    if (jstart_id.eq.cible_jstart) then
      call glob_allocate(frvpar_glob_tmp,0,dom_r-1,0, &
        Nbproc_r*(Nvpar+1),'frvpar_glob_tmp')
    else
      call glob_allocate(frvpar_glob_tmp,0,1,0, &
        1,'frvpar_glob_tmp')
    end if
    !*** for f(theta,vpar) cross-sections ***
    call glob_allocate(fthetavpar_loc_tmp,0,dom_theta-1,0,Nvpar, &
      'fthetavpar_loc_tmp')
    if (istart_id.eq.cible_istart) then
      call glob_allocate(fthetavpar_glob_tmp,0,dom_theta-1,0, &
        Nbproc_theta*(Nvpar+1),'fthetavpar_glob_tmp')
    else
      call glob_allocate(fthetavpar_glob_tmp,0,1,0, &
        1,'fthetavpar_glob_tmp')
    end if
!emodif
      
    !*** allocation of arrays for output saving ***
    !-> for Phi2D saving 
    call glob_allocate(Phirth_diag,0,Nr,0,Ntheta,'Phirth_diag')
    call glob_allocate(Phirphi_diag,0,Nr,0,Nphi,'Phirphi_diag')
    call glob_allocate(Phithphi_diag,0,Ntheta,0,Nphi, &
      'Phithphi_diag')
    !-> for f2D saving 
    call glob_allocate(frtheta_v0_diag,0,Nr,0,Ntheta, &
      'frtheta_v0_diag')
    call glob_allocate(frtheta_vmax_diag,0,Nr,0,Ntheta, &
      'frtheta_vmax_diag')
    call glob_allocate(fphivpar_mu0_diag,0,Nphi,0,Nvpar, &
      'fphivpar_mu0_diag')
    call glob_allocate(fphivpar_mumax_diag,0,Nphi,0,Nvpar, &
      'fphivpar_mumax_diag')
    call glob_allocate(frvpar_mu0_diag,0,Nr,0,Nvpar, &
      'frvpar_mu0_diag')
    call glob_allocate(frvpar_mumax_diag,0,Nr,0,Nvpar, &
      'frvpar_mumax_diag')
    call glob_allocate(fthetavpar_mu0_diag,0,Ntheta,0,Nvpar, &
      'fthetavpar_mu0_diag')
    call glob_allocate(fthetavpar_mumax_diag,0,Ntheta,0,Nvpar, &
      'fthetavpar_mumax_diag')
    call glob_allocate(fvparmu_diag,0,Nvpar,0,Nmu,'fvparmu_diag')
  end subroutine new_Pcross_section
      
  !-------------------------------------------------
  ! Destructor 
  !-------------------------------------------------
  subroutine delete_Pcross_section
    use globals, only : plocal_id, cible
      
    call glob_deallocate(frtheta_loc_tmp)
    call glob_deallocate(frtheta_glob_tmp)
    call glob_deallocate(frvpar_loc_tmp) 
    call glob_deallocate(frvpar_glob_tmp)
    call glob_deallocate(fthetavpar_loc_tmp) 
    call glob_deallocate(fthetavpar_glob_tmp)
      
    !*** deallocation of arrays for output saving ***
    !-> for Phi2D saving 
    call glob_deallocate(Phirth_diag)
    call glob_deallocate(Phirphi_diag)
    call glob_deallocate(Phithphi_diag)
    !->  for f2D saving 
    call glob_deallocate(frtheta_v0_diag)
    call glob_deallocate(frtheta_vmax_diag)
    call glob_deallocate(fphivpar_mu0_diag)
    call glob_deallocate(fphivpar_mumax_diag)
    call glob_deallocate(frvpar_mu0_diag)
    call glob_deallocate(frvpar_mumax_diag)
    call glob_deallocate(fthetavpar_mu0_diag)
    call glob_deallocate(fthetavpar_mumax_diag)
    call glob_deallocate(fvparmu_diag)
  end subroutine delete_Pcross_section
      
  !---------------------------------------------- 
  ! computation of the 2D cross-sections of Phi
  !---------------------------------------------- 
  subroutine compute_Phi2DCS(geom, poiss)
    use globals, only : ir_Phi, itheta_Phi, iphi_Phi
    use poisson_class
    type(geometry), intent(in)    :: geom
    type(poisson) , intent(inout) :: poiss
      
    integer :: ir, itheta, iphi
      
    !*** cross-sections of the electric potential ***
    ! placer dans une routine compute_Phi2DCS
    ir_Phi      = int(geom%Nr/2)
    itheta_Phi  = 0
    iphi_Phi    = 0
    do itheta = 0,geom%Ntheta
      do ir = 0,geom%Nr
        Phirth_diag(ir,itheta) = poiss%Phi(ir,itheta,iphi_Phi)
      end do
    end do
    do iphi = 0,geom%Nphi
      do ir = 0,geom%Nr
        Phirphi_diag(ir,iphi) = poiss%Phi(ir,itheta_Phi,iphi)
      end do
    end do
    do iphi = 0,geom%Nphi
      do itheta = 0,geom%Ntheta
        Phithphi_diag(itheta,iphi) = poiss%Phi(ir_Phi,itheta,iphi)
      end do
    end do
  end subroutine compute_Phi2DCS
      
  !------------------------------------------------------------- 
  ! distribution function computation in r and theta
  !  directions f for fixed values of phi and vpar
  !  ( i.e for iphi=phig(iphi_fix) and ivpar=vparg(ivpar_fix) )
  !------------------------------------------------------------- 
  subroutine compute_frtheta(f,iphi_fix,ivpar_fix,frtheta)
    use globals, only : cible_mu, istart, iend, jstart, jend, &
      Nbproc_r, Nbproc_theta, dom_r, dom_theta, Nr, Ntheta, &
      plocal_id, mpi_comm_mu
    use MPIutils_module
    use interpolation_module
    type(fdistribu5d)            , intent(in)    :: f
    integer                      , intent(in)    :: iphi_fix
    integer                      , intent(in)    :: ivpar_fix
    real(RKIND), dimension(0:,0:), intent(inout) :: frtheta
    
    integer :: ierr
    integer :: ir, itheta
    integer :: p_r, p_theta
    integer :: nbelts
      
!bmodif
!VG!    do itheta = f%jstart,f%jend
!VG!      do ir = f%istart,f%iend
!VG!        frtheta_loc_tmp(ir,itheta) = &
!VG!          f%values(ir,itheta,iphi_fix,ivpar_fix)
!VG!      enddo
!VG!    enddo
!VG!    call MPI_GATHER(frtheta_loc_tmp,1, &
!VG!      mpitype_globbloc,frtheta_glob_tmp,1, &
!VG!      mpitype_globbloc,cible_mu,mpi_comm_mu,ierr)
      
    ! -> solution to avoid temporary the problem 
    !   of MPI_GATHER instruction
    do itheta = f%jstart,f%jend
      do ir = f%istart,f%iend
        frtheta_loc_tmp(ir-f%istart,itheta-f%jstart) = &
          f%values(ir,itheta,iphi_fix,ivpar_fix)
      enddo
    enddo
    nbelts = dom_r * dom_theta
    call comm_gather2D(frtheta_loc_tmp,nbelts,0,dom_theta, &
      cible_mu,frtheta_glob_tmp,mpi_comm_mu) 
!emodif
      
    if (plocal_id .eq. cible_mu) then
      do p_theta = 0,Nbproc_theta-1
        do p_r = 0,Nbproc_r-1
          do itheta = 0,dom_theta -1
            do ir = 0,dom_r-1
              frtheta(p_r*dom_r+ir,p_theta*dom_theta+itheta) = &
                frtheta_glob_tmp(ir, &
                (p_theta*Nbproc_r+p_r)*dom_theta+itheta)
            enddo
          enddo
        enddo
      enddo
      !*** periodic condition in theta ***
      do ir = 0,f%n1
        frtheta(ir,f%n2) = frtheta(ir,0)
      end do
    end if
  end subroutine compute_frtheta
  
  
  !----------------------------------------------------------- 
  ! distribution function computation in phi, v parallel 
  !  directions for fixed values of r and theta
  !  (i.e. for r=rg(ir_fix), theta=thetag(itheta_fix))
  ! Rk : it supposes that the processor contains the 
  !      point (ir_fix,itheta_fix)
  !-----------------------------------------------------------
  subroutine compute_fphivpar(f,ir_fix,itheta_fix,fphivpar)
    type(fdistribu5d)            , intent(in)    :: f
    integer                      , intent(in)    :: ir_fix
    integer                      , intent(in)    :: itheta_fix
    real(RKIND), dimension(0:,0:), intent(inout) :: fphivpar
      
    integer :: iphi, ivpar
      
    fphivpar = 0._RKIND
    do ivpar = 0,f%n4
      do iphi = 0,f%n3-1
        fphivpar(iphi,ivpar) = &
          f%values(ir_fix,itheta_fix,iphi,ivpar)
      enddo
    enddo
    fphivpar(f%n3,0:f%n4) = fphivpar(0,0:f%n4)
  end subroutine compute_fphivpar
      
  !-----------------------------------------------------------
  ! distribution function computation in v parallel, mu
  !  directions for fixed values of r, theta and phi
  !  (i.e. for r=rg(ir_fix), theta=thetag(itheta_fix) and
  !   phi=phig(iphi_fix))
  ! Rk : it supposes that the processor contains the 
  !      point (ir_fix,itheta_fix)
  !-----------------------------------------------------------
  subroutine compute_fvparmu(f,ir_fix,itheta_fix,iphi_fix,fvparmu)
    use globals, only : istart, iend, jstart, jend, &
      cible, mu_id, Nvpar, Nmu
    include "mpiperso.h"
    type(fdistribu5d)            , intent(in)    :: f
    integer                      , intent(in)    :: ir_fix
    integer                      , intent(in)    :: itheta_fix
    integer                      , intent(in)    :: iphi_fix
    real(RKIND), dimension(0:,0:), intent(inout) :: fvparmu
    
    integer :: ierr
    integer :: ivpar, nb_send, imu, nb_reqr, itag_vparmu
      
    real(RKIND), dimension(0:Nvpar)               :: fvpar_tmp
    integer    , dimension(MPI_STATUS_SIZE,0:Nmu) :: status
    integer    , dimension(0:Nmu)                 :: reqr
    
    nb_send = Nvpar+1
    nb_reqr = Nmu+1
    if (pglobal_id .eq. cible) then
      do imu = 0,Nmu
        itag_vparmu = 8000 + imu
        call MPI_IRECV(fvparmu(0,imu),nb_send, &
          MPI_REAL8,MPI_ANY_SOURCE,itag_vparmu, &
          MPI_COMM_WORLD,reqr(imu),ierr)
      enddo
    endif
      
    if ( (ir_fix.ge.istart).and.(ir_fix.le.iend).and. &
      (itheta_fix.ge.jstart).and.(itheta_fix.le.jend) ) then
      do ivpar = 0,Nvpar
        fvpar_tmp(ivpar) = &
          f%values(ir_fix,itheta_fix,iphi_fix,ivpar)
      end do
      itag_vparmu = 8000 + mu_id
      call MPI_SSEND(fvpar_tmp,nb_send,MPI_REAL8,cible, &
        itag_vparmu,MPI_COMM_WORLD,ierr)
    end if
    if (pglobal_id.eq.cible) &
      call MPI_WAITALL(Nmu+1,reqr(0),status(1,0),ierr)
  end subroutine compute_fvparmu
      
  !----------------------------------------------------------- 
  ! distribution function computation in r, v parallel 
  !  directions for fixed values of theta and phi
  !  (i.e theta=thetag(itheta_fix), phi=phig(iphi_fix)
  ! Rk : it supposes that the processor contains the 
  !      point ir_fix
  !-----------------------------------------------------------
  subroutine compute_frvpar(f,itheta_fix,iphi_fix,frvpar)
    use globals, only : istart, iend, jstart, jend, &
      jstart_id, cible_jstart, Nbproc_r, mpi_comm_column, dom_r
    use MPIutils_module
    include "mpiperso.h"
    type(fdistribu5d)            , intent(in)    :: f
    integer                      , intent(in)    :: itheta_fix
    integer                      , intent(in)    :: iphi_fix
    real(RKIND), dimension(0:,0:), intent(inout) :: frvpar
      
    integer :: ierr
    integer :: ir, ivpar, p_r
    integer :: nbelts
      
!bmodif
!VG!    do ivpar = 0,f%n4
!VG!      do ir = f%istart,f%iend
!VG!        frvpar_loc_tmp(ir,ivpar) = &
!VG!          f%values(ir,itheta_fix,iphi_fix,ivpar)
!VG!      end do
!VG!    end do
!VG!    call MPI_GATHER(frvpar_loc_tmp,1, &
!VG!      mpitype_globbloc_rvpar,frvpar_glob_tmp,1, &
!VG!      mpitype_globbloc_rvpar,cible_jstart,mpi_comm_column,ierr)
    ! -> solution to avoid temporary the problem 
    !of MPI_GATHER instruction
    do ivpar = 0,f%n4
      do ir = f%istart,f%iend
        frvpar_loc_tmp(ir-f%istart,ivpar) = &
          f%values(ir,itheta_fix,iphi_fix,ivpar)
      end do
    end do
    nbelts = dom_r*(f%n4+1)
    call comm_gather2D(frvpar_loc_tmp,nbelts,0,f%n4+1, &
      cible_jstart,frvpar_glob_tmp,mpi_comm_column) 
!emodif
      
    if (jstart_id .eq. cible_jstart) then
      do ivpar = 0,f%n4
        do p_r = 0,Nbproc_r-1
          do ir = 0,dom_r-1
            frvpar(p_r*dom_r+ir,ivpar) = &
              frvpar_glob_tmp(ir,p_r*(f%n4+1)+ivpar)
          enddo
        enddo
      enddo
    end if
  end subroutine compute_frvpar
      
  !----------------------------------------------------------- 
  ! distribution function computation in theta, v parallel 
  !  directions for fixed values of r and phi
  !  (i.e r=rg(ir_fix), phi=phig(iphi_fix))
  ! Rk : it supposes that the processor contains the 
  !      point itheta_fix
  !-----------------------------------------------------------
  subroutine compute_fthetavpar(f,ir_fix,iphi_fix,fthetavpar)
    use globals, only : istart, iend, jstart, jend, &
      istart_id, cible_istart, Nbproc_theta, &
      mpi_comm_row, dom_theta
    use MPIutils_module
    include "mpiperso.h"
    type(fdistribu5d)            , intent(in)    :: f
    integer                      , intent(in)    :: ir_fix
    integer                      , intent(in)    :: iphi_fix
    real(RKIND), dimension(0:,0:), intent(inout) :: fthetavpar
      
    integer :: ierr
    integer :: itheta, ivpar, p_theta
    integer :: nbelts
      
    do ivpar = 0,f%n4
      do itheta = f%jstart,f%jend
        fthetavpar_loc_tmp(itheta-f%jstart,ivpar) = &
          f%values(ir_fix,itheta,iphi_fix,ivpar)
      end do
    end do
    nbelts = dom_theta*(f%n4+1)
    call comm_gather2D(fthetavpar_loc_tmp,nbelts,0,f%n4+1, &
      cible_istart,fthetavpar_glob_tmp,mpi_comm_row) 
      
    if (istart_id .eq. cible_istart) then
      do ivpar = 0,f%n4
        do p_theta = 0,Nbproc_theta-1
          do itheta = 0,dom_theta-1
            fthetavpar(p_theta*dom_theta+itheta,ivpar) = &
              fthetavpar_glob_tmp(itheta,p_theta*(f%n4+1)+ivpar)
          enddo
        enddo
      enddo
    end if
  end subroutine compute_fthetavpar
      
  !------------------------------------------------------------- 
  ! Saving of the three cross-sections of f:
  !  - f(r,theta)  for iphi, ivpar, imu fixed
  !  - f(phi,vpar) for ir, itheta, imu fixed
  !  - f(r,vpar)   for itheta, iphi, imu fixed
  !-------------------------------------------------------------
  subroutine compute_cross_section(f,geom, &
    ir_diag,itheta_diag,iphi_diag, &
    ivpar_diag,imu_diag,frtheta,fphivpar,frvpar,fthetavpar)
    use globals, only : Nbproc_r, Nbproc_theta, Nbproc_mu, &
      pglobal_id, plocal_id, mu_id, istart_id, jstart_id, &
      cible, cible_mu, cible_istart, cible_jstart, &
      istart, iend, jstart, jend, &
      Rarray1_NrNtheta, Rarray_NthetaNvpar, &
      Rarray_NphiNvpar, Rarray_NrNvpar
    include "mpiperso.h"
    type(fdistribu5d)             , intent(in) :: f
    type(geometry)                , intent(in) :: geom
    integer                       , intent(in) :: ir_diag
    integer                       , intent(in) :: itheta_diag
    integer                       , intent(in) :: iphi_diag
    integer                       , intent(in) :: ivpar_diag
    integer                       , intent(in) :: imu_diag
    real(RKIND), dimension(0:,0:), intent(out) :: frtheta
    real(RKIND), dimension(0:,0:), intent(out) :: fphivpar
    real(RKIND), dimension(0:,0:), intent(out) :: frvpar
    real(RKIND), dimension(0:,0:), intent(out) :: fthetavpar
      
    real(RKIND), dimension(:,:), pointer :: frtheta_proc
    real(RKIND), dimension(:,:), pointer :: fphivpar_proc
    real(RKIND), dimension(:,:), pointer :: frvpar_proc
    real(RKIND), dimension(:,:), pointer :: fthetavpar_proc
      
    integer                                 :: ierr
    integer                                 :: nb_send1
    integer                                 :: nb_send2
    integer                                 :: nb_send3
    integer                                 :: nb_send4
    integer                                 :: itag_rtheta
    integer                                 :: itag_rvpar
    integer                                 :: itag_thetavpar
    integer                                 :: itag_phivpar
    integer, dimension(MPI_STATUS_SIZE,1:4) :: status
    integer, dimension(1:4)                 :: reqr
      
    frtheta_proc    => Rarray1_NrNtheta
    nb_send1        = (geom%Nr+1)*(geom%Ntheta+1)
    fphivpar_proc   => Rarray_NphiNvpar
    nb_send2        = (geom%Nphi+1)*(geom%Nvpar+1)
    frvpar_proc     => Rarray_NrNvpar
    nb_send3        = (geom%Nr+1)*(geom%Nvpar+1)
    fthetavpar_proc => Rarray_NthetaNvpar
    nb_send4        = (geom%Ntheta+1)*(geom%Nvpar+1)
      
    !*** initialization of the tags ***
    itag_rtheta    = 11
    itag_phivpar   = 12
    itag_rvpar     = 13
    itag_thetavpar = 14
    reqr           = 0
      
    if ((mu_id.eq.imu_diag)) then
      call compute_frtheta(f,iphi_diag,ivpar_diag,frtheta_proc)
      if ( (ir_diag.ge.istart).and.(ir_diag.le.iend).and. &
        (itheta_diag.ge.jstart).and.(itheta_diag.le.jend) ) then
        call compute_fphivpar(f,ir_diag,itheta_diag,fphivpar_proc)
      endif
      if ((itheta_diag.ge.jstart).and.(itheta_diag.le.jend)) then
        call compute_frvpar(f,itheta_diag,iphi_diag,frvpar_proc)
      endif
      if ((ir_diag.ge.istart).and.(ir_diag.le.iend)) then
        call compute_fthetavpar(f,ir_diag,iphi_diag,fthetavpar_proc)
      endif
    endif
      
    if (pglobal_id.eq.cible) then
      call MPI_IRECV(frtheta(0,0),nb_send1, &
        MPI_REAL8,MPI_ANY_SOURCE,itag_rtheta, &
        MPI_COMM_WORLD,reqr(1),ierr)
      call MPI_IRECV(fphivpar(0,0),nb_send2, &
        MPI_REAL8,MPI_ANY_SOURCE,itag_phivpar, &
        MPI_COMM_WORLD,reqr(2),ierr)
      call MPI_IRECV(frvpar(0,0),nb_send3, &
        MPI_REAL8,MPI_ANY_SOURCE,itag_rvpar, &
        MPI_COMM_WORLD,reqr(3),ierr)
      call MPI_IRECV(fthetavpar(0,0),nb_send4, &
        MPI_REAL8,MPI_ANY_SOURCE,itag_thetavpar, &
        MPI_COMM_WORLD,reqr(4),ierr)
    end if
      
    !*** cross-sections for mu=imu_diag ***
    if ((mu_id.eq.imu_diag)) then
       ! -> f(r,theta) for iphi, ivpar, imu fixed
       if (plocal_id.eq.cible_mu) then
         call MPI_SSEND(frtheta_proc(0,0),nb_send1, &
           MPI_REAL8,cible,itag_rtheta,MPI_COMM_WORLD,ierr)
       end if
       ! -> f(phi,vpar) for ir, itheta, imu fixed
       if ( (ir_diag.ge.istart).and.(ir_diag.le.iend).and. &
            (itheta_diag.ge.jstart).and.(itheta_diag.le.jend) ) then
         call MPI_SSEND(fphivpar_proc(0,0),nb_send2, &
           MPI_REAL8,cible,itag_phivpar,MPI_COMM_WORLD,ierr)
       end if
       ! -> f(r,vpar) for itheta, iphi, imu fixed
       if ((itheta_diag.ge.jstart).and.(itheta_diag.le.jend)) then
         if (jstart_id.eq.cible_jstart) then
           call MPI_SSEND(frvpar_proc(0,0),nb_send3, &
             MPI_REAL8,cible,itag_rvpar,MPI_COMM_WORLD,ierr)
         end if
       end if
       ! -> f(theta,vpar) for ir, iphi, imu fixed
       if ((ir_diag.ge.istart).and.(ir_diag.le.iend)) then
         if (istart_id.eq.cible_istart) then
           call MPI_SSEND(fthetavpar_proc(0,0),nb_send4, &
             MPI_REAL8,cible,itag_thetavpar,MPI_COMM_WORLD,ierr)
         end if
       end if
    end if
      
    if (pglobal_id.eq.cible) &
      call MPI_WAITALL(4,reqr(1),status(1,1),ierr)
  end subroutine compute_cross_section
      
  !------------------------------------------------------------- 
  ! Saving of the cross-sections associated to imu=Nmu/4,
  ! corresponding to trapped particles.
  !-------------------------------------------------------------
  subroutine compute_ftrapped(f,geom)
    use globals, only : ir_f2D_trapped, itheta_f2D_trapped, &
      iphi_f2D_trapped, ivpar_f2D_trapped, imu_f2D_trapped
      
    type(fdistribu5d), intent(in) :: f
    type(geometry)   , intent(in) :: geom
      
    !*** cross-sections for mu=max(mu)/2 ***
    ir_f2D_trapped     = int(geom%Nr/2)
    itheta_f2D_trapped = 0
    iphi_f2D_trapped   = int(geom%Nphi/3)
    ivpar_f2D_trapped  = geom%ivpar0
    imu_f2D_trapped    = int(geom%Nmu/2)
    call compute_cross_section(f,geom,ir_f2D_trapped, &
      itheta_f2D_trapped,iphi_f2D_trapped,ivpar_f2D_trapped, &
      imu_f2D_trapped,frtheta_v0_diag, &
      fphivpar_mumax_diag,frvpar_mumax_diag,fthetavpar_mumax_diag)
      
    !*** (v parallel, mu) cross-sections for theta=0 ***
    call compute_fvparmu(f,ir_f2D_trapped,itheta_f2D_trapped, &
      iphi_f2D_trapped,fvparmu_diag)
  end subroutine compute_ftrapped
      
  !------------------------------------------------------------- 
  ! Saving of the cross-sections associated to mu=0.,
  ! corresponding to passing particles.
  !-------------------------------------------------------------
  subroutine compute_fpassing(f,geom)
    use globals     , only :     nb_vth0, &
      ir_f2D_passing, itheta_f2D_passing, iphi_f2D_passing, &
      ivpar_f2D_passing, imu_f2D_passing
    use utils_module, only : locate, locate_nonequidistant
    type(fdistribu5d), intent(in) :: f
    type(geometry)   , intent(in) :: geom
      
    real(RKIND) :: vpar_diag, mu_diag
      
    !*** cross-sections for mu=0. ***
    ir_f2D_passing     = int(geom%Nr/2)
    itheta_f2D_passing = 0
    iphi_f2D_passing   = int(geom%Nphi/3)
    vpar_diag          = 3._RKIND*geom%vparg(geom%Nvpar)/ &
      real(nb_vth0)    
    mu_diag            = max(0._RKIND,geom%mug(0))
    call locate(vpar_diag,geom%vparg,geom%Nvpar, &
      geom%dvpar," compute_fpassing "//char(0),ivpar_f2D_passing)
!baoter (DMU)
!VG!    call locate(mu_diag,geom%mug,geom%Nmu,geom%dmu,imu_f2D_passing)
    call locate_nonequidistant(mu_diag,geom%mug,geom%Nmu, &
      " compute_fpassing "//char(0),imu_f2D_passing)
!eaoter
    call compute_cross_section(f,geom,ir_f2D_passing,itheta_f2D_passing, &
      iphi_f2D_passing,ivpar_f2D_passing,imu_f2D_passing,frtheta_vmax_diag, &
      fphivpar_mu0_diag,frvpar_mu0_diag,fthetavpar_mu0_diag)
  end subroutine compute_fpassing
end module Pcross_section_module
