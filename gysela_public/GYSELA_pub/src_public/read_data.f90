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
      
!----------------------------------------
! file : read_data.f90
! date : 16/09/2005
! used for the reading of :
! - mesh datas
! - the equilibrium datas
! - restart file (if in restart mode)
!----------------------------------------
      
!------------------------------------
! parallel data reading
!------------------------------------
subroutine parallel_data
  use prec_const
  use globals, only :  pglobal_id, &
    Nbproc_r, Nbproc_theta, Nbproc_mu, bloc_phi, &
    Nbthread, transpose4D, uout_res, outputproc
  use MPIutils_module
  implicit none
  include "mpiperso.h"
  
  namelist /PARALLEL/ Nbproc_r, Nbproc_theta, Nbproc_mu, &
    bloc_phi, Nbthread, transpose4D
      
  integer            :: ierr
  integer            :: position
  integer, parameter :: lbuffer = 1024
  character*1        :: buffer(lbuffer)
      
  if (pglobal_id.eq.outputproc) then    
    open(UNIT=21,FILE='DATA') 
    read(21,PARALLEL)
    write(uout_res,PARALLEL)
    position = 0
    CALL MPI_PACK(Nbproc_r,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(Nbproc_theta,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(Nbproc_mu,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(bloc_phi,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(Nbthread,1,MPI_INTEGER,buffer,lbuffer,position, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(transpose4D,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(buffer,lbuffer,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
    close(21)
  else
    CALL MPI_BCAST(buffer,lbuffer,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
    position = 0
    CALL MPI_UNPACK(buffer,lbuffer,position,Nbproc_r,1, &
      MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,Nbproc_theta,1, &
      MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,Nbproc_mu,1, &
      MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,bloc_phi,1, &
      MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,Nbthread,1, &
      MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,transpose4D,1, &
      MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  end if
      
  ! For optimization of the OpenMP parallelization of 2D advection,
  ! bloc_phi must be at least equal to the number of threads.
  if (bloc_phi .lt. Nbthread) then
    bloc_phi = Nbthread
    if (pglobal_id.eq.0) then
      write(6,*) 'ATTENTION replacement : bloc_phi = ',bloc_phi
    end if
  end if
      
  ! If the number of processors in theta direction is greater
  ! than 1, the transposition is required
  if (Nbproc_theta.gt.1) then
    if (.not.transpose4D) then
      transpose4D = .true.
      if (pglobal_id.eq.0) &
        write(6,*) 'ATTENTION transpose4D is forced to .true.'
    end if
  end if
end subroutine parallel_data
      
!------------------------------------
! mesh data reading
!------------------------------------
subroutine mesh_data
  use prec_const
  use globals, only :  pglobal_id, Nr, Ntheta, Nphi, Nvpar, Nmu, &
    a, rhomin, rhomax, Ltheta, Lphi, aspect_ratio, &
    nb_vth0, Lmu, mumin, integral_vperp, uout_res, outputproc
  implicit none
  include "mpiperso.h"
  
  namelist /MESH/ Nr, Ntheta, Nphi, Nvpar, Nmu, a, rhomin, rhomax, &
    Ltheta, Lphi, aspect_ratio, nb_vth0, Lmu, mumin, integral_vperp
      
  integer            :: ierr
  integer            :: position
  integer, parameter :: lbuffer = 1024
  character*1        :: buffer(lbuffer)
      
  if (pglobal_id.eq.outputproc) then
    open(UNIT=21,FILE='DATA') 
    read(21,MESH)
    write(uout_res,MESH)
    position = 0
    CALL MPI_PACK(Nr,1,MPI_INTEGER,buffer,lbuffer,position, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(Ntheta,1,MPI_INTEGER,buffer,lbuffer,position, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(Nphi,1,MPI_INTEGER,buffer,lbuffer,position, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(Nvpar,1,MPI_INTEGER,buffer,lbuffer,position, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(Nmu,1,MPI_INTEGER,buffer,lbuffer,position, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(a,1,MPI_REAL8,buffer,lbuffer,position, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(rhomin,1,MPI_REAL8,buffer,lbuffer,position, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(rhomax,1,MPI_REAL8,buffer,lbuffer,position, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(Ltheta,1,MPI_REAL8,buffer,lbuffer,position, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(Lphi,1,MPI_REAL8,buffer,lbuffer,position, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(aspect_ratio,1,MPI_REAL8,buffer, &
      lbuffer,position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(nb_vth0,1,MPI_REAL8,buffer,lbuffer,position, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(Lmu,1,MPI_REAL8,buffer,lbuffer,position, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(mumin,1,MPI_REAL8,buffer,lbuffer,position, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(integral_vperp,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(buffer,lbuffer,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
    close(21)
  else
    CALL MPI_BCAST(buffer,lbuffer,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
    position = 0
    CALL MPI_UNPACK(buffer,lbuffer,position,Nr,1,MPI_INTEGER, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,Ntheta,1,MPI_INTEGER, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,Nphi,1,MPI_INTEGER, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,Nvpar,1,MPI_INTEGER, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,Nmu,1,MPI_INTEGER, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,a,1,MPI_REAL8, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,rhomin,1,MPI_REAL8, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,rhomax,1,MPI_REAL8, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,Ltheta,1,MPI_REAL8, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,Lphi,1,MPI_REAL8, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,aspect_ratio,1, &
      MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,nb_vth0,1, &
      MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,Lmu,1, &
      MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,mumin,1, &
      MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,integral_vperp,1, &
      MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  end if
end subroutine mesh_data
      
!----------------------------------------
! verification of the mesh number
!----------------------------------------
subroutine verif_mesh_number
  use prec_const
  use globals, only : pglobal_id, Nr, Ntheta, Nphi, Nvpar, Nmu, &
    Nbproc_r, Nbproc_theta, Lmu
  use MPIutils_module
  implicit none
      
  integer     :: pow	
  real(RKIND) :: rem1, rem2, divr, divtheta, divphi, divvpar	  
      
  if ((Nr.eq.0).or.(Ntheta.eq.0) &
    .or.(Nphi.eq.0).or.(Nvpar.eq.0)) then
    print*,'Problem with one of the dimension which is equal to 0'
    stop
  end if
      
  !------------------------------------------------------
  ! verify the number of points in (r,theta) direction :  
  !  - Nr+1 must be a power of 2                       
  !     (for the parallel transposition)               
  !  -  Ntheta must be a power of 2                    
  !     (for Fast Fourier Transform)                   
  !------------------------------------------------------
  !*** -> (Nr+1) verification ***
  rem1     = modulo(real(log(real(Nr+1))),real(log(2._RKIND))) 
  divr     = log(real(Nr+1))/log(2._RKIND)
  if (rem1.ne.0) then
    pow = int(divr)
    Nr  = max(1,2**pow-1)
    if (pglobal_id.eq.0) then
      write(6,*) 'ATTENTION changment : Nr = ',Nr
    endif
  endif
      
  !*** -> Ntheta verification  ***
  rem2     = modulo(real(log(real(Ntheta))),real(log(2._RKIND))) 
  divtheta = log(real(Ntheta))/log(2._RKIND)
  if (rem2.ne.0) then
    pow    = int(divtheta)
    Ntheta = max(1,2**pow)
    if (pglobal_id.eq.0) then
      write(6,*) 'ATTENTION changment : Ntheta = ',Ntheta
    endif
  end if
      
  !-------------------------------------------------------------
  !  verify the number of points in (phi,vpar) direction : 
  !   -  Nphi must be a power of 2                         
  !      (for Fast Fourier Transform)                      
  !-------------------------------------------------------------
  !*** -> Nphi verification ***
  rem1    = modulo(real(log(real(Nphi))),real(log(2._RKIND))) 
  divphi  = log(real(Nphi))/log(2._RKIND)
  if (rem1.ne.0) then
    pow   = int(divphi)
    Nphi  = max(1,2**pow)
    if (pglobal_id.eq.0) then
      write(6,*) 'ATTENTION replacement : Nphi = ',Nphi
    end if
  end if
      
  if (((Nr+1)/Nbproc_r).lt.32) then
    if (pglobal_id.eq.0) then
      print*, "Attention (Nr+1)/Nbproc_r ", &
        "must be equal to or greater than 32 "
    end if
    call ppexit()
    call exit(-1)
  end if
  if ((Ntheta/Nbproc_theta).lt.32) then
    if (pglobal_id.eq.0) then
      print*, "Attention Ntheta/Nbproc_theta ", &
        "must be equal to or greater than 32 "
    end if
    call ppexit()
    call exit(-1)
  end if
      
  !-------------------------------------------------------------
  !  verify the number of points in mu direction
  !-------------------------------------------------------------
  if ((Lmu.eq.0).and.(Nmu.ne.0)) then
    if (pglobal_id.eq.0) then
      print*, 'Lmu=0 and Nmu not equal to 0 is not possible'
    end if
    call ppexit()
    call exit(-1)
  end if
end subroutine verif_mesh_number
      
!----------------------------------------
! equilibrium data reading
!----------------------------------------
subroutine equil_data
  use prec_const
  use MPIutils_module, only : ppexit
  use globals, only : pglobal_id, Ntheta, Nphi, &
    plasma_current, rbar_in_feq, use_psibar, Zi, tau0, &
    read_n0, read_Te, read_Ti, read_q, profile_choice, rpeak, &
    kappan, kappaTi, kappaTe, deltarn, deltarTi, deltarTe, &
    epsilon, n, m, init_perturb, &
    magnetic_shear, q0, deltaq, alphaq, &
    reversed_shear, qmin, qa, rho_qmin, &
    B_curvature, zonal_flow, &
    filter_choice, filter_deltam, K_curv, &
    ExB_shear, max_shear, rho_shear, deltar_shear, &
    uout_res, outputproc
  implicit none
  include "mpiperso.h"
  
  namelist /EQUIL/ plasma_current, rbar_in_feq, &
    use_psibar, Zi, tau0, &
    read_n0, read_Te, read_Ti, read_q, profile_choice,rpeak, &
    kappan, kappaTi, kappaTe, deltarn, deltarTi, deltarTe, &
    magnetic_shear, q0, deltaq, alphaq, &
    reversed_shear, qmin, qa, rho_qmin, &
    B_curvature, epsilon, n, m, init_perturb, &
    zonal_flow, filter_choice, filter_deltam, &
    ExB_shear, max_shear, rho_shear, deltar_shear
      
  integer            :: ierr
  integer            :: position
  integer, parameter :: lbuffer = 1024
  character*1        :: buffer(lbuffer)
      
  if (pglobal_id.eq.0) then
    open(UNIT=21,FILE='DATA') 
    read(21,EQUIL)
    write(uout_res,EQUIL)
    position = 0
    CALL MPI_PACK(plasma_current,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(rbar_in_feq,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(use_psibar,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(Zi,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(tau0,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(read_n0,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(read_Te,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(read_Ti,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(read_q,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(profile_choice,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(rpeak,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(kappan,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(kappaTi,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(kappaTe,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(deltarn,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(deltarTi,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(deltarTe,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(magnetic_shear,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(q0,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(deltaq,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(alphaq,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(reversed_shear,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(qmin,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(qa,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(rho_qmin,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(B_curvature,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(epsilon,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(m,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(n,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(init_perturb,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(zonal_flow,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(filter_choice,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(filter_deltam,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(ExB_shear,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(max_shear,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(rho_shear,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(deltar_shear,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(buffer,lbuffer,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
    close(21)
  else
    CALL MPI_BCAST(buffer,lbuffer,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
    position = 0
    CALL MPI_UNPACK(buffer,lbuffer,position,plasma_current, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,rbar_in_feq, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,use_psibar, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,Zi, &
      1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,tau0, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,read_n0, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,read_Te, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,read_Ti, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,read_q, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,profile_choice, &
      1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,rpeak, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,kappan, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,kappaTi, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,kappaTe, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,deltarn, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,deltarTi, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,deltarTe, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,magnetic_shear, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,q0, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,deltaq, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,alphaq, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,reversed_shear, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,qmin, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,qa, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,rho_qmin, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,B_curvature, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,epsilon, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,m,1, &
      MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,n,1, &
      MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,init_perturb,1, &
      MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,zonal_flow,1, &
      MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,filter_choice,1, &
      MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,filter_deltam, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,ExB_shear,1, &
      MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,max_shear, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,rho_shear, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,deltar_shear, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      
  !*** Definition of the coefficient according to the case treated
  K_curv = 0._RKIND
  if (magnetic_shear) then
    if (B_curvature) then
      K_curv = 1._RKIND
    end if
  else
    !-> if no magnetic shear, rbar = r
    if (rbar_in_feq) then
      if (pglobal_id.eq.0) &
        print*,' The case magnetic_shear = false', &
          ' and rbar_in_feq = true is not possible.'
      rbar_in_feq = .false.
      if (pglobal_id.eq.0) &
        print*,'BECAREFUL: rbar_in_feq has been put equal to', &
        rbar_in_feq
    end if
    !-> if no magnetic shear, no B curvature can be taken 
    !->   into account
    if (B_curvature) then
      if (pglobal_id.eq.0) &
        print*,' The magnetic_shear = false', &
        ' and  case B_curvature = true is not possible'
      call ppexit()
      call exit(-1)
    end if
    !-> if no magnetic shear, the diagonal filter 
    !->  (filter_choice=1) is not possible
    if ((filter_choice.eq.1).or.(filter_choice.eq.2)) then
      if (pglobal_id.eq.0) &
        print*,' The case magnetic_shear = false', &
          ' and filter_choice = 1 is not possible.'
      call ppexit()
      call exit(-1)
    end if
  end if
      
  !*** Test coherence of (m,n) modes chosen
  if ((m.ge.Ntheta/2).or.(n.ge.Nphi/2)) then
    if (pglobal_id.eq.0) &
      print*,'BECAREFUL: initial m or n non consistent with', &
      ' Ntheta or Nphi', &
      ' changing to m =',2*Ntheta/3,'and n =',2*Nphi/3
    m = 2*Ntheta/3
    n = 2*Nphi/3
!AS!    call ppexit()
!AS!    call exit(-1)
  end if
end subroutine equil_data
      
!----------------------------------------
! algorithm data reading
!----------------------------------------
subroutine algorithm_data
  use prec_const
  use globals, only : pglobal_id, transpose4D, &
    LF_average, deltat, nbiter, &
    diag_nbstep, refresh_nbstep, &
    leapfrog, advec2Ddt, &
    hffiltertype, hffilterfreq, hffilterfref, &
    reducedbegin, reducedend, reduceddt, &
    uout_res, outputproc
  use MPIutils_module, only : ppexit
  implicit none
  include "mpiperso.h"
  
  namelist /ALGORITHM/ leapfrog, LF_average, &
    advec2Ddt, hffiltertype, hffilterfreq, &
    hffilterfref, reducedbegin, reducedend, reduceddt, &
    deltat, nbiter, diag_nbstep, &
    refresh_nbstep
      
  integer            :: ierr
  integer            :: position
  integer, parameter :: lbuffer = 1024
  character*1        :: buffer(lbuffer)
  
  if (pglobal_id.eq.0) then
    open(UNIT=21,FILE='DATA') 
    read(21,ALGORITHM)
    write(uout_res,ALGORITHM)
    position = 0
    CALL MPI_PACK(leapfrog,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(LF_average,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(advec2Ddt,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(hffiltertype,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(hffilterfreq,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(hffilterfref,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(reduceddt,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(reducedbegin,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(reducedend,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(deltat,1,MPI_REAL8,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(nbiter,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(diag_nbstep,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(refresh_nbstep,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(buffer,lbuffer,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
    close(21)
  else
    CALL MPI_BCAST(buffer,lbuffer,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
    position = 0
    CALL MPI_UNPACK(buffer,lbuffer,position,leapfrog, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr) 
    CALL MPI_UNPACK(buffer,lbuffer,position,LF_average, &
      1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,advec2Ddt, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,hffiltertype, &
      1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,hffilterfreq, &
      1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,hffilterfref, &
       1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,reduceddt, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,reducedbegin, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,reducedend, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,deltat, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,nbiter, &
      1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,diag_nbstep, &
      1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,refresh_nbstep, &
      1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  end if
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)   
      
end subroutine algorithm_data
      
!----------------------------------------
! test data reading 
!----------------------------------------
subroutine test_data
  use prec_const
  use globals, only : pglobal_id, memory_test, &
    gyroaverage, slab_geometry, Rosenbluth_Hinton, &
    m0n0_eq0, single_m, single_n, &
    solve_poiss, RHS_only, &
    coeff_pol, feq_choice, &
    vpar0_option, init_perturb, Phi00BCr_Neumann, &
    hvpar_in_fperturb, uout_res, outputproc
  use MPIutils_module, only : ppexit
  implicit none
  include "mpiperso.h"
  
  namelist /TEST/ memory_test, gyroaverage, &
    slab_geometry, Rosenbluth_Hinton, &
    m0n0_eq0, single_m, single_n, solve_poiss, RHS_only, &
    coeff_pol, feq_choice, &
    vpar0_option, Phi00BCr_Neumann, hvpar_in_fperturb
      
  integer            :: ierr
  integer            :: position
  integer, parameter :: lbuffer = 1024
  character*1        :: buffer(lbuffer)
  
  if (pglobal_id.eq.outputproc) then
    open(UNIT=21,FILE='DATA') 
    read(21,TEST)
    write(uout_res,TEST)
    position = 0
    CALL MPI_PACK(memory_test,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(gyroaverage,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(slab_geometry,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(Rosenbluth_Hinton,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(m0n0_eq0,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(single_m,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(single_n,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(solve_poiss,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(coeff_pol,1,MPI_REAL8,buffer,lbuffer,position, &
      MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(RHS_only,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(feq_choice,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(vpar0_option,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(Phi00BCr_Neumann,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(hvpar_in_fperturb,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(buffer,lbuffer,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
    close(21)
  else
    CALL MPI_BCAST(buffer,lbuffer,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
    position = 0
    CALL MPI_UNPACK(buffer,lbuffer,position,memory_test, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,gyroaverage, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,slab_geometry, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,Rosenbluth_Hinton, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,m0n0_eq0, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,single_m, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,single_n, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,solve_poiss, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,coeff_pol, &
      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,RHS_only,1,MPI_LOGICAL,&
      MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,feq_choice, &
      1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,vpar0_option, &
      1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,Phi00BCr_Neumann, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,hvpar_in_fperturb, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  end if
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  
  !*** check the coherence of the datas for the ***
  !***   Rosenbluth-Hinton test                 ***
  if (Rosenbluth_Hinton) then
    if (m0n0_eq0) then
      if (pglobal_id.eq.0) &
        print*,"The mode (0,0) cannot be forced to 0", &
        " in the case of the Rosenbluth-Hinton test"
      call ppexit()
      call exit(-1)
    end if
  end if
  !*** force "init_perturb" to case 3 if "single_n" is true  ***
  !***                              4 if "single_m" is true  ***
  !***                              1 if "single_m" and      ***
  !***                                  "single_m" are true  ***
  !***  => avoids bad loops on "n" and "m" in subroutine     ***
  !***    init_fperturb                                      ***
  if (single_n) then
    init_perturb = 3
  end if
  if (single_m) then
    init_perturb = 4
  end if  
  if (single_n.and.single_m) then
    init_perturb = 1
  end if
end subroutine test_data
      
!----------------------------------------
! output data reading 
!----------------------------------------
subroutine output_data
  use prec_const
  use globals, only : pglobal_id, integration_CS, diag_level,  &
    Phi3D_saving, f5D_saving, FMoments3D_saving, CalviExp_saving, &
    Phi3D_nbstep, f5D_nbstep, FMoments3D_nbstep, CalviExp_nbstep, &
    rst_saving, uout_res, outputproc
  implicit none
  include "mpiperso.h"
  
  namelist /OUTPUT/ integration_CS, diag_level, &
    Phi3D_saving, FMoments3D_saving, f5D_saving, CalviExp_saving, &
    Phi3D_nbstep, FMoments3D_nbstep, f5D_nbstep, CalviExp_nbstep, &
    rst_saving
      
  integer            :: ierr
  integer            :: position
  integer, parameter :: lbuffer = 1024
  character*1        :: buffer(lbuffer)
  
  if (pglobal_id.eq.0) then
    open(UNIT=21,FILE='DATA') 
    read(21,OUTPUT)
    write(uout_res,OUTPUT)
    position = 0
    CALL MPI_PACK(integration_CS,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(diag_level,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(Phi3D_saving,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(FMoments3D_saving,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(f5D_saving,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(CalviExp_saving,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(Phi3D_nbstep,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(FMoments3D_nbstep,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(f5D_nbstep,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(CalviExp_nbstep,1,MPI_INTEGER,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_PACK(rst_saving,1,MPI_LOGICAL,buffer,lbuffer, &
      position,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(buffer,lbuffer,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
    close(21)
  else
    CALL MPI_BCAST(buffer,lbuffer,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
    position = 0
    CALL MPI_UNPACK(buffer,lbuffer,position,integration_CS, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,diag_level, &
      1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,Phi3D_saving, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,FMoments3D_saving, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,f5D_saving, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,CalviExp_saving, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,Phi3D_nbstep, &
      1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,FMoments3D_nbstep, &
      1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,f5D_nbstep, &
      1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,CalviExp_nbstep, & 
      1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(buffer,lbuffer,position,rst_saving, &
      1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
  end if
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
end subroutine output_data
      
!-----------------------------------------------------
! Verification of the compatibility between the number
! of processors and the number of points
!-----------------------------------------------------
subroutine verif_parallel_data
  use globals, only : Nr, Ntheta, Nphi, Nmu, pglobal_id,&
    Nbproc_tot, Nbproc_r, Nbproc_theta, Nbproc_mu, bloc_phi
  use MPIutils_module, only : ppexit
  implicit none
      
  !*** tests on the number of processors in each direction ***
  !-> Verification of Nbproc_tot
  if (Nbproc_tot.ne.(Nbproc_r*Nbproc_theta*Nbproc_mu)) then
    if (pglobal_id.eq.0) then
      print*,' Nbproc_tot not equal to ', &
        'Nbproc_r*Nbproc_theta*Nbproc_mu'
      print*,' => Nbproc_tot   = ',Nbproc_tot
      print*,'    Nbproc_r     = ',Nbproc_r
      print*,'    Nbproc_theta = ',Nbproc_theta
      print*,'    Nbproc_mu    = ',Nbproc_mu
    end if
    call ppexit()
    call exit(-1)
  end if
      
  !-> Verification of Nbproc_r
  if (Nbproc_r.gt.(Nr+1)) then
    if (pglobal_id.eq.0) then
      print*,' Nbproc_r greater than Nr+1 '
      print*,' => Nbproc_r = ',Nbproc_r
      print*,'    Nr       = ',Nr
    end if
    call ppexit()
    call exit(-1)
  end if
  if (mod((Nr+1),Nbproc_r) .ne. 0) then
    if (pglobal_id.eq.0) &
      print *,"Nbproc_r does not divide (Nr+1)"
    call ppexit()
    call exit(-1)
  endif
  
  !-> Verification of Nbproc_theta
  if (Nbproc_theta.gt.(Ntheta+1)) then
    if (pglobal_id.eq.0) then
      print*,' Nbproc_theta greater than Ntheta+1 '
      print*,' => Nbproc_theta = ',Nbproc_theta
      print*,'    Ntheta       = ',Ntheta
    end if
    call ppexit()
    call exit(-1)
  end if
  if (mod((Ntheta),Nbproc_theta) .ne. 0) then 
    if (pglobal_id.eq.0) &
      print *,"Nbproc_theta does not divide Ntheta"
    call ppexit()
    call exit(-1)
  endif
      
  !-> Verification of bloc_phi
  if (mod((Nphi),bloc_phi) .ne. 0) then 
    if (pglobal_id.eq.0) &
      print *,"bloc_phi does not divide (Nphi)"
    call ppexit()
    call exit(-1)
  endif
      
  !-> Verification of Nbproc_mu
  if (Nbproc_mu.gt.(Nmu+1)) then
    if (pglobal_id.eq.0) then
      print*,' Nbproc_mu greater than Nmu+1 '
      print*,' => Nbproc_mu = ',Nbproc_mu
      print*,'    Nmu       = ',Nmu
    end if
    call ppexit()
    call exit(-1)
  end if
end subroutine verif_parallel_data
      
!----------------------------------------
!  reading of the input file .dat
!----------------------------------------
subroutine read_input
  use globals
  use MPIutils_module, only : ppexit
  implicit none
      
  real(RKIND) :: max_memory_estimated
      
  call parallel_data
  call mesh_data
  call equil_data
  call algorithm_data
  call test_data
  call output_data
      
  !*** Verification of the number of points in the mesh ***
  call verif_mesh_number  
      
  !*** Verification of the compatibility between the number ***
  !*** of processors and the number of points               ***
  call verif_parallel_data
      
end subroutine read_input
