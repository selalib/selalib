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
      
!-------------------------------------------------------
! file : gys_alloc.f90
! date : 24/11/2010
!  - allocation of all the arrays required for a 
!  GYSELA run
!-------------------------------------------------------
module gys_alloc_module
  implicit none
      
  !******************************
  contains
  !******************************       
      
  !-----------------------------------------------
  ! temporary global arrays allocation
  !-----------------------------------------------
  subroutine allocate_temporary_arrays(geom)
    use globals
    use geometry_class
    use mem_alloc_module
    use MPIutils_module
    implicit none
    type(geometry), intent(in)  :: geom
      
    integer :: nbmoments, i
      
    !*** arrays for parallel QN solver ***
    nbmoments = 2
    call glob_allocate(Rarray_PNrPNthetaNphi_nbM2,0,(nbmoments-1), &
      istart,iend,jstart,jend,0,Nphi-1, &
      'Rarray_PNrPNthetaNphi_nbM2')
    call glob_allocate(Rarray_NrNthetamuphi_nbM2, &
      0,nbmoments*dom_r-1,0,dom_theta*Nbproc_loc*(Nmu+1)*dom_mapphi, &
      'Rarray_NrNthetamuphi_nbM2')
    !*** arrays for fluid moment computations ***
    nbmoments = 8
    call glob_allocate(Rarray_PNrPNthetaNphi_nbM8,0,(nbmoments-1), &
      istart,iend,jstart,jend,0,Nphi-1, &
      'Rarray_PNrPNthetaNphi_nbM8')
    call glob_allocate(Rarray_NrNthetamuphi_nbM8, &
      0,nbmoments*dom_r-1,0,dom_theta*Nbproc_loc*(Nmu+1)*dom_mapphi, &
      'Rarray_NrNthetamuphi_nbM8')
      
    !*** array allocations for gyroaveraged    ***
    !*** real array allocations ***
    !*** -> array(0:Nr) ***
    call glob_allocate(Rarray1_Nr,0,geom%Nr,'Rarray1_Nr')
    call glob_allocate(Rarray2_Nr,0,geom%Nr,'Rarray2_Nr')
    call glob_allocate(Rarray3_Nr,0,geom%Nr,'Rarray3_Nr')
    call glob_allocate(Rarray4_Nr,0,geom%Nr,'Rarray4_Nr')
    !*** -> array(0:Ntheta) ***
    call glob_allocate(Rarray1_Ntheta,0,geom%Ntheta,'Rarray1_Ntheta')
    call glob_allocate(Rarray2_Ntheta,0,geom%Ntheta,'Rarray2_Ntheta')
    !*** -> array(0:Nphi) ***
    call glob_allocate(Rarray1_Nphi,0,geom%Nphi,'Rarray1_Nphi')
    call glob_allocate(Rarray2_Nphi,0,geom%Nphi,'Rarray2_Nphi')
    !*** -> array(0:Nvpar) ***
    call glob_allocate(Rarray1_Nvpar,0,geom%Nvpar,'Rarray1_Nvpar')
    call glob_allocate(Rarray2_Nvpar,0,geom%Nvpar,'Rarray2_Nvpar')
    call glob_allocate(Rarray3_Nvpar,0,geom%Nvpar,'Rarray3_Nvpar')
    call glob_allocate(Rarray4_Nvpar,0,geom%Nvpar,'Rarray4_Nvpar')
    call glob_allocate(Rarray5_Nvpar,0,geom%Nvpar,'Rarray5_Nvpar')
    call glob_allocate(Rarray6_Nvpar,0,geom%Nvpar,'Rarray6_Nvpar')
    call glob_allocate(Rarray7_Nvpar,0,geom%Nvpar,'Rarray7_Nvpar')
    !*** -> array(-1:Nphi+1) ***
    call glob_allocate(Rarray_m1Nphip1,-1,geom%Nphi+1, &
      'Rarray_m1Nphip1')
    !*** -> array(-1:Nvpar+1) ***
    call glob_allocate(Rarray_m1Nvparp1,-1,geom%Nvpar+1, &
      'Rarray_m1Nvparp1')
    !*** -> array(0:1,0:Ntheta) ***
    call glob_allocate(Rarray_1Ntheta,0,1,0,geom%Ntheta, &
      'Rarray_1Ntheta')
    !*** -> array(0:Nr,0:Ntheta) ***
    call glob_allocate(Rarray1_NrNtheta,0,geom%Nr,0,geom%Ntheta, &
      'Rarray1_NrNtheta')
    call glob_allocate(Rarray2_NrNtheta,0,geom%Nr,0,geom%Ntheta, &
      'Rarray2_NrNtheta')
    call glob_allocate(Rarray3_NrNtheta,0,geom%Nr,0,geom%Ntheta, &
      'Rarray3_NrNtheta')
    call glob_allocate(Rarray4_NrNtheta,0,geom%Nr,0,geom%Ntheta, &
      'Rarray4_NrNtheta')
    call glob_allocate(Rarray5_NrNtheta,0,geom%Nr,0,geom%Ntheta, &
      'Rarray5_NrNtheta')
    call glob_allocate(Rarray6_NrNtheta,0,geom%Nr,0,geom%Ntheta, &
      'Rarray6_NrNtheta')
    call glob_allocate(Rarray7_NrNtheta,0,geom%Nr,0,geom%Ntheta, &
      'Rarray7_NrNtheta')
    !*** -> array(0:Nr,0:Nphi) ***
    call glob_allocate(Rarray_NrNphi,0,geom%Nr,0,geom%Nphi, &
      'Rarray_NrNphi')
    call glob_allocate(Rarray2_NrNphi,0,geom%Nr,0,geom%Nphi, &
      'Rarray2_NrNphi')
    !*** -> array(0:Nr,0:Nvpar) ***
    call glob_allocate(Rarray_NrNvpar,0,geom%Nr,0,geom%Nvpar, &
      'Rarray_NrNvpar')
    !*** -> array(0:Ntheta,0:Nvpar) ***
    call glob_allocate(Rarray_NthetaNvpar,0,geom%Ntheta,&
      0,geom%Nvpar,'Rarray_NthetaNvpar')
    !*** -> array(0:Nphi,0:Nvpar) ***
    call glob_allocate(Rarray_NphiNvpar,0,geom%Nphi,0,geom%Nvpar, &
      'Rarray_NphiNvpar')
    !*** -> array(0:Nvpar,0:Nmu) ***
    call glob_allocate(Rarray_NvparNmu,0,geom%Nvpar,0,geom%Nmu, &
      'Rarray_NvparNmu')
    !*** -> array(-1:Nr+1,-1:Ntheta+1) ***
    call glob_allocate(Rarray1_m1Nrp1m1Nthetap1,-1,geom%Nr+1, &
      -1,geom%Ntheta+1,'Rarray1_m1Nrp1m1Nthetap1')
    call glob_allocate(Rarray2_m1Nrp1m1Nthetap1,-1,geom%Nr+1, &
      -1,geom%Ntheta+1,'Rarray2_m1Nrp1m1Nthetap1')
    call glob_allocate(Rarray3_m1Nrp1m1Nthetap1,-1,geom%Nr+1, &
      -1,geom%Ntheta+1,'Rarray3_m1Nrp1m1Nthetap1')
    !*** -> array(-1:Nphi+1,-1:Nvpar+1) ***
    call glob_allocate(Rarray_m1Nphip1m1Nvparp1,-1,geom%Nphi+1, &
      -1,geom%Nvpar+1,'Rarray_m1Nphip1m1Nvparp1')
    !*** -> array(0:Nr,0:Ntheta,iphistart:iphistart+dom_mapphi-1) ***
    call glob_allocate(Rarray1_NrNthetaDomphi,0,geom%Nr, &
      0,geom%Ntheta,iphistart,iphistart+dom_mapphi-1, &
      'Rarray1_NrNthetaDomphi')
    call glob_allocate(Rarray2_NrNthetaDomphi,0,geom%Nr, &
      0,geom%Ntheta,iphistart,iphistart+dom_mapphi-1, &
      'Rarray2_NrNthetaDomphi')
    call glob_allocate(Rarray3_NrNthetaDomphi,0,geom%Nr, &
      0,geom%Ntheta,iphistart,iphistart+dom_mapphi-1, &
      'Rarray3_NrNthetaDomphi')
    call glob_allocate(Rarray4_NrNthetaDomphi,0,geom%Nr, &
      0,geom%Ntheta,iphistart,iphistart+dom_mapphi-1, &
      'Rarray4_NrNthetaDomphi')
    call glob_allocate(Rarray5_NrNthetaDomphi,0,geom%Nr, &
      0,geom%Ntheta,iphistart,iphistart+dom_mapphi-1, &
      'Rarray5_NrNthetaDomphi')
    !*** -> array(0:Nr,0:Ntheta,0:Nphi) ***
    call glob_allocate(Rarray1_NrNthetaNphi,0,geom%Nr, &
      0,geom%Ntheta,0,geom%Nphi,'Rarray1_NrNthetaNphi')
    call glob_allocate(Rarray2_NrNthetaNphi,0,geom%Nr, &
      0,geom%Ntheta,0,geom%Nphi,'Rarray2_NrNthetaNphi')
    !*** -> array(istart:iend,jstart:jend,0:Nphi) ***
    call glob_allocate(Rarray1_PNrPNthetaNphi, &
      istart,iend,jstart,jend,0,geom%Nphi, &
      'Rarray1_PNrPNthetaNphi')
    call glob_allocate(Rarray2_PNrPNthetaNphi, &
      istart,iend,jstart,jend,0,geom%Nphi, &
      'Rarray2_PNrPNthetaNphi')
    call glob_allocate(Rarray3_PNrPNthetaNphi, &
      istart,iend,jstart,jend,0,geom%Nphi, &
      'Rarray3_PNrPNthetaNphi')
    !*** complex array allocation ***
    !*** -> array(1:Nr+1,1:Ntheta,1:Nphi) ***
      call glob_allocate(Carray_1Nrp11Ntheta1Nphi,1,geom%Nr+1, &
        1,geom%Ntheta,1,geom%Nphi,'Carray_1Nrp11Ntheta1Nphi')
    !*** array for counting the number of particles out of domain ***
    call glob_allocate(nbpart_rleft_iter,0,Nmu, &
      'nbpart_rleft_iter')
    call glob_allocate(nbpart_rright_iter,0,Nmu, &
      'nbpart_rright_iter')
    call glob_allocate(nbpart_thleft_iter,0,Nmu, &
      'nbpart_thleft_iter')
    call glob_allocate(nbpart_thright_iter,0,Nmu, &
      'nbpart_thright_iter')
    if (plocal_id.eq.0) then
      call glob_allocate(nbpart_rleft_muid,0,nbiter, &
        'nbpart_rleft_muid')
      call glob_allocate(nbpart_rright_muid,0,nbiter, &
        'nbpart_rright_muid')
      call glob_allocate(nbpart_thleft_muid,0,nbiter, &
        'nbpart_thleft_muid')
      call glob_allocate(nbpart_thright_muid,0,nbiter, &
        'nbpart_thright_muid')
    end if
  end subroutine allocate_temporary_arrays
      
  !-----------------------------------------------
  ! temporary global arrays deallocation
  !-----------------------------------------------
  subroutine deallocate_temporary_arrays
    use globals
    use mem_alloc_module
    implicit none
      
    integer :: i
      
    !*** array deallocation for parallel QN solver ***
    call glob_deallocate(Rarray_PNrPNthetaNphi_nbM2)
    call glob_deallocate(Rarray_NrNthetamuphi_nbM2)
    !*** array deallocation for fluid moments ***
    call glob_deallocate(Rarray_PNrPNthetaNphi_nbM8)
    call glob_deallocate(Rarray_NrNthetamuphi_nbM8)
    !*** real array deallocations ***
    call glob_deallocate(Rarray1_Nr)
    call glob_deallocate(Rarray2_Nr)
    call glob_deallocate(Rarray3_Nr)
    call glob_deallocate(Rarray4_Nr)
    call glob_deallocate(Rarray1_Ntheta)
    call glob_deallocate(Rarray2_Ntheta)
    call glob_deallocate(Rarray1_Nphi)
    call glob_deallocate(Rarray2_Nphi)
    call glob_deallocate(Rarray1_Nvpar)
    call glob_deallocate(Rarray2_Nvpar)
    call glob_deallocate(Rarray3_Nvpar)
    call glob_deallocate(Rarray4_Nvpar)
    call glob_deallocate(Rarray5_Nvpar)
    call glob_deallocate(Rarray6_Nvpar)
    call glob_deallocate(Rarray7_Nvpar)
    call glob_deallocate(Rarray_m1Nphip1)
    call glob_deallocate(Rarray_m1Nvparp1)
    call glob_deallocate(Rarray_1Ntheta)
    call glob_deallocate(Rarray1_NrNtheta)
    call glob_deallocate(Rarray2_NrNtheta)
    call glob_deallocate(Rarray3_NrNtheta)
    call glob_deallocate(Rarray4_NrNtheta)
    call glob_deallocate(Rarray5_NrNtheta)
    call glob_deallocate(Rarray6_NrNtheta)
    call glob_deallocate(Rarray7_NrNtheta)
    call glob_deallocate(Rarray_NrNphi)
    call glob_deallocate(Rarray2_NrNphi)
    call glob_deallocate(Rarray_NrNvpar)
    call glob_deallocate(Rarray_NthetaNvpar)
    call glob_deallocate(Rarray_NphiNvpar)
    call glob_deallocate(Rarray_NvparNmu)
    call glob_deallocate(Rarray1_m1Nrp1m1Nthetap1)
    call glob_deallocate(Rarray2_m1Nrp1m1Nthetap1)
    call glob_deallocate(Rarray3_m1Nrp1m1Nthetap1)
    call glob_deallocate(Rarray_m1Nphip1m1Nvparp1)
    call glob_deallocate(Rarray1_NrNthetaDomphi)
    call glob_deallocate(Rarray2_NrNthetaDomphi)
    call glob_deallocate(Rarray3_NrNthetaDomphi)
    call glob_deallocate(Rarray4_NrNthetaDomphi)
    call glob_deallocate(Rarray5_NrNthetaDomphi)
    call glob_deallocate(Rarray1_NrNthetaNphi)
    call glob_deallocate(Rarray2_NrNthetaNphi)
    call glob_deallocate(Rarray1_PNrPNthetaNphi)
    call glob_deallocate(Rarray2_PNrPNthetaNphi)
    call glob_deallocate(Rarray3_PNrPNthetaNphi)
    !*** complex array deallocations ***
      call glob_deallocate(Carray_1Nrp11Ntheta1Nphi)
    !*** array for counting the number of particles out of domain ***
    call glob_deallocate(nbpart_rleft_iter)
    call glob_deallocate(nbpart_rright_iter)
    call glob_deallocate(nbpart_thleft_iter)
    call glob_deallocate(nbpart_thright_iter)
    if (plocal_id.eq.0) then
      call glob_deallocate(nbpart_rleft_muid)
      call glob_deallocate(nbpart_rright_muid)
      call glob_deallocate(nbpart_thleft_muid)
      call glob_deallocate(nbpart_thright_muid)
    end if
  end subroutine deallocate_temporary_arrays
end module gys_alloc_module
