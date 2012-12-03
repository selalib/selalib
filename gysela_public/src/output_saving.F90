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
! file : output_saving.f90
! date : 16/09/2005
!  used for saving the physic results
!------------------------------------------------
module output_saving_module
  use geometry_class
  use init_profile_class
  use init_magnetic_class
  use globals
  use mem_alloc_module
#ifndef NOHDF5
  use HDF5
  use HDF5_io_module
#endif
      
  implicit none
  
  !******************************
  contains
  !******************************
      
      
  !***********************************************************
  !   SAVING DIRECTLY IN HDF5 FORMAT
  !***********************************************************
  !---------------------------------------------- 
  ! Output saving of the coordinate system
  !  . covariant components of the metric tensor
  !  . jacobian in space 
  !  -> directly in HDF5 binary format
  !----------------------------------------------   
  subroutine HDF5_coordsys_saving(geom,coord_sys)
    use coord_system_class
    type(geometry)     , intent(in) :: geom
    type(coord_system) , intent(in) :: coord_sys
#ifndef NOHDF5
    integer           :: Nbr, Nbc, icr
    character(LEN=30) :: coordsys_file
    integer(HID_T)    :: file_id 
    integer           :: ierr
      
    coordsys_file = "coord_system.h5"
    call HDF5_create(coordsys_file,file_id,ierr)
    if (ierr.ne.0) then
      print*,'pglobal_id = ',pglobal_id, &
        ' ==> error for opening of ',coordsys_file
    end if
      
    Nbr = geom%Nr+1
    Nbc = geom%Ntheta+1
    !*** covariant components of the metric tensor ***
    call HDF5_array2D_saving(file_id,coord_sys%g11,Nbr,Nbc, &
      'g11'//char(0))
    call HDF5_array2D_saving(file_id,coord_sys%g12,Nbr,Nbc, &
      'g12'//char(0))
    call HDF5_array2D_saving(file_id,coord_sys%g13,Nbr,Nbc, &
      'g13'//char(0))
    call HDF5_array2D_saving(file_id,coord_sys%g21,Nbr,Nbc, &
      'g21'//char(0))
    call HDF5_array2D_saving(file_id,coord_sys%g22,Nbr,Nbc, &
      'g22'//char(0))
    call HDF5_array2D_saving(file_id,coord_sys%g23,Nbr,Nbc, &
      'g23'//char(0))
    call HDF5_array2D_saving(file_id,coord_sys%g31,Nbr,Nbc, &
      'g31'//char(0))
    call HDF5_array2D_saving(file_id,coord_sys%g32,Nbr,Nbc, &
      'g32'//char(0))
    call HDF5_array2D_saving(file_id,coord_sys%g33,Nbr,Nbc, &
      'g33'//char(0))
      
    !*** jacobian in space ***
    call HDF5_array2D_saving(file_id,jacobian_space,Nbr,Nbc, &
      'jacob_space'//char(0))
      
    !***  \int dtheta dphi Js with Js the jacobian in space ***
    call HDF5_array1D_saving(file_id,intdthetadphi_Js(0:),Nbr, &
      'intdthetadphi_Js'//char(0))
    call HDF5_close(file_id)
#endif
  end subroutine HDF5_coordsys_saving
      
  !---------------------------------------------- 
  ! Output saving of the initial state 
  !  -> directly in HDF5 binary format
  !----------------------------------------------   
  subroutine HDF5_init_saving(geom,init_prof,init_magnet &
    ,init_curr &
    )
    use globals, only : numfmt_rst, nb_restart_diag
    use coord_system_class, only : R
    type(geometry)     , intent(in) :: geom
    type(init_profile) , intent(in) :: init_prof
    type(init_magnetic), intent(in) :: init_magnet
    type(init_current) , intent(in) :: init_curr
#ifndef NOHDF5
    integer           :: Nbr, Nbc, icr
    character(LEN=30) :: init_state_file
    integer(HID_T)    :: file_id   
    integer           :: ierr
      
    real(RKIND) :: rhostar
      
    write(init_state_file,'(A,'//numfmt_rst//',A)') &
      "init_state", nb_restart_diag, ".h5"
    call HDF5_create(trim(init_state_file),file_id,ierr)
    if (ierr.ne.0) then
      print*,'pglobal_id = ',pglobal_id, &
        ' ==> error for opening of ',init_state_file
    end if
      
    !*** initial geometry ****
    rhostar = 1._RKIND/a
    call HDF5_integer_saving(file_id,geom%Nr,'Nr'//char(0))
    call HDF5_integer_saving(file_id,geom%Ntheta,'Ntheta'//char(0))
    call HDF5_integer_saving(file_id,geom%Nphi,'Nphi'//char(0))
    call HDF5_integer_saving(file_id,geom%Nvpar,'Nvpar'//char(0))
    call HDF5_integer_saving(file_id,geom%Nmu,'Nmu'//char(0))
    call HDF5_real_saving(file_id,rhostar,'rhostar'//char(0))
    call HDF5_real_saving(file_id,R0,'R0'//char(0))
    Nbc = geom%Nr+1
    call HDF5_array1D_saving(file_id,geom%rg(0:),Nbc,'rg'//char(0))
    Nbc = geom%Ntheta+1
    call HDF5_array1D_saving(file_id,geom%thetag(0:),Nbc, &
      'thetag'//char(0))
    Nbc = geom%Nphi+1 
    call HDF5_array1D_saving(file_id,geom%phig(0:),Nbc, &
      'phig'//char(0))
    Nbc = geom%Nvpar+1
    call HDF5_array1D_saving(file_id,geom%vparg(0:),Nbc, &
      'vparg'//char(0))
    Nbc = geom%Nmu+1
    call HDF5_array1D_saving(file_id,geom%mug(0:),Nbc, &
      'mug'//char(0))
      
    !*** initial profiles ****
    Nbc = geom%Nr+1
    call HDF5_array1D_saving(file_id,init_prof%Te(0:),Nbc, &
      'Te'//char(0))
    call HDF5_array1D_saving(file_id,init_prof%Ti,Nbc, &
      'Ti'//char(0))
    call HDF5_array1D_saving(file_id,init_prof%n0,Nbc, &
      'n0'//char(0))
    call HDF5_array1D_saving(file_id,init_prof%dTedr,Nbc, &
      'dTedr'//char(0))
    call HDF5_array1D_saving(file_id,init_prof%dn0dr,Nbc, &
      'dn0dr'//char(0))
    call HDF5_array1D_saving(file_id,init_prof%dlogTidr,Nbc, &
      'dlogTidr'//char(0))
    call HDF5_array1D_saving(file_id,init_prof%iota,Nbc, &
      'iota'//char(0))
    call HDF5_array1D_saving(file_id,init_prof%shear,Nbc, &
      'shear'//char(0))
    call HDF5_array1D_saving(file_id,init_prof%psi,Nbc, &
      'psi'//char(0))
      
    !*** Dissipation coefficient profiles ****
    Nbr = 1
    Nbc = geom%Nr+1
    call HDF5_array1D_saving(file_id,init_prof%fric_coeff(0:), &
      Nbc,'coefNu'//char(0))
    call HDF5_array1D_saving(file_id,init_prof%diff_coeffDr(0:), &
      Nbc,'coefDr'//char(0))
      
    !*** initial profile of R(r,theta) ***
    Nbr = geom%Nr+1
    Nbc = geom%Ntheta+1
    call HDF5_array2D_saving(file_id,R,Nbr,Nbc,'R'//char(0))
      
    !*** initial profile of B(r,theta) and its derivatives ***
    Nbr = geom%Nr+1
    Nbc = geom%Ntheta+1
    call HDF5_array2D_saving(file_id,init_magnet%B_norm, &
      Nbr,Nbc,'B'//char(0))
    call HDF5_array2D_saving(file_id,init_magnet%dBdr, &
      Nbr,Nbc,'dBdr'//char(0))
    call HDF5_array2D_saving(file_id,init_magnet%dBdtheta, &
      Nbr,Nbc,'dBdtheta'//char(0))
    
    !*** saving of the n_eq(r) and dneq_dr(r) profiles ***
    Nbr = 1
    Nbc = geom%Nr+1
    call HDF5_array1D_saving(file_id,neq_r,Nbc,'neq'//char(0))
    call HDF5_array1D_saving(file_id,dneq_dr,Nbc,'dneq_dr'//char(0))
    !*** saving of the Ti_eq(r) profiles ***
    Nbr = 1
    Nbc = geom%Nr+1
    call HDF5_array1D_saving(file_id,Tieq_r,Nbc,'Tieq'//char(0))
      
    !*** initial profile of mu0J_phi(r,theta) ***
    Nbr = geom%Nr+1
    Nbc = geom%Ntheta+1
    call HDF5_array2D_saving(file_id,init_curr%mu0_Jphi, &
      Nbr,Nbc,'mu0_Jphi'//char(0))
      
    call HDF5_close(file_id)
#endif
  end subroutine HDF5_init_saving
      
  !---------------------------------------------- 
  ! Output saving of the initial state 
  !  -> in ASCII format
  !----------------------------------------------   
  subroutine ASCII_init1D_saving(geom, init_prof &
    ,init_magnet &
    )
    use globals, only : numfmt_rst, nb_restart_diag        
    use coord_system_class, only : R
    type(geometry)     , intent(in) :: geom
    type(init_profile) , intent(in) :: init_prof
    type(init_magnetic), intent(in) :: init_magnet
    integer            :: ir
    integer, parameter :: uinit = 30
    character(LEN=30)  :: init_state_file
      
    write(init_state_file,'(A,'//numfmt_rst//',A)') &
      "init_state1D", nb_restart_diag, ".dat"
    open(uinit, file = init_state_file, status = 'UNKNOWN', &
      form = 'FORMATTED')
      write(uinit,'(A12,16A20)') &
        'rg','Te','Ti','n0','dTedr', &
        'dlogTidr','dn0dr','neq','dneq_dr', &
        'Tieq','iota','shear','psi','coefNu','coefDr'
    do ir = 0,geom%Nr 
        write(uinit,'(1pe12.3,14(1pe20.12))') &
          geom%rg(ir), init_prof%Te(ir), &
          init_prof%Ti(ir), init_prof%n0(ir), &
          init_prof%dTedr(ir), init_prof%dlogTidr(ir), &
          init_prof%dn0dr(ir), neq_r(ir), &
          dneq_dr(ir), Tieq_r(ir), init_prof%iota(ir), &
          init_prof%shear(ir), init_prof%psi(ir), &
          init_prof%fric_coeff(ir), init_prof%diff_coeffDr(ir)
    end do
    close(uinit)
      
!VG!    !*** initial profile of R(r,theta) ***
!VG!    Nbr = geom%Nr+1
!VG!    Nbc = geom%Ntheta+1
!VG!    call HDF5_array2D_saving(file_id,R,Nbr,Nbc,'R'
!VG!
!VG!    !*** initial profile of B(r,theta) and its derivatives ***
!VG!    Nbr = geom%Nr+1
!VG!    Nbc = geom%Ntheta+1
!VG!    call HDF5_array2D_saving(file_id,init_magnet%B_norm, &
!VG!      Nbr,Nbc,'B'
!VG!    call HDF5_array2D_saving(file_id,init_magnet%dBdr, &
!VG!      Nbr,Nbc,'dBdr'
!VG!    call HDF5_array2D_saving(file_id,init_magnet%dBdtheta, &
!VG!      Nbr,Nbc,'dBdtheta'
  end subroutine ASCII_init1D_saving
      
  !---------------------------------------------- 
  ! Output saving of conservation laws 
  !  -> directly in HDF5 binary format
  !----------------------------------------------   
  subroutine HDF5_CL_saving(geom,idiag_num)
    use globals       , only : nbelectrons_diag, nbions_diag
    use physics_module, only : L2norm_diag, &
      entropy_diag, Enkin_diag, Enpot_diag
    type(geometry), intent(in) :: geom
    integer       , intent(in) :: idiag_num
     
#ifndef NOHDF5
    character(LEN=50) :: conservation_laws_file    
    integer(HID_T)    :: file_id   
    integer           :: Nbr, Nbc, icr
    integer           :: ierr
      
    real(RKIND)       :: rdiag_level
      
    if (idiag_num.ne.-1) then
      write(conservation_laws_file,'(A,'//numfmt//',A)') &
        "conservation_laws", idiag_num, ".h5"
    else
      conservation_laws_file = "conservation_laws_init.h5"
    end if
    call HDF5_create(trim(conservation_laws_file),file_id,ierr)
    if (ierr.ne.0) then
      print*,'pglobal_id = ',pglobal_id, &
        ' ==> error for opening of ',conservation_laws_file
    end if
      
    Nbr = 1
    Nbc = 1
    ! -> save : 'time_diag   ' 
    call HDF5_real_saving(file_id,time_diag,'time_diag'//char(0))
    ! -> save : 'deltat   ' 
    call HDF5_real_saving(file_id,deltat,'deltat'//char(0))
    ! -> save : 'diag_level   ' 
    rdiag_level = MIN(4.,float(diag_level))
    call HDF5_real_saving(file_id,rdiag_level,'diag_level'//char(0))
    ! -> save : 'nbions_diag  '
    call HDF5_real_saving(file_id,nbions_diag,'nbions'//char(0))
    ! -> save : 'nbelectrons_diag '
    call HDF5_real_saving(file_id,nbelectrons_diag, &
      'nbelec'//char(0))
    if (diag_level.gt.1) then
      ! -> save : 'Enkin_diag   '
    call HDF5_real_saving(file_id,Enkin_diag,'Enkin'//char(0))
      ! -> save : 'Enpot_diag   '
    call HDF5_real_saving(file_id,Enpot_diag,'Enpot'//char(0))
      ! -> save : 'entropy_diag   '
    call HDF5_real_saving(file_id,entropy_diag,'entropy'//char(0))
      ! -> save : 'L2norm_diag '
    call HDF5_real_saving(file_id,L2norm_diag,'L2norm'//char(0))
    end if
    call HDF5_close(file_id)
#else
    call ASCII_CL_saving(geom,idiag_num)
#endif
  end subroutine HDF5_CL_saving
      
  !---------------------------------------------- 
  ! Output saving of conservation laws 
  !  -> in ASCII format
  !----------------------------------------------   
  subroutine ASCII_CL_saving(geom,idiag_num)  
    use globals, only : nbelectrons_diag, nbions_diag
    use physics_module, only : L2norm_diag, &
      entropy_diag, Enkin_diag, Enpot_diag
    type(geometry), intent(in) :: geom
    integer       , intent(in) :: idiag_num
      
    integer, parameter :: ucl = 31     
    character(LEN=30)  :: conservation_laws_file    
      
    if (idiag_num.ne.-1) then
      write(conservation_laws_file,'(A,'//numfmt//',A)') &
        "conservation_laws", idiag_num, ".dat"
    else
      conservation_laws_file = "conservation_laws_init.dat"
    end if
    open(ucl, file = conservation_laws_file, status = 'UNKNOWN', &
      form = 'FORMATTED')
    write(ucl,'(A12,6A20)') &
      'time_diag','nbions','nbelectrons','Enkin', &
      'Enpot','entropy','entropy','L2norm'
    write(ucl,'(1pe12.3,6(1pe20.12))') &
      time_diag, nbions_diag, nbelectrons_diag, Enkin_diag, &
      Enpot_diag, entropy_diag, L2norm_diag
    close(ucl)
  end subroutine ASCII_CL_saving
      
  !---------------------------------------------- 
  ! Output saving of radial profiles linked to 
  !  guiding-centers
  !  -> directly in HDF5 binary format
  !----------------------------------------------   
  subroutine HDF5_rprof_GC_saving(geom,idiag_num)  
    use physics_module, only : nvpoloGCr_turb_diag, &
      nvpoloGCr_neo_diag, nvpoloGCr_vpar_diag, &
      LphiGCr_diag, TpolGCr_diag, GammaGCr_turb_diag, &
      GammaGCr_neo_diag, RSphiGCr_turb_diag, &
      RSphiGCr_neo_diag, QGCr_turb_diag, QGCr_neo_diag, &
      niGCr_diag, pressGCr_diag, pressGCr_perp_diag, &
      pressGCr_par_diag, dpressGCr_dt_diag, dLphiGCr_dt_diag, &
      QGCr_perp_turb_diag, QGCr_perp_neo_diag, &
      QGCr_par_turb_diag, QGCr_par_neo_diag, &
      QGCr_dtvpar_diag, &
      dpGCr_dt_perp_diag, dpGCr_dt_par_diag
    type(geometry), intent(in) :: geom
    integer       , intent(in) :: idiag_num
     
#ifndef NOHDF5
    character(LEN=30) :: rprof_GC_file    
    integer(HID_T)    :: file_id   
    integer           :: Nbr, Nbc, icr
    real(RKIND)       :: rdiag_level
    integer           :: ierr
      
    if (idiag_num.ne.-1) then
      write(rprof_GC_file,'(A,'//numfmt//',A)') &
        "rprof_GC", idiag_num, ".h5"
    else
      rprof_GC_file = "rprof_GC_init.h5"
    end if
    call HDF5_create(trim(rprof_GC_file),file_id,ierr)
    if (ierr.ne.0) then
      print*,'pglobal_id = ',pglobal_id, &
        ' ==> error for opening of ',rprof_GC_file
    end if
      
    ! -> save : 'time_diag   ' 
    Nbr   = 1
    Nbc   = 1
    call HDF5_real_saving(file_id,time_diag,'time_diag'//char(0))
    ! -> save : 'deltat   ' 
    call HDF5_real_saving(file_id,deltat,'deltat'//char(0))
    ! -> save : 'diag_level   ' 
    rdiag_level = MIN(4.,float(diag_level))
    call HDF5_real_saving(file_id,rdiag_level,'diag_level'//char(0))
    Nbr   = 1
    Nbc   = geom%Nr+1
    ! -> save : 'rg    '
    call HDF5_array1D_saving(file_id,geom%rg,Nbc,'rg'//char(0))
    ! -> save : 'nvpoloGCr_turb_diag '
    call HDF5_array1D_saving(file_id,nvpoloGCr_turb_diag,Nbc, &
      'nvpoloGCr_turb'//char(0))
    ! -> save : 'nvpoloGCr_neo_diag '
    call HDF5_array1D_saving(file_id,nvpoloGCr_neo_diag,Nbc, &
      'nvpoloGCr_neo'//char(0))
    ! -> save : 'nvpoloGCr_vpar_diag '
    call HDF5_array1D_saving(file_id,nvpoloGCr_vpar_diag,Nbc,&
      'nvpoloGCr_vpar'//char(0))
    ! -> save : 'LphiGCr_diag   '
    call HDF5_array1D_saving(file_id,LphiGCr_diag,Nbc, &
      'LphiGCr'//char(0))
    ! -> save : 'TpolGCr   '
    call HDF5_array1D_saving(file_id,TpolGCr_diag,Nbc, &
      'TpolGCr'//char(0))
    ! -> save : 'GammaGCr_turb_diag   '
    call HDF5_array1D_saving(file_id,GammaGCr_turb_diag,Nbc, &
      'GammaGCr_turb'//char(0))
    ! -> save : 'GammaGCr_neo_diag   '
    call HDF5_array1D_saving(file_id,GammaGCr_neo_diag,Nbc, &
      'GammaGCr_neo'//char(0))
    ! -> save : 'RSphiGCr_turb_diag   '
    call HDF5_array1D_saving(file_id,RSphiGCr_turb_diag,Nbc, &
      'RSphiGCr_turb'//char(0))
    ! -> save : 'RSphiGCr_neo_diag   '
    call HDF5_array1D_saving(file_id,RSphiGCr_neo_diag,Nbc, &
      'RSphiGCr_neo'//char(0))
    ! -> save : 'QGCr_turb_diag   '
    call HDF5_array1D_saving(file_id,QGCr_turb_diag,Nbc, &
      'QGCr_turb'//char(0))
    ! -> save : 'QGCr_neo_diag   '
    call HDF5_array1D_saving(file_id,QGCr_neo_diag,Nbc, &
      'QGCr_neo'//char(0))
    ! -> save : 'niGCr_diag   '
    call HDF5_array1D_saving(file_id,niGCr_diag,Nbc, &
      'niGCr'//char(0))
    ! -> save : 'pressGCr_diag   '
    call HDF5_array1D_saving(file_id,pressGCr_diag,Nbc, &
      'pressGCr'//char(0))
    ! -> save : 'pressGCr_perp_diag   '
    call HDF5_array1D_saving(file_id,pressGCr_perp_diag,Nbc, &
      'pressGCr_perp'//char(0))
    ! -> save : 'pressGCr_par_diag   '
    call HDF5_array1D_saving(file_id,pressGCr_par_diag,Nbc, &
      'pressGCr_par'//char(0))
    ! -> save : 'dpressGCr_dt_diag   '
    call HDF5_array1D_saving(file_id,dpressGCr_dt_diag,Nbc, &
      'dpressGCr_dt'//char(0))
    ! -> save : 'dLphiGCr_dt_diag   '
    call HDF5_array1D_saving(file_id,dLphiGCr_dt_diag,Nbc, &
      'dLphiGCr_dt'//char(0))
    ! -> save : 'QGCr_perp_turb_diag   '
    call HDF5_array1D_saving(file_id,QGCr_perp_turb_diag, &
      Nbc,'QGCr_perp_turb'//char(0))
    ! -> save : 'QGCr_perp_neo_diag   '
    call HDF5_array1D_saving(file_id,QGCr_perp_neo_diag, &
      Nbc,'QGCr_perp_neo'//char(0))
    ! -> save : 'QGCr_par_turb_diag   '
    call HDF5_array1D_saving(file_id,QGCr_par_turb_diag,Nbc, &
      'QGCr_par_turb'//char(0))
    ! -> save : 'QGCr_par_neo_diag   '
    call HDF5_array1D_saving(file_id,QGCr_par_neo_diag,Nbc, &
      'QGCr_par_neo'//char(0))
    ! -> save : 'QGCr_dtvpar_diag   '
    call HDF5_array1D_saving(file_id,QGCr_dtvpar_diag,Nbc, &
      'QGCr_dtvpar'//char(0))
    ! -> save : 'dpGCr_dt_perp_diag   '
    call HDF5_array1D_saving(file_id,dpGCr_dt_perp_diag,Nbc, &
      'dpGCr_dt_perp'//char(0))
    ! -> save : 'dpGCr_dt_par_diag   '
    call HDF5_array1D_saving(file_id,dpGCr_dt_par_diag,Nbc, &
      'dpGCr_dt_par'//char(0))
    call HDF5_close(file_id)
#else
    call ASCII_rprof_GC_saving(geom,idiag_num)
#endif
  end subroutine HDF5_rprof_GC_saving
      
  !---------------------------------------------- 
  ! Output saving of radial profiles linked to 
  !  particles
  !  -> directly in HDF5 binary format
  !----------------------------------------------   
  subroutine HDF5_rprof_part_saving(geom,idiag_num)  
    use physics_module, only : nir_diag, nuparr_diag, &
      pressr_diag, nvpolor_mag_diag, Gammar_turb_diag, &
      Gammar_neo_diag, Qr_turb_diag, Qr_neo_diag, &
      RSthetar_diag, RSparr_diag
    use poisson_class , only : Phi00_diag
    type(geometry), intent(in) :: geom
    integer       , intent(in) :: idiag_num
     
#ifndef NOHDF5
    character(LEN=30) :: rprof_part_file    
    integer(HID_T)    :: file_id   
    integer           :: Nbr, Nbc, icr
    real(RKIND)       :: rdiag_level
    integer           :: ierr
      
    if (idiag_num.ne.-1) then
      write(rprof_part_file,'(A,'//numfmt//',A)') &
        "rprof_part", idiag_num, ".h5"
    else
      rprof_part_file = "rprof_part_init.h5"
    end if
    call HDF5_create(trim(rprof_part_file),file_id,ierr)
    if (ierr.ne.0) then
      print*,'pglobal_id = ',pglobal_id, &
        ' ==> error for opening of ',rprof_part_file
    end if
      
    ! -> save : 'time_diag   ' 
    Nbr   = 1
    Nbc   = 1
    call HDF5_real_saving(file_id,time_diag,'time_diag'//char(0))
    ! -> save : 'deltat   ' 
    call HDF5_real_saving(file_id,deltat,'deltat'//char(0))
    ! -> save : 'diag_level   ' 
    rdiag_level = MIN(4.,float(diag_level))
    call HDF5_real_saving(file_id,rdiag_level,'diag_level'//char(0))
    Nbr   = 1
    Nbc   = geom%Nr+1
    ! -> save : 'rg    '
    call HDF5_array1D_saving(file_id,geom%rg,Nbc,'rg'//char(0))
    ! -> save : 'nir_diag    '
    call HDF5_array1D_saving(file_id,nir_diag,Nbc,'nir'//char(0))
    ! -> save : 'nuparr_diag    '
    call HDF5_array1D_saving(file_id,nuparr_diag,Nbc, &
      'nuparr'//char(0))
    ! -> save : 'pressr_diag '
    call HDF5_array1D_saving(file_id,pressr_diag,Nbc, &
      'pressr'//char(0))
    ! -> save : 'nvpolor_mag_diag '
    call HDF5_array1D_saving(file_id,nvpolor_mag_diag,Nbc, &
      'nvpolor_mag'//char(0))
    ! -> save : 'Phi00_diag '
    call HDF5_array1D_saving(file_id,Phi00_diag,Nbc, &
      'Phi00'//char(0))
    ! -> save : 'Gammar_turb_diag   '
    call HDF5_array1D_saving(file_id,Gammar_turb_diag,Nbc, &
      'Gammar_turb'//char(0))
    ! -> save : 'Gammar_neo_diag   '
    call HDF5_array1D_saving(file_id,Gammar_neo_diag,Nbc, &
      'Gammar_neo'//char(0))
    ! -> save : 'Qr_turb_diag   '
    call HDF5_array1D_saving(file_id,Qr_turb_diag,Nbc, &
      'Qr_turb'//char(0))
    ! -> save : 'Qr_neo_diag   '
    call HDF5_array1D_saving(file_id,Qr_neo_diag,Nbc, &
      'Qr_neo'//char(0))
    ! -> save : 'RSthetar_diag   '
    call HDF5_array1D_saving(file_id,RSthetar_diag,Nbc, &
      'RSthetar'//char(0))
    ! -> save : 'RSparr_diag   '
    call HDF5_array1D_saving(file_id,RSparr_diag,Nbc, &
      'RSparr'//char(0))
    call HDF5_close(file_id)
#else
    call ASCII_rprof_part_saving(geom,idiag_num)
#endif
  end subroutine HDF5_rprof_part_saving
      
  !---------------------------------------------- 
  ! Output saving of radial profiles linked 
  !  to guiding-centers
  !  -> in ASCII format
  !----------------------------------------------   
  subroutine ASCII_rprof_GC_saving(geom,idiag_num)  
    use physics_module, only : nvpoloGCr_turb_diag, &
      nvpoloGCr_neo_diag, nvpoloGCr_vpar_diag, &
      LphiGCr_diag, TpolGCr_diag, GammaGCr_turb_diag, &
      GammaGCr_neo_diag, RSphiGCr_turb_diag, &
      RSphiGCr_neo_diag
      
    type(geometry), intent(in) :: geom
    integer       , intent(in) :: idiag_num
      
    integer, parameter :: uprof = 32
    character(LEN=30)  :: rprof_GC_file    
    integer            :: ir
      
    if (idiag_num.ne.-1) then
      write(rprof_GC_file,'(A,'//numfmt//',A)') &
        "rprof_GC", idiag_num, ".dat"
    else
      rprof_GC_file = "rprof_GC_init.dat"
    end if
    open(uprof, file = rprof_GC_file, status = 'UNKNOWN', &
      form = 'FORMATTED')
    write(uprof,'(A12,9A20)') &
      'rg', 'nvpoloGCr_turb', 'nvpoloGCr_neo', 'nvpoloGCr_vpar', &
      'LphiGCr', 'TpolGCr', 'GammaGCr_turb', 'GammaGCr_neo', &
      'RSphiGCr_turb', 'RSphiGCr_neo'
    do ir = 0,geom%Nr
      write(uprof,'(1pe12.3,9(1pe20.12))') &
        geom%rg(ir), nvpoloGCr_turb_diag(ir), &
        nvpoloGCr_neo_diag(ir), nvpoloGCr_vpar_diag(ir), &
        LphiGCr_diag(ir), TpolGCr_diag(ir), &
        GammaGCr_turb_diag(ir), GammaGCr_neo_diag(ir), &    
        RSphiGCr_turb_diag(ir), RSphiGCr_neo_diag(ir) 
    end do
    close(uprof)
  end subroutine ASCII_rprof_GC_saving
      
  !---------------------------------------------- 
  ! Output saving of radial profiles linked 
  !  to particles
  !  -> in ASCII format
  !----------------------------------------------   
  subroutine ASCII_rprof_part_saving(geom,idiag_num)  
    use physics_module, only : nir_diag, nuparr_diag, &
      pressr_diag, nvpolor_mag_diag, Gammar_turb_diag, &
      Gammar_neo_diag, Qr_turb_diag, Qr_neo_diag, &
      RSthetar_diag, RSparr_diag
    use poisson_class        , only : Phi00_diag
      
    type(geometry), intent(in) :: geom
    integer       , intent(in) :: idiag_num
      
    integer, parameter :: uprof = 32
    character(LEN=30)  :: rprof_part_file    
    integer            :: ir
      
    if (idiag_num.ne.-1) then
      write(rprof_part_file,'(A,'//numfmt//',A)') &
        "rprof_part", idiag_num, ".dat"
    else
      rprof_part_file = "rprof_part_init.dat"
    end if
    open(uprof, file = rprof_part_file, status = 'UNKNOWN', &
      form = 'FORMATTED')
    write(uprof,'(A12,11A20)') &
      'rg', 'Phi00', 'nir', 'nuparr', 'pressr', &
      'nvpolor_mag', 'Qr_turb', 'Qr_neo', 'Gammar_turb', &
      'Gammar_neo', 'RSthetar', 'RSparr'
    do ir = 0,geom%Nr
      write(uprof,'(1pe12.3,11(1pe20.12))') &
        geom%rg(ir), Phi00_diag(ir), &
        nir_diag(ir), nuparr_diag(ir), &
        pressr_diag(ir), nvpolor_mag_diag(ir), &
        Qr_turb_diag(ir), Qr_neo_diag(ir), &
        Gammar_turb_diag(ir), Gammar_neo_diag(ir), &
        RSthetar_diag(ir), RSparr_diag(ir)
    end do
    close(uprof)
  end subroutine ASCII_rprof_part_saving
      
  !---------------------------------------------- 
  ! Output saving of 2D cross-sections of Phi
  !  -> directly in HDF5 binary format
  !----------------------------------------------   
  subroutine HDF5_Phi2D_saving(geom,idiag_num)  
    use Pcross_section_module, only : Phirth_diag, &
      Phirphi_diag, Phithphi_diag
    type(geometry), intent(in) :: geom
    integer       , intent(in) :: idiag_num
     
#ifndef NOHDF5
    character(LEN=30) :: Phi2D_file    
    integer(HID_T)    :: file_id        
    integer           :: Nbr, Nbc, icr
    integer           :: ierr
      
    if (idiag_num.ne.-1) then
      write(Phi2D_file,'(A,'//numfmt//',A)') &
        "Phi2D", idiag_num, ".h5"
    else
      Phi2D_file = "Phi2D_init.h5"
    end if
    call HDF5_create(trim(Phi2D_file),file_id,ierr)
    if (ierr.ne.0) then
      print*,'pglobal_id = ',pglobal_id, &
        ' ==> error for opening of ',Phi2D_file
    end if
      
    ! -> save ir_Phi, itheta_Phi, iphi_Phi
    call HDF5_integer_saving(file_id,ir_Phi,'ir_Phi'//char(0))
    call HDF5_integer_saving(file_id,itheta_Phi, &
      'itheta_Phi'//char(0))
    call HDF5_integer_saving(file_id,iphi_Phi,'iphi_Phi'//char(0))
      
    ! -> save : 'time_diag   ' 
    call HDF5_real_saving(file_id,time_diag,'time_diag'//char(0))
    ! -> save : 'deltat   ' 
    call HDF5_real_saving(file_id,deltat,'deltat'//char(0))
    ! -> save : 'Phirth_diag '
    Nbr    = geom%Nr+1
    Nbc    = geom%Ntheta+1
    call HDF5_array2D_saving(file_id,Phirth_diag(0:,0:),Nbr,Nbc, &
      'Phirth'//char(0))
    ! -> save : 'Phirphi_diag   '
    Nbr    = geom%Nr+1
    Nbc    = geom%Nphi+1
    call HDF5_array2D_saving(file_id,Phirphi_diag(0:,0:),Nbr,Nbc, &
      'Phirphi'//char(0))
    ! -> save : 'Phithphi_diag  '
    Nbr    = geom%Ntheta+1
    Nbc    = geom%Nphi+1
    call HDF5_array2D_saving(file_id,Phithphi_diag(0:,0:),Nbr,Nbc, &
      'Phithphi'//char(0))
      
    call HDF5_close(file_id)
#endif
  end subroutine HDF5_Phi2D_saving
      
  !---------------------------------------------- 
  ! Output saving of 2D cross-sections of f
  !  -> directly in HDF5 binary format
  !----------------------------------------------   
  subroutine HDF5_f2D_saving(geom,idiag_num)
    use globals, only : ir_f2D_trapped, itheta_f2D_trapped, &
      iphi_f2D_trapped, ivpar_f2D_trapped, imu_f2D_trapped, &
      ir_f2D_passing, itheta_f2D_passing, iphi_f2D_passing, &
      ivpar_f2D_passing, imu_f2D_passing
    use Pcross_section_module, only : frtheta_v0_diag, &
      frtheta_vmax_diag, fphivpar_mu0_diag, fphivpar_mumax_diag, &
      frvpar_mu0_diag, frvpar_mumax_diag, &
      fthetavpar_mu0_diag, fthetavpar_mumax_diag, fvparmu_diag
    type(geometry), intent(in) :: geom
    integer       , intent(in) :: idiag_num
     
#ifndef NOHDF5
    character(LEN=30) :: f2D_file    
    integer(HID_T)    :: file_id        
    integer           :: Nbr, Nbc, icr
    integer           :: ierr
      
    if (idiag_num.ne.-1) then
      write(f2D_file,'(A,'//numfmt//',A)') &
        "f2D", idiag_num, ".h5"
    else
      f2D_file = "f2D_init.h5"
    end if
    call HDF5_create(trim(f2D_file),file_id,ierr)
    if (ierr.ne.0) then
      print*,'pglobal_id = ',pglobal_id, &
        ' ==> error for opening of ',f2D_file
    end if
      
    ! -> save : ir_diag, itheta_diag, iphi_diag, ivpar_diag, &
    !    imu_diag for trapped and passing particules
    call HDF5_integer_saving(file_id,ir_f2D_trapped, &
      'ir_f2D_trapped'//char(0))
    call HDF5_integer_saving(file_id,itheta_f2D_trapped, &
      'itheta_f2D_trapped'//char(0))
    call HDF5_integer_saving(file_id,iphi_f2D_trapped, &
      'iphi_f2D_trapped'//char(0))
    call HDF5_integer_saving(file_id,ivpar_f2D_trapped, &
      'ivpar_f2D_trapped'//char(0))
    call HDF5_integer_saving(file_id,imu_f2D_trapped, &
      'imu_f2D_trapped'//char(0))
    call HDF5_integer_saving(file_id,ir_f2D_passing, &
      'ir_f2D_passing'//char(0))
    call HDF5_integer_saving(file_id,itheta_f2D_passing, &
      'itheta_f2D_passing'//char(0))
    call HDF5_integer_saving(file_id,iphi_f2D_passing, &
      'iphi_f2D_passing'//char(0))
    call HDF5_integer_saving(file_id,ivpar_f2D_passing, &
      'ivpar_f2D_passing'//char(0))
    call HDF5_integer_saving(file_id,imu_f2D_passing, &
      'imu_f2D_passing'//char(0))
    ! -> save : 'time_diag   ' 
    call HDF5_real_saving(file_id,time_diag,'time_diag'//char(0))
    ! -> save : 'deltat   ' 
    call HDF5_real_saving(file_id,deltat,'deltat'//char(0))
    ! -> save : 'frth_v0_diag  '
    Nbr    = geom%Nr+1
    Nbc    = geom%Ntheta+1
    call HDF5_array2D_saving(file_id,frtheta_v0_diag(0:,0:), &
      Nbr,Nbc,'frthv0'//char(0))
    ! -> save : 'frth_vm_diag  '
    Nbr    = geom%Nr+1
    Nbc    = geom%Ntheta+1
    call HDF5_array2D_saving(file_id,frtheta_vmax_diag(0:,0:), &
      Nbr,Nbc,'frthvm'//char(0))
    ! -> save : 'fphivpar_mu0_diag '
    Nbr    = geom%Nphi+1
    Nbc    = geom%Nvpar+1
    call HDF5_array2D_saving(file_id,fphivpar_mu0_diag(0:,0:), &
      Nbr,Nbc,'fphivpar_mu0'//char(0))
    ! -> save : 'fphivpar_mumax_diag '
    Nbr    = geom%Nphi+1
    Nbc    = geom%Nvpar+1
    call HDF5_array2D_saving(file_id,fphivpar_mumax_diag(0:,0:), &
      Nbr,Nbc,'fphivpar_mumax'//char(0))
    ! -> save : 'frvpar_mu0_diag '
    Nbr    = geom%Nr+1
    Nbc    = geom%Nvpar+1
    call HDF5_array2D_saving(file_id,frvpar_mu0_diag(0:,0:), &
      Nbr,Nbc,'frvpar_mu0'//char(0))
    ! -> save : 'frvpar_mumax_diag '
    Nbr    = geom%Nr+1
    Nbc    = geom%Nvpar+1
    call HDF5_array2D_saving(file_id,frvpar_mumax_diag(0:,0:), &
      Nbr,Nbc,'frvpar_mumax'//char(0))
    ! -> save : 'fthetavpar_mu0_diag '
    Nbr    = geom%Ntheta+1
    Nbc    = geom%Nvpar+1
    call HDF5_array2D_saving(file_id,fthetavpar_mu0_diag(0:,0:), &
      Nbr,Nbc,'fthvpar_mu0'//char(0))
    ! -> save : 'fthetavpar_mumax_diag '
    Nbr    = geom%Ntheta+1
    Nbc    = geom%Nvpar+1
    call HDF5_array2D_saving(file_id,fthetavpar_mumax_diag(0:,0:), &
      Nbr,Nbc,'fthvpar_mumax'//char(0))
    ! -> save : 'fvparmu_diag '
    Nbr    = geom%Nvpar+1
    Nbc    = geom%Nmu+1
    call HDF5_array2D_saving(file_id,fvparmu_diag(0:,0:), &
      Nbr,Nbc,'fvparmu'//char(0))
      
    call HDF5_close(file_id)
#endif
  end subroutine HDF5_f2D_saving
      
  !-----------------------------------------------
  ! - Initialisation of the global variables for
  !  the iteration, diagnostics ...
  ! - Initialisation of the output file
  !  and saving of the mesh 
  ! - Allocation for the output array saving
  !-----------------------------------------------
  subroutine output_first_writing(geom)
    use globals
    use geometry_class
    implicit none
    type(geometry), intent(in)  :: geom
      
    if (.not. restart) then
      nb_diag = 0
    else
      nb_diag = -1
    end if
      
    ir_save     = int(geom%Nr/2)
    itheta_save = int(geom%Ntheta/3)
    iphi_save   = int(geom%Nphi/3)
  end subroutine output_first_writing
end module output_saving_module
