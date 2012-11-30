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
      
!-----------------------------------------------------
! file : resu3D_saving.f90
! date : 25/09/2008
!  Saving of the 3D results in HDF5 format, like :
!   - the electrostatic potential Phi(r,theta,phi) 
!   - the fluid moments : J0M0(r,theta,phi) and
!       J0M2(r,theta,phi)
!  Rk : The 3D arrays are supposed parallelized 
!       and each processor saves its local part
!-----------------------------------------------------
module resu3D_saving_module
  use prec_const
#ifndef NOHDF5
  use HDF5
  use HDF5_io_module
#endif
      
  implicit none
      
  real(RKIND), dimension(:,:,:), pointer, private :: resu3D_loc_HDF5
      
  !******************************
  contains
  !******************************
      
  !-------------------------------------------------
  ! Constructor 
  !-------------------------------------------------
  subroutine new_resu3D_saving()
    use globals, only : istart, iend, inumber, &
      jstart, jend, jnumber, Nphi
      
    call glob_allocate(resu3D_loc_HDF5,0,inumber-1, &
      0,jnumber-1,0,Nphi-1,'resu3D_loc_HDF5')
  end subroutine new_resu3D_saving
      
  !-------------------------------------------------
  ! Destructor 
  !-------------------------------------------------
  subroutine delete_resu3D_saving
    call glob_deallocate(resu3D_loc_HDF5)
  end subroutine delete_resu3D_saving
      
  !************************************************************
  ! Reading and Writing when the 3D array is not distributed
  !  on several processors
  !************************************************************
  !------------------------------------------------------------- 
  ! saves 3D results associated to the diagnostic "diag_num" 
  !  (where the 3D array is on one unique processor)
  !-------------------------------------------------------------
  subroutine HDF5_resu3D_saving(var3D_name,diag_num,geom,resu3D)
    use globals, only : R0, time_diag
    use geometry_class
    character*(*)                , intent(in) :: var3D_name
    integer                      , intent(in) :: diag_num
    type(geometry)               , intent(in) :: geom
    real(RKIND), dimension(:,:,:), pointer    :: resu3D
      
#ifndef NOHDF5
    integer           :: i, j, iphi
    character(LEN=50) :: file_name
    character(LEN=50) :: var3D_name_tmp
    integer(HID_T)    :: file_id   ! file identifier
      
    !*** saving of the local part of resu3D ***
    var3D_name_tmp = trim(var3D_name)//'_3D'
    file_name      = create_filename3D(var3D_name_tmp,diag_num)
    call HDF5_create(trim(file_name),file_id)
    call HDF5_real_saving(file_id,time_diag,"time_diag")
    call HDF5_integer_saving(file_id,geom%Nr,"Nr")
    call HDF5_integer_saving(file_id,geom%Ntheta,"Ntheta")
    call HDF5_integer_saving(file_id,geom%Nphi,"Nphi")
    call HDF5_real_saving(file_id,R0,"R0")
    call HDF5_array1D_saving(file_id,geom%rg(0:),geom%Nr+1,"rg")
    call HDF5_array1D_saving(file_id,geom%thetag(0:), &
      geom%Ntheta+1,"thetag")
    call HDF5_array1D_saving(file_id,geom%phig(0:), &
      geom%Nphi+1,"phig")
    call HDF5_array3D_saving(file_id,resu3D,geom%Nr+1, &
      geom%Ntheta+1,geom%Nphi+1,trim(var3D_name_tmp))
    call HDF5_close(file_id)
#endif
  end subroutine HDF5_resu3D_saving
      
  !************************************************************
  ! Reading and Writing when the 3D array is distributed
  !  on several processors
  !************************************************************
  !----------------------------------------------------
  ! saving of the HDF5 master file for a specific
  !  3D result, with the name '<var3D_name>3D_master.h5'
  !----------------------------------------------------
  subroutine HDF5_master_resu3D(var3D_name,diag_num,geom,nb_proc)
    use globals, only : R0, time_diag
    use geometry_class
    character*(*) , intent(in) :: var3D_name
    integer       , intent(in) :: diag_num
    type(geometry), intent(in) :: geom
    integer       , intent(in) :: nb_proc
      
#ifndef NOHDF5
    integer           :: len_varname
    integer           :: iproc, error       
    integer(SIZE_T)   :: isize
    integer(HID_T)    :: file_id   
    integer(HID_T)    :: group_id, group_proc_id
    character(LEN=50) :: master_filename  
    character(LEN=30) :: diag_rep 
    character(LEN=30) :: proc_name 
    character(LEN=60) :: proc_rep  
      
    !*** Open HDF5 file ***
    call H5open_f(error)
      
    !*** Create a new file using default properties ***
    len_varname     = len(trim(var3D_name))
    master_filename = trim(var3D_name)// &
      "3D_master_d    .h5"
    write(master_filename(len_varname+12:len_varname+15), &
      '(i4.4)') diag_num    
    call H5Fcreate_f(trim(master_filename), &
      H5F_ACC_TRUNC_F,file_id, error)     
      
    !*** Saving of the global geometry ***
    call HDF5_real_saving(file_id,time_diag,"time_diag")
    call HDF5_integer_saving(file_id,geom%Nr,"Nr")
    call HDF5_integer_saving(file_id,geom%Ntheta-1,"Ntheta")
    call HDF5_integer_saving(file_id,geom%Nphi-1,"Nphi")
    call HDF5_integer_saving(file_id,nb_proc,"nb_proc")
    call HDF5_real_saving(file_id,R0,"R0")
    call HDF5_array1D_saving(file_id,geom%rg(0:),geom%Nr+1,"rg")
    call HDF5_array1D_saving(file_id,geom%thetag(0:), &
      geom%Ntheta,"thetag")
    call HDF5_array1D_saving(file_id,geom%phig(0:), &
      geom%Nphi,"phig")
      
    !*** Creating of the group : diag_num ***
    diag_rep = "/diag_    "
    write(diag_rep(7:10),'(i4.4)') diag_num
    !-> 10*nbproc because nb_proc*/proc_<num>
    isize = 10*nb_proc
    call H5Gcreate_f(file_id,trim(diag_rep),group_id, &
      error,size_hint=isize)
    isize = 0
    do iproc = 0,nb_proc-1
      proc_name = "/proc_    "
      write(proc_name(7:10),'(i4.4)') iproc
      proc_rep  = trim(diag_rep)//trim(proc_name)
      call H5Gcreate_f(group_id,trim(proc_rep), &
        group_proc_id,error,size_hint=isize)
      call H5Gset_comment_f(group_proc_id, &
        trim(proc_rep),"processor number",error)
      call H5Gclose_f(group_proc_id,error)
    end do
    call H5Gclose_f(group_id,error)    
    call H5Fclose_f(file_id,error)
#endif
  end subroutine HDF5_master_resu3D
      
  !------------------------------------------------------------- 
  ! saves 3D results associated to the diagnostic "diag_num" 
  !  (where the 3D array is distributed on several processors)
  !-------------------------------------------------------------
  subroutine HDF5_r3D_saving(var3D_name,diag_num,geom,resu3D)
    use globals
    use geometry_class
    character*(*)                , intent(in) :: var3D_name
    integer                      , intent(in) :: diag_num
    type(geometry)               , intent(in) :: geom
    real(RKIND), dimension(:,:,:), pointer    :: resu3D
      
#ifndef NOHDF5
    integer           :: i, j, iphi
    character(LEN=50) :: file_name
    character(LEN=50) :: var3D_name_loc, var3D_name_tmp
    integer(HID_T)    :: file_id   ! file identifier
      
    real(RKIND), &
      dimension(0:inumber-1,0:jnumber-1,0:Nphi-1) :: resu3D_tmp
      
    !*** saving of the local part of resu3D ***
    if (mu_id.eq.0) then
      do iphi = 0,Nphi-1
        do j = jstart,jend
          do i = istart,iend
            resu3D_tmp(i-istart,j-jstart,iphi) = resu3D(i,j,iphi)
          end do
        end do
      end do
      var3D_name_tmp = trim(var3D_name)//'3D'
      var3D_name_loc = trim(var3D_name)//'3D_loc'
      file_name      = create_filename3D(var3D_name_tmp, &
        diag_num,plocal_id)
      call HDF5_create(trim(file_name),file_id)
      call HDF5_integer_saving(file_id,istart,"istart")
      call HDF5_integer_saving(file_id,iend,"iend")
      call HDF5_integer_saving(file_id,jstart,"jstart")
      call HDF5_integer_saving(file_id,jend,"jend")
      call HDF5_integer_saving(file_id,Nphi-1,"Nphi")
      call HDF5_array1D_saving(file_id,geom%rg(istart:), &
        inumber,"rg_loc")
      call HDF5_array1D_saving(file_id,geom%thetag(jstart:), &
        jnumber,"thetag_loc")
      call HDF5_array1D_saving(file_id,geom%phig(0:),Nphi,"phig")
      call HDF5_array3D_saving(file_id,resu3D_tmp, &
        inumber,jnumber,Nphi,trim(var3D_name_loc))
      call HDF5_close(file_id)
    end if
      
    !*** writting of the resu3D_master.h5 file ***
    if (pglobal_id.eq.0) then
      call HDF5_master_resu3D(trim(var3D_name),diag_num, &
        geom,Nbproc_r*Nbproc_theta)   
    end if
#endif
  end subroutine HDF5_r3D_saving
      
  !----------------------------------------------------------- 
  ! read 3D results associated to the diagnostic "diag_num"
  !  which are stored in "nb_proc" parallel files 
  !  "<var3D_name>3D_d<diag_num>_p<proc_num>.h5" by using 
  !  "<var3D_name>3D_master.h5 file" 
  !-----------------------------------------------------------
  subroutine HDF5_r3D_reading(var3D_name,diag_num,nb_proc)
    use mem_alloc_module
    character*(*), intent(in) :: var3D_name
    integer      , intent(in) :: diag_num
    integer      , intent(in) :: nb_proc
      
#ifndef NOHDF5
    integer(HID_T), &
      dimension(0:nb_proc-1)    :: proc_file_id   ! file identifier
    integer(HID_T)              :: file_id        ! file identifier
    integer                     :: error          ! error flag
    character(LEN=50)           :: master_filename
    character(LEN=30)           :: diag_rep, proc_rep
    character(LEN=30)           :: mount_rep
    character(LEN=50)           :: proc_file_name
    character(LEN=50)           :: var_name
    character(LEN=50)           :: var3D_name_loc, var3D_name_tmp
    integer                     :: iproc, len_varname
                                
    integer                     :: ir, itheta, iphi
    integer                     :: i1, i2, i3
    integer                     :: Nr_tmp, Ntheta_tmp, Nphi_tmp
    integer                     :: istart_tmp, iend_tmp
    integer                     :: jstart_tmp, jend_tmp
    integer                     :: inumber_tmp, jnumber_tmp
    real(RKIND), &
      dimension(:)    , pointer :: rg_tmp, thetag_tmp, phig_tmp
    real(RKIND), &
      dimension(:)    , pointer :: rg_loc, thetag_loc
    real(RKIND), &
      dimension(:)    , pointer :: rg_glob, thetag_glob
    real(RKIND), &
      dimension(:,:,:), pointer :: resu3D_tmp, resu3D_loc
      
    diag_rep = "/diag_    "
    write(diag_rep(7:10),'(i4.4)') diag_num
      
    !*** Open <var_name>3D_d<diag_num>_master.h5 file ***
    !***  and mount all the files                     ***
    !***  <var_name>3D_d<diag_num>_p<proc_num>.h5     ***
    !***  associated                                  ***
    len_varname     = len(trim(var3D_name))
    master_filename = &
      trim(var3D_name)//"3D_master_d    .h5"
    write(master_filename(len_varname+12:len_varname+15), &
      '(i4.4)') diag_num
    call H5Fopen_f(trim(master_filename), &
      H5F_ACC_RDONLY_F,file_id,error)
      
    var3D_name_tmp = trim(var3D_name)//'3D'
    var3D_name_loc = trim(var3D_name)//'3D_loc'
    do iproc = 0,nb_proc-1
       proc_file_name = &
         create_filename3D(trim(var3D_name_tmp),diag_num,iproc)
       call H5Fopen_f(trim(proc_file_name), &
         H5F_ACC_RDONLY_F,proc_file_id(iproc),error)
       proc_rep  = "/proc_    "
       write(proc_rep(7:10),'(i4.4)') iproc
       mount_rep = trim(diag_rep)//trim(proc_rep)
       call H5Fmount_f(file_id,trim(mount_rep), &
         proc_file_id(iproc),error)
    end do
    
    !*** reading of the global geometry data for array allocating ***
    call HDF5_integer_reading(file_id,Nr_tmp,"Nr")
    call HDF5_integer_reading(file_id,Ntheta_tmp,"Ntheta")
    call HDF5_integer_reading(file_id,Nphi_tmp,"Nphi")
    call temp_allocate(rg_tmp,0,Nr_tmp,'rg_tmp')
    call temp_allocate(thetag_tmp,0,Ntheta_tmp,'thetag_tmp')
    call temp_allocate(phig_tmp,0,Nphi_tmp,'phig_tmp')
    call temp_allocate(resu3D_tmp,0,Nr_tmp,0,Ntheta_tmp, &
      0,Nphi_tmp,'resu3D_tmp')
    call HDF5_array1D_reading(file_id,rg_tmp,"rg")
    call HDF5_array1D_reading(file_id,thetag_tmp,"thetag")
    call HDF5_array1D_reading(file_id,phig_tmp,"phig")
      
    !*** global array allocation (used to check the geometry) ***
    call temp_allocate(rg_glob,0,Nr_tmp,'rg_glob')
    call temp_allocate(thetag_glob,0,Ntheta_tmp,'thetag_glob')
      
    !*** reading of each part in the nb_proc files ***
    do iproc = 0,nb_proc-1
      proc_rep  = "/proc_    "
      write(proc_rep(7:10),'(i4.4)') iproc
      mount_rep = trim(diag_rep)//trim(proc_rep)
      ! -> reading of istart
      var_name = create_variable_name(mount_rep,"istart")
      call HDF5_integer_reading(file_id,istart_tmp,trim(var_name))
      ! -> reading of iend
      var_name = create_variable_name(mount_rep,"iend")
      call HDF5_integer_reading(file_id,iend_tmp,trim(var_name))
      ! -> reading of jstart
      var_name = create_variable_name(mount_rep,"jstart")
      call HDF5_integer_reading(file_id,jstart_tmp,trim(var_name))
      ! -> reading of jend
      var_name = create_variable_name(mount_rep,"jend")
      call HDF5_integer_reading(file_id,jend_tmp,trim(var_name))
      ! -> local array allocation
      if (iproc.eq.0) then
        inumber_tmp = iend_tmp-istart_tmp+1
        jnumber_tmp = jend_tmp-jstart_tmp+1
        call temp_allocate(resu3D_loc,0,inumber_tmp-1, &
          0,jnumber_tmp-1,0,Nphi_tmp,'resu3D_loc')
      end if
      ! -> reading of rg(istart:iend)
      var_name = create_variable_name(mount_rep,"rg_loc")
      rg_loc => rg_glob(istart_tmp:)
      call HDF5_array1D_reading(file_id,rg_loc,trim(var_name))
      ! -> reading of thetag(jstart:jend)
      var_name = create_variable_name(mount_rep,"thetag_loc")
      thetag_loc => thetag_glob(jstart_tmp:)
      call HDF5_array1D_reading(file_id,thetag_loc,trim(var_name))
      ! -> reading of resu3D(istart:iend,jstart:jend,0:Nphi-1)
      var_name = create_variable_name(mount_rep, &
        trim(var3D_name_loc))
      call HDF5_array3D_reading(file_id,resu3D_loc,trim(var_name))
      i3 = 0
      do iphi = 0,Nphi_tmp
        i2 = 0
        do itheta = jstart_tmp,jend_tmp
          i1 = 0
          do ir = istart_tmp,iend_tmp
            resu3D_tmp(ir,itheta,iphi) = resu3D_loc(i1,i2,i3)
            i1 = i1 + 1
          end do
          i2 = i2 + 1
        end do
        i3 = i3 + 1
      end do
    end do
    !*** periodic boundary conditions in theta ***
    thetag_glob(Ntheta_tmp)                    = &
      thetag_tmp(Ntheta_tmp)
    resu3D_tmp(0:Nr_tmp,Ntheta_tmp,0:Nphi_tmp) = &
      resu3D_tmp(0:Nr_tmp,0,0:Nphi_tmp)
      
    !*** Closing of all the files ***
    do iproc = 0,nb_proc-1
       call H5Fclose_f(proc_file_id(iproc),error)    
    end do
    call H5Fclose_f(file_id,error)
    
    !*** Array deallocation ***
    call temp_deallocate(rg_tmp)
    call temp_deallocate(thetag_tmp)
    call temp_deallocate(phig_tmp)
    call temp_deallocate(resu3D_tmp)
    call temp_deallocate(rg_glob)
    call temp_deallocate(thetag_glob)
    call temp_deallocate(resu3D_loc)
#endif
  end subroutine HDF5_r3D_reading
end module resu3D_saving_module
