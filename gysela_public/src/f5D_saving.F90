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
      
!---------------------------------------------------------
! file : f5D_saving.f90
! date : 25/01/2006
!  Saving in HDF5 format of the distribution function
!  f(r,theta,phi,vpar) parametrized by a value of mu 
!  -> f is parallelized and each processor saves is 
!     local part
!---------------------------------------------------------
module f5D_saving_module
  use clock_module
  use prec_const
#ifndef NOHDF5
  use HDF5
  use HDF5_io_module
#endif
      
  implicit none
      
  real(RKIND), dimension(:,:,:,:,:), &
    pointer, public :: Rarray5D_loc_HDF5
  
  !******************************
  contains
  !******************************
      
  !-------------------------------------------------
  ! Constructor 
  !-------------------------------------------------
  subroutine new_f5D_saving()
    use globals, only : istart, iend, inumber, &
      jstart, jend, jnumber, Nphi, Nvpar, mu_id
      
    ! -> create a virtual 5D dimension 
    !   (for simplified used of ExtractSlice tools)
    call glob_allocate(Rarray5D_loc_HDF5,0,inumber-1,0,jnumber-1, &
      0,Nphi-1,0,Nvpar,0,0,'Rarray5D_loc_HDF5')
  end subroutine new_f5D_saving
      
  !-------------------------------------------------
  ! Destructor 
  !-------------------------------------------------
  subroutine delete_f5D_saving
    call glob_deallocate(Rarray5D_loc_HDF5)
  end subroutine delete_f5D_saving
      
  !--------------------------------------------- 
  ! saving of the HDF5 f5D_master.h5 file
  !---------------------------------------------
  subroutine HDF5_f5Dmaster_saving(diag_num,geom)
    use globals, only : Nbproc_r, Nbproc_theta, Nbproc_mu, &
      time_diag, R0
    use geometry_class
    integer       , intent(in) :: diag_num
    type(geometry), intent(in) :: geom
      
#ifndef NOHDF5
    integer                         :: error
    integer(HID_T)                  :: file_id   
    character(LEN=50)               :: master_filename
    integer          , dimension(5) :: dimSizeTot
    integer          , dimension(5) :: nbBloc
    character(LEN=30)               :: num, name 
    real(4)          , dimension(5) :: dimMin, dimMax
      
    dimMin(1) = geom%rg(0)
    dimMin(2) = geom%thetag(0)
    dimMin(3) = geom%phig(0)
    dimMin(4) = geom%vparg(0)
    dimMin(5) = geom%mug(0)
      
    dimMax(1) = geom%rg(geom%Nr)
    dimMax(2) = geom%thetag(geom%Ntheta)
    dimMax(3) = geom%phig(geom%Nphi)
    dimMax(4) = geom%vparg(geom%Nvpar)
    dimMax(5) = geom%mug(geom%Nmu)
      
    ! As we are reading dataset from C program, 
    !  dimension are inverted
    dimSizeTot(1) = geom%Nr+1
    dimSizeTot(2) = geom%Ntheta
    dimSizeTot(3) = geom%Nphi
    dimSizeTot(4) = geom%Nvpar+1
    dimSizeTot(5) = geom%Nmu+1
      
    nbBloc(1) = Nbproc_r
    nbBloc(2) = Nbproc_theta
    nbBloc(3) = 1
    nbBloc(4) = 1
    nbBloc(5) = Nbproc_mu
    !*** Open HDF5 file ***
    call H5open_f(error)
      
    !*** Create a new file using default properties ***
    write(num,'(i9)') diag_num
    name="f5D_"//adjustl(num)
    master_filename=trim(name)//".master.h5"
    call H5Fcreate_f(trim(master_filename), &
      H5F_ACC_TRUNC_F,file_id, error)     
      
    !*** Saving of the global geometry ***
    call HDF5_array1D_saving_int(file_id,dimSizeTot,5,"dimSizeTot")
    call HDF5_array1D_saving_int(file_id,nbBloc,5,"nbBloc")
    call HDF5_array1D_saving_r4(file_id,dimMin,5,"dimMin")
    call HDF5_array1D_saving_r4(file_id,dimMax,5,"dimMax")
    call HDF5_real_saving(file_id,time_diag,"tTime")
    call HDF5_real_saving(file_id,R0,"R0")
      
    call H5Fclose_f(file_id,error)
#endif
  end subroutine HDF5_f5Dmaster_saving
      
  !----------------------------------------------------------- 
  ! saves f5D values associated to the diagnostic "diag_num" 
  !-----------------------------------------------------------
  subroutine HDF5_f5D_saving(diag_num,geom,f4D_mu)
    use globals, only : istart, iend, inumber, &
      jstart, jend, jnumber, mu_id ,&
      Nphi, Nvpar, pglobal_id
    use geometry_class
    integer                        , intent(in) :: diag_num
    type(geometry)                 , intent(in) :: geom
    real(RKIND), dimension(:,:,:,:), pointer    :: f4D_mu
      
#ifndef NOHDF5
    character(LEN=50) :: file_name
    integer(HID_T)    :: file_id   ! file identifier
    real(RKIND), dimension(:,:,:,:,:), pointer :: f5D_tmp
      
    integer :: i, j, iphi, ivpar
      
    call clck_time(bclock_HDF5)
    !*** saving of the local part of f4D_mu ***
    f5D_tmp => Rarray5D_loc_HDF5
    do ivpar = 0,Nvpar
      do iphi = 0,Nphi-1
        do j = jstart,jend
          do i = istart,iend
            f5D_tmp(i-istart,j-jstart,iphi,ivpar,0) = &
              f4D_mu(i,j,iphi,ivpar)
          end do
        end do
      end do
    end do
      
    file_name = create_filename5D("f5D",diag_num,pglobal_id)
    call HDF5_create(trim(file_name),file_id)
    call HDF5_array5D_saving(file_id,f5D_tmp,inumber,jnumber,&
      Nphi,Nvpar+1,1,"f5D_loc")
    call HDF5_close(file_id)
      
    !*** writting of the f5D_master.h5 file ***
    if (pglobal_id.eq.0) then
      call HDF5_f5Dmaster_saving(diag_num,geom)    
    end if
#endif
    call clck_time(eclock_HDF5)
    call clck_diff(bclock_HDF5,eclock_HDF5,global_time_5D_saving)
  end subroutine HDF5_f5D_saving
end module f5D_saving_module
