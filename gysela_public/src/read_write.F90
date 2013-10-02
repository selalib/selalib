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
! file : read_write.f90
! date : 16/09/2005
! used for the reading of :
! - mesh datas
! - the equilibrium datas
! - restart file (if in restart mode)
!----------------------------------------
module read_write_module
#ifndef NOHDF5
  use HDF5
  use HDF5_io_module
#endif
  implicit none
      
  !-------------------------------------
  !  PUBLIC VARIABLES
  !-------------------------------------  
  ! Output file names
  character(LEN=100), public       :: parameter_file
  character(LEN=100), public       :: file_list
  integer           , public, save :: upar       = 24
  integer           , public, save :: ulist      = 25
  integer           , public, save :: uout_CL    = 26
  integer           , public, save :: uout_rprof = 27
  integer           , public, save :: uout_r3    = 28
!R3 logical           , public       :: file_r3_open !R3
      
  !******************************
  contains
  !******************************
  !------------------------------------------
  ! write parameters in the text result file
  !------------------------------------------
  subroutine Write_parameters()
    use globals  
    implicit none
      
    namelist /PARALLEL/ Nbproc_r, Nbproc_theta, Nbproc_mu, &
      bloc_phi, Nbthread, transpose4D
    namelist /MESH/ Nr, Ntheta, Nphi, Nvpar, Nmu, &
      a, rhomin, rhomax, Ltheta, Lphi, &
      aspect_ratio, nb_vth0, Lmu, mumin, integral_vperp
    namelist /EQUIL/ restart, plasma_current, &
      rbar_in_feq, use_psibar, &
      Zi, tau0, read_n0, read_Te, read_Ti, read_q, &
      profile_choice, rpeak, kappan, kappaTi, kappaTe, & 
      deltarn, deltarTi, deltarTe, &
      magnetic_shear, q0, deltaq, alphaq, &
      reversed_shear, qmin, qa, rho_qmin, &
      B_curvature, epsilon, n, m, init_perturb, &
      zonal_flow, filter_choice, filter_deltam, &
      ExB_shear, max_shear, rho_shear, deltar_shear
     namelist /ALGORITHM/ leapfrog, LF_average, &
       advec2Ddt, hffiltertype, hffilterfreq, hffilterfref, &
       reducedbegin, reducedend, reduceddt, &
       deltat, nbiter, diag_nbstep, refresh_nbstep
    namelist /TEST/ memory_test, gyroaverage, &
      slab_geometry, Rosenbluth_Hinton, &
      m0n0_eq0, single_m, single_n, solve_poiss, coeff_pol, &
      RHS_only, &
      feq_choice, vpar0_option, Phi00BCr_Neumann, &
      hvpar_in_fperturb
    namelist /OUTPUT/ integration_CS, diag_level, &
      Phi3D_saving, FMoments3D_saving, f5D_saving, &
      CalviExp_saving, Phi3D_nbstep, FMoments3D_nbstep, &
      f5D_nbstep, CalviExp_nbstep, rst_saving
      
    if (pglobal_id.eq.0) then
      write(parameter_file,'(A,'//numfmt_rst//',A)') &
        "gysela_data", nb_restart_diag, ".out"
      open(upar,file = parameter_file,status = 'REPLACE', &
        form = 'FORMATTED')
      write(upar,PARALLEL)
      write(upar,MESH)
      write(upar,EQUIL)
      write(upar,ALGORITHM)
      write(upar,TEST)
      write(upar,OUTPUT)
      close(upar)
    end if
  end subroutine Write_parameters
      
  !------------------------------------------
  ! Read the restart files
  !------------------------------------------
  subroutine read_f5d_rst(filename,initial_time,geom,init_prof, &
    argname,f5d)
    use globals, only : deltat, nbiter_prev, nb_restart, nb_diag, &
      Phi3D_nb_diag, FMoments3D_nb_diag, f5D_nb_diag, mu_id, &
      num_restart
    use prec_const
    use geometry_class
    use init_profile_class
    use fdistribu5d_class
    implicit none
      
    character(LEN=*)  , intent(in)    :: filename
    real(RKIND)       , intent(inout) :: initial_time
    type(geometry)    , intent(in)    :: geom
    type(init_profile), intent(in)    :: init_prof
    character(LEN=*)  , intent(in)    :: argname
    type(fdistribu5d) , intent(inout) :: f5d
      
    character(LEN=100) :: structname
      
#ifndef NOHDF5
    logical            :: fexist
    character(LEN=100) :: filename_tmp
    integer(HID_T)     :: file_id  
    integer            :: ierr
      
    filename_tmp = filename
    inquire(file=filename_tmp,exist=fexist)
    if (.not.(fexist)) then
      write(6,'(3A)') ' Error : the file ',filename_tmp, &
        ' does not exist'
    end if
      
    call HDF5_open(trim(filename_tmp),file_id)
    call HDF5_real_reading(file_id,deltat,"deltat")
    call HDF5_integer_reading(file_id,nbiter_prev,"nbiter_prev")
    call HDF5_integer_reading(file_id,nb_restart,"nb_restart")
    call HDF5_integer_reading(file_id,nb_diag,"nb_diag")
    call HDF5_integer_reading(file_id,Phi3D_nb_diag,"Phi3D_nb_diag")
    call HDF5_integer_reading(file_id,FMoments3D_nb_diag, &
      "FMoments3D_nb_diag")
    call HDF5_integer_reading(file_id,f5D_nb_diag,"f5D_nb_diag")
    call HDF5_real_reading(file_id,initial_time,"initial_time")
    write(structname,'(A,A)') argname,"%mu_value"
    call HDF5_real_reading(file_id,f5d%mu_value, &
      structname(1:len_trim(structname)))
    write(structname,'(A,A)') argname,"%values"
    call HDF5_array4D_reading(file_id, &
      f5d%values,structname(1:len_trim(structname)),ierr)
    call HDF5_close(file_id)
    if (ierr.ne.0) then
      print*,'pglobal_id = ',pglobal_id, &
        ' ==> error in reading of the 4D structure '
      call exit(-1)
    end if
#endif
  end subroutine read_f5d_rst
      
  !------------------------------------------
  ! Write the restart file
  !------------------------------------------
  subroutine write_f5d_rst(initial_time,argname,f5d)
    use globals, only : rst_saving, &
      deltat, nbiter_prev, nb_restart, nb_diag, &
      Phi3D_nb_diag, FMoments3D_nb_diag, f5D_nb_diag , &
      num_restart
    use prec_const
    use clock_module
    use fdistribu5d_class
      
#ifdef __INTEL_COMPILER
    implicit none
    include "mpif.h"
#else
    use MPI
    implicit none
#endif
      
!R3 #include "r3_info.h" !R3
    real(RKIND)      , intent(in) :: initial_time
    character(LEN=*) , intent(in) :: argname
    type(fdistribu5d), pointer    :: f5d
      
#ifndef NOHDF5
    integer        :: inumber_buf, jnumber_buf
    integer(HID_T) :: file_id   ! HDF5 file identifier
    integer        :: ierr, icount_open
#endif
    character(LEN=100) :: structname
    character(LEN=100) :: filename
      
    if (rst_saving) then                     
!R3      call r3_info_begin(r3_info_index_0, &   !R3
!R3        'CPU_write_f5d_res')                  !R3
      call clck_time(bclock_writerst)
    end if
      
#ifndef NOHDF5
    if (rst_saving) then
      write(filename,'(A,I1,A,A,A,I5.5,A)') &
        "gysela.rst.n",num_restart,".",trim(argname), &
        ".p",pglobal_id,".h5"
      
      !*** writting of the restart files ***
      ierr        = 1
      icount_open = 0
      do while ( (ierr.ne.0).and.(icount_open.lt.5) ) 
        call HDF5_create(trim(filename),file_id,ierr)
        icount_open = icount_open+1
      end do
      if (ierr.ne.0) then
        print*,'pglobal_id = ',pglobal_id, &
          ' ==> error for opening of ',filename
      end if
      call HDF5_real_saving(file_id,deltat,"deltat")
      call HDF5_integer_saving(file_id,nbiter_prev,"nbiter_prev")
      call HDF5_integer_saving(file_id,nb_restart,"nb_restart")
      call HDF5_integer_saving(file_id,nb_diag,"nb_diag")
      call HDF5_integer_saving(file_id, &
        Phi3D_nb_diag,"Phi3D_nb_diag")
      call HDF5_integer_saving(file_id,FMoments3D_nb_diag, &
        "FMoments3D_nb_diag")
      call HDF5_integer_saving(file_id,f5D_nb_diag,"f5D_nb_diag")
      call HDF5_real_saving(file_id,initial_time,"initial_time")
      write(structname,'(A,A)') argname,"%mu_value"
      call HDF5_real_saving(file_id, &
        f5d%mu_value,structname(1:len_trim(structname)))
      inumber_buf = f5d%iend_buf-f5d%istart_buf+1
      jnumber_buf = f5d%jend_buf-f5d%jstart_buf+1
      write(structname,'(A,A)') argname,"%values"
      call HDF5_array4D_saving(file_id, &
        f5d%values(f5d%istart_buf:,f5d%jstart_buf:,0:,0:), &
        inumber_buf,jnumber_buf,f5d%n3+1, &
        f5d%n4+1,structname(1:len_trim(structname)))
      call HDF5_close(file_id)
    end if
#endif
      
    if (rst_saving) then
!R3      call r3_info_end (r3_info_index_0)  !R3
      call clck_time(eclock_writerst)
      call clck_diff(bclock_writerst,eclock_writerst, &
        global_time_writerst)
    end if
  end subroutine write_f5d_rst
      
  !-----------------------------------------------
  ! Initialisation of the restart
  !-----------------------------------------------
  subroutine init_restart
    use globals          , only : restart, &
      pglobal_id, num_restart, nb_restart_diag
    implicit none
    include "mpiperso.h"
      
    integer :: ios, ierr
      
    if (pglobal_id.eq.0) then
      write(file_list,'(A)')"file_list.out"
      open(ulist,file = file_list, status="old", iostat=ios)
      close(ulist)
      
      if (ios.ne.0) then 
        restart         = .false.
        num_restart     = 1
        nb_restart_diag = 0
      else
        restart = .true.
        !*** nb_restart file reading ***  
        open(ulist,file = file_list, status="old", iostat=ios)
        if (ios.eq.0) then 
          read(ulist,*) nb_restart_diag
          nb_restart_diag = nb_restart_diag+1
        end if
        close(ulist)
        !*** num_restart file reading ***  
        write(file_list,'(A)')"num_restart.out"
        open(ulist,file = file_list, status="old", iostat=ios)
        if (ios.eq.0) then 
          read(ulist,'(I1)') num_restart
        end if
        close(ulist)
      end if
    end if
    call MPI_BCAST(restart,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(num_restart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      
#ifdef NOHDF5
    if (restart) then
      if (pglobal_id.eq.0) then
        print*,'RESTART NOT POSSIBLE : ', &
          'The restart files are not saved in the mode NOHDF5'
      end if
      stop
    end if
#endif
  end subroutine init_restart
      
  !-----------------------------------------------
  ! end of the restart
  !-----------------------------------------------
  subroutine end_restart()
    use globals          , only : pglobal_id, &
      num_restart, nb_restart
    implicit none
    include "mpiperso.h"
      
    if (pglobal_id.eq.0) then
      write(file_list,'(A)') "num_restart.out"
      open(ulist,file = file_list, status="replace", &
        form = 'FORMATTED')
      write(ulist,'(I1)') num_restart
      close(ulist)
      
      write(file_list,'(A)')"file_list.out"
      open(ulist,file = file_list, status="replace", &
        form = 'FORMATTED')
      write(ulist,*) nb_restart
      close(ulist)
    end if
  end subroutine end_restart
end module read_write_module
