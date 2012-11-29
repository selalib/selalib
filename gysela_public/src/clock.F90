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
      
!---------------------------------------------
! file : clock.f90
! date : 31/03/2006
!  variables for CPU time computation
!   -> used by G. Latu for local spline
!      development 
!---------------------------------------------
module clock_module
  use prec_const
  use globals, only : pglobal_id, Nbproc_tot
  implicit none   
      
  integer :: timeprec
  
#ifdef PWR
  PARAMETER (TIMEPREC=8)
  integer(8), external :: irtc
#else
  PARAMETER (TIMEPREC=8)
#endif
  !-> number of period per seconds
  integer(TIMEPREC) :: nb_periods_sec
  !-> maximum value of the clock counter   
  integer(TIMEPREC) :: nb_periods_max
      
  integer(TIMEPREC), public :: bclock_init
  integer(TIMEPREC), public :: eclock_init
  integer(TIMEPREC), public :: bclock_advec1Dvpar
  integer(TIMEPREC), public :: eclock_advec1Dvpar
  integer(TIMEPREC), public :: bclock_advec1Dphi
  integer(TIMEPREC), public :: eclock_advec1Dphi
  integer(TIMEPREC), public :: bclock_advec2D
  integer(TIMEPREC), public :: eclock_advec2D
  integer(TIMEPREC), public :: bclock_advec1Dr
  integer(TIMEPREC), public :: eclock_advec1Dr
  integer(TIMEPREC), public :: bclock_advec1Dtheta
  integer(TIMEPREC), public :: eclock_advec1Dtheta
  integer(TIMEPREC), public :: bclock_advance_diffr        !PRIV_DIFF
  integer(TIMEPREC), public :: eclock_advance_diffr        !PRIV_DIFF
  integer(TIMEPREC), public :: bclock_advance_Krook        !PRIV_DIFF
  integer(TIMEPREC), public :: eclock_advance_Krook        !PRIV_DIFF
  integer(TIMEPREC), public :: bclock_meanTandV
  integer(TIMEPREC), public :: eclock_meanTandV
  integer(TIMEPREC), public :: bclock_velocities
  integer(TIMEPREC), public :: eclock_velocities
  integer(TIMEPREC), public :: bclock_poisson
  integer(TIMEPREC), public :: eclock_poisson
  integer(TIMEPREC), public :: bclock_compJ0
  integer(TIMEPREC), public :: eclock_compJ0
  integer(TIMEPREC), public :: bclock_fncopy
  integer(TIMEPREC), public :: eclock_fncopy
  integer(TIMEPREC), public :: bclock_diagnostics
  integer(TIMEPREC), public :: eclock_diagnostics
  integer(TIMEPREC), public :: bclock_diag3D
  integer(TIMEPREC), public :: eclock_diag3D  
  integer(TIMEPREC), public :: bclock_HDF5_diag
  integer(TIMEPREC), public :: eclock_HDF5_diag
  integer(TIMEPREC), public :: bclock_HDF5
  integer(TIMEPREC), public :: eclock_HDF5
  integer(TIMEPREC), public :: bclock_physics
  integer(TIMEPREC), public :: eclock_physics
  integer(TIMEPREC), public :: bclock_CrossSection
  integer(TIMEPREC), public :: eclock_CrossSection
  integer(TIMEPREC), public :: bclock_writerst
  integer(TIMEPREC), public :: eclock_writerst
      
  real(RKIND) :: global_time_init
  real(RKIND) :: global_time_advec1Dvpar
  real(RKIND) :: global_time_advec1Dphi
  real(RKIND) :: global_time_advec2D
  real(RKIND) :: global_time_advec1Dr
  real(RKIND) :: global_time_advec1Dtheta
  real(RKIND) :: global_time_advance_diffr                 !PRIV_DIFF
  real(RKIND) :: global_time_advance_Krook                 !PRIV_DIFF
  real(RKIND) :: global_time_meanTandV
  real(RKIND) :: global_time_coll_refresh
  real(RKIND) :: global_time_poisson
  real(RKIND) :: global_time_compJ0
  real(RKIND) :: global_time_fncopy
  real(RKIND) :: global_time_diag3D
  real(RKIND) :: global_time_diagnostics
  real(RKIND) :: global_time_1D2D_diag
  real(RKIND) :: global_time_5D_saving
  real(RKIND) :: global_time_physics
  real(RKIND) :: global_time_CrossSection
  real(RKIND) :: global_time_writerst
  real(RKIND) :: global_time
      
  !******************************
  contains
  !******************************
      
#ifdef PWR
  !--------------------------------- 
  ! bip used to compute the CPU time 
  !---------------------------------
  subroutine clck_time(p_time)
    integer(TIMEPREC), intent(out) :: p_time
      
    p_time=irtc()
  end subroutine clck_time
#else
#ifdef CTIMER 
  !--------------------------------- 
  ! bip used to compute the CPU time 
  !---------------------------------
  subroutine clck_time(p_time)
    integer(TIMEPREC), intent(out) :: p_time
    call clockget(p_time)
  end subroutine clck_time
      
#else
  !--------------------------------- 
  ! bip used to compute the CPU time 
  !---------------------------------
  subroutine clck_time(p_time)
    integer(TIMEPREC), intent(out) :: p_time
    
    call system_clock(count=p_time)
  end subroutine clck_time
#endif
#endif
      
  !--------------------------------- 
  ! compute difference between 2 bips 
  !---------------------------------
  subroutine clck_diff(p_time1,p_time2,p_difftime)
    integer(TIMEPREC), intent(in)    :: p_time1, p_time2
    real(RKIND)      , intent(inout) :: p_difftime
      
    integer(TIMEPREC) :: nb_periods
    
    nb_periods = p_time2-p_time1
    if (p_time2 < p_time1) then
       nb_periods = nb_periods + nb_periods_max
    endif
    p_difftime = p_difftime+real(nb_periods)/real(nb_periods_sec)
  end subroutine clck_diff
      
  !--------------------------------- 
  ! compute difference between 2 bips 
  !---------------------------------
  subroutine clck_ldiff(p_time1,p_time2,p_difftime)
    integer(TIMEPREC), intent(in)    :: p_time1, p_time2
    real(RKIND)      , intent(inout) :: p_difftime
      
    integer(TIMEPREC) :: nb_periods
    
    nb_periods = p_time2-p_time1
    if (p_time2 < p_time1) then
       nb_periods = nb_periods + nb_periods_max
    endif
    p_difftime = real(nb_periods)/real(nb_periods_sec)
  end subroutine clck_ldiff
      
  !---------------------------------------
  ! initialization of the time variables
  !--------------------------------------- 
  subroutine initialization_time
#ifdef CTIMER
    call clockrate(nb_periods_sec)
    call clockmaxperiod(nb_periods_max)
#else
    call system_clock(count_rate=nb_periods_sec, &
      count_max=nb_periods_max)
#endif
    
    if (pglobal_id.eq.0) then
      global_time_init            = 0._RKIND
      global_time_advec1Dvpar     = 0._RKIND
      global_time_advance_diffr   = 0._RKIND               !PRIV_DIFF
      global_time_advance_Krook   = 0._RKIND               !PRIV_DIFF
      global_time_meanTandV       = 0._RKIND
      global_time_coll_refresh    = 0._RKIND
      global_time_advec1Dphi      = 0._RKIND
      global_time_advec2D         = 0._RKIND
      global_time_poisson         = 0._RKIND
      global_time_fncopy          = 0._RKIND
      global_time_compJ0          = 0._RKIND
      global_time_diag3D          = 0._RKIND
      global_time_physics         = 0._RKIND
      global_time_diagnostics     = 0._RKIND
      global_time_1D2D_diag       = 0._RKIND
      global_time_5D_saving       = 0._RKIND
      global_time_CrossSection    = 0._RKIND
      global_time_writerst        = 0._RKIND
      global_time                 = 0._RKIND
    end if
  end subroutine initialization_time
      
  !------------------------------
  ! CPU time writing
  !------------------------------ 
  subroutine write_time()
!CP!    use globals, only : rst_saving, PSM, uout_res, outputproc
    use globals, only : rst_saving, uout_res, outputproc
    call clck_diff(bclock_init,eclock_init, &
      global_time_init)
      
    if (pglobal_id.eq.outputproc) then
      write(uout_res,*) ' '
      write(uout_res,*) & 
        ' **************************************', &
        '************************** '
      write(uout_res,*) ' ****             CPU TIME COMPUTATION', &
        '                       **** '
      write(uout_res,*) ' **************************************', &
        '************************** '
      write(uout_res,*) &
        ' time for initialization              = ', &
        global_time_init
      write(uout_res,*) &
        ' time for advec1D in vpar             = ', &
        global_time_advec1Dvpar
      write(uout_res,*) &
        ' time for advec1D in phi              = ', &
        global_time_advec1Dphi
        write(uout_res,*) &
          ' time for advec2D in (r,theta)        = ', &
          global_time_advec2D
      write(uout_res,*) &
        ' time for Tmean and Vmean computation = ', &
        global_time_meanTandV
      write(uout_res,*) &
        ' time for solving QN equation         = ', &
        global_time_poisson
      write(uout_res,*) &
        ' time for copy of fn fnm1 fnp1        = ', &
        global_time_fncopy
      write(uout_res,*) &
        ' time for gyro3D of Phi               = ', &
        global_time_compJ0
      write(uout_res,*) &
        ' time for all diagnostics             = ', &
        global_time_diagnostics
      write(uout_res,*) &
        ' time for diagnostics 3D              = ', &
        global_time_diag3D
      write(uout_res,*) &
        ' time for physical value computation  = ', &
        global_time_physics
      write(uout_res,*) &
        ' time for f cross-sections            = ', &
        global_time_CrossSection
      write(uout_res,*) &
        ' time for 1D/2D diagnostic saving     = ', &
        global_time_1D2D_diag
      write(uout_res,*) &
        ' time for HDF5 f5D file saving        = ', &
        global_time_5D_saving
      
      global_time = global_time_init + &
        global_time_advec1Dvpar + global_time_advec1Dphi + &
        global_time_advance_diffr + &                      !PRIV_DIFF
        global_time_advance_Krook + &                      !PRIV_DIFF
        global_time_meanTandV + global_time_poisson + &
        global_time_fncopy + global_time_compJ0 + &
        global_time_diagnostics + global_time_coll_refresh
        global_time = global_time + global_time_advec2D
      if (rst_saving) then
        write(uout_res,*) &
          ' time for restart file writing        = ', &
          global_time_writerst
        global_time = global_time + global_time_writerst
      end if
      
      write(uout_res,*) ' **************************************', &
        '************************** '
      write(uout_res,*) &
        ' Total time (without file saving)     = ', &
        global_time - global_time_1D2D_diag - &
        global_time_5D_saving - global_time_writerst - &
        global_time_diag3D
      
      write(uout_res,*) ' **************************************', &
        '************************** '
      write(uout_res,*) &
        ' Total time                           = ', &
        global_time
      write(uout_res,*) ' **************************************', &
        '************************** '
    end if    
  end subroutine write_time
end module clock_module
