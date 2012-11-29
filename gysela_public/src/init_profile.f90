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
      
!------------------------------------------------
! file : init_profile.f90
! date : 21/06/2000
!  Initialization of the density and temperature
!  profiles
!------------------------------------------------
module init_profile_class
  use globals, only : profile_choice, magnetic_shear
  use prec_const
  use mem_alloc_module
  use geometry_class
  use spline1d_class
      
  implicit none
      
  !*** init_profile type definition ***
  type :: init_profile
    integer                                 :: BCr_left, BCr_right
    real(RKIND)                             :: tau0
    real(RKIND)                             :: rpeak
    real(RKIND)                             :: deltarn
    real(RKIND)                             :: deltarTe
    real(RKIND)                             :: deltarTi
    real(RKIND)                             :: invLn
    real(RKIND)                             :: invLTe
    real(RKIND)                             :: invLTi
    real(RKIND)   , dimension(:)  , pointer :: func_tanh
    real(RKIND)   , dimension(:)  , pointer :: n0, Te, Ti
    real(RKIND)   , dimension(:)  , pointer :: n0_coef1d
    real(RKIND)   , dimension(:)  , pointer :: Te_coef1d
    real(RKIND)   , dimension(:)  , pointer :: Ti_coef1d
    real(RKIND)   , dimension(:)  , pointer :: dn0dr
    real(RKIND)   , dimension(:)  , pointer :: dTidr
    real(RKIND)   , dimension(:)  , pointer :: dTedr
    real(RKIND)   , dimension(:)  , pointer :: dlogTidr
    real(RKIND)   , dimension(:)  , pointer :: d2n0dr2
    real(RKIND)   , dimension(:)  , pointer :: d2Tidr2
    real(RKIND)                             :: Te0, n0norm
    real(RKIND)   , dimension(:)  , pointer :: iota, shear, psi
    complex(CKIND), dimension(:)  , pointer :: cphi00_src
    real(RKIND)   , dimension(:)  , pointer :: fric_coeff
    real(RKIND)   , dimension(:)  , pointer :: diff_coeffDr
    real(RKIND)   , dimension(:)  , pointer :: diff_coeffrDr_hf
    real(RKIND)   , dimension(:)  , pointer :: diff_A
    real(RKIND)   , dimension(:)  , pointer :: diff_B
    real(RKIND)   , dimension(:)  , pointer :: diff_C
    real(RKIND)                             :: q_rp
    real(RKIND)                             :: psi_rp
    real(RKIND)                             :: Ti_rp
    real(RKIND)                             :: n0_rp
    real(RKIND)                             :: psimin, psimax
    integer                                 :: ipsimin, ipsimax
    type(nspline1d)                         :: nspline1d_r
  end type init_profile
      
  !******************************
  contains
  !******************************
    
  !---------------------------------------- 
  ! sech = cosh^-1 definition
  !---------------------------------------- 
  function sech(x)
    real(RKIND), intent(in) :: x   
    real(RKIND) :: sech
    
    sech = 1._RKIND/cosh(x)
  end function sech  
      
     
  !---------------------------------------------------- 
  ! Constructor : Memory allocation for density 
  ! and temperature profiles
  !----------------------------------------------------
  subroutine new_init_profile(ithis,geom)
    use globals, only : memory_test, BC_Hermite, &
      tau0, rpeak, deltarn, &
      deltarTe, deltarTi, kappan, kappaTe, kappaTi, R0
    type(init_profile), intent(out) :: ithis
    type(geometry)    , intent(in)  :: geom
    
    integer     :: Nr_loc, Ntheta_loc, Nvpar_loc
    real(RKIND) :: dr_loc, a_norm
    
    Nr_loc     = geom%Nr
    Ntheta_loc = geom%Ntheta
    Nvpar_loc  = geom%Nvpar
    dr_loc     = geom%dr
    if (.not.memory_test) then
      a_norm     = geom%rg(Nr_loc)-geom%rg(0)
      
      !*** variable initialization ***
      ithis%tau0     = tau0
      ithis%rpeak    = geom%rg(0)+rpeak*a_norm
      ithis%deltarn  = deltarn*a_norm
      ithis%deltarTe = deltarTe*a_norm
      ithis%deltarTi = deltarTi*a_norm
      ! -> kappa = R/L 
      ithis%invLn    = kappan/R0        
      ithis%invLTe   = kappaTe/R0
      ithis%invLTi   = kappaTi/R0
    end if
      
    !*** memory allocation ***
    call glob_allocate(ithis%func_tanh,0,Nr_loc,'func_tanh')
    call glob_allocate(ithis%Te,0,Nr_loc,'Te')
    call glob_allocate(ithis%Ti,0,Nr_loc,'Ti')       
    call glob_allocate(ithis%n0,0,Nr_loc,'n0')
    call glob_allocate(ithis%dTedr,0,Nr_loc,'dTedr')
    call glob_allocate(ithis%dTidr,0,Nr_loc,'dTidr')
    call glob_allocate(ithis%dn0dr,0,Nr_loc,'dn0dr')
    call glob_allocate(ithis%dlogTidr,0,Nr_loc,'dlogTidr')        
    call glob_allocate(ithis%d2n0dr2,0,Nr_loc,'d2n0dr2')
    call glob_allocate(ithis%d2Tidr2,0,Nr_loc,'d2Tidr2')
    call glob_allocate(ithis%Te_coef1d,-1,Nr_loc+1,'Te_coef1d')
    call glob_allocate(ithis%Ti_coef1d,-1,Nr_loc+1,'Ti_coef1d')
    call glob_allocate(ithis%n0_coef1d,-1,Nr_loc+1,'n0_coef1d')  
    call glob_allocate(ithis%iota,0,Nr_loc,'iota')
    call glob_allocate(ithis%shear,0,Nr_loc,'shear')
    call glob_allocate(ithis%psi,0,Nr_loc,'psi')
    call glob_allocate(ithis%cphi00_src,0,Nr_loc,'cphi00_src')
    call glob_allocate(ithis%fric_coeff,0,Nr_loc,'fric_coeff')
    call glob_allocate(ithis%diff_coeffDr,0,Nr_loc,'diff_coeffDr')
    call glob_allocate(ithis%diff_coeffrDr_hf,0,Nr_loc+1,&
      'diff_coeffrDr_hf')
    call glob_allocate(ithis%diff_A,0,Nr_loc,'diff_A')
    call glob_allocate(ithis%diff_B,0,Nr_loc,'diff_B')
    call glob_allocate(ithis%diff_C,0,Nr_loc,'diff_C')
      
    !*** initialization of the boundary conditions ***
    ithis%BCr_left  = BC_Hermite
    ithis%BCr_right = BC_Hermite
      
    !*** initialization of the spline in r direction ***
    call new_spline1d_natural(ithis%nspline1d_r,Nr_loc,dr_loc)
  end subroutine new_init_profile
  
      
  !---------------------------------------------------- 
  ! destructor
  !---------------------------------------------------- 
  subroutine del_init_profile(ithis)
    type(init_profile), intent(out) :: ithis
      
    call glob_deallocate(ithis%func_tanh)
    call glob_deallocate(ithis%Te)         
    call glob_deallocate(ithis%dTedr) 
    call glob_deallocate(ithis%Ti)    
    call glob_deallocate(ithis%dTidr)      
    call glob_deallocate(ithis%dlogTidr)
    call glob_deallocate(ithis%n0)        
    call glob_deallocate(ithis%dn0dr)    
    call glob_deallocate(ithis%d2n0dr2)
    call glob_deallocate(ithis%d2Tidr2)
    call glob_deallocate(ithis%Te_coef1d) 
    call glob_deallocate(ithis%Ti_coef1d) 
    call glob_deallocate(ithis%n0_coef1d)
    call glob_deallocate(ithis%iota)
    call glob_deallocate(ithis%shear)
    call glob_deallocate(ithis%psi)
    call glob_deallocate(ithis%cphi00_src)
    call glob_deallocate(ithis%fric_coeff)
    call glob_deallocate(ithis%diff_coeffDr)
    call glob_deallocate(ithis%diff_coeffrDr_hf)
    call glob_deallocate(ithis%diff_A)
    call glob_deallocate(ithis%diff_B)
    call glob_deallocate(ithis%diff_C)
    call del_spline1d_natural(ithis%nspline1d_r)
  end subroutine del_init_profile
      
  !**************************************************************
  !  CHOICE OF INITIAL PROFILES AS FUNCTION OF HYPERBOLIC
  !   COSINUS :
  !  -> EX: 1/T(r) dT(r)/dr = -(1/LT)*func(cosh^-2)
  !**************************************************************
  !---------------------------------------------------- 
  !  n0 density profile initialization
  !  this profile is given by the density gradient 
  !  profile, defined by 3 parameters : 
  !    - kappan, rpeak and deltar
  !  and
  !    1/n0(r) dn0(r)/dr = -(1/Ln)*cosh^-2(r-rpeak/deltar)
  !  with (1/Ln) = kappan/R
  !---------------------------------------------------- 
  subroutine init_n0_1(ithis,geom)
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
      
    real(RKIND) :: r, dr_loc, rth, tmp
    integer     :: ir, Nr_loc
    real(RKIND) :: n0norm_tmp
    
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
    
    !*** compute n0 solution of :                            ***
    !***  2/(n0(r)+n0(r-1))*(n0(r)-n0(r-1))/dr               ***
    !***                  = -(1/LTi)*cosh^-2(r-rpeak/deltar) ***
    ithis%n0(0) = 10._RKIND**19
    do ir=1,Nr_loc
      r            = geom%rg(ir-1)
      rth          = r + dr_loc*0.5_RKIND
      tmp          = -ithis%invLn * &
        sech((rth-ithis%rpeak)/ithis%deltarn)**2
      tmp          = 0.5_RKIND*dr_loc*tmp
      ithis%n0(ir) = (1._RKIND+tmp)/(1._RKIND-tmp)*ithis%n0(ir-1)
    enddo
    
    !*** normalisation of the density at int(n0(r)rdr)/int(rdr) ***
    ! -> computation of int(n0(r)rdr)
    n0norm_tmp = 0._RKIND
    do ir = 1,Nr_loc-1
      n0norm_tmp = n0norm_tmp + ithis%n0(ir)*geom%rg(ir)
    enddo
    n0norm_tmp = n0norm_tmp + 0.5_RKIND* &
      (ithis%n0(0)*geom%rg(0) + ithis%n0(Nr_loc)*geom%rg(Nr_loc))
    ! -> division by int(rdr)
    n0norm_tmp = n0norm_tmp*2._RKIND*dr_loc/ & 
      (geom%rg(Nr_loc)**2-geom%rg(0)**2)
      
    ithis%n0norm       = n0norm_tmp
    ithis%n0(0:Nr_loc) = ithis%n0(0:Nr_loc)/ithis%n0norm
  end subroutine init_n0_1
      
  !-------------------------------------------------------- 
  ! derivative of ionic density profile
  !  dn0(r)/dr = -n0(r)*(1/Ln)*cosh^-2(r-rpeak/deltar)
  !--------------------------------------------------------     
  subroutine init_dn0dr_1(ithis,geom)
    type(init_profile), intent(inout) :: ithis    
    type(geometry)    , intent(in)    :: geom
      
    real(RKIND) :: r     
    integer     :: ir, Nr_loc
        
    !*** local variable initialization ***
    Nr_loc = geom%Nr
        
    !*** temperature gradient profile for electron ***
    do ir = 0,Nr_loc
      r               = geom%rg(ir)
      ithis%dn0dr(ir) = -ithis%n0(ir)* &
        ithis%invLn*sech((r-ithis%rpeak)/ithis%deltarn)**2
    end do
  end subroutine init_dn0dr_1
      
  !---------------------------------------------------- 
  ! electronic temperature profile initialization
  !  this profile is given by the temperature gradient 
  !  profile, defined by 3 parameters : 
  !    - kappaTe, rpeak and deltar
  !  and
  !    1/Te(r) dTe(r)/dr = -invLTe*cosh^-2(r-rpeak/deltar)
  !  with (1/LTe) = kappaTe/R
  !----------------------------------------------------     
  subroutine init_Te_1(ithis,geom)
    type(init_profile), intent(inout) :: ithis    
    type(geometry)    , intent(in)    :: geom
      
    real(RKIND) :: r, dr_loc, rth, tmp, w0, w1
    integer     :: ir, Nr_loc
        
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
            
    !*** compute Te solution of :                            ***
    !***  2/(Te(r)+Te(r-1))*(Te(r)-Te(r-1))/dr               ***
    !***                  = -(1/LTe)*cosh^-2(r-rpeak/deltar) ***
    ithis%Te(0) = 1.3_RKIND	!keV
    do ir = 1,Nr_loc
      r            = geom%rg(ir-1)
      rth          = r+ dr_loc*0.5_RKIND
      tmp          = -0.5_RKIND*dr_loc*ithis%invLTe &
        *sech((rth-ithis%rpeak)/ithis%deltarTe)**2
      ithis%Te(ir) = (1._RKIND+tmp)/(1._RKIND-tmp)*ithis%Te(ir-1)
    enddo
    
    !*** normalisation of the temperature to 1 at r=rpeak ***
    ir                 = int((ithis%rpeak-geom%rg(0))/dr_loc)
    w1                 = (ithis%rpeak-geom%rg(ir))/dr_loc
    w0                 = 1._RKIND-w1
    ithis%Te0          = w0*ithis%Te(ir)+w1*ithis%Te(ir+1)
    ithis%Te(0:Nr_loc) = ithis%Te(0:Nr_loc)/ithis%Te0
  end subroutine init_Te_1
      
  !-------------------------------------------------------- 
  ! derivative of electronic temperature profile 
  !  dTe(r)/dr = -Te(r)*(1/LTe)*cosh^-2(r-rpeak/deltar)
  !--------------------------------------------------------     
  subroutine init_dTedr_1(ithis,geom)
    type(init_profile), intent(inout) :: ithis    
    type(geometry)    , intent(in)    :: geom
      
    real(RKIND) :: r     
    integer     :: ir, Nr_loc
        
    !*** local variable initialization ***
    Nr_loc = geom%Nr
        
    !*** temperature gradient profile for electron ***
    do ir = 0,Nr_loc
      r               = geom%rg(ir)
      ithis%dTedr(ir) = -ithis%Te(ir)* &
        ithis%invLTe*sech((r-ithis%rpeak)/ithis%deltarTe)**2
    end do
  end subroutine init_dTedr_1
       
       
  !---------------------------------------------------- 
  ! ionic temperature profile initialization
  !  this profile is given by the temperature gradient 
  !  profile, defined by 3 parameters : 
  !    - kappaTi, rpeak and deltar
  !  and
  !    1/Ti(r) dTi(r)/dr = -kappaTi*cosh^-2(r-rpeak/deltar)
  !  with (1/LTi) = kappaTi/R
  !----------------------------------------------------   
  subroutine init_Ti_1(ithis,geom)
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
    
    real(RKIND) :: r, dr_loc, rth, tmp, w0, w1, Ti0
    integer     :: ir, Nr_loc
    
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
        
    !*** compute Ti solution of :                            ***
    !***  2/(Ti(r)+Ti(r-1))*(Ti(r)-Ti(r-1))/dr               ***
    !***                  = -(1/LTi)*cosh^-2(r-rpeak/deltar) ***
    ithis%Ti(0) = ithis%Te(0)
    do ir = 1,Nr_loc
      r            = geom%rg(ir-1)
      rth          = r+ dr_loc*0.5_RKIND
      tmp          = -0.5_RKIND*dr_loc*ithis%invLTi &
        *sech((rth-ithis%rpeak)/ithis%deltarTi)**2
      ithis%Ti(ir) = (1._RKIND+tmp)/(1._RKIND-tmp)*ithis%Ti(ir-1)
    enddo
    
    !*** normalisation of the temperature to ***
    !*** Te0 (Te0 = Ti0/tau0)                ***    
    ir                 = int((ithis%rpeak-geom%rg(0))/dr_loc)
    w1                 = (ithis%rpeak-geom%rg(ir))/dr_loc
    w0                 = 1._RKIND-w1
    Ti0                = w0*ithis%Ti(ir)+w1*ithis%Ti(ir+1)
    ithis%Ti(0:Nr_loc) = (ithis%Ti(0:Nr_loc)/Ti0)*ithis%tau0
  end subroutine init_Ti_1  
  
  
  !----------------------------------------------------------- 
  ! ionic temperature gradient profile initialization
  !  this profile is by 3 parameters : 
  !    - kappaTi, rpeak and deltar
  !  and
  !    1/Ti(r) dTi(r)/dr = -(1/LTi)*cosh^-2(r-rpeak/deltar)
  !  with (1/LTi) = kappaTi/R
  !-----------------------------------------------------------   
  subroutine init_dlogTidr_1(ithis,geom)
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
    
    real(RKIND) :: r, dr_loc 
    integer     :: ir, Nr_loc
    
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
    
    !*** temperature gradient profile for ion ***
    do ir = 0,Nr_loc
      r                  = geom%rg(ir)
      ithis%dlogTidr(ir) = -ithis%invLTi* &
        sech((r-ithis%rpeak)/ithis%deltarTi)**2
    end do
  end subroutine init_dlogTidr_1
      
  !**************************************************************
  !  CHOICE OF INITIAL PROFILES SUCH THAT DN/DT 
  !   AND D(NT)/DT ARE CONSTANT
  !**************************************************************
  !-------------------------------------------------------- 
  ! computation of func(tanh) = [tanh(0.5*(rg-rmin-rbuff)) 
  !                - tanh(0.5*(rg-rmax+rbuff))]/max(hr) 
  !-------------------------------------------------------- 
  subroutine init_func_tanh_1(ithis,geom)
    type(init_profile), intent(inout) :: ithis    
    type(geometry)    , intent(in)    :: geom
      
    integer     :: ir, Nr_loc
    real(RKIND) :: ri, dr_loc
    real(RKIND) :: rbuff, rmin, rmax, max_func
    real(RKIND) :: func_tanh_tmp
        
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
    
    !*** computation of func(tanh) and its maximum ***
    rbuff    = 0.25_RKIND*Nr_loc*dr_loc
    rmin     = geom%rg(0)
    rmax     = geom%rg(Nr_loc)
    max_func = 0._RKIND
    do ir = 0,Nr_loc
      ri = geom%rg(ir)
      func_tanh_tmp       = tanh(0.5_RKIND*(ri-rmin-rbuff)) - &
        tanh(0.5_RKIND*(ri-rmax+rbuff))
      ithis%func_tanh(ir) = func_tanh_tmp
      max_func            = max(func_tanh_tmp,max_func)
    end do
    ithis%func_tanh = ithis%func_tanh/max_func
  end subroutine init_func_tanh_1
      
  !---------------------------------------------------- 
  !  n0 density profile initialization such that
  !    dn0(r)/dr = -(1/Ln)*func(tanh(r))
  !  so :
  !   n0(ir) = n0(ir-1) - delta_r*(1/Ln)*func(tanh(r))
  !---------------------------------------------------- 
  subroutine init_n0_2(ithis,geom)
    use globals, only : a, Rarray1_Nr
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
      
    integer     :: ir, Nr_loc
    real(RKIND) :: dr_loc
    real(RKIND) :: tmp
      
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
      
    !*** computation of n0 ***
    ithis%n0(0) = 1._RKIND
    do ir = 1,Nr_loc
      tmp          = -ithis%invLn*ithis%func_tanh(ir)
      ithis%n0(ir) = ithis%n0(ir-1)+tmp*dr_loc
    enddo
    ithis%n0 = ithis%n0 - ithis%n0(Nr_loc) + 1._RKIND
  end subroutine init_n0_2
      
  !-------------------------------------------------------- 
  ! derivative of ionic density profile equal to a constant 
  !    dn0(r)/dr = -(1/Ln)*func(tanh(r))
  !--------------------------------------------------------     
  subroutine init_dn0dr_2(ithis,geom)
    type(init_profile), intent(inout) :: ithis    
    type(geometry)    , intent(in)    :: geom
      
    integer     :: ir, Nr_loc
    real(RKIND) :: dr_loc
        
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
        
    !*** computation of dn0/dr ***
    do ir = 0,Nr_loc
      ithis%dn0dr(ir) = -ithis%invLn*ithis%func_tanh(ir)
    end do
  end subroutine init_dn0dr_2
 
 
  !---------------------------------------------------- 
  !  electronic temperature profile initialization
  !   such that the derivative of the pressure profile
  !   is constant :
  !  d(n0*Te)/dt = -(1/LP)*func(tanh(r))
  !   with 1/LP = (1/Ln)+(1/LTe)
  !   so :
  !  n0(ir)*Te(ir) = n0(ir-1)*Te(ir-1) 
  !                  - delta_r*(1/LP)*func(tanh(ri)) 
  !---------------------------------------------------- 
  subroutine init_Te_2(ithis,geom)
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
      
    integer     :: ir, Nr_loc
    real(RKIND) :: dr_loc
    real(RKIND) :: tmp
      
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
      
    ithis%Te(0) = 1.3_RKIND	!keV
    !*** compute of Te ***
    do ir = 1,Nr_loc
      tmp         = -(ithis%invLn+ithis%invLTe)* &
        ithis%func_tanh(ir)
      ithis%Te(ir) = (ithis%n0(ir-1)/ithis%n0(ir))* &
        ithis%Te(ir-1) + tmp*dr_loc/ithis%n0(ir)
    enddo
    ithis%Te = ithis%Te - ithis%Te(Nr_loc) + 1._RKIND
  end subroutine init_Te_2
      
  !---------------------------------------------------- 
  ! ionic temperature gradient profile initialization
  !   such that the derivative of the pressure profile
  !   is constant :
  !  d(n0*Te)/dt = -(1/LP)*func(tanh(r))
  !   with (1/LP) = (1/Ln)+(1/LTe)
  !   so :
  !  dTe/dr = -(1/LP)*func(tanh(r))/n0-(Te/n0)*dn0/dr
  !---------------------------------------------------- 
  subroutine init_dTedr_2(ithis,geom)
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
    
    integer     :: ir, Nr_loc
    real(RKIND) :: dr_loc
      
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
    
    !*** temperature gradient profile for ion ***
    do ir = 0,Nr_loc
      ithis%dTedr(ir) = -(ithis%invLn+ithis%invLTe) * &
        ithis%func_tanh(ir)/ithis%n0(ir) - &
        (ithis%Te(ir)/ithis%n0(ir))*ithis%dn0dr(ir)
    end do
  end subroutine init_dTedr_2
      
  !---------------------------------------------------- 
  !  ionic temperature profile initialization
  !   such that the derivative of the pressure profile
  !   is constant :
  !  d(n0*Ti)/dt = -(1/LP)*func(tanh(r))
  !   with 1/LP = (1/Ln)+(1/LTi)
  !   so :
  !  n0(ir)*Ti(ir) = n0(ir-1)*Ti(ir-1) 
  !                  - delta_r*(1/LP)*func(tanh(ri))
  !---------------------------------------------------- 
  subroutine init_Ti_2(ithis,geom)
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
      
    integer     :: ir, Nr_loc
    real(RKIND) :: dr_loc
    real(RKIND) :: tmp
      
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
      
    !*** compute of Ti ***
    ithis%Ti(0) = ithis%Te(0)
    do ir = 1,Nr_loc
      tmp         = -(ithis%invLn+ithis%invLTi)* &
        ithis%func_tanh(ir)
      ithis%Ti(ir) = (ithis%n0(ir-1)/ithis%n0(ir))* &
        ithis%Ti(ir-1) + tmp*dr_loc/ithis%n0(ir)
    enddo
    ithis%Ti = ithis%Ti - ithis%Ti(Nr_loc) + 1._RKIND
  end subroutine init_Ti_2
      
  !---------------------------------------------------- 
  ! ionic temperature gradient profile initialization
  !   such that the derivative of the pressure profile
  !   is constant :
  !  d(n0*Ti)/dt = -(1/LP)*func(tanh(r))
  !   with 1/LP = (1/Ln)+(1/LTi)
  !   so :
  ! (1/Ti)dTi/dr = -(1/LP)*func(tanh(r))/(n0*Ti)
  !              -(1/n0)*dn0/dr
  !---------------------------------------------------- 
  subroutine init_dlogTidr_2(ithis,geom)
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
    
    integer     :: ir, Nr_loc
    real(RKIND) :: dr_loc
      
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
    
    !*** temperature gradient profile for ion ***
    do ir = 0,Nr_loc
      ithis%dlogTidr(ir) = -(ithis%invLn+ithis%invLTi) * &
        ithis%func_tanh(ir)/(ithis%n0(ir)*ithis%Ti(ir)) - &
        ithis%dn0dr(ir)/ithis%n0(ir)
    end do
  end subroutine init_dlogTidr_2
      
  !---------------------------------------------------- 
  !  electronic temperature profile initialization
  !   such that the derivative of the temperature 
  !   profile is constant :
  !  d(Te)/dt = -(1/LTe)*func(tanh(r))
  !   so :
  !  Te(ir) = Te(ir-1) - delta_r*(1/LTe)*func(tanh(ri)) 
  !---------------------------------------------------- 
  subroutine init_Te_3(ithis,geom)
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
      
    integer     :: ir, Nr_loc
    real(RKIND) :: dr_loc
    real(RKIND) :: tmp
      
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
      
    ithis%Te(0) = 1.3_RKIND	!keV
    !*** compute of Te ***
    do ir = 1,Nr_loc
      tmp         = - ithis%invLTe*ithis%func_tanh(ir)
      ithis%Te(ir) = ithis%Te(ir-1) + tmp*dr_loc
    enddo
    ithis%Te = ithis%Te - ithis%Te(Nr_loc) + 1._RKIND
  end subroutine init_Te_3
      
  !---------------------------------------------------- 
  ! ionic temperature gradient profile initialization
  !   such that the derivative of the temperature
  !    profile is constant :
  !  dTe/dt = -(1/LTe)*func(tanh(r))
  !---------------------------------------------------- 
  subroutine init_dTedr_3(ithis,geom)
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
    
    integer     :: ir, Nr_loc
    real(RKIND) :: dr_loc
      
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
    
    !*** temperature gradient profile for ion ***
    do ir = 0,Nr_loc
      ithis%dTedr(ir) = -ithis%invLTe * &
        ithis%func_tanh(ir)
    end do
  end subroutine init_dTedr_3
      
  !---------------------------------------------------- 
  !  ionic temperature profile initialization
  !   such that the derivative of the temperature 
  !   profile is constant :
  !  dTi/dt = -(1/Ti)*func(tanh(r))
  !---------------------------------------------------- 
  subroutine init_Ti_3(ithis,geom)
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
      
    integer     :: ir, Nr_loc
    real(RKIND) :: dr_loc
    real(RKIND) :: tmp
      
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
      
    !*** compute of Ti ***
    ithis%Ti(0) = ithis%Te(0)
    do ir = 1,Nr_loc
      tmp         = -ithis%invLTi*ithis%func_tanh(ir)
      ithis%Ti(ir) = ithis%Ti(ir-1) + tmp*dr_loc
    enddo
    ithis%Ti = ithis%Ti - ithis%Ti(Nr_loc) + 1._RKIND
  end subroutine init_Ti_3
      
  !---------------------------------------------------- 
  ! ionic temperature gradient profile initialization
  !   such that the derivative of the temperature 
  !   profile is constant :
  !     dTi/dt = -(1/LTi)*func(tanh(r))
  !  so  (1/Ti)dTi/dt = -(1/LTi)*func(tanh(r))/Ti
  !---------------------------------------------------- 
  subroutine init_dlogTidr_3(ithis,geom)
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
    
    integer     :: ir, Nr_loc
    real(RKIND) :: dr_loc
      
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
    
    !*** temperature gradient profile for ion ***
    do ir = 0,Nr_loc
      ithis%dlogTidr(ir) = -ithis%invLTi * &
        ithis%func_tanh(ir)/ithis%Ti(ir)
    end do
  end subroutine init_dlogTidr_3
      
  !**************************************************************
  !  CHOICE OF INITIAL PROFILES AS FUNCTION OF HYPERBOLIC
  !   COSINUS FOR BENCHMARK WITH ORB5
  !  -> EX: 1/T(r) dT(r)/dr = -(1/LT)*
  !       (-1+cosh^-2((r-rmin)/deltar)+cosh^-2((r-rmax)/deltar))
  !**************************************************************
  !-------------------------------------------------------- 
  ! computation of func(tanh) = [tanh(0.5*(rg-rmin-rbuff1)) 
  !                - tanh(0.5*(rg-rmax+rbuff2))]/max(hr) 
  !-------------------------------------------------------- 
  subroutine init_func_tanh_4(ithis,geom)
    type(init_profile), intent(inout) :: ithis    
    type(geometry)    , intent(in)    :: geom
      
    integer     :: ir, Nr_loc
    real(RKIND) :: ri, Lr_loc
    real(RKIND) :: rmin, rmax, rp
    real(RKIND) :: rbuff1, rbuff2
    real(RKIND) :: func_tanh_tmp, max_func
        
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    Lr_loc = geom%Lr
    
    !*** computation of func(tanh) and its maximum ***
    rbuff1   = 0.2_RKIND*Lr_loc 
    rbuff2   = 0.2_RKIND*Lr_loc 
    rmin     = geom%rg(0)
    rmax     = geom%rg(Nr_loc)
    rp       = geom%rg(Nr_loc/2)
    max_func = 0._RKIND
    do ir = 0,Nr_loc
      ri = geom%rg(ir)
      func_tanh_tmp       = tanh(0.5_RKIND*(ri-rmin-rbuff1)) - &
        tanh(0.5_RKIND*(ri-rmax+rbuff2))
      ithis%func_tanh(ir) = func_tanh_tmp
      max_func            = max(func_tanh_tmp,max_func)
    end do
    ithis%func_tanh = ithis%func_tanh/(max_func)
  end subroutine init_func_tanh_4
      
  !-------------------------------------------------------------- 
  !  n0 density profile initialization
  !  this profile is given by the density gradient 
  !  profile, defined by 3 parameters : 
  !    - kappan, rpeak and deltar
  !  and
  !    1/n0(r) dn0(r)/dr = -(1/Ln)*func_sech(r)*func_tanh(r)
  !  with (1/Ln) = kappan/R
  !  and where 
  !   func_sech(r) = -1+cosh^-2((r-rmin)/deltar) + 
  !      cosh^-2((r-rmax)/deltar)
  !-------------------------------------------------------------- 
      
  subroutine init_n0_4(ithis,geom)
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
      
    real(RKIND) :: r, dr_loc, rth, tmp
    integer     :: ir, Nr_loc
    real(RKIND) :: n0norm_tmp
    real(RKIND) :: Lr_loc, rmin, rmax, func_sech
    
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
    Lr_loc = geom%Lr
    rmin   = geom%rg(0)     !+0.1_RKIND*Lr_loc
    rmax   = geom%rg(Nr_loc)!-0.1_RKIND*Lr_loc
    
    !*** compute n0 solution of :                ***
    !***  2/(n0(r)+n0(r-1))*(n0(r)-n0(r-1))/dr   ***
    !***                  = -(1/Ln)*func_sech(r) ***
    ithis%n0(0) = 10._RKIND**19
    do ir =1,Nr_loc
      r            = geom%rg(ir-1)
      rth          = r + dr_loc*0.5_RKIND
      func_sech    = -1._RKIND + &
        sech((rth-rmin)/ithis%deltarn)**2 + &
        sech((rth-rmax)/ithis%deltarn)**2
      tmp          = 0.5_RKIND*dr_loc * &
        ithis%invLn*func_sech*ithis%func_tanh(ir)
      ithis%n0(ir) = (1._RKIND+tmp)/(1._RKIND-tmp)*ithis%n0(ir-1)
    enddo
    
    !*** normalisation of the density at int(n0(r)rdr)/int(rdr) ***
    ! -> computation of int(n0(r)rdr)
    n0norm_tmp = 0._RKIND
    do ir = 1,Nr_loc-1
      n0norm_tmp = n0norm_tmp + ithis%n0(ir)*geom%rg(ir)
    enddo
    n0norm_tmp = n0norm_tmp + 0.5_RKIND* &
      (ithis%n0(0)*geom%rg(0) + ithis%n0(Nr_loc)*geom%rg(Nr_loc))
    ! -> division by int(rdr)
    n0norm_tmp = n0norm_tmp*2._RKIND*dr_loc/ & 
      (geom%rg(Nr_loc)**2-geom%rg(0)**2)
      
    ithis%n0norm       = n0norm_tmp
    ithis%n0(0:Nr_loc) = ithis%n0(0:Nr_loc)/ithis%n0norm
  end subroutine init_n0_4
      
  !---------------------------------------------------------------
  ! derivative of ionic density profile
  !  dn0(r)/dr = -n0(r)*(1/Ln)*
  !       (-1+cosh^-2((r-rmin)/deltar)+cosh^-2((r-rmax)/deltar))
  !---------------------------------------------------------------
  subroutine init_dn0dr_4(ithis,geom)
    type(init_profile), intent(inout) :: ithis    
    type(geometry)    , intent(in)    :: geom
      
    integer     :: ir, Nr_loc
    real(RKIND) :: Lr_loc, r, rmin, rmax, func_sech
        
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    Lr_loc = geom%Lr
    rmin   = geom%rg(0)     !+0.1_RKIND*Lr_loc
    rmax   = geom%rg(Nr_loc)!-0.1_RKIND*Lr_loc
        
    !*** temperature gradient profile for electron ***
    do ir = 0,Nr_loc
      r               = geom%rg(ir)
      func_sech       = -1._RKIND + &
        sech((r-rmin)/ithis%deltarn)**2 + &
        sech((r-rmax)/ithis%deltarn)**2      
      ithis%dn0dr(ir) = ithis%n0(ir)* &
        ithis%invLn*func_sech*ithis%func_tanh(ir)
    end do
  end subroutine init_dn0dr_4
      
  !------------------------------------------------------------
  ! electronic temperature profile initialization
  !  this profile is given by the temperature gradient 
  !  profile, defined by 3 parameters : 
  !    - kappaTe, rpeak and deltar
  !  and
  !    1/Te(r) dTe(r)/dr = -invLTe*func_sech(r)
  !  with (1/LTe) = kappaTe/R
  !  and where
  !   func_sech(r) = (-1+cosh^-2((r-rmin)/deltar) + 
  !    cosh^-2((r-rmax)/deltar))
  !------------------------------------------------------------
  subroutine init_Te_4(ithis,geom)
    type(init_profile), intent(inout) :: ithis    
    type(geometry)    , intent(in)    :: geom
      
    integer     :: ir, Nr_loc
    real(RKIND) :: r, dr_loc, rth, tmp, w0, w1
    real(RKIND) :: Lr_loc, rmin, rmax, func_sech
        
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
    Lr_loc = geom%Lr
    rmin   = geom%rg(0)     !+0.1_RKIND*Lr_loc
    rmax   = geom%rg(Nr_loc)!-0.1_RKIND*Lr_loc
            
    !*** compute Te solution of :                 ***
    !***  2/(Te(r)+Te(r-1))*(Te(r)-Te(r-1))/dr    ***
    !***                  = -(1/LTe)*func_sech(r) ***
    ithis%Te(0) = 2._RKIND	!keV
    do ir = 1,Nr_loc
      r            = geom%rg(ir-1)
      rth          = r+ dr_loc*0.5_RKIND
      func_sech    = -1._RKIND + &
        sech((rth-rmin)/ithis%deltarTe)**2 + &
        sech((rth-rmax)/ithis%deltarTe)**2
      tmp          = 0.5_RKIND*dr_loc * &
        ithis%invLTe*func_sech*ithis%func_tanh(ir)
      ithis%Te(ir) = (1._RKIND+tmp)/(1._RKIND-tmp)*ithis%Te(ir-1)
    enddo
    
    !*** normalisation of the temperature to 1 at r=rpeak ***
    ir                 = int((ithis%rpeak-geom%rg(0))/dr_loc)
    w1                 = (ithis%rpeak-geom%rg(ir))/dr_loc
    w0                 = 1._RKIND-w1
    ithis%Te0          = w0*ithis%Te(ir)+w1*ithis%Te(ir+1)
    ithis%Te(0:Nr_loc) = ithis%Te(0:Nr_loc)/ithis%Te0
  end subroutine init_Te_4
      
  !--------------------------------------------------------------
  ! derivative of electronic temperature profile 
  !  dTe(r)/dr = -Te(r)*(1/LTe)*
  !       (-1+cosh^-2((r-rmin)/deltar)+cosh^-2((r-rmax)/deltar))
  !--------------------------------------------------------------
  subroutine init_dTedr_4(ithis,geom)
    type(init_profile), intent(inout) :: ithis    
    type(geometry)    , intent(in)    :: geom
      
    real(RKIND) :: Lr_loc, r, rmin, rmax, func_sech
    integer     :: ir, Nr_loc
        
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    Lr_loc = geom%Lr
    rmin   = geom%rg(0)     !+0.1_RKIND*Lr_loc
    rmax   = geom%rg(Nr_loc)!-0.1_RKIND*Lr_loc
        
    !*** temperature gradient profile for electron ***
    do ir = 0,Nr_loc
      r               = geom%rg(ir)
      func_sech       = -1._RKIND + &
        sech((r-rmin)/ithis%deltarTe)**2 + &
        sech((r-rmax)/ithis%deltarTe)**2
      ithis%dTedr(ir) = ithis%Te(ir) * &
        ithis%invLTe*func_sech*ithis%func_tanh(ir)
    end do
  end subroutine init_dTedr_4
       
       
  !-------------------------------------------------------------
  ! ionic temperature profile initialization
  !  this profile is given by the temperature gradient 
  !  profile, defined by 3 parameters : 
  !    - kappaTi, rpeak and deltar
  !  and
  !    1/Ti(r) dTi(r)/dr = -invLTi*func_sech(r)*func_tanh(r)
  !  with (1/LTi) = kappaTi/R
  !  and where
  !   func_sech(r) = (-1+cosh^-2((r-rmin)/deltar) + 
  !      cosh^-2((r-rmax)/deltar))
  !-------------------------------------------------------------
  subroutine init_Ti_4(ithis,geom)
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
    
    integer     :: ir, Nr_loc
    real(RKIND) :: Lr_loc, r, dr_loc
    real(RKIND) :: rth, tmp, w0, w1, Ti0
    real(RKIND) :: rmin, rmax, func_sech
    
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
    Lr_loc = geom%Lr
    rmin   = geom%rg(0)     !+0.1_RKIND*Lr_loc
    rmax   = geom%rg(Nr_loc)!-0.1_RKIND*Lr_loc
      
    !*** compute Ti solution of :                            ***
    !***  2/(Ti(r)+Ti(r-1))*(Ti(r)-Ti(r-1))/dr               ***
    !***                  = -(1/LTi)*cosh^-2(r-rpeak/deltar) ***
    ithis%Ti(0) = ithis%Te(0)
    do ir = 1,Nr_loc
      r            = geom%rg(ir-1)
      rth          = r+ dr_loc*0.5_RKIND
      func_sech    = -1._RKIND + &
        sech((rth-rmin)/(ithis%deltarTi))**2 + &
        sech((rth-rmax)/ithis%deltarTi)**2      
      tmp          = 0.5_RKIND*dr_loc * &
        ithis%invLTi*func_sech*ithis%func_tanh(ir)
      ithis%Ti(ir) = (1._RKIND+tmp)/(1._RKIND-tmp)*ithis%Ti(ir-1)
    enddo
    
    !*** normalisation of the temperature to ***
    !***   Te0 (Te0 = Ti0/tau0)              ***    
    ir                 = int((ithis%rpeak-geom%rg(0))/dr_loc)
    w1                 = (ithis%rpeak-geom%rg(ir))/dr_loc
    w0                 = 1._RKIND-w1
    Ti0                = w0*ithis%Ti(ir)+w1*ithis%Ti(ir+1)
    ithis%Ti(0:Nr_loc) = (ithis%Ti(0:Nr_loc)/Ti0)*ithis%tau0        
  end subroutine init_Ti_4
  
  
  !-------------------------------------------------------------
  ! ionic temperature gradient profile initialization
  !  this profile is by 3 parameters : 
  !    - kappaTi, rpeak and deltar
  !  and
  !    1/Ti(r) dTi(r)/dr = -(1/LTi)*func_sech(r)*func_tanh(r)
  !  with (1/LTi) = kappaTi/R
  !  and where
  !   func_sech(r) = (-1+cosh^-2((r-rmin)/deltar) + 
  !     cosh^-2((r-rmax)/deltar))
  !-------------------------------------------------------------
  subroutine init_dlogTidr_4(ithis,geom)
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
    
    integer     :: ir, Nr_loc
    real(RKIND) :: Lr_loc, r, rmin, rmax, func_sech
      
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    Lr_loc = geom%Lr
    rmin   = geom%rg(0)     !+0.1_RKIND*Lr_loc
    rmax   = geom%rg(Nr_loc)!-0.1_RKIND*Lr_loc
      
    do ir = 0,Nr_loc
      r                  = geom%rg(ir)
      func_sech          = -1._RKIND + &
        sech((r-rmin)/(ithis%deltarTi))**2 + &
        sech((r-rmax)/ithis%deltarTi)**2
      ithis%dlogTidr(ir) = ithis%invLTi * &
        func_sech*ithis%func_tanh(ir) 
    end do
  end subroutine init_dlogTidr_4
      
  !**************************************************************
  !  For TEST: initialisation with simple analytical profiles
  !   for the collisions
  !  -> n0 is linear
  !  -> Te = Ti = cte
  !**************************************************************
  !-------------------------------------------------------- 
  ! computation of func(tanh) = [tanh(0.5*(rg-rmin-rbuff1)) 
  !                - tanh(0.5*(rg-rmax+rbuff2))]/max(hr) 
  !-------------------------------------------------------- 
  subroutine init_func_tanh_5(ithis,geom)
    type(init_profile), intent(inout) :: ithis    
    type(geometry)    , intent(in)    :: geom
      
    integer     :: ir, Nr_loc
    real(RKIND) :: ri, Lr_loc
    real(RKIND) :: rmin, rmax, rp
    real(RKIND) :: rbuff1, rbuff2
    real(RKIND) :: func_tanh_tmp, max_func
        
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    Lr_loc = geom%Lr
    
    !*** computation of func(tanh) and its maximum ***
    rbuff1   = 0.2_RKIND*Lr_loc 
    rbuff2   = 0.3_RKIND*Lr_loc 
    rmin     = geom%rg(0)
    rmax     = geom%rg(Nr_loc)
    rp       = geom%rg(Nr_loc/2)
    max_func = 0._RKIND
    do ir = 0,Nr_loc
      ri = geom%rg(ir)
      func_tanh_tmp       = tanh(0.5_RKIND*(ri-rmin-rbuff1)) - &
        tanh(0.5_RKIND*(ri-rmax+rbuff2))
      ithis%func_tanh(ir) = func_tanh_tmp
      max_func            = max(func_tanh_tmp,max_func)
    end do
    ithis%func_tanh = ithis%func_tanh/(max_func)
  end subroutine init_func_tanh_5
      
  !-------------------------------------------------------------- 
  !  n0 density profile initialization
  !  this profile is linear, analytic, independant 
  !  of the 3 parameters : 
  !    - kappan, rpeak and deltar
  !  and
  !    n0(r) = -2*Dn*r/(rMax-rMin) + 1/(rMax - rMin) * 
  !          ( np*(rMax - rMin) + Dn*(rMax - rMin) )
  !-------------------------------------------------------------- 
  subroutine init_n0_5(ithis,geom)
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
      
    real(RKIND) :: r, dr_loc, tmp
    integer     :: ir, Nr_loc
    real(RKIND) :: n0norm_tmp
    real(RKIND) :: Lr_loc, rmin, rmax
    real(RKIND) :: np, Dn, AAA, BBB
    
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
    Lr_loc = geom%Lr
    rmin   = geom%rg(0)     !+0.1_RKIND*Lr_loc
    rmax   = geom%rg(Nr_loc)!-0.1_RKIND*Lr_loc
    np     = 1._RKIND
    Dn     = 0.05_RKIND
    
    !*** compute n0 which is also "nGi_eq"          ***
    !*** CAREFUL: modify also in "fequil4d.f90",    ***
    !*** subroutine "compute_nGieq" for consistency *** 
    do ir =0,Nr_loc
      r            = geom%rg(ir)
      AAA          = -2*Dn/(rmax-rmin)
      BBB          = ( Dn*(rmax+rmin) + np*(rmax-rmin) ) / &
        (rmax-rmin)
      ithis%n0(ir) = AAA*r + BBB
    enddo
      
    !*** normalisation of the density at int(n0(r)rdr)/int(rdr) ***
    ! -> computation of int(n0(r)rdr)
    n0norm_tmp = 0._RKIND
    do ir = 1,Nr_loc-1
      n0norm_tmp = n0norm_tmp + ithis%n0(ir)*geom%rg(ir)
    enddo
    n0norm_tmp = n0norm_tmp + 0.5_RKIND* &
      (ithis%n0(0)*geom%rg(0) + ithis%n0(Nr_loc)*geom%rg(Nr_loc))
    ! -> division by int(rdr)
    n0norm_tmp = n0norm_tmp*2._RKIND*dr_loc/ & 
      (geom%rg(Nr_loc)**2-geom%rg(0)**2)
      
    ithis%n0norm       = n0norm_tmp
    ithis%n0(0:Nr_loc) = ithis%n0(0:Nr_loc)/ithis%n0norm
  end subroutine init_n0_5
      
  !---------------------------------------------------------------
  ! derivative of ionic density profile
  !  dn0(r)/dr = -2*Dn/(rMax-rMin)
  !---------------------------------------------------------------
  subroutine init_dn0dr_5(ithis,geom)
    type(init_profile), intent(inout) :: ithis    
    type(geometry)    , intent(in)    :: geom
      
    integer     :: ir, Nr_loc
    real(RKIND) :: Lr_loc, r, rmin, rmax, func_sech
    real(RKIND) :: np, Dn
        
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    Lr_loc = geom%Lr
    rmin   = geom%rg(0)     !+0.1_RKIND*Lr_loc
    rmax   = geom%rg(Nr_loc)!-0.1_RKIND*Lr_loc
    np     = 1._RKIND
    Dn     = 0.05_RKIND
        
    !*** temperature gradient profile for electron ***
    do ir = 0,Nr_loc
      ithis%dn0dr(ir) = -2*Dn/(rmax-rmin)
    end do
  end subroutine init_dn0dr_5
      
  !------------------------------------------------------------
  ! electronic temperature profile initialization
  !  this profile is taken constant for test purposes
  ! The 3 parameters : 
  !    - kappaTe, rpeak and deltar
  !  do not intervene
  !------------------------------------------------------------
  subroutine init_Te_5(ithis,geom)
    type(init_profile), intent(inout) :: ithis    
    type(geometry)    , intent(in)    :: geom
      
    integer     :: ir, Nr_loc
    real(RKIND) :: r, dr_loc, rth, tmp, w0, w1
    real(RKIND) :: Lr_loc, rmin, rmax, func_sech
        
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
    Lr_loc = geom%Lr
    rmin   = geom%rg(0)     !+0.1_RKIND*Lr_loc
    rmax   = geom%rg(Nr_loc)!-0.1_RKIND*Lr_loc
            
    !*** compute Te = 1 everywhere ***
    do ir = 0,Nr_loc
      ithis%Te(ir) = 1._RKIND
    enddo
  end subroutine init_Te_5
      
  !-----------------------------------------------------
  ! derivative of electronic temperature profile 
  !  dTe(r)/dr = 0
  !-----------------------------------------------------
  subroutine init_dTedr_5(ithis,geom)
    type(init_profile), intent(inout) :: ithis    
    type(geometry)    , intent(in)    :: geom
      
    real(RKIND) :: Lr_loc, r, rmin, rmax, func_sech
    integer     :: ir, Nr_loc
        
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    Lr_loc = geom%Lr
    rmin   = geom%rg(0)     !+0.1_RKIND*Lr_loc
    rmax   = geom%rg(Nr_loc)!-0.1_RKIND*Lr_loc
        
    !*** temperature gradient profile for electron ***
    do ir = 0,Nr_loc
      ithis%dTedr(ir) = 10._RKIND**(-30)
    end do
  end subroutine init_dTedr_5
       
       
  !------------------------------------------------------
  ! ionic temperature profile initialization
  !  this profile is taken constant for test purposes
  ! The 3 parameters : 
  !    - kappaTe, rpeak and deltar
  !  do not intervene
  !------------------------------------------------------
  subroutine init_Ti_5(ithis,geom)
    type(init_profile), intent(inout) :: ithis    
    type(geometry)    , intent(in)    :: geom
      
    integer     :: ir, Nr_loc
    real(RKIND) :: r, dr_loc, rth, tmp, w0, w1
    real(RKIND) :: Lr_loc, rmin, rmax, func_sech
        
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
    Lr_loc = geom%Lr
    rmin   = geom%rg(0)     !+0.1_RKIND*Lr_loc
    rmax   = geom%rg(Nr_loc)!-0.1_RKIND*Lr_loc
            
    !*** compute Ti = 1 everywhere ***
    do ir = 0,Nr_loc
      ithis%Ti(ir) = 1._RKIND
    enddo
  end subroutine init_Ti_5
      
  !----------------------------------------------------
  ! ionic temperature gradient profile initialization
  !  this profile is null
  !----------------------------------------------------
  subroutine init_dlogTidr_5(ithis,geom)
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
    
    integer     :: ir, Nr_loc
    real(RKIND) :: Lr_loc, r, rmin, rmax, func_sech
      
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    Lr_loc = geom%Lr
    rmin   = geom%rg(0)     !+0.1_RKIND*Lr_loc
    rmax   = geom%rg(Nr_loc)!-0.1_RKIND*Lr_loc
      
    do ir = 0,Nr_loc
      ithis%dlogTidr(ir) = 10._RKIND**(-30)
    end do
  end subroutine init_dlogTidr_5
      
  !-----------------------------------------------------------
  ! Initialisation of the hyperbolic tangent 
  !   used for profile definition
  !-----------------------------------------------------------
  subroutine init_func_tanh(ithis,geom)
    use globals, only : profile_choice
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
      
    select case (profile_choice)
    case (1,2,3)
      call init_func_tanh_1(ithis,geom)
    case (4)
      call init_func_tanh_4(ithis,geom)
    case (5)
      call init_func_tanh_5(ithis,geom)
    case default
      if (pglobal_id.eq.0) then
        print*, 'profile_choice = ', profile_choice, &
          ' is not a proper choice'
        stop
      end if
    end select
  end subroutine init_func_tanh
      
  !-----------------------------------------------------------
  ! Initialisation of ion density profile
  !  + computation of the first and second derivatives
  !-----------------------------------------------------------
  subroutine init_n0(ithis,geom)
    use globals            , only : read_n0, profile_choice
    use utils_module       , only : deriv1
    use read_profile_module, only : read_ascii_profile
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
      
    if (.not.read_n0) then
      select case (profile_choice)
      case (1)
        call init_n0_1(ithis,geom)
        call init_dn0dr_1(ithis,geom)
      case (2,3)
        call init_n0_2(ithis,geom)
        call init_dn0dr_2(ithis,geom)
      case (4)
        call init_n0_4(ithis,geom)
        call init_dn0dr_4(ithis,geom)
      case (5)
        call init_n0_5(ithis,geom)
        call init_dn0dr_5(ithis,geom)
      end select
    else
      call read_ascii_profile('n0.dat',geom,ithis%n0)
      call deriv1(ithis%n0,ithis%dn0dr,geom%Nr,geom%dr,0)
    end if
    !*** computation of the second derivative of n0 ***
    call deriv1(ithis%dn0dr,ithis%d2n0dr2,geom%Nr,geom%dr,0)
  end subroutine init_n0
      
  !-----------------------------------------------------------
  ! Initialisation of electron temperature profile
  !  + computation of the first and second derivatives
  !-----------------------------------------------------------
  subroutine init_Te(ithis,geom)
    use globals            , only : read_Te, profile_choice
    use utils_module       , only : deriv1
    use read_profile_module, only : read_ascii_profile
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
      
    if (.not.read_Te) then
      select case (profile_choice)
      case (1)
        call init_Te_1(ithis,geom)
        call init_dTedr_1(ithis,geom)
      case (2)
        call init_Te_2(ithis,geom)
        call init_dTedr_2(ithis,geom)
      case (3)
        call init_Te_3(ithis,geom)
        call init_dTedr_3(ithis,geom)
      case (4)
        call init_Te_4(ithis,geom)
        call init_dTedr_4(ithis,geom)
      case (5)
        call init_Te_5(ithis,geom)
        call init_dTedr_5(ithis,geom)
      end select
    else
      call read_ascii_profile('Te.dat',geom,ithis%Te)
      call deriv1(ithis%Te,ithis%dTedr,geom%Nr,geom%dr,0)
    end if
  end subroutine init_Te
      
  !-----------------------------------------------------------
  ! Initialisation of ion temperature profile
  !  + computation of the first and second derivatives
  !-----------------------------------------------------------
  subroutine init_Ti(ithis,geom)
    use globals            , only : read_Ti, profile_choice
    use utils_module       , only : deriv1
    use read_profile_module, only : read_ascii_profile
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
      
    integer :: ir_loc
      
    if (.not.read_Ti) then
      select case (profile_choice)
      case (1)
        call init_Ti_1(ithis,geom)
        call init_dlogTidr_1(ithis,geom)
      case (2)
        call init_Ti_2(ithis,geom)
        call init_dlogTidr_2(ithis,geom)
      case (3)
        call init_Ti_3(ithis,geom)
        call init_dlogTidr_3(ithis,geom)
      case (4)
        call init_Ti_4(ithis,geom)
        call init_dlogTidr_4(ithis,geom)
      case (5)
        call init_Ti_5(ithis,geom)
        call init_dlogTidr_5(ithis,geom)
      end select
      !-> computation of dTi/dr
      do ir_loc = 0,geom%Nr
        ithis%dTidr(ir_loc) = ithis%dlogTidr(ir_loc) * &
          ithis%Ti(ir_loc)
      end do
    else
      call read_ascii_profile('Ti.dat',geom,ithis%Ti)
      call deriv1(ithis%Ti,ithis%dTidr,geom%Nr,geom%dr,0)
      !-> computation of dlogTi/dr
      do ir_loc = 0,geom%Nr
        ithis%dlogTidr(ir_loc) = ithis%dTidr(ir_loc) / &
          ithis%Ti(ir_loc)
      end do
    end if
    !*** computation of the second derivative of Ti ***
    call deriv1(ithis%dTidr,ithis%d2Tidr2,geom%Nr,geom%dr,0)
  end subroutine init_Ti
      
  !---------------------------------------------- 
  ! compute the 1D-spline coefficient of Te
  ! ( used for cubic spline interpolation )
  ! ( for example in energy computation )
  !---------------------------------------------- 
  subroutine compute_Tecoef1d(ithis,geom)
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
    
    real(RKIND)                 :: dr_loc 
    integer                     :: Nr_loc
    real(RKIND), dimension(0:1) :: deriv_Te
      
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
      
    !*** compute the cubic spline coefficient ***
    deriv_Te(0) = ithis%dTedr(Nr_loc)
    deriv_Te(1) = ithis%dTedr(0)
    call natural_spline_coef(ithis%nspline1d_r,ithis%Te,&
      ithis%BCr_left,ithis%BCr_right,deriv_Te)
    ithis%Te_coef1d = ithis%nspline1d_r%scoef
  end subroutine compute_Tecoef1d
      
  !---------------------------------------------- 
  ! compute the 1D-spline coefficient of Te
  ! ( used for cubic spline interpolation )
  ! ( for example in energy computation )
  !---------------------------------------------- 
  subroutine compute_Ticoef1d(ithis,geom)
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
    
    real(RKIND)                 :: dr_loc 
    integer                     :: Nr_loc
    real(RKIND), dimension(0:1) :: deriv_Ti
      
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
      
    !*** compute the cubic spline coefficient ***
    deriv_Ti(0) = ithis%dlogTidr(Nr_loc)*ithis%Ti(Nr_loc)
    deriv_Ti(1) = ithis%dlogTidr(0)*ithis%Ti(0)
    call natural_spline_coef(ithis%nspline1d_r,ithis%Ti,&
      ithis%BCr_left,ithis%BCr_right,deriv_Ti)
    ithis%Ti_coef1d = ithis%nspline1d_r%scoef
  end subroutine compute_Ticoef1d
      
  !---------------------------------------------- 
  ! compute the 1D-spline coefficient of Te
  ! ( used for cubic spline interpolation )
  ! ( for example in energy computation )
  !---------------------------------------------- 
  subroutine compute_n0coef1d(ithis,geom)
    type(init_profile), intent(inout) :: ithis  
    type(geometry)    , intent(in)    :: geom
    
    real(RKIND)                 :: dr_loc 
    integer                     :: Nr_loc
    real(RKIND), dimension(0:1) :: deriv_n0
    
    !*** local variable initialization ***
    Nr_loc = geom%Nr
    dr_loc = geom%dr
      
    !*** compute the cubic spline coefficient ***
    deriv_n0(0) = ithis%dn0dr(Nr_loc)
    deriv_n0(1) = ithis%dn0dr(0) 
    call natural_spline_coef(ithis%nspline1d_r,ithis%n0,&
      ithis%BCr_left,ithis%BCr_right,deriv_n0)
    ithis%n0_coef1d = ithis%nspline1d_r%scoef
  end subroutine compute_n0coef1d
      
  !-------------------------------------------------------- 
  ! Compute iota the inverse of the safety factor
  !  iota(r) = Kshear * 1/q(r) with
  !     q(r) = q0 + deltaq*(r/a)^alphaq
  !  and Kshear = 1 if the magnetic shear is taken into 
  !   account and 0 otherwise
  !-------------------------------------------------------- 
  subroutine compute_iota(ithis,geom)
    use globals, only : magnetic_shear, read_q, &
      a, q0, deltaq, alphaq, reversed_shear, qmin, qa, rho_qmin
    use read_profile_module, only : read_ascii_profile
    type(init_profile), intent(out) :: ithis
    type(geometry)    , intent(in)  :: geom 
      
    integer     :: ir, Nr_loc    
    real(RKIND) :: ri, qr
    real(RKIND) :: rhomin2, rhomin4, rhomin6, &
      var, rho, rho2, C2, C3
      
    !*** local variable initialization ***
    Nr_loc = geom%Nr
      
    !*** safety factor ***
    if (magnetic_shear) then
      if (.not.read_q) then
        if (reversed_shear) then
          rhomin2 = rho_qmin*rho_qmin
          rhomin4 = rhomin2*rhomin2
          rhomin6 = rhomin4*rhomin2
          var     = rhomin4*(1-rhomin2)*(1-rhomin2)
          C2      = (rhomin6*(qa-qmin) + &
            (1-rhomin2)**3*(q0-qmin))/var
          C3      = (rhomin4*(qa-qmin) - &
            (1-rhomin2)**2*(q0-qmin))/var
          do ir = 0,Nr_loc
            ri   = geom%rg(ir)
            rho  = ri/a
            rho2 = rho*rho
            qr   = qmin + C2*(rho2-rhomin2)**2 + &
              C3*(rho2-rhomin2)**3
            ithis%iota(ir) = 1._RKIND/qr
          end do
        else
          do ir = 0,Nr_loc
            ri = geom%rg(ir)
            qr = q0      
            if (ri.ne.0._RKIND) &
              qr = qr + deltaq*exp(alphaq*log(ri/a))
            ithis%iota(ir) = 1._RKIND/qr
          end do
        end if
      else
        !-> reading of the q profile
        call read_ascii_profile('safety_factor.dat',geom,ithis%iota)
        !-> computation of 1/q
        do ir = 0,Nr_loc
          ithis%iota(ir) = 1._RKIND/ithis%iota(ir)
        end do
      end if
    else
      !-> iota = 0._RKIND 
      do ir = 0,Nr_loc
        ithis%iota(ir) = 0._RKIND
      end do
    end if
  end subroutine compute_iota
      
  !------------------------------------------------------------ 
  ! Compute the magnetic shear
  !   s(r) = r/q dq/dr
  ! as q(r) = q0 + deltaq*(r/a)^alphaq
  !   s(r) = deltaq*alphaq/q*(r/a)^alphaq
  !
  ! Rk : if .not.magnetic_shear, then iota = 0. so s(r) = 0
  !------------------------------------------------------------ 
  subroutine compute_magnetic_shear(ithis,geom)
    use globals, only : a, q0, deltaq, alphaq, &
      reversed_shear, qmin, qa, rho_qmin
    type(init_profile), intent(out) :: ithis
    type(geometry)    , intent(in)  :: geom 
      
    integer     :: ir, Nr_loc    
    real(RKIND) :: ri
    real(RKIND) :: rhomin2, rhomin4, rhomin6, &
                   var, rho, rho2, C2, C3, dqdr
      
    !*** local variable initialization ***
    Nr_loc = geom%Nr
      
    !*** computation of s profile ***
    if (reversed_shear) then
      rhomin2 = rho_qmin*rho_qmin
      rhomin4 = rhomin2*rhomin2
      rhomin6 = rhomin4*rhomin2
      var     = rhomin4*(1-rhomin2)*(1-rhomin2)
      C2      = (rhomin6*(qa-qmin)+(1-rhomin2)**3*(q0-qmin))/var
      C3      = (rhomin4*(qa-qmin)-(1-rhomin2)**2*(q0-qmin))/var
      do ir = 0,Nr_loc
        ri   = geom%rg(ir)
        rho  = ri/a
        rho2 = rho*rho
        dqdr = 4._RKIND*C2*(rho2-rhomin2)*rho/a + &
               6._RKIND*C3*(rho2-rhomin2)**2*rho/a
        ithis%shear(ir) = ri*ithis%iota(ir) * dqdr
      end do
    else
      do ir = 0,Nr_loc
        ri              = geom%rg(ir)
        ithis%shear(ir) = deltaq*alphaq*ithis%iota(ir) * &
          exp(alphaq*log(ri/a))
      end do
    end if
  end subroutine compute_magnetic_shear
      
  !---------------------------------------------------------------
  ! Computation of psi(r) = -\int_0^r r/q(r) dr using the formula:
  !    psi(r) = -\int_0^r0 r/q(r) dr - \int_r0^r r/q(r) dr
  ! Rk : if .not.magnetic_shear, then iota = 0. so psi(r) = 0
  !---------------------------------------------------------------
  subroutine compute_psi(ithis,geom)
    use mem_alloc_module
    type(init_profile), intent(out) :: ithis
    type(geometry)    , intent(in)  :: geom
    
    integer                           :: ir
    real(RKIND)                       :: ri, iota_ri
    real(RKIND)                       :: psi_r
    real(RKIND), dimension(0:geom%Nr) :: H_r
      
    ! -> H(r) = -r/q(r)
    H_r = -geom%rg*ithis%iota
    ! -> psi(r0) = -\int_0^r0 H(r) dr = [H(r0)+H(0)]*(r0-0)/2
    ithis%psi(0) = H_r(0)*geom%rg(0)*0.5_RKIND 
      
    ! -> computation of psi(r) = -\int_r0^r H(r) dr
    !  with psi(ri) = psi(ri-1) + \int_ri-1^ri H(r) dr
    do ir = 1,geom%Nr
      ithis%psi(ir) = ithis%psi(ir-1) + &
        0.5_RKIND*(H_r(ir)+H_r(ir-1))*geom%dr
    end do
      
    !-> computation of psimin = min(psi) and psimax = max(psi)
    !->  and the corresponding indices ipsimin and ipsimax
    if ( (maxval(ithis%psi)+minval(ithis%psi)) .ne. 0._RKIND ) then 
      ithis%psimin = 1.e+15_RKIND
      ithis%psimax = -1.e+15_RKIND
      do ir = 0,geom%Nr
        ithis%psimin = min(ithis%psimin,ithis%psi(ir))
        ithis%psimax = max(ithis%psimax,ithis%psi(ir))
        if (ithis%psi(ir).eq.ithis%psimin) ithis%ipsimin = ir
        if (ithis%psi(ir).eq.ithis%psimax) ithis%ipsimax = ir
      end do
    else
      ithis%ipsimin = 0
      ithis%ipsimax = 0
    end if
  end subroutine compute_psi
      
  !-----------------------------------------------------------
  ! Definition of phi00_src(r) used in poisson.f90:
  !    phi00 = phi00 + phi00_src for the (m,n)=(0,0) mode
  ! This prescribed shape leads to an additional ExB shear
  !-----------------------------------------------------------
  subroutine compute_phi00_src(ithis,geom)
    use globals, only : a, ExB_shear, max_shear, &
      rho_shear, deltar_shear
    type(init_profile), intent(out) :: ithis
    type(geometry)    , intent(in)  :: geom
      
    integer     :: ir
    real(RKIND) :: rho, norm, lncosh0, lncosh, phi00_src
      
    if (ExB_shear) then
      norm    = max_shear*(a*deltar_shear)**2*float(geom%Ntheta)
      lncosh0 = log(cosh((geom%rg(0)/a-rho_shear)/deltar_shear))
      do ir = 0,geom%Nr
        rho       = geom%rg(ir)/a
        lncosh    = log(cosh((rho-rho_shear)/deltar_shear))
        phi00_src = (lncosh - lncosh0)*norm
        ithis%cphi00_src(ir) = phi00_src
      end do
    else
      do ir = 0,geom%Nr
        ithis%cphi00_src(ir) = (0._RKIND,0._RKIND)
      end do
    end if
  end subroutine compute_phi00_src
      
  !---------------------------------------------------------- 
  ! Computation of a profile at the point r = rpeak
  !---------------------------------------------------------
  subroutine compute_value_atrpeak(ithis,geom,fct,fct_atrpeak)
    use globals     , only : Nr
    use utils_module, only : locate
    type(init_profile)             , intent(inout) :: ithis
    type(geometry)                 , intent(in)    :: geom
    real(RKIND)   , dimension(0:Nr), intent(in)    :: fct
    real(RKIND)                    , intent(out)   :: fct_atrpeak
      
    integer                         :: i, ipos
    real(RKIND)   , dimension(-1:2) :: sbase_fct
      
    !*** computation of the cubic spline coefficients ***
    call natural_spline_coef(ithis%nspline1d_r,fct, &
      ithis%BCr_left,ithis%BCr_right)
    !*** interpolation at r=rpeak ***
    call locate(ithis%rpeak,geom%rg,geom%Nr,geom%dr, &
      " compute_value_atrpeak "//char(0),ipos)
    call spline_basis(geom%rg(ipos),ithis%rpeak, &
      geom%rg(ipos+1),geom%dr,sbase_fct)
    fct_atrpeak = 0._RKIND
    do i = -1,2
      fct_atrpeak = fct_atrpeak + &
        ithis%nspline1d_r%scoef(ipos+i)*sbase_fct(i)
    enddo
  end subroutine compute_value_atrpeak
      
  !----------------------------------------- 
  ! Initialization of all the profiles  
  !-----------------------------------------
  subroutine init_all_profile(ithis,geom)
    use globals, only : pglobal_id
    use read_profile_module
    type(init_profile), intent(out) :: ithis
    type(geometry)    , intent(in)  :: geom 
    
    real(RKIND) :: iota_rp, dt
      
    !*** profile object construction ***
    call new_init_profile(ithis,geom)
      
    if (.not.memory_test) then
      !*** initialisation of the hyperbolic tangent used ***
      !***   for profile definition                      ***
      call init_func_tanh(ithis,geom)
      
      !*** initialisation of the ion density profile ***
      call init_n0(ithis,geom)
      
      !*** initialisation of the electron temperature profile ***
      call init_Te(ithis,geom)
      
      !*** initialisation of the ion temperature profile ***
      call init_Ti(ithis,geom)
      
      !*** cubic spline coefficient for Te ***
      call compute_Tecoef1d(ithis,geom) 
      
      !*** cubic spline coefficient for Ti ***
      call compute_Ticoef1d(ithis,geom) 
      
      !*** cubic spline coefficient for n0 ***
      call compute_n0coef1d(ithis,geom) 
      
      !*** iota profile initialization ***
      call compute_iota(ithis,geom)
      call compute_magnetic_shear(ithis,geom)
      call compute_psi(ithis,geom)
      
      !*** ExB shear initialization ***
      call compute_phi00_src(ithis,geom)
      
      !*** initialization of the values at r=rpeak ***
      ! -> psi(rp) with rp = r(peak)
      call compute_value_atrpeak(ithis,geom,ithis%psi,ithis%psi_rp)
      ! -> q(rp) with rp = r(peak)
      call compute_value_atrpeak(ithis,geom,ithis%iota,iota_rp)
      if (iota_rp.ne.0._RKIND) then
        ithis%q_rp = 1._RKIND/iota_rp
      else 
        ithis%q_rp = 0._RKIND
      end if
      call compute_value_atrpeak(ithis,geom,ithis%Ti,ithis%Ti_rp)
      call compute_value_atrpeak(ithis,geom,ithis%n0,ithis%n0_rp)
    end if
  end subroutine init_all_profile
end module init_profile_class
