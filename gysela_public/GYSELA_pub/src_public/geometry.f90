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
! file : geometry.f90
! date : 16/09/2005
! 5d geometry definition (r,theta,phi,vparallel,mu)
!--------------------------------------------------
module geometry_class
  use prec_const
  use mem_alloc_module
  implicit none
  
  !*** geometry type definition ***
  type :: geometry
    !->  origin coordinates
    real(RKIND) :: rmin, thetamin, phimin, vparmin, mumin
    !-> grid step size
    real(RKIND) :: dr, dtheta, dphi, dvpar
    !-> grid length            
    real(RKIND) :: Lr, Ltheta, Lphi, Lvpar, Lmu
    !-> number of grid points
    integer     :: Nr, Ntheta, Nphi, Nvpar, Nmu
    ! coordinates of the points in the 5 directions
    real(RKIND), dimension(:), pointer :: rg
    real(RKIND), dimension(:), pointer :: thetag
    real(RKIND), dimension(:), pointer :: phig
    real(RKIND), dimension(:), pointer :: vparg
    real(RKIND), dimension(:), pointer :: mug
    ! coefficients for the integral in all directions 
    !   with trapeze method
    real(RKIND), dimension(:), pointer :: coeff_intdr
    real(RKIND), dimension(:), pointer :: coeff_intdtheta
    real(RKIND), dimension(:), pointer :: coeff_intdphi
    real(RKIND), dimension(:), pointer :: coeff_intdvpar
    real(RKIND), dimension(:), pointer :: coeff_intdmu
    ! array for cos(theta) and sin(theta)
    real(RKIND), dimension(:), pointer :: cos_theta, sin_theta
    ! index in vparg of vparg = 0.
    integer :: ivpar0
  end type geometry
      
  integer, parameter :: integration_scheme = 2
    
  !******************************
  contains
  !******************************
      
  !----------------------------------------------------
  ! initialisation of the space coordinates (r,theta,phi)
  !----------------------------------------------------  
  subroutine new_space_coord(gthis,rmin,thetamin,phimin, &
    Lr,Ltheta,Lphi,Nr,Ntheta,Nphi)
    use globals, only : memory_test
    type(geometry) , intent(out)   :: gthis
    real(RKIND)    , intent(in)    :: rmin, thetamin, phimin 
    real(RKIND)    , intent(in)    :: Lr, Ltheta, Lphi       
    integer        , intent(in)    :: Nr, Ntheta, Nphi       
    
    integer :: i1, i2, i3
    
    !*** variables initialisation ***
    !-> minimal values
    gthis%rmin     = rmin
    gthis%thetamin = thetamin
    gthis%phimin   = phimin
    !-> length
    gthis%Lr       = Lr
    gthis%Ltheta   = Ltheta
    gthis%Lphi     = Lphi
    !-> number of points
    gthis%Nr       = Nr
    gthis%Ntheta   = Ntheta
    gthis%Nphi     = Nphi
    !-> grid step size
    gthis%dr       = Lr/float(Nr)
    gthis%dtheta   = Ltheta/float(Ntheta)
    gthis%dphi     = Lphi/float(Nphi)
    
    !*** coordinate array allocation ***
    call glob_allocate(gthis%rg,0,Nr,'rg')
    call glob_allocate(gthis%thetag,0,Ntheta,'thetag')
    call glob_allocate(gthis%phig,0,Nphi,'phig')
    call glob_allocate(gthis%cos_theta,0,Ntheta,'cos_theta')
    call glob_allocate(gthis%sin_theta,0,Ntheta,'sin_theta')
            
    if (.not.memory_test) then
      !*** coordinate array initialisation ***    
      do i1 = 0,Nr
        gthis%rg(i1) = gthis%rmin + i1*gthis%dr
      enddo      
      
      do i2 = 0,Ntheta
        gthis%thetag(i2) = gthis%thetamin + i2*gthis%dtheta
      enddo
      do i3 = 0,Nphi
        gthis%phig(i3) = gthis%phimin + i3*gthis%dphi
      enddo
      
      !*** initialization of the array for cosinus and sinus ***
      do i2 = 0,Ntheta
        gthis%cos_theta(i2) = cos(gthis%thetag(i2))
        gthis%sin_theta(i2) = sin(gthis%thetag(i2))
      end do
    end if
      
    !*** initialization of the coefficients for the ***
    !*** integral computation in r, theta, phi      ***
    call compute_coeff_intspace(gthis)
  end subroutine new_space_coord
  
      
  !-------------------------------------------------------
  ! initialisation of the velocity coordinates 
  !   -> parallel velocity
  !   -> perpendicular velocity
  !-------------------------------------------------------     
  subroutine new_velocity_coord(gthis,vparmin,mumin, &
    Lvpar,Lmu,Nvpar,Nmu)
    use globals, only : memory_test, nb_vth0, integral_vperp
    use utils_module, only : locate
    type(geometry), intent(inout) :: gthis
    real(RKIND)   , intent(in)    :: vparmin, mumin  
    real(RKIND)   , intent(in)    :: Lvpar, Lmu      
    integer       , intent(in)	  :: Nvpar, Nmu      
      
    real(RKIND) :: dmu
    real(RKIND) :: vperp, dvperp, vperp_min, vperp_max
    integer     :: i4, i5
    
    !*** variables initialization ***
    !-> minimal values
    gthis%vparmin = vparmin
    gthis%mumin   = mumin
    !-> length
    gthis%Lvpar   = Lvpar
    gthis%Lmu     = Lmu
    !-> number of points
    gthis%Nvpar   = Nvpar
    gthis%Nmu     = Nmu
    !-> grid step size
    gthis%dvpar   = Lvpar/float(Nvpar)
    
    !*** coordinate array allocation ***
    call glob_allocate(gthis%vparg,0,Nvpar,'vparg')   
    call glob_allocate(gthis%mug,0,Nmu,'mug')
         
    if (.not.memory_test) then
      !*** parallel velocity array initialization ***  
      do i4 = 0,Nvpar
        gthis%vparg(i4) = vparmin+i4*gthis%dvpar
      enddo
      if (Nmu.ne.0) then
        if (integral_vperp) then
          vperp_min = sqrt(2._RKIND*mumin)
          vperp_max = sqrt(2._RKIND*(Lmu+mumin))
          dvperp    = (vperp_max-vperp_min)/float(Nmu)
          do i5 = 0,Nmu
            vperp         = vperp_min+i5*dvperp
            gthis%mug(i5) = vperp*vperp/2._RKIND
          end do
        else
          dmu = gthis%Lmu/float(gthis%Nmu)
          do i5 = 0,Nmu
            gthis%mug(i5) = gthis%mumin+i5*dmu
          enddo
        end if
      else
        if (integral_vperp) then
          vperp_min    = sqrt(2._RKIND*mumin)
          gthis%mug(0) = vperp_min*vperp_min/2._RKIND
        else
          gthis%mug(0) = gthis%mumin
        end if
      end if
      !*** find the index ivpar0 corresponding to vpar = 0 ***
      call locate(0._RKIND,gthis%vparg,gthis%Nvpar, &
        gthis%dvpar," new_velocity_coord "//char(0), &
        gthis%ivpar0)
    end if
      
    !*** initialization of the coefficients for the ***
    !*** integral computation in vpar and mu        ***
    call compute_coeff_intvelocity(gthis)
  end subroutine new_velocity_coord
      
  
  !-----------------------------------------------------------
  ! Computation of the coefficients required for a
  !  numerical integration in the three spatial directions 
  !  (r,theta,phi), using one of the following methods :
  !    1) trapeze method
  !    2) Simpson method
  !    3) Villarceau-Boole method
  !    4) Hardy method
  !  and storage in the arrays "coeff_intdr", "coeff_intdtheta"
  !   and "coeff_intdphi"
  !------------------------------------------------------------
  subroutine compute_coeff_intspace(gthis)
    use globals, only : memory_test
    use utils_module
    type(geometry), intent(inout) :: gthis
      
    !*** array allocation ***
    call glob_allocate(gthis%coeff_intdr,0,gthis%Nr, &
      'coeff_intdr')
    call glob_allocate(gthis%coeff_intdtheta,0,gthis%Ntheta, &
      'coeff_intdtheta')
    call glob_allocate(gthis%coeff_intdphi,0,gthis%Nphi, &
      'coeff_intdphi')
      
    if (.not.memory_test) then
      select case (integration_scheme)
      case(1)
        call compute_trapeze_coeff(gthis%Nr,gthis%dr, &
          .false.,gthis%coeff_intdr)
        call compute_trapeze_coeff(gthis%Ntheta,gthis%dtheta, &
          .true.,gthis%coeff_intdtheta)
        call compute_trapeze_coeff(gthis%Nphi,gthis%dphi, &
          .true.,gthis%coeff_intdphi)
      case(2)
        call compute_simpson_coeff(gthis%Nr,gthis%dr, &
          .false.,gthis%coeff_intdr)
        call compute_simpson_coeff(gthis%Ntheta,gthis%dtheta, &
          .true.,gthis%coeff_intdtheta)
        call compute_simpson_coeff(gthis%Nphi,gthis%dphi, &
          .true.,gthis%coeff_intdphi)
      case(3)
        call compute_villarceau_coeff(gthis%Nr,gthis%dr, &
          gthis%coeff_intdr)
        call compute_villarceau_coeff(gthis%Ntheta,gthis%dtheta, &
          gthis%coeff_intdtheta)
        call compute_villarceau_coeff(gthis%Nphi,gthis%dphi, &
          gthis%coeff_intdphi)
      case(4)
        call compute_hardy_coeff(gthis%Nr,gthis%dr, &
          gthis%coeff_intdr)
        call compute_hardy_coeff(gthis%Ntheta,gthis%dtheta, &
          gthis%coeff_intdtheta)
        call compute_hardy_coeff(gthis%Nphi,gthis%dphi, &
          gthis%coeff_intdphi)
      end select
    end if
  end subroutine compute_coeff_intspace
      
  !------------------------------------------------------------
  ! Computation of the coefficients required for a
  !  numerical integration in the both velocity 
  !  directions (vpar,mu) using one of the following methods :
  !    1) trapeze method
  !    2) Simpson method
  !    3) Villarceau-Boole method
  !    4) Hardy method
  !  and storage in the arrays "coeff_intdvpar", "coeff_intdmu"
  !------------------------------------------------------------
  subroutine compute_coeff_intvelocity(gthis)
    use globals, only : memory_test, integral_vperp
    use utils_module
    type(geometry), intent(inout) :: gthis
      
    real(RKIND), dimension(0:gthis%Nmu) :: integr_coeff
    real(RKIND) :: dmu
    real(RKIND) :: vperp, dvperp, vperp_min, vperp_max
    integer     :: imu
      
    !*** array allocation ***
    call glob_allocate(gthis%coeff_intdvpar,0,gthis%Nvpar, &
      'coeff_intdvpar')
    call glob_allocate(gthis%coeff_intdmu,0,gthis%Nmu,&
      'coeff_intdmu')
      
    if (.not.memory_test) then
      select case (integration_scheme)
      case(1)
        call compute_trapeze_coeff(gthis%Nvpar,gthis%dvpar, &
          .false.,gthis%coeff_intdvpar)
      case(2)
        call compute_simpson_coeff(gthis%Nvpar,gthis%dvpar, &
          .false.,gthis%coeff_intdvpar)
      case(3)
        call compute_villarceau_coeff(gthis%Nvpar,gthis%dvpar, &
          gthis%coeff_intdvpar)
      case(4)
        call compute_hardy_coeff(gthis%Nvpar,gthis%dvpar, &
          gthis%coeff_intdvpar)
      end select
      !*** weights for integral in mu ***
      if (gthis%Nmu.eq.0) then
        gthis%coeff_intdmu(0) = 1._RKIND
      else
        if (integral_vperp) then
          vperp_min = sqrt(2._RKIND*gthis%mumin)
          vperp_max = sqrt(2._RKIND*(gthis%Lmu+gthis%mumin))
          dvperp    = (vperp_max-vperp_min)/float(gthis%Nmu)
          if (integration_scheme.eq.1) then
            ! trapeze
            integr_coeff = 0._RKIND
            do imu = 0,gthis%Nmu-1
              integr_coeff(imu)   = integr_coeff(imu)   + 1._RKIND
              integr_coeff(imu+1) = integr_coeff(imu+1) + 1._RKIND
            end do
            integr_coeff = dvperp*integr_coeff/2._RKIND
          else
            ! simpson
            integr_coeff(0)         = 1._RKIND
            integr_coeff(gthis%Nmu) = 1._RKIND
            do imu = 0,gthis%Nmu-1
              if (mod(imu,2).eq.0) then
                integr_coeff(imu) = 2._RKIND
              else
                integr_coeff(imu) = 4._RKIND
              end if
            end do
            integr_coeff = dvperp*integr_coeff/3._RKIND
          end if
          do imu = 0,gthis%Nmu
            vperp                   = sqrt(2._RKIND*gthis%mug(imu))
            gthis%coeff_intdmu(imu) = vperp*integr_coeff(imu)
          end do
        else
          dmu = gthis%Lmu/float(gthis%Nmu)
          call compute_simpson_coeff(gthis%Nmu,dmu, &
            .false.,gthis%coeff_intdmu)
        end if
      end if
    end if
  end subroutine compute_coeff_intvelocity
      
  !----------------------------------------------------- 
  ! geometry destructor
  !-----------------------------------------------------   
  subroutine del_geometry(gthis)
    type(geometry), intent(inout) :: gthis
      
    ! -> mesh grid
    call glob_deallocate(gthis%rg)
    call glob_deallocate(gthis%thetag)
    call glob_deallocate(gthis%phig)
    call glob_deallocate(gthis%vparg)
    call glob_deallocate(gthis%mug)
    ! -> integration coefficients
    call glob_deallocate(gthis%coeff_intdr)
    call glob_deallocate(gthis%coeff_intdtheta)
    call glob_deallocate(gthis%coeff_intdphi)
    call glob_deallocate(gthis%coeff_intdvpar)
    call glob_deallocate(gthis%coeff_intdmu)
    ! -> cos(theta) and sin(theta)
    call glob_deallocate(gthis%cos_theta)
    call glob_deallocate(gthis%sin_theta)
  end subroutine del_geometry
      
end module geometry_class
