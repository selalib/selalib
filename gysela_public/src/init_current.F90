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
      
!--------------------------------------------------------------
! file : init_current.f90
! date : 07/17/2009
!  Initialisation of the plasma current J, i.e :
!    - covariant coordinates
!    - contravariant coordinates
!--------------------------------------------------------------
module init_current_class
  use prec_const
  use globals, only : plasma_current, R0
  use mem_alloc_module
  use geometry_class
  use coord_system_class
  use init_profile_class
      
  implicit none
      
  !*** plasma current definition ***
  type :: init_current
    !-> covariant components of J
    real(RKIND), dimension(:,:), pointer :: mu0_Jr
    real(RKIND), dimension(:,:), pointer :: mu0_Jtheta
    real(RKIND), dimension(:,:), pointer :: mu0_Jphi
    !-> contravariant components of J
    real(RKIND), dimension(:,:), pointer :: mu0_J_gradr
    real(RKIND), dimension(:,:), pointer :: mu0_J_gradtheta
    real(RKIND), dimension(:,:), pointer :: mu0_J_gradphi
 end type init_current
    
  !******************************
  contains
  !******************************
     
  !---------------------------------------------------- 
  ! Constructor 
  !----------------------------------------------------
  subroutine new_init_current(pcthis,geom)
    type(init_current), intent(out) :: pcthis
    type(geometry)      , intent(in)  :: geom
    
    !-> allocation for covariant components
    call glob_allocate(pcthis%mu0_Jr,0,geom%Nr, &
      0,geom%Ntheta,'mu0_Jr')
    call glob_allocate(pcthis%mu0_Jtheta,0,geom%Nr, &
      0,geom%Ntheta,'mu0_Jtheta')
    call glob_allocate(pcthis%mu0_Jphi,0,geom%Nr, &
      0,geom%Ntheta,'mu0_Jphi')
    !-> allocation for contravariant components
    call glob_allocate(pcthis%mu0_J_gradr,0,geom%Nr, &
      0,geom%Ntheta,'mu0_J_gradr')
    call glob_allocate(pcthis%mu0_J_gradtheta,0,geom%Nr, &
      0,geom%Ntheta,'mu0_J_gradtheta')
    call glob_allocate(pcthis%mu0_J_gradphi,0,geom%Nr, &
      0,geom%Ntheta,'mu0_J_gradphi')
  end subroutine new_init_current
  
      
  !---------------------------------------------------- 
  ! destructor
  !---------------------------------------------------- 
  subroutine del_init_current(pcthis)
    type(init_current), intent(inout) :: pcthis
      
    !-> deallocation for covariant components
    call glob_deallocate(pcthis%mu0_Jr)
    call glob_deallocate(pcthis%mu0_Jtheta)
    call glob_deallocate(pcthis%mu0_Jphi)
    !-> deallocation for contravariant components
    call glob_deallocate(pcthis%mu0_J_gradr)
    call glob_deallocate(pcthis%mu0_J_gradtheta)
    call glob_deallocate(pcthis%mu0_J_gradphi)
  end subroutine del_init_current
      
  !------------------------------------------------------- 
  ! Initialisation of the covariant coordinates of 
  !  the plasma current, as following:
  !  if Bstar = false,
  !    Jr = 0, Jtheta = 0 and Jphi = 0
  !  else the current is of the form
  !    vec_J = J_T R grad_phi
  !  with
  !    mu0 J_T = B0 R0/R * zeta/r * 
  !       (1 + r/zeta*dzeta/dr - r/R*cos(theta))
  !            = B0 R0/R * (zeta/r + dzeta/dr - 
  !                         zeta/R*cos(theta))
  !  where zeta = r/(qR0)
  !  i.e. mu0 Jphi = mu0*J_T*R ; mu0 Jr = mu0 Jtheta = 0
  !-------------------------------------------------------
  subroutine init_plasma_current(pcthis,geom, &
    coord_sys,init_prof,zeta)
    type(init_current)        , intent(inout) :: pcthis
    type(geometry)            , intent(in)    :: geom
    type(coord_system)        , intent(in)    :: coord_sys 
    type(init_profile)        , intent(in)    :: init_prof
    real(RKIND), dimension(0:), intent(in)    :: zeta
    
    integer     :: ir, itheta
    real(RKIND) :: ri, s_i, iota_i, zeta_i, dzetadr_i
    real(RKIND) :: cos_thetaj, inv_Rij
      
    !*** allocation of the arrays ***
    call new_init_current(pcthis,geom)
      
    !*** initialization of the covariant coordinates ***
    if (.not.memory_test) then
      if ( (.not.plasma_current).or.(.not.magnetic_shear) ) then
        do itheta = 0,geom%Ntheta
          do ir = 0,geom%Nr
            pcthis%mu0_Jr(ir,itheta)     = 0._RKIND 
            pcthis%mu0_Jtheta(ir,itheta) = 0._RKIND
            pcthis%mu0_Jphi(ir,itheta)   = 0._RKIND
          end do
        end do
      else
        do itheta = 0,geom%Ntheta
          cos_thetaj = geom%cos_theta(itheta)
          do ir = 0,geom%Nr
            pcthis%mu0_Jr(ir,itheta)     = 0._RKIND
            pcthis%mu0_Jtheta(ir,itheta) = 0._RKIND
            !-> mu0 Jphi = mu0*J_T*R = B0 R0 * (zeta/r + 
            !     dzeta/dr - zeta/R*cos(theta))
            ri      = geom%rg(ir)
            inv_Rij = inv_R(ir,itheta)
            s_i     = init_prof%shear(ir)
            iota_i  = init_prof%iota(ir)
            zeta_i  = zeta(ir)
            !---> dzeta/dr = 1/(qR0)*(1-r/q*dq/dr) = 1/(qR0)*(1-s)
            dzetadr_i = (1._RKIND-s_i) * iota_i/R0
            !---> mu0 Jphi
            pcthis%mu0_Jphi(ir,itheta) = R0 * (zeta_i/ri + &
              dzetadr_i-zeta_i*inv_Rij*cos_thetaj)
          end do
        end do
      end if
      !*** computation of the contravariant coordinates ***
      call compute_contravariant_vector(coord_sys, &
        pcthis%mu0_Jr,pcthis%mu0_Jtheta,pcthis%mu0_Jphi, &
        pcthis%mu0_J_gradr,pcthis%mu0_J_gradtheta, &
        pcthis%mu0_J_gradphi)
    end if
  end subroutine init_plasma_current
end module init_current_class
