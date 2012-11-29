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
!------------------------------------------------------------------
! file : advec2D_BSL.f90
! date : 27/06/2000
!   2D advection solving, i.e. solving of :
!      df/dt+vgc.grad(f)= 0 
!   using a Backward Semi-Lagrangian method ;
!   where the projections in (r,theta) of the guiding-center
!   velocity are defined as :
!    . vgc_gradr = vpar*bstar_gradr + vExB_gradr + vD_gradr 
!    . vgc_gradtheta = vpar*bstar_gradtheta + 
!                      vExB_gradtheta + vD_gradtheta
!   with vExB the ExB drift velocity and 
!   vD the curvature drift velocity.
!
! This 2D advection is performed in the case of the
!  Backward Semi-Lagrangian (BSL) scheme, i.e:
! Let dt_star and the distance dist be defined as:
!  . dt_star = t_{n+1} - t_star
!  . dist = X(t_{n+1}) - X(t_star)
!  (knowing that X(t_{n+1}) is a node of the mesh)
! Let U be the advection field.
! Then,
!  . [X(t_{n+1})-X(t_star)]/dt_star = U(X_middle)
!     with X_middle = [X(t_{n+1})+X(t_star)]/2
!  which is equivalent to :
!   -> dist = dt_star*U(X(t_{n+1})-dist/2)  [*]
!  which gives using a Taylor expansion at second order
!   -> dist = dt_star* { U(X(t_{n+1}))- 
!             dist/2 partial_x U(X(t_{n+1})) + O(x^2)}
!  and replacing dist by its previous expression [*] and
!  keeping the terms of first order, then:
!  . dist = dt_star*U(X(t_{n+1})) -
!           0.5*(dt_star)^2*U(X(t_{n+1}))*partial_x U(X(t_{n+1}))
!
! Rk : The distance 'dist' is computed in the 'Taylor' subroutines
!      while the position of the feet of the characteristic
!      X(t_star) is computed in the 'advec2D' subroutines, 
!      by using the formula:
!       . X(t_star) = X(t_{n+1}) - dist 
!------------------------------------------------------------------
module advec2D_BSL_module
  use clock_module
  use globals, only : Nr, Ntheta, Nphi, Nvpar
  use prec_const
  use geometry_class
  use spline1d_class
  use interpolation_module
  use utils_module
      
  implicit none
        
  include "remap_inline.h"
  include "Bstar_inline.h"
  include "velocities_inline.h"
  !******************************
  contains
  !******************************
#include "remap_inline.f90"
#include "Bstar_inline.f90"
#include "velocities_inline.f90"
      
  !------------------------------------------------------------
  !   computation of the derivates of u
  !     by interpolation by cubic splines of u for a knot of 
  !     the mesh (r(ir),theta(itheta))
  !------------------------------------------------------------
  subroutine advecfield(ir,itheta,geom,ux_scoef,uy_scoef,du)
    integer                        , intent(in)  :: ir
    integer                        , intent(in)  :: itheta
    type(geometry)                 , intent(in)  :: geom
    real(RKIND), dimension(-1:,-1:), intent(in)  :: ux_scoef
    real(RKIND), dimension(-1:,-1:), intent(in)  :: uy_scoef
    real(RKIND), dimension(4)      , intent(out) :: du
      
    real(RKIND)                  :: rp, thetap
    real(RKIND)                  :: duxg_dr
    real(RKIND)                  :: duxg_dth
    real(RKIND)                  :: duyg_dr
    real(RKIND)                  :: duyg_dth
    real(RKIND), dimension(-1:2) :: srbase, stbase
    real(RKIND), dimension(-1:2) :: srbase_prime, stbase_prime
    integer                      :: i1, i2
      
    rp     = geom%rg(ir)
    thetap = geom%thetag(itheta)
      
    !*** calculate spline basis ***
    call spline_basis(ir,geom%Nr,srbase)
    call spline_basis(itheta,geom%Ntheta,stbase)
      
    !*** Calculate dux_dr(rc,thetac), dux_dth(rc,thetac)  ***
    !***   duy_dr(rc,thetac) and duy_dth(rc,thetac) ***
    call spline_basisderiv(ir,geom%Nr,geom%dr,srbase_prime)
    call spline_basisderiv(itheta,geom%Ntheta, &
      geom%dtheta,stbase_prime)
    duxg_dr  = 0._RKIND
    duxg_dth = 0._RKIND
    duyg_dr  = 0._RKIND
    duyg_dth = 0._RKIND
    do i2 = -1,2
      do i1 = -1,2
        duxg_dr = duxg_dr + ux_scoef(ir+i1,itheta+i2)*  & 
          srbase_prime(i1)*stbase(i2)
        duxg_dth = duxg_dth + ux_scoef(ir+i1,itheta+i2)*  & 
          srbase(i1)*stbase_prime(i2)
        duyg_dr = duyg_dr + uy_scoef(ir+i1,itheta+i2)*  & 
          srbase_prime(i1)*stbase(i2)
        duyg_dth = duyg_dth + uy_scoef(ir+i1,itheta+i2)*  & 
          srbase(i1)*stbase_prime(i2)
      end do
    end do
    du(1) = duxg_dr
    du(2) = duxg_dth
    du(3) = duyg_dr
    du(4) = duyg_dth
  end subroutine advecfield
      
  !*******************************************************
  ! GLOBAL 2D ADVECTION IN CARTESIAN COORDINATES BY 
  !  USING A TAYLOR ALGORITHM
  !*******************************************************
  !------------------------------------------------------------
  ! - Solving by a Taylor algorithm of :
  !    . alpha = dt*ux(xi-alpha,yj-beta)
  !    . beta  = dt*uy(xi-alpha,yj-beta)
  !  with the cubic spline interpolation of the electric field
  !
  ! - this algorithm is used for characteristic 
  !    feet computation (in cartesian coordinates)
  ! - Solution: 
  !    alpha = dt*ux(xi,yj) - 
  !            0.5*dt*dt*{ (dux/dx)(xi,yj)*ux(xi,yj) 
  !                        +(dux/dy)(xi,yj)*uy(xi,yj) }
  !    beta  = dt*uy(xi,yj) - 
  !            0.5*dt*dt*{ (duy/dx)(xi,yj)*ux(xi,yj) 
  !                        +(duy/dy)(xi,yj)*uy(xi,yj) }
  !------------------------------------------------------------
  subroutine taylor_xy(geom,dt,ux,uy,ux_scoef,uy_scoef, &
    alpha_ij,beta_ij)
    use globals, only : istart, iend, jstart, jend
    type(geometry)                 , intent(in)    :: geom
    real(RKIND)                    , intent(in)    :: dt
    real(RKIND), dimension(0:,0:)  , intent(in)    :: ux, uy
    real(RKIND), dimension(-1:,-1:), intent(in)    :: ux_scoef
    real(RKIND), dimension(-1:,-1:), intent(in)    :: uy_scoef
    real(RKIND), dimension(0:,0:)  , intent(inout) :: alpha_ij
    real(RKIND), dimension(0:,0:)  , intent(inout) :: beta_ij
      
    integer     :: ir, itheta
    real(RKIND) :: rp, invrp, cos_thetap, sin_thetap
    real(RKIND) :: uxg, uyg
    !-> u derivates at (ri,thetaj)
    real(RKIND), dimension(4) :: du
    real(RKIND) :: duxg_dr, duxg_dth, duxg_dx, duxg_dy
    real(RKIND) :: duyg_dr, duyg_dth, duyg_dx, duyg_dy
      
    do itheta = 0,Ntheta
      do ir = 0,Nr
        alpha_ij(ir,itheta) = 0._RKIND
        beta_ij(ir,itheta)  = 0._RKIND
      end do
    end do
    do itheta = 0,Ntheta-1
      cos_thetap = geom%cos_theta(itheta)
      sin_thetap = geom%sin_theta(itheta)
      do ir = 1,Nr-1
        rp    = geom%rg(ir)
        invrp = (1._RKIND/rp)
        uxg   = ux(ir,itheta)
        uyg   = uy(ir,itheta)
        !*** computation of the first derivatives of ux and uy ***
        !*** by interpolation in the (r,theta) space           ***
        call advecfield(ir,itheta,geom,ux_scoef,uy_scoef,du)
        duxg_dr  = du(1)
        duxg_dth = du(2)
        duyg_dr  = du(3)
        duyg_dth = du(4)
      
        !*** computation of dux(rp,thetap)/dx           ***
        !*** and dux(rp,thetap)/dy by :                 ***
        !***  . dux(rp,thetap)/dx = cos(thetap)*dux/dr  ***
        !***       -(1/r)*sin(thetap)*dux/dtheta        ***
        !***  . dux(rp,thetap)/dy = sin(thetap)*dux/dr  ***
        !***       +(1/r)*cos(thetap)*dux/dtheta        ***
        if (rp.ne.0._RKIND) then
          duxg_dx = cos_thetap*duxg_dr - &
            invrp*sin_thetap*duxg_dth
          duxg_dy = sin_thetap*duxg_dr + & 
            invrp*cos_thetap*duxg_dth
          duyg_dx = cos_thetap*duyg_dr - &
            invrp*sin_thetap*duyg_dth
          duyg_dy = sin_thetap*duyg_dr + &
            invrp*cos_thetap*duyg_dth 
        else
          duxg_dx = cos_thetap*duxg_dr
          duxg_dy = sin_thetap*duxg_dr
          duyg_dx = cos_thetap*duyg_dr
          duyg_dy = sin_thetap*duyg_dr
        endif
        
        !*** computation of the displacement ***
        alpha_ij(ir,itheta) = dt*uxg - &
          dt*dt*(duxg_dx*uxg+duxg_dy*uyg)
        beta_ij(ir,itheta)  = dt*uyg - &
          dt*dt*(duyg_dx*uxg+duyg_dy*uyg)
      enddo	! end do ir
    enddo	! end do itheta
    !*** boundary conditions in theta ***
    do ir = 0,Nr
      alpha_ij(ir,Ntheta) = alpha_ij(ir,0)
      beta_ij(ir,Ntheta)  = beta_ij(ir,0)
    end do
  end subroutine taylor_xy
      
  !-----------------------------------------------------
  ! 2D advection solving (in cartesian coordinates)
  !  - compute the characteristic feet and 
  !  - interpole the distribution function at this point
  ! (using the explicit taylor development)
  !----------------------------------------------------- 
  subroutine advec2D_xy(geom,dt,init_prof, &
    init_magnet,init_curr,f)
    use globals
    use fdistribu5d_class
    use init_profile_class
    use init_magnetic_class
    use init_current_class
    use utils_module
    use clock_module
    type(geometry)     , intent(in)    :: geom
    real(RKIND)        , intent(in)    :: dt
    type(init_profile) , intent(in)    :: init_prof
    type(init_magnetic), intent(in)    :: init_magnet
    type(init_current) , intent(in)    :: init_curr
    type(fdistribu5d)  , intent(inout) :: f
  
    integer     :: ir, itheta, iphi, ivpar
    real(RKIND) :: ri, vparl, rmin, rmax
    real(RKIND) :: cos_thetaj, sin_thetaj
    real(RKIND) :: drdt_tmp, dthetadt_tmp
    real(RKIND) :: dxdr_tmp, dxdtheta_tmp
    real(RKIND) :: dydr_tmp, dydtheta_tmp
    real(RKIND) :: bstar_gradr_tmp, bstar_gradtheta_tmp 
    real(RKIND) :: vExB_gradr_tmp, vExB_gradtheta_tmp
    real(RKIND) :: vD_gradtheta_tmp, vD_gradr_tmp
    real(RKIND) :: rfeet, tfeet, xfeet, yfeet
    real(RKIND) :: finterpol
    logical     :: bound_r
      
    real(RKIND), dimension(:,:), pointer :: scoef_rtheta
    real(RKIND), dimension(:,:), pointer :: rhs
    real(RKIND), dimension(:,:), pointer :: ux, uy
    real(RKIND), dimension(:,:), pointer :: ux_scoef2D, uy_scoef2D
    real(RKIND), dimension(:,:), pointer :: alpha_ij
    real(RKIND), dimension(:,:), pointer :: beta_ij
      
    !-> components of ExB drift velocity  (in (r,theta,phi))
    real(RKIND), dimension(:), pointer :: vExB_gradr_1D
    real(RKIND), dimension(:), pointer :: vExB_gradtheta_1D
      
    !-> components of the curvature drift velocity
    real(RKIND), dimension(:), pointer :: vD_gradr_1D
    real(RKIND), dimension(:), pointer :: vD_gradtheta_1D
      
    call clck_time(bclock_advec2D)
      
    rhs           => Rarray1_NrNtheta
    scoef_rtheta  => Rarray1_m1Nrp1m1Nthetap1
    ux            => Rarray2_NrNtheta
    uy            => Rarray3_NrNtheta
    ux_scoef2D    => Rarray2_m1Nrp1m1Nthetap1
    uy_scoef2D    => Rarray3_m1Nrp1m1Nthetap1
    alpha_ij      => Rarray4_NrNtheta
    beta_ij       => Rarray5_NrNtheta
    
    vExB_gradr_1D     => Rarray1_Nr
    vExB_gradtheta_1D => Rarray2_Nr
    vD_gradr_1D       => Rarray3_Nr
    vD_gradtheta_1D   => Rarray4_Nr
      
    rmin = geom%rg(0)
    rmax = geom%rg(geom%Nr)
      
    do ivpar = 0,f%n4
      vparl = geom%vparg(ivpar)
      do iphi = 0,f%n3-1
        !*** cubic spline computation of the  ***
        !***  distribution function           ***
        do itheta = 0,Ntheta-1
          do ir = 0,Nr
            rhs(ir,itheta) = f%values(ir,itheta,iphi,ivpar)
          end do
        end do
        !***  -> periodic condition in theta direction ***
        do ir = 0,Nr
          rhs(ir,Ntheta) = rhs(ir,0) 
        end do
        !***  -> boundary derivate computation ***
        call compute_spline_rtheta(f,rhs,scoef_rtheta)
        !**************************************************
        !*** Computation of the 2D advection in the     *** 
        !***  cartesian coordinates                     ***
        !***   dx/dt = ux and  dy/dt = uy               ***
        !*** where :                                    ***
        !*** . ux = dx/dr*dr/dt + dx/dtheta*dtheta/dt   ***
        !*** . uy = dy/dr*dr/dt + dy/dtheta*dtheta/dt   ***
        !*** with dx/dr     = cos(theta) ;              ***
        !***      dx/dtheta = -r*sin(theta) and         ***
        !***      dy/dr     = sin(theta) ;              ***
        !***      dy/dtheta = r*cos(theta)              ***
        !*** and  . dr/dt     = vpar*bstar_gradr +      ***
        !***                    vExB_gradr + vD_gradr   ***
        !***      . dtheta/dt = vpar*bstar_gradtheta +  ***
        !***            vExB_gradtheta + vD_gradtheta   ***
        !**************************************************
        do itheta = 0,Ntheta
          cos_thetaj = geom%cos_theta(itheta)
          sin_thetaj = geom%sin_theta(itheta)
          !-> compute vExB.gradr and vExB.gradtheta
          call compute_vExBgradr(geom,init_magnet,init_curr, &
            dJ0Phidtheta,dJ0Phidphi,0,Nr,itheta,itheta, &
            iphi,iphi,ivpar,ivpar,vExB_gradr_1D)
          call compute_vExBgradtheta(geom,init_magnet,init_curr, &
            dJ0Phidr,dJ0Phidphi,0,Nr,itheta,itheta, &
            iphi,iphi,ivpar,ivpar,vExB_gradtheta_1D)
          !-> compute vD.gradr and vD.gradtheta
          call compute_vDgradr(geom,init_magnet,init_curr, &
            0,Nr,itheta,itheta,ivpar,ivpar,vD_gradr_1D)
          call compute_vDgradtheta(geom,init_magnet,init_curr, &
            0,Nr,itheta,itheta,ivpar,ivpar,vD_gradtheta_1D)
          do ir = 0,Nr
            !-> compute bstar.gradr and bstar.gradtheta
#ifdef NOPRECOMPUTE 
              call compute_bstar_contravariant(geom, &
                init_magnet,init_curr,ir,itheta,ivpar, &
                Sbstar_gradx1=bstar_gradr_tmp, &
                Sbstar_gradx2=bstar_gradtheta_tmp)
#else
              call precomputed_bstar_gradr(ir,itheta,ivpar, &
                bstar_gradr_tmp)
              call precomputed_bstar_gradtheta(ir,itheta,ivpar, &
                bstar_gradtheta_tmp)
#endif
            !-> initialisation of the velocity projections
            vExB_gradr_tmp       = vExB_gradr_1D(ir)
            vExB_gradtheta_tmp   = vExB_gradtheta_1D(ir)
            vD_gradr_tmp         = vD_gradr_1D(ir)
            vD_gradtheta_tmp     = vD_gradtheta_1D(ir)
            !-> dx/dr, dx/dtheta, dy/dr and dy/dtheta
            ri           = geom%rg(ir)
            dxdr_tmp     = cos_thetaj
            dxdtheta_tmp = -ri*sin_thetaj
            dydr_tmp     = sin_thetaj
            dydtheta_tmp = ri*cos_thetaj
            !-> dr/dt and dtheta/dt
            drdt_tmp     = vparl*bstar_gradr_tmp + &
              vExB_gradr_tmp + vD_gradr_tmp 
            dthetadt_tmp = vparl*bstar_gradtheta_tmp + &
              vExB_gradtheta_tmp + vD_gradtheta_tmp
            !-> computation of ux and uy
            ux(ir,itheta) = dxdr_tmp*drdt_tmp + &
              dxdtheta_tmp*dthetadt_tmp
            uy(ir,itheta) = dydr_tmp*drdt_tmp + &
              dydtheta_tmp*dthetadt_tmp
          end do
        end do
        call compute_scoef2D(geom,f%nspline1d_r, &
          f%BCr_left,f%BCr_right, &
          f%pspline1d_theta,ux,ux_scoef2D)
        call compute_scoef2D(geom,f%nspline1d_r, &
          f%BCr_left,f%BCr_right, &
          f%pspline1d_theta,uy,uy_scoef2D)
        call taylor_xy(geom,dt,ux,uy, &
          ux_scoef2D,uy_scoef2D,alpha_ij,beta_ij)  
        do itheta = 0,Ntheta-1
          do ir = 0,Nr
            call ignore_r_boundary(ir,bound_r)
            if (.not.bound_r) then
              !*** computes the feet characteristic       ***
              rfeet = geom%rg(ir)
              tfeet = geom%thetag(itheta)
              !*** transformation of (r,theta) into (x,y) ***
              call cyl2cart(rfeet,tfeet,xfeet,yfeet)
              !*** 2D advection                           ***
              xfeet = xfeet - alpha_ij(ir,itheta)
              yfeet = yfeet - beta_ij(ir,itheta)
              !*** transformation of (x,y) into (r,theta) ***
              call cart2cyl(xfeet,yfeet,rfeet,tfeet)    
              !*** Check that r and theta are on the      ***
              !***   subdomain treated                    ***
              call r_verif(rmin,rmax,rfeet)
              !*** f interpolation                        ***
              call interpol2d_rtheta(geom,scoef_rtheta,geom%rg, &
                geom%thetag,rfeet,tfeet,finterpol)
              f%values(ir,itheta,iphi,ivpar) = finterpol
            end if
          end do
        end do
        !*** periodic conditions in theta ***
        f%values(0:Nr,Ntheta,iphi,ivpar) = &
          f%values(0:Nr,0,iphi,ivpar)
      end do
    end do
    call clck_time(eclock_advec2D)
    call clck_diff(bclock_advec2D,eclock_advec2D, &
      global_time_advec2D)
  end subroutine advec2D_xy
      
  !************************************************************
  ! LOCAL ADVECTIONS BY USING LOCAL CUBIC SPLINES
  !************************************************************
  !-------------------------------------------------------
  !   computation of ux(xg,yg) and uy(xg,yg) and 
  !     the derivates of u by interpolation by cubic 
  !     splines of u for a knot of the mesh 
  !     (r(ir),theta(itheta))
  ! -> developped by G. Latu for local splines
  !-------------------------------------------------------
  subroutine local_advecfield(ir,itheta,geom,ux_scoef,uy_scoef,du)
    use globals, only : istart, iend, jstart, jend
    integer             , intent(in)  :: ir, itheta
    type(geometry)      , intent(in)  :: geom
    real(RKIND), &
      dimension(-2:,-2:), intent(in)  :: ux_scoef, uy_scoef
    real(RKIND), &
      dimension(4)      , intent(out) :: du
      
    real(RKIND)                  :: duxg_dr
    real(RKIND)                  :: duxg_dth
    real(RKIND)                  :: duyg_dr
    real(RKIND)                  :: duyg_dth
    real(RKIND), dimension(-1:2) :: srbase
    real(RKIND), dimension(-1:2) :: stbase
    real(RKIND), dimension(-1:2) :: srbase_prime
    real(RKIND), dimension(-1:2) :: stbase_prime
    integer                      :: i1, i2
      
    !*** Calculate dux_dr(rc,thetac), dux_dth(rc,thetac)  ***
    !***   duy_dr(rc,thetac) and duy_dth(rc,thetac) ***
    srbase_prime(-1) = -1._RKIND/(2*geom%dr) 
    srbase_prime( 0) = 0._RKIND                
    srbase_prime( 1) = 1._RKIND/(2*geom%dr)
    srbase(-1)       = 1._RKIND/6._RKIND
    srbase( 0)       = 4._RKIND/6._RKIND
    srbase( 1)       = 1._RKIND/6._RKIND
    stbase_prime(-1) = -1._RKIND/(2*geom%dtheta)
    stbase_prime( 0) = 0._RKIND
    stbase_prime( 1) = 1._RKIND/(2*geom%dtheta)
    stbase(-1)       = 1._RKIND/6._RKIND
    stbase( 0)       = 4._RKIND/6._RKIND
    stbase( 1)       = 1._RKIND/6._RKIND
      
    duxg_dr  = 0._RKIND
    duxg_dth = 0._RKIND
    duyg_dr  = 0._RKIND
    duyg_dth = 0._RKIND
    do i2 = -1,1
      do i1 = -1,1
        duxg_dr = duxg_dr + &
          ux_scoef(ir-istart+i1,itheta-jstart+i2) * & 
          srbase_prime(i1)*stbase(i2)
        duxg_dth = duxg_dth + &
          ux_scoef(ir-istart+i1,itheta-jstart+i2) * & 
          srbase(i1)*stbase_prime(i2)
        duyg_dr = duyg_dr + &
          uy_scoef(ir-istart+i1,itheta-jstart+i2) * & 
          srbase_prime(i1)*stbase(i2)
        duyg_dth = duyg_dth + &
          uy_scoef(ir-istart+i1,itheta-jstart+i2) * & 
          srbase(i1)*stbase_prime(i2)
      end do
    end do
    du(1) = duxg_dr
    du(2) = duxg_dth
    du(3) = duyg_dr
    du(4) = duyg_dth
  end subroutine local_advecfield
      
  !-------------------------------------------------------------
  ! - Solving by a Taylor algorithm of :
  !    . alpha = dt*ux(xi-alpha,yj-beta)
  !    . beta  = dt*uy(xi-alpha,yj-beta)
  !  with the cubic spline interpolation of the electric field
  !
  ! - this algorithm is used for characteristic 
  !    feet computation (in cartesian coordinates)
  ! - Solution: 
  !    alpha = dt*ux(xi,yj) - 
  !            0.5*dt*dt*{ (dux/dx)(xi,yj)*ux(xi,yj) 
  !                        +(dux/dy)(xi,yj)*uy(xi,yj) }
  !    beta  = dt*uy(xi,yj) - 
  !            0.5*dt*dt*{ (duy/dx)(xi,yj)*ux(xi,yj) 
  !                        +(duy/dy)(xi,yj)*uy(xi,yj) }
  ! -> developped by G. Latu for local splines
  !------------------------------------------------------------
  subroutine local_taylor(geom,dt,ux_scoef,uy_scoef, &
    alpha_ij,beta_ij,rhsux,rhsuy)
    use globals, only : istart, iend, jstart, jend, &
      istart_buf, jstart_buf, stencil
    type(geometry)                     , intent(in)    :: geom
    real(RKIND)                        , intent(in)    :: dt
    real(RKIND), dimension(-2:,-2:)    , intent(in)    :: ux_scoef
    real(RKIND), dimension(-2:,-2:)    , intent(in)    :: uy_scoef
    real(RKIND), dimension(0:,0:)      , intent(inout) :: alpha_ij
    real(RKIND), dimension(0:,0:)      , intent(inout) :: beta_ij
    real(RKIND), &
      dimension(istart_buf:,jstart_buf:), intent(in)   :: rhsux
    real(RKIND), &
      dimension(istart_buf:,jstart_buf:), intent(in)   :: rhsuy
      
    integer     :: ir, itheta  
    real(RKIND) :: rp, thetap, costhetap, sinthetap, invrp, dt2
    !* ux(rp,thetap) and uy(rp,thetap) *
    real(RKIND) :: uxg, uyg, gux, guy
    !* u derivates at (rp,thetap) *			
    real(RKIND), dimension(4) :: du
    real(RKIND) :: duxg_dr, duxg_dth, duxg_dx, duxg_dy
    real(RKIND) :: duyg_dr, duyg_dth, duyg_dx, duyg_dy
      
    do itheta = jstart,jend
      do ir = istart,iend
        alpha_ij(ir,itheta) = 0._RKIND
        beta_ij(ir,itheta)  = 0._RKIND
      end do
    end do
    dt2 = dt*dt
    do itheta = jstart,jend
      sinthetap = geom%sin_theta(itheta)
      costhetap = geom%cos_theta(itheta)
      do ir = istart,iend
        rp    = geom%rg(ir)
        invrp = (1._RKIND/rp)
        !*** computation of ux(rp,thetap) and uy(rp,thetap)  ***
        !***  by interpolation in the (r,theta) space ***
        call local_advecfield(ir,itheta,geom,ux_scoef,uy_scoef,du)
        duxg_dr  = du(1)
        duxg_dth = du(2)
        duyg_dr  = du(3)
        duyg_dth = du(4)
        uxg      = rhsux(ir, itheta)
        uyg      = rhsuy(ir, itheta)
      
        !*** computation of dux(rp,thetap)/dx           ***
        !*** and dux(rp,thetap)/dy by :                 ***
        !***  . dux(rp,thetap)/dx = cos(thetap)*dux/dr  ***
        !***       -(1/r)*sin(thetap)*dux/dtheta        ***
        !***  . dux(rp,thetap)/dy = sin(thetap)*dux/dr  ***
        !***       +(1/r)*cos(thetap)*dux/dtheta        ***
        if (rp.ne.0._RKIND) then
          duxg_dx = costhetap*duxg_dr - &
            invrp*sinthetap*duxg_dth
          duxg_dy = sinthetap*duxg_dr + & 
            invrp*costhetap*duxg_dth
          duyg_dx = costhetap*duyg_dr - &
            invrp*sinthetap*duyg_dth
          duyg_dy = sinthetap*duyg_dr + &
            invrp*costhetap*duyg_dth 
        else
          duxg_dx = costhetap*duxg_dr
          duxg_dy = sinthetap*duxg_dr
          duyg_dx = costhetap*duyg_dr
          duyg_dy = sinthetap*duyg_dr
        endif
      
        !*** computation of the displacement ***
        alpha_ij(ir,itheta) = dt*uxg - &
          0.5_RKIND*dt2*(duxg_dx*uxg+duxg_dy*uyg)
        beta_ij(ir,itheta)  = dt*uyg - &
          0.5_RKIND*dt2*(duyg_dx*uxg+duyg_dy*uyg)
      enddo	! end do ir
    enddo	! end do itheta
  end subroutine local_taylor
      
  !-----------------------------------------------------
  ! Copy of the receive buffers 
  !-----------------------------------------------------
  subroutine comm_copyrecvbuffers(d_r,d_theta,zz,geom, &
    init_prof,init_magnet,init_curr,f,ivpar) 
    use globals
#ifdef _OPENMP
    use OMPutils_module, only : omp_get_thread_num
#endif
    use OMPutils_module, only : Romp_vect_mpstencil, &
      Romp_fillde_112, Romp1_0Nr, Romp2_0Nr, Romp3_0Nr, &
      Romp4_0Nr
    use MPIutils_module
    use local_spline_module
    use fdistribu5d_class
    use init_profile_class
    use init_magnetic_class
    use init_current_class
    use clock_module
    type(geometry)     , intent(in) :: geom
    real(RKIND)        , intent(in) :: d_r
    real(RKIND)        , intent(in) :: d_theta
    integer            , intent(in) :: zz
    integer            , intent(in) :: ivpar
    type(init_profile) , intent(in) :: init_prof
    type(init_magnetic), intent(in) :: init_magnet
    type(init_current) , intent(in) :: init_curr
    type(fdistribu5d)  , intent(in) :: f
      
    ! -> variables for OpenMP parallelization
    integer     :: b_phi, ir, itheta, iphi, i, j, k, zifi, bifi
    integer     :: tid
    integer     :: istart_clamp,iend_clamp
    real(RKIND) :: NWcorner, NEcorner, SWcorner, SEcorner
    real(RKIND) :: derivtmp
    real(RKIND) :: rder_dr, lder_dr
    real(RKIND), dimension(:), pointer :: vect, fillde
      
    ! -> variables for 2D advection
    integer     :: q
    real(RKIND) :: ri, vparl
    real(RKIND) :: cos_thetaj, sin_thetaj
    real(RKIND) :: drdt_tmp, dthetadt_tmp
    real(RKIND) :: dxdr_tmp, dxdtheta_tmp
    real(RKIND) :: dydr_tmp, dydtheta_tmp
    real(RKIND) :: bstar_gradr_tmp, bstar_gradtheta_tmp
    real(RKIND) :: vExB_gradr_tmp, vExB_gradtheta_tmp
    real(RKIND) :: vD_gradr_tmp, vD_gradtheta_tmp
    real(RKIND) :: rfeet, tfeet, xfeet, yfeet
    real(RKIND) :: finterpol
      
    !-> components of ExB drift velocity
    real(RKIND), dimension(:), pointer :: vExB_gradr_1D
    real(RKIND), dimension(:), pointer :: vExB_gradtheta_1D
      
    !-> components of curvature drift velocity
    real(RKIND), dimension(:), pointer :: vD_gradr_1D
    real(RKIND), dimension(:), pointer :: vD_gradtheta_1D
      
    istart_clamp = max(0,istart_buf)
    iend_clamp   = min(Nr,iend_buf)
      
#ifdef _OPENMP
!$OMP PARALLEL private(tid,vect,fillde,iphi, &
!$OMP lder_dr,rder_dr,b_phi,bifi,zifi,ir,i,k, &
!$OMP NWcorner,NEcorner,SWcorner,SEcorner,j, &
!$OMP itheta,derivtmp,q,finterpol,ri,vparl, &
!$OMP cos_thetaj,sin_thetaj, &
!$OMP drdt_tmp,dthetadt_tmp,dxdr_tmp,dxdtheta_tmp, &
!$OMP dydr_tmp,dydtheta_tmp, &
!$OMP bstar_gradr_tmp,bstar_gradtheta_tmp, &
!$OMP vExB_gradr_tmp,vExB_gradtheta_tmp, &
!$OMP vD_gradr_tmp,vD_gradtheta_tmp, &
!$OMP vExB_gradr_1D,vExB_gradtheta_1D, &
!$OMP vD_gradr_1D,vD_gradtheta_1D, &
!$OMP rfeet,tfeet,xfeet,yfeet) default(shared) 
!$OMP BARRIER
    tid = 1+omp_get_thread_num()
#else
    tid = 1
#endif
    vect   => Romp_vect_mpstencil(tid)%val
    fillde => Romp_fillde_112(tid)%val
 
    vExB_gradr_1D     => Romp1_0Nr(tid)%val
    vExB_gradtheta_1D => Romp2_0Nr(tid)%val
    vD_gradr_1D       => Romp3_0Nr(tid)%val
    vD_gradtheta_1D   => Romp4_0Nr(tid)%val
      
    if (zz .ne. 0) then
      lder_dr = leftconst / d_r
      rder_dr = rightconst / d_r
      b_phi   = mod(zz-bloc_phi,2*bloc_phi)
!$OMP DO SCHEDULE(static) 
      do iphi = 0, 3*bloc_phi-1
        bifi = iphi
        zifi = 3*b_phi+bifi
        !*** f values ***
        do ir = istart, iend
          i = ir - istart
          tmprhs2(ir,jstart-2,zifi) = rbufWW(i,bifi) 
          tmprhs2(ir,jend+1,zifi)   = rbufEE(i,bifi)  
          tmprhs2(ir,jend+3,zifi)   = rbufEE(2*dom_r+i,bifi)  
          tmprhs2(ir,jend+2,zifi)   = rbufEE(dom_r+i,bifi) + &
            leftderive10( &
            tmprhs2(ir,jend+1-stencil:jend,zifi),d_theta) 
          tmprhs2(ir,jstart-1,zifi) = rbufWW(dom_r+i,bifi) + &
            rightderive10( &
            tmprhs2(ir,jstart+1:jstart+stencil,zifi),d_theta)
        end do
        !*** partial difference in r and theta of f ***
        !***   in the corners of f                  ***
        do k = 1,stencil
          vect(k) = rightderive10( &
            tmprhs2(istart+1:istart+stencil,jstart+k,zifi),d_r)
        enddo
        NWcorner = rightderive10(vect(1:stencil),d_theta)
        do k = -stencil,-1
          vect(k) = rightderive10( &
            tmprhs2(istart+1:istart+stencil,jend+1+k,zifi),d_r)
        enddo
        NEcorner = leftderive10(vect(-stencil:-1),d_theta)
        do k = 1,stencil
          vect(k) = leftderive10( &
            tmprhs2(iend+1-stencil:iend,jstart+k,zifi),d_r)
        enddo
        SWcorner = rightderive10(vect(1:stencil),d_theta)
        do k = -stencil,-1
          vect(k) = leftderive10( &
            tmprhs2(iend+1-stencil:iend,jend+1+k,zifi),d_r)
        enddo
        SEcorner = leftderive10(vect(-stencil:-1),d_theta)
      
        if (iend .ne. Nr) then
          do itheta = jstart,jend
            j = itheta - jstart
            tmprhs2(iend+2,itheta,zifi) = &
              rbufSS(dom_theta+j,bifi) + leftderive10( &
              tmprhs2(iend+1-stencil:iend, itheta , zifi), d_r)
            tmprhs2(iend+1,itheta,zifi) = rbufSS(j,bifi)
            tmprhs2(iend+3,itheta,zifi) = rbufSS(2*dom_theta+j,bifi)
          end do
          tmprhs2(iend+1,jstart-2,zifi) = rbufSW(0,bifi)
          tmprhs2(iend+3,jstart-2,zifi) = rbufSW(1,bifi)
          tmprhs2(iend+1,jend+1,zifi)   = rbufSE(0,bifi)
          tmprhs2(iend+3,jend+1,zifi)   = rbufSE(1,bifi)
          tmprhs2(iend+1,jend+3,zifi)   = rbufSE(2,bifi)
          tmprhs2(iend+3,jend+3,zifi)   = rbufSE(3,bifi)
          tmprhs2(iend+1,jstart-1,zifi) = &
            rbufSS(3*dom_theta+0,bifi) + rbufSW(4,bifi) !$2 10
          tmprhs2(iend+3,jstart-1,zifi) = rbufSW(5,bifi) + &
            rbufSS(3*dom_theta+1,bifi) !$11 3
          tmprhs2(iend+2,jstart-2,zifi) = &
            rbufWW(2*dom_r+1,bifi) + rbufSW(3,bifi) !$9 8
          tmprhs2(iend+2,jend+1,zifi)   = rbufSE(5,bifi) + &
            rbufEE(3*dom_r+2,bifi) !$1 5
          tmprhs2(iend+1,jend+2,zifi)   = rbufSE(6,bifi) + &
            rbufSS(3*dom_theta+2,bifi) !$2 10
          tmprhs2(iend+3,jend+2,zifi)   = rbufSE(7,bifi) + &
            rbufSS(3*dom_theta+3,bifi) !$3 11
          tmprhs2(iend+2,jend+3,zifi)   = &
            rbufEE(3*dom_r+3,bifi) + rbufSE(8,bifi) !$7 6
          tmprhs2(iend+2,jstart-1,zifi) = &
            SWcorner + rbufSS(3*dom_theta+4,bifi) + &
            rbufWW(2*dom_r+3,bifi) + rbufSW(2,bifi) !*7 8 9
          tmprhs2(iend+2,jend+2,zifi)   = SEcorner + &
            rbufSS(3*dom_theta+5,bifi) + &
            rbufEE(3*dom_r+5,bifi) + rbufSE(4,bifi) !*11 12 13
        else
          !-> compact domain in r direction, f is constant 
          !->   outside the domain
          do itheta = jstart,jend
            j = itheta - jstart
            tmprhs2(iend+2,itheta,zifi) = &
              rder_dr*tmprhs2(iend,itheta,zifi) + leftderive10( &
              tmprhs2(iend+1-stencil:iend,itheta,zifi),d_r)
            tmprhs2(iend+1,itheta,zifi) = tmprhs2(iend,itheta,zifi)
            tmprhs2(iend+3,itheta,zifi) = tmprhs2(iend,itheta,zifi)
          end do
          !-> compact domain in r direction, f is constant 
          !->   outside the domain
          tmprhs2(iend+1,jstart-2,zifi) = &
            tmprhs2(iend,jstart-2,zifi)
          tmprhs2(iend+3,jstart-2,zifi) = &
            tmprhs2(iend,jstart-2,zifi)
          tmprhs2(iend+1,jend+1,zifi)   = tmprhs2(iend,jend+1,zifi)
          tmprhs2(iend+3,jend+1,zifi)   = tmprhs2(iend,jend+1,zifi)
          tmprhs2(iend+1,jend+3,zifi)   = tmprhs2(iend,jend+3,zifi)
          tmprhs2(iend+3,jend+3,zifi)   = tmprhs2(iend,jend+3,zifi)
          !-> compact domain in r direction, f is constant 
          !->   outside the domain
          tmprhs2(iend+1,jstart-1,zifi) = &
            tmprhs2(iend,jstart-1,zifi)
          tmprhs2(iend+3,jstart-1,zifi) = &
            tmprhs2(iend,jstart-1,zifi)
          tmprhs2(iend+2,jstart-2,zifi) = rbufWW(2*dom_r+1,bifi) + &
            rder_dr * tmprhs2(iend,jstart-2,zifi)!$9 8
          tmprhs2(iend+2,jend+1,zifi)   = &
            tmprhs2(iend,jend+1,zifi)*rder_dr + &
            rbufEE(3*dom_r+2,bifi) !$1 5
          tmprhs2(iend+1,jend+2,zifi)   = tmprhs2(iend,jend+2,zifi)
          tmprhs2(iend+3,jend+2,zifi)   = tmprhs2(iend,jend+2,zifi)
          tmprhs2(iend+2,jend+3,zifi)   = rbufEE(3*dom_r+3,bifi) + &
            rder_dr * tmprhs2(iend,jend+3,zifi) !$7 6
          !-> compact domain in r direction, f is constant 
          !->   outside the domain
          derivtmp = rder_dr * tmprhs2(iend,jstart-1,zifi)
          tmprhs2(iend+2,jstart-1,zifi) = SWcorner + derivtmp &
            + rbufWW(2*dom_r+3,bifi) 
          derivtmp = rder_dr * tmprhs2(iend,jend+2,zifi)
          tmprhs2(iend+2,jend+2,zifi) = SEcorner + derivtmp &
            + rbufEE(3*dom_r+5,bifi) 
        end if
      
        if (istart .ne. 0) then
          do itheta = jstart,jend
            j = itheta - jstart
            tmprhs2(istart-1,itheta,zifi) = &
              rbufNN(dom_theta+j,bifi) + rightderive10( &
              tmprhs2(istart+1:istart+stencil,itheta,zifi),d_r) 
            tmprhs2(istart-2,itheta,zifi) = rbufNN(j,bifi) 
          end do
          tmprhs2(istart-2,jstart-2,zifi) = rbufNW(0,bifi)
          tmprhs2(istart-2,jend+1,zifi)   = rbufNE(0,bifi)
          tmprhs2(istart-2,jend+3,zifi)   = rbufNE(1,bifi)
          tmprhs2(istart-1,jstart-2,zifi) = rbufNW(2,bifi) + &
            rbufWW(2*dom_r+0,bifi) !$9 8      
          tmprhs2(istart-2,jstart-1,zifi) = rbufNW(3,bifi) + &
            rbufNN(2*dom_theta+1,bifi) !$12 4     
          tmprhs2(istart-2,jend+2,zifi)   = rbufNE(3,bifi) + &
            rbufNN(2*dom_theta+0,bifi) !$4 12
          tmprhs2(istart-1,jend+1,zifi)   = &
            rbufEE(3*dom_r+0,bifi) + rbufNE(4,bifi) !$1 5
          tmprhs2(istart-1,jend+3,zifi)   = rbufNE(5,bifi) + &
            rbufEE(3*dom_r+1,bifi) !$ 7 6
          tmprhs2(istart-1,jstart-1,zifi) = NWcorner + &
            rbufNN(2*dom_theta+2,bifi) + &
            rbufWW(2*dom_r+2,bifi) + rbufNW(1,bifi) !*1 2 3
          tmprhs2(istart-1,jend+2,zifi) = NEcorner + &
            rbufNN(2*dom_theta+3,bifi) + &
            rbufEE(3*dom_r+4,bifi) + rbufNE(2,bifi) !*4 5 6
        else
          !-> compact domain in r direction, f is constant 
          !->   outside the domain
          do itheta = jstart,jend
            j = itheta - jstart
            tmprhs2(istart-1,itheta,zifi) = &
              lder_dr*tmprhs2(istart,itheta,zifi) + &
              rightderive10( &
              tmprhs2(istart+1:istart+stencil,itheta ,zifi),d_r) 
            tmprhs2(istart-2,itheta,zifi) = &
              tmprhs2(istart,itheta,zifi)
          end do
          !-> compact domain in r direction, f is constant 
          !->   outside the domain
          tmprhs2(istart-2,jstart-2,zifi) = &
            tmprhs2(istart,jstart-2,zifi)
          tmprhs2(istart-2,jend+1,zifi)   = &
            tmprhs2(istart,jend+1,zifi)
          tmprhs2(istart-2,jend+3,zifi)   = &
            tmprhs2(istart,jend+3,zifi)
          !-> compact domain in r direction, f is constant 
          !->   outside the domain
          tmprhs2(istart-1,jstart-2,zifi) = &
            lder_dr*tmprhs2(istart,jstart-2,zifi) + &
            rbufWW(2*dom_r+0,bifi) !$9 8      
          tmprhs2(istart-2,jstart-1,zifi) = &
            tmprhs2(istart,jstart-1,zifi)
          tmprhs2(istart-2,jend+2,zifi)   = &
            tmprhs2(istart,jend+2,zifi)
          tmprhs2(istart-1,jend+1,zifi)   = &
            rbufEE(3*dom_r+0,bifi) + &
            lder_dr * tmprhs2(istart,jend+1,zifi) !$1 5
          tmprhs2(istart-1,jend+3,zifi)   = &
            lder_dr * tmprhs2(istart,jend+3,zifi) + &
            rbufEE(3*dom_r+1,bifi) !$ 7 6
          !-> compact domain in r direction, f is constant 
          !->   outside the domain
          derivtmp = lder_dr * tmprhs2(istart,jstart-1,zifi)
          tmprhs2(istart-1,jstart-1,zifi) = NWcorner + derivtmp &
            + rbufWW(2*dom_r+2,bifi) 
          derivtmp = lder_dr * tmprhs2(istart,jend+2,zifi) 
          tmprhs2(istart-1,jend+2,zifi) = NEcorner + derivtmp &
            + rbufEE(3*dom_r+4,bifi) 
        end if
      enddo
!$OMP END DO
    end if
      
    if (zz .ne. Nphi) then
      vparl = geom%vparg(ivpar)
!$OMP DO SCHEDULE(static)
      do iphi = zz,zz+bloc_phi-1
        do k = jstart_buf,jend_buf
          !**************************************************
          !*** Computation of the 2D advection in the     *** 
          !***  cartesian coordinates                     ***
          !***   dx/dt = ux and  dy/dt = uy               ***
          !*** where :                                    ***
          !*** . ux = dx/dr*dr/dt + dx/dtheta*dtheta/dt   ***
          !*** . uy = dy/dr*dr/dt + dy/dtheta*dtheta/dt   ***
          !*** with dx/dr     = cos(theta) ;              ***
          !***      dx/dtheta = -r*sin(theta) and         ***
          !***      dy/dr     = sin(theta) ;              ***
          !***      dy/dtheta = r*cos(theta)              ***
          !*** and  . dr/dt     = vpar*bstar_gradr +      ***
          !***                    vExB_gradr + vD_gradr   ***
          !***      . dtheta/dt = vpar*bstar_gradtheta +  ***
          !***            vExB_gradtheta + vD_gradtheta   ***
          !**************************************************
          itheta     = mod(k+2*Ntheta,Ntheta)
          cos_thetaj = geom%cos_theta(itheta)
          sin_thetaj = geom%sin_theta(itheta)
          !-> compute vExB.gradr and vExB.gradtheta
          call compute_vExBgradr(geom,init_magnet,init_curr, &
            dJ0Phidtheta,dJ0Phidphi,istart_clamp,iend_clamp, &
            itheta,itheta,iphi,iphi,ivpar,ivpar,vExB_gradr_1D)
          call compute_vExBgradtheta(geom,init_magnet,init_curr, &
            dJ0Phidr,dJ0Phidphi,istart_clamp,iend_clamp, &
            itheta,itheta,iphi,iphi,ivpar,ivpar,vExB_gradtheta_1D)
          !-> compute vD.gradr and vD.gradtheta
          call compute_vDgradr(geom,init_magnet,init_curr, &
            istart_clamp,iend_clamp,itheta,itheta, &
            ivpar,ivpar,vD_gradr_1D)
          call compute_vDgradtheta(geom,init_magnet,init_curr, &
            istart_clamp,iend_clamp,itheta,itheta, &
            ivpar,ivpar,vD_gradtheta_1D)
          do q = istart_buf,iend_buf
            if (q.lt.0) then
              ir = 0
            else if (q.gt.Nr) then
              ir = Nr
            else
              ir = q
            endif
            !-> compute bstar.gradr and bstar.gradtheta
#ifdef NOPRECOMPUTE 
              call compute_bstar_contravariant(geom, &
                init_magnet,init_curr,ir,itheta,ivpar, &
                Sbstar_gradx1=bstar_gradr_tmp, &
                Sbstar_gradx2=bstar_gradtheta_tmp)
#else
              call precomputed_bstar_gradr(ir,itheta,ivpar, &
                bstar_gradr_tmp)
              call precomputed_bstar_gradtheta(ir,itheta,ivpar, &
                bstar_gradtheta_tmp)
#endif
            !-> initialisation of the velocity projections
            vExB_gradr_tmp     = vExB_gradr_1D(ir-istart_clamp)
            vExB_gradtheta_tmp = vExB_gradtheta_1D(ir-istart_clamp)
            vD_gradr_tmp       = vD_gradr_1D(ir-istart_clamp)
            vD_gradtheta_tmp   = vD_gradtheta_1D(ir-istart_clamp)
            !-> dx/dr, dx/dtheta, dy/dr and dy/dtheta
            ri           = geom%rg(ir)
            dxdr_tmp     = cos_thetaj
            dxdtheta_tmp = -ri*sin_thetaj
            dydr_tmp     = sin_thetaj
            dydtheta_tmp = ri*cos_thetaj
            !-> dr/dt and dtheta/dt
            drdt_tmp     = vparl*bstar_gradr_tmp + &
              vExB_gradr_tmp + vD_gradr_tmp 
            dthetadt_tmp = vparl*bstar_gradtheta_tmp + &
              vExB_gradtheta_tmp + vD_gradtheta_tmp
            !-> computation of ux and uy
            b_phi                = 3*mod(iphi,2*bloc_phi)
            tmprhs2(q,k,b_phi)   = f%values(ir,k,iphi,ivpar)
            tmprhs2(q,k,b_phi+1) = dxdr_tmp*drdt_tmp + &
              dxdtheta_tmp*dthetadt_tmp
            tmprhs2(q,k,b_phi+2) = dydr_tmp*drdt_tmp + &
              dydtheta_tmp*dthetadt_tmp
          end do
        end do
      enddo
!$OMP END DO
    endif
      
    if (zz .ne. Nphi) then
      b_phi = mod(zz,2*bloc_phi)
!$OMP DO SCHEDULE(static) 
      do iphi = 0,3*bloc_phi-1
        bifi       = iphi
        zifi       = 3*b_phi+bifi
        fillde(1)  = rightderive10( &
          tmprhs2(istart+1:istart+stencil,jstart,zifi),d_r) 
        fillde(2)  = rightderive10( &
          tmprhs2(istart,jstart+1:jstart+stencil,zifi),d_theta)
        fillde(3)  = rightderive10( &
          tmprhs2(istart+1,jstart+1:jstart+stencil,zifi), &
          d_theta)
        fillde(4)  = rightderive10( &
          tmprhs2(iend,jstart+1:jstart+stencil,zifi),d_theta)
        fillde(5)  = leftderive10( &
          tmprhs2(iend+1-stencil:iend,jstart,zifi),d_r) 
        fillde(6)  = rightderive10( &
          tmprhs2(istart+1:istart+stencil,jstart+1,zifi),d_r) 
        fillde(7)  = leftderive10( &
          tmprhs2(iend+1-stencil:iend,jstart+1,zifi),d_r) 
        fillde(8)  = rightderive10( &
          tmprhs2(istart+1:istart+stencil,jend,zifi),d_r) 
        fillde(9)  = leftderive10( &
          tmprhs2(iend+1-stencil:iend,jend,zifi),d_r) 
        fillde(10) = leftderive10( &
          tmprhs2(istart,jend+1-stencil:jend,zifi),d_theta)
        fillde(11) = leftderive10( &
          tmprhs2(istart+1,jend+1-stencil:jend,zifi),d_theta)
        fillde(12) = leftderive10( &
          tmprhs2(iend,jend+1-stencil:jend,zifi),d_theta)
        vect       = 0
        do k = 1,stencil
          vect(k) = rightderive10( &
            tmprhs2(istart+1:istart+stencil,jstart+k,zifi),d_r)
        enddo
        NWcorner = rightderive10(vect(1:stencil),d_theta)
        vect     = 0
        do k = -stencil,-1
          vect(k) = rightderive10( &
            tmprhs2(istart+1:istart+stencil,jend+1+k,zifi),d_r)
        enddo
        NEcorner = leftderive10(vect(-stencil:-1),d_theta)
        vect     = 0
        do k = 1,stencil
          vect(k) = leftderive10( &
            tmprhs2(iend+1-stencil:iend,jstart+k,zifi),d_r)
        enddo
        SWcorner = rightderive10(vect(1:stencil),d_theta)
        vect     = 0
        do k = -stencil,-1
          vect(k) = leftderive10( &
            tmprhs2(iend+1-stencil:iend,jend+1+k,zifi),d_r)
        enddo
        SEcorner = leftderive10(vect(-stencil:-1),d_theta)
      
        !--------------------------------------------------------
      
        do ir = istart,iend
          i = ir - istart
          sbufWW(i,bifi)         = tmprhs2(ir,jstart,zifi)
          sbufWW(dom_r+i,bifi)   = rightderive10( &
            tmprhs2(ir,jstart+1:jstart+stencil,zifi),d_theta)
          sbufWW(2*dom_r+i,bifi) = tmprhs2(ir,jstart+1,zifi) ! DYS  
          sbufEE(i,bifi)         = tmprhs2(ir,jend,zifi)
          sbufEE(dom_r+i,bifi)   = leftderive10( &
            tmprhs2(ir,jend+1-stencil:jend,zifi),d_theta)
        enddo
        sbufWW(3*dom_r+0,bifi) = fillde(1) !$1
        sbufWW(3*dom_r+1,bifi) = fillde(6) !$6
        sbufWW(3*dom_r+2,bifi) = fillde(5) !$5
        sbufWW(3*dom_r+3,bifi) = fillde(7) !$7
        sbufWW(3*dom_r+4,bifi) = NWcorner  !*5
        sbufWW(3*dom_r+5,bifi) = SWcorner  !*11
      
        sbufEE(2*dom_r+0,bifi) = fillde(8) !$8
        sbufEE(2*dom_r+1,bifi) = fillde(9) !$9
        sbufEE(2*dom_r+2,bifi) = NEcorner  !*2
        sbufEE(2*dom_r+3,bifi) = SEcorner  !*8
      
        !--------------------------------------------------------
      
        do itheta = jstart,jend
          j = itheta - jstart
          sbufNN(j,bifi)             = tmprhs2(istart,itheta,zifi)
          sbufNN(dom_theta+j,bifi)   = rightderive10( &
            tmprhs2(istart+1:istart+stencil,itheta,zifi),d_r) 
          sbufNN(2*dom_theta+j,bifi) = &
            tmprhs2(istart+1,itheta,zifi) ! DYS
      
          sbufSS(j,bifi)             = tmprhs2(iend,itheta,zifi)
          sbufSS(dom_theta+j,bifi)   = leftderive10( &
            tmprhs2(iend+1-stencil:iend,itheta,zifi),d_r)
        enddo
      
        sbufNN(3*dom_theta+0,bifi) = fillde(2)  !$2
        sbufNN(3*dom_theta+1,bifi) = fillde(3)  !$3
        sbufNN(3*dom_theta+2,bifi) = fillde(10) !$10
        sbufNN(3*dom_theta+3,bifi) = fillde(11) !$11
        sbufNN(3*dom_theta+4,bifi) = NWcorner   !*7
        sbufNN(3*dom_theta+5,bifi) = NEcorner   !*10
      
        sbufSS(2*dom_theta+0,bifi) = fillde(12) !$12
        sbufSS(2*dom_theta+1,bifi) = fillde(4)  !$4
        sbufSS(2*dom_theta+2,bifi) = SWcorner   !*1
        sbufSS(2*dom_theta+3,bifi) = SEcorner   !*4
      
        !----------
      
        sbufNW(0,bifi) = tmprhs2(istart,jstart,zifi)
        sbufNW(1,bifi) = tmprhs2(istart+1,jstart,zifi)   ! DYS
        sbufNW(2,bifi) = tmprhs2(istart,jstart+1,zifi)   ! DYS
        sbufNW(3,bifi) = tmprhs2(istart+1,jstart+1,zifi) ! DYS
        sbufNW(4,bifi) = NWcorner  !*12
        sbufNW(5,bifi) = fillde(1) !$1
        sbufNW(6,bifi) = fillde(2) !$2
        sbufNW(7,bifi) = fillde(3) !$3
        sbufNW(8,bifi) = fillde(6) !$6
      
        sbufNE(0,bifi) = tmprhs2(istart,jend,zifi)
        sbufNE(1,bifi) = tmprhs2(istart+1,jend,zifi) ! DYS
        sbufNE(2,bifi) = NEcorner   !*9
        sbufNE(3,bifi) = fillde(8)  !$8
        sbufNE(4,bifi) = fillde(10) !$10
        sbufNE(5,bifi) = fillde(11) !$11
      
        sbufSW(0,bifi) = tmprhs2(iend,jstart,zifi)
        sbufSW(1,bifi) = tmprhs2(iend,jstart+1,zifi) ! DYS
        sbufSW(2,bifi) = SWcorner  !*6
        sbufSW(3,bifi) = fillde(4) !$4
        sbufSW(4,bifi) = fillde(5) !$5
        sbufSW(5,bifi) = fillde(7) !$7
      
        sbufSE(0,bifi) = tmprhs2(iend,jend,zifi)
        sbufSE(1,bifi) = SEcorner   !*3
        sbufSE(2,bifi) = fillde(9)  !$9
        sbufSE(3,bifi) = fillde(12) !$12
      enddo
!$OMP END DO
    endif
!$OMP BARRIER
!$OMP END PARALLEL
  end subroutine comm_copyrecvbuffers
      
  !-----------------------------------------------------
  ! 2D advection solving (in cartesian coordinates)
  !  - compute the characteristic feet and 
  !  - interpole the distribution function at this point
  ! (using the explicit taylor development)
  ! -> developped by G. Latu for local splines
  !----------------------------------------------------- 
  subroutine local_advec2D_xy(geom,dt, &
    init_prof,init_magnet,init_curr,f,maxalp_xy)
    use globals
    use MPIutils_module
#ifdef _OPENMP
    use OMPutils_module, only : omp_get_thread_num
#endif
    use OMPutils_module, only : Romp1_0Nr_0Ntheta, &
      Romp2_0Nr_0Ntheta
    use fdistribu5d_class
    use init_profile_class
    use init_magnetic_class
    use init_current_class
    use utils_module
    use clock_module
    implicit none
!R3 #include "r3_info.h" !R3
    type(geometry)     , intent(in)    :: geom
    real(RKIND)        , intent(in)    :: dt
    type(init_profile) , intent(in)    :: init_prof
    type(init_magnetic), intent(in)    :: init_magnet
    type(init_current) , intent(in)    :: init_curr
    type(fdistribu5d)  , intent(inout) :: f
    real(RKIND)        , intent(out)   :: maxalp_xy
      
    integer     :: ir, itheta, b_phi, iphi, ivpar, zz, k, q
    integer     :: iter_save
    real(RKIND) :: rfeet, tfeet, xfeet, yfeet
    real(RKIND) :: finterpol, outofdom_percent
      
    real(RKIND), dimension(:,:), pointer :: alpha_ij
    real(RKIND), dimension(:,:), pointer :: beta_ij
      
    integer           :: tid
    integer(TIMEPREC) :: tdeb, tfin
    real(RKIND)       :: tcomm2, tcomm3, ttmprhs, tadvect
      
!R3 call r3_info_begin(r3_info_index_0,'CPU_local_advec2d_xy') !R3
      
    call clck_time(bclock_advec2D)
      
    maxalp_xy  = 0._RKIND
      
    !*** counter for particles out of domain for each iter ***
    nbpart_rleft_iter(mu_id)   = 0
    nbpart_rright_iter(mu_id)  = 0
    nbpart_thleft_iter(mu_id)  = 0
    nbpart_thright_iter(mu_id) = 0
      
    ttmprhs = 0._RKIND
    tcomm2  = 0._RKIND
    tcomm3  = 0._RKIND
    tadvect = 0._RKIND
      
    call clck_time(tdeb)
    call clck_time(tfin)
    call clck_diff(tdeb,tfin,tcomm3)
      
    do ivpar = 0,Nvpar
      !*** cubic spline computation of the ***
      !***   distribution function         ***
      do zz = 0,Nphi-1+bloc_phi,bloc_phi
        call clck_time(tdeb)
        ! Derivative computation on previous bloc of 'z'
        call comm_copyrecvbuffers(geom%dr,geom%dtheta,zz,geom, &
          init_prof,init_magnet,init_curr,f,ivpar) 
        call clck_time(tfin)
        call clck_diff(tdeb,tfin,ttmprhs)
      
        if (zz .ne. Nphi) then
          ! Derivative computation on actual bloc of 'z' and send
          call clck_time(tdeb)
          call comm_irecv
          call comm_ssend
          call comm_wait_recv
          call clck_time(tfin)
          call clck_diff(tdeb,tfin,tcomm2)
        endif
      
        if (zz .ne. 0) then
          b_phi = mod(zz-bloc_phi,2*bloc_phi)
          call clck_time(tdeb)
      
#ifdef _OPENMP
!$OMP PARALLEL private(alpha_ij,beta_ij,tid,iphi,itheta,ir,rfeet, &
!$OMP tfeet,xfeet,yfeet,finterpol) reduction(max:maxalp_xy) default(shared)
!$OMP BARRIER
          tid = 1+omp_get_thread_num()
#else
          tid = 1
#endif
          alpha_ij => Romp1_0Nr_0Ntheta(tid)%val
          beta_ij  => Romp2_0Nr_0Ntheta(tid)%val
      
!$OMP DO SCHEDULE(static,1)
          do iphi = zz-bloc_phi,zz-1
            call hermite(f%hhspline2D_rtheta(tid),&
              tmprhs2(:,:,3*mod(iphi,2*bloc_phi)))
            call hermite(f%uxspline2D_rtheta(tid),&
              tmprhs2(:,:,3*mod(iphi,2*bloc_phi)+1))
            call hermite(f%uyspline2D_rtheta(tid),&
              tmprhs2(:,:,3*mod(iphi,2*bloc_phi)+2))
            !**************************************************
            !*** Computation of the 2D advection in the     *** 
            !***  cartesian coordinates                     ***
            !***   dx/dt = ux and  dy/dt = uy               ***
            !*** where :                                    ***
            !*** . ux = dx/dr*dr/dt + dx/dtheta*dtheta/dt   ***
            !*** . uy = dy/dr*dr/dt + dy/dtheta*dtheta/dt   ***
            !*** with dx/dr     = cos(theta) ;              ***
            !***      dx/dtheta = -r*sin(theta) and         ***
            !***      dy/dr     = sin(theta) ;              ***
            !***      dy/dtheta = r*cos(theta)              ***
            !*** and  . dr/dt     = vpar*bstar_gradr +      ***
            !***            vExB_gradr + vD_gradr           ***
            !***      . dtheta/dt = vpar*bstar_gradtheta +  ***
            !***            vExB_gradtheta + vD_gradtheta   ***
            !**************************************************
            call local_taylor(geom,dt, &
              f%uxspline2d_rtheta(tid)%coef, &
              f%uyspline2d_rtheta(tid)%coef, alpha_ij,beta_ij, &
              tmprhs2(:,:,3*mod(iphi,2*bloc_phi)+1), &
              tmprhs2(:,:,3*mod(iphi,2*bloc_phi)+2))
      
            !*** computation of the maximum displacement ***
            !***  normalized by the space step           ***     
            maxalp_xy = max(maxalp_xy,maxval(abs( &
              alpha_ij(f%istart_modif:f%iend_modif, &
              jstart:jend)))/geom%dr)
            maxalp_xy = max(maxalp_xy,maxval(abs( &
              beta_ij(f%istart_modif:f%iend_modif, &
              jstart:jend)))/geom%dtheta)
      
            do itheta = jstart,jend
              do ir = f%istart_modif,f%iend_modif
                !*** computes the feet characteristic       ***
                rfeet = geom%rg(ir)
                tfeet = geom%thetag(itheta)
                !*** transformation of (r,theta) into (x,y) ***
                call cyl2cart(rfeet,tfeet,xfeet,yfeet)
                !*** 2D advection                           ***
                xfeet = xfeet - alpha_ij(ir,itheta)
                yfeet = yfeet - beta_ij(ir,itheta)
                !*** transformation of (x,y) into (r,theta) ***
                call cart2cyl(xfeet,yfeet,rfeet,tfeet)    
                !*** Check that r and theta are on the      ***
                !***  subdomain treated                     ***
                call local_r_verif(f%hhspline2D_rtheta(tid),rfeet)
                if (Nbproc_theta.ne.1) then
                  call local_theta_verif(f%hhspline2D_rtheta(tid), &
                    geom,tfeet)  
                end if
                !*** f interpolation ***
                call interpol2d(f%hhspline2D_rtheta(tid), &
                  rfeet,tfeet,finterpol)
                f%values(ir,itheta,iphi,ivpar) = finterpol
              end do
            end do
          end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
          call clck_time(tfin)
          call clck_diff(tdeb,tfin,tadvect)
        end if
      end do
    end do
      
#ifdef TIMER    
    write(6,'(I4,A,F13.7,A,F13.7,A,F13.7,A,F13.7,A,F13.7)') &
      pglobal_id, " Temps advection2D,  init rhs ", &
      ttmprhs," comm2 ",tcomm2," comm3 ",tcomm3," advec2D ",tadvect
#endif
      
    !*** counter for particles out of domain ***
    if (plocal_id.eq.0) then
      iter_save                      = max(iter_run,0)          
      nbpart_rleft_muid(iter_save)   = nbpart_rleft_iter(mu_id)
      nbpart_rright_muid(iter_save)  = nbpart_rright_iter(mu_id)
      nbpart_thleft_muid(iter_save)  = nbpart_thleft_iter(mu_id)
      nbpart_thright_muid(iter_save) = nbpart_thright_iter(mu_id)
      if ((nbpart_rleft_iter(mu_id) + &
        nbpart_rright_iter(mu_id)).ne.0) then
        outofdom_percent = &
          ((nbpart_rright_iter(mu_id)+nbpart_rleft_iter(mu_id))/ &
          float(Nbpoints_muid))*100._RKIND
        write(6,'(A,I4,A,f6.2,A)') &
          '      -> WARNING : particles out ' // &
          'of domain in r for mu_id     (',mu_id,') = ', &
          outofdom_percent,' %'
      end if
      if ((nbpart_thleft_iter(mu_id) + &
        nbpart_thright_iter(mu_id)).ne.0) then
        outofdom_percent = &
          ((nbpart_thright_iter(mu_id)+nbpart_thleft_iter(mu_id))/ &
          float(Nbpoints_muid))*100._RKIND
        write(6,'(A,I4,A,f6.2,A)') &
          '      -> WARNING : particles out ' // &
          'of domain in theta for mu_id (',mu_id,') = ', &
          outofdom_percent,' %'
      end if
    end if
    call clck_time(eclock_advec2D)
    call clck_diff(bclock_advec2D,eclock_advec2D, &
      global_time_advec2D)
      
!R3 call r3_info_end (r3_info_index_0) !R3
  end subroutine local_advec2D_xy
      
  !************************************************************
  ! ADVECTIONS BY USING A TRANSPOSITION OF THE
  !  4D DISTRIBUTION FUNCTION AND THEN BY USING
  !  GLOBAL CUBIC SPLINES
  !
  !  (Rk : For the moment the global cubic splines correspond 
  !   to the local cubic splines on one unique patch)  
  !************************************************************
  !--------------------------------------------------------------
  ! - Solving by a Taylor algorithm of :
  !    . alpha = dt*ux(xi-alpha,yj-beta)
  !    . beta  = dt*uy(xi-alpha,yj-beta)
  !  with the cubic spline interpolation of the electric field
  !
  ! - this algorithm is used for characteristic 
  !    feet computation (in cartesian coordinates)
  ! - Solution: 
  !    alpha = dt*ux(xi,yj)- 
  !            0.5*dt*dt*{ (dux/dx)(xi,yj)*ux(xi,yj) 
  !                                +(dux/dy)(xi,yj)*uy(xi,yj))
  !    beta  = dt*uy(xi,yj) - 
  !            0.5*dt*dt*{ (duy/dx)(xi,yj)*ux(xi,yj) 
  !                                +(duy/dy)(xi,yj)*uy(xi,yj))
  ! -> developped by G. Latu and used in the case of
  !    4D transposition of the distribution function 
  !--------------------------------------------------------------
  subroutine transpose_taylor(geom,dt, &
    ux_scoef,uy_scoef,alpha_ij,beta_ij,rhsux,rhsuy)
    use globals, only : stencil, bufsize, Nr, Ntheta
    type(geometry)                 , intent(in)    :: geom
    real(RKIND)                    , intent(in)    :: dt
    real(RKIND), dimension(-2:,-2:), intent(in)    :: ux_scoef
    real(RKIND), dimension(-2:,-2:), intent(in)    :: uy_scoef
    real(RKIND), dimension(0:,0:)  , intent(inout) :: alpha_ij
    real(RKIND), dimension(0:,0:)  , intent(inout) :: beta_ij
    real(RKIND), &
      dimension(-bufsize:,-bufsize:), intent(in)   :: rhsux, rhsuy
      
    integer     :: ir, itheta  
    real(RKIND) :: rp, thetap, costhetap, sinthetap, invrp, dt2
    !* ux(rp,thetap) and uy(rp,thetap) *
    real(RKIND) :: uxg, uyg, gux, guy
    !* u derivates at (rp,thetap) *			
    real(RKIND), dimension(4) :: du
    real(RKIND) :: duxg_dr, duxg_dth, duxg_dx, duxg_dy
    real(RKIND) :: duyg_dr, duyg_dth, duyg_dx, duyg_dy
    integer     :: igstart, igend, igstart_buf, igend_buf, igsize
    integer     :: jgstart, jgend, jgstart_buf, jgend_buf, jgsize
      
    igstart     = 0
    igend       = Nr
    jgstart     = 0
    jgend       = Ntheta-1
    igsize      = igend - igstart + 1
    jgsize      = jgend - jgstart + 1
    igstart_buf = igstart - bufsize
    jgstart_buf = jgstart - bufsize
    igend_buf   = igend + bufsize
    jgend_buf   = jgend + bufsize
      
    do itheta = jgstart,jgend
      do ir = igstart,igend
        alpha_ij(ir,itheta) = 0._RKIND
        beta_ij(ir,itheta)  = 0._RKIND
      end do
    end do
    dt2 = dt*dt
    do itheta = jgstart,jgend
      sinthetap = geom%sin_theta(itheta)
      costhetap = geom%cos_theta(itheta)
      do ir = igstart,igend
        rp    = geom%rg(ir)
        invrp = (1._RKIND/rp)
        !*** computation of ux(rp,thetap) and uy(rp,thetap)  ***
        !***  by interpolation in the (r,theta) space ***
        call transpose_advecfield(ir,itheta,geom, &
          ux_scoef,uy_scoef,du)
        duxg_dr  = du(1)
        duxg_dth = du(2)
        duyg_dr  = du(3)
        duyg_dth = du(4)
        uxg      = rhsux(ir, itheta)
        uyg      = rhsuy(ir, itheta)
      
        !*** computation of dux(rp,thetap)/dx           ***
        !*** and dux(rp,thetap)/dy by :                 ***
        !***  . dux(rp,thetap)/dx = cos(thetap)*dux/dr  ***
        !***       -(1/r)*sin(thetap)*dux/dtheta        ***
        !***  . dux(rp,thetap)/dy = sin(thetap)*dux/dr  ***
        !***       +(1/r)*cos(thetap)*dux/dtheta        ***
        if (rp.ne.0._RKIND) then
          duxg_dx = costhetap*duxg_dr - &
            invrp*sinthetap*duxg_dth
          duxg_dy = sinthetap*duxg_dr + & 
            invrp*costhetap*duxg_dth
          duyg_dx = costhetap*duyg_dr - &
            invrp*sinthetap*duyg_dth
          duyg_dy = sinthetap*duyg_dr + &
            invrp*costhetap*duyg_dth 
        else
          duxg_dx = costhetap*duxg_dr
          duxg_dy = sinthetap*duxg_dr
          duyg_dx = costhetap*duyg_dr
          duyg_dy = sinthetap*duyg_dr
        endif
      
        !*** computation of the displacement ***
        alpha_ij(ir,itheta) = dt*uxg - &
          0.5_RKIND*dt2*(duxg_dx*uxg+duxg_dy*uyg)
        beta_ij(ir,itheta)  = dt*uyg - &
          0.5_RKIND*dt2*(duyg_dx*uxg+duyg_dy*uyg)
      enddo	! end do ir
    enddo	! end do itheta
  end subroutine transpose_taylor
      
  !-----------------------------------------------------
  ! Copy of the receive buffers 
  !-----------------------------------------------------
  subroutine transpose_fillbuffers(d_r,d_theta,zz, &
    geom,init_prof,init_magnet,init_curr,f,ivpar) 
    use globals
#ifdef _OPENMP
    use OMPutils_module, only : omp_get_thread_num
#endif
    use OMPutils_module, only : Romp1_0Nr, Romp2_0Nr, &
      Romp3_0Nr, Romp4_0Nr, Romp_vect_mpstencil, &
      Romp_fillde_112
    use MPIutils_module
    use local_spline_module
    use fdistribu5d_class
    use init_profile_class
    use init_magnetic_class
    use init_current_class
    use clock_module
    type(geometry)     , intent(in) :: geom
    real(RKIND)        , intent(in) :: d_theta, d_r
    integer            , intent(in) :: zz, ivpar
    type(init_profile) , intent(in) :: init_prof
    type(init_magnetic), intent(in) :: init_magnet
    type(init_current) , intent(in) :: init_curr
    type(fdistribu5d)  , intent(in) :: f
      
    ! -> variables for OpenMP parallelization
    integer     :: ir, itheta, iphi, i, j, k, zifi, bifi
    integer     :: igstart, igend, igstart_buf, igend_buf, igsize
    integer     :: jgstart, jgend, jgstart_buf, jgend_buf, jgsize
    integer     :: tid
    real(RKIND) :: NWcorner, NEcorner, SWcorner, SEcorner
    real(RKIND) :: derivtmp
    real(RKIND) :: rder_dr, lder_dr
    real(RKIND), dimension(:), pointer :: vect, fillde
      
    ! -> variables for 2D advection
    integer     :: q
    real(RKIND) :: ri, vparl
    real(RKIND) :: cos_thetaj, sin_thetaj
    real(RKIND) :: drdt_tmp, dthetadt_tmp
    real(RKIND) :: dxdr_tmp, dxdtheta_tmp
    real(RKIND) :: dydr_tmp, dydtheta_tmp
    real(RKIND) :: bstar_gradr_tmp, bstar_gradtheta_tmp
    real(RKIND) :: vExB_gradr_tmp, vExB_gradtheta_tmp
    real(RKIND) :: vD_gradr_tmp, vD_gradtheta_tmp
    real(RKIND) :: rfeet, tfeet, xfeet, yfeet
    real(RKIND) :: finterpol
      
    !-> components of ExB drift velocity  (in (r,theta,phi))
    real(RKIND), dimension(:), pointer :: vExB_gradr_1D
    real(RKIND), dimension(:), pointer :: vExB_gradtheta_1D
      
    !-> components of the curvature drift velocity
    real(RKIND), dimension(:), pointer :: vD_gradr_1D
    real(RKIND), dimension(:), pointer :: vD_gradtheta_1D
      
#ifdef _OPENMP
!$OMP PARALLEL private(tid,vect,fillde,iphi, &
!$OMP lder_dr,rder_dr,bifi,zifi,ir,i,k, &
!$OMP NWcorner,NEcorner,SWcorner,SEcorner, &
!$OMP j,itheta,derivtmp,q,finterpol,ri,vparl, &
!$OMP cos_thetaj,sin_thetaj, &
!$OMP drdt_tmp,dthetadt_tmp,dxdr_tmp,dxdtheta_tmp, &
!$OMP dydr_tmp,dydtheta_tmp, &
!$OMP bstar_gradr_tmp,bstar_gradtheta_tmp, &
!$OMP vExB_gradr_tmp,vExB_gradtheta_tmp, &
!$OMP vD_gradr_tmp,vD_gradtheta_tmp, &
!$OMP vExB_gradr_1D,vExB_gradtheta_1D, &
!$OMP vD_gradr_1D,vD_gradtheta_1D, &
!$OMP rfeet,tfeet,xfeet,yfeet,igstart_buf,igend_buf, &
!$OMP igstart,igend,jgstart_buf,jgend_buf,jgstart,jgend, &
!$OMP igsize,jgsize) default(shared) 
!$OMP BARRIER
    tid = 1+omp_get_thread_num()
#else
    tid = 1
#endif
    igstart     = 0
    igend       = Nr
    jgstart     = 0
    jgend       = Ntheta-1
    igsize      = igend - igstart + 1
    jgsize      = jgend - jgstart + 1
    igstart_buf = igstart - bufsize
    jgstart_buf = jgstart - bufsize
    igend_buf   = igend + bufsize
    jgend_buf   = jgend + bufsize
      
    vect   => Romp_vect_mpstencil(tid)%val
    fillde => Romp_fillde_112(tid)%val
 
    vExB_gradr_1D     => Romp1_0Nr(tid)%val
    vExB_gradtheta_1D => Romp2_0Nr(tid)%val
    vD_gradr_1D       => Romp3_0Nr(tid)%val
    vD_gradtheta_1D   => Romp4_0Nr(tid)%val
      
    vparl = geom%vparg(ivpar)
!$OMP DO SCHEDULE(static)
    do iphi = zz, zz+bloc_phi-1
      do k = jgstart_buf,jgend_buf
        !**************************************************
        !*** Computation of the 2D advection in the     *** 
        !***  cartesian coordinates                     ***
        !***   dx/dt = ux and  dy/dt = uy               ***
        !*** where :                                    ***
        !*** . ux = dx/dr*dr/dt + dx/dtheta*dtheta/dt   ***
        !*** . uy = dy/dr*dr/dt + dy/dtheta*dtheta/dt   ***
        !*** with dx/dr     = cos(theta) ;              ***
        !***      dx/dtheta = -r*sin(theta) and         ***
        !***      dy/dr     = sin(theta) ;              ***
        !***      dy/dtheta = r*cos(theta)              ***
        !*** and  . dr/dt     = vpar*bstar_gradr +      ***
        !***                    vExB_gradr + vD_gradr   ***
        !***      . dtheta/dt = vpar*bstar_gradtheta +  ***
        !***            vExB_gradtheta + vD_gradtheta   ***
        !**************************************************
        itheta     = mod(k+8*Ntheta,Ntheta)
        cos_thetaj = geom%cos_theta(itheta)
        sin_thetaj = geom%sin_theta(itheta)
        !-> compute vExB.gradr and vExB.gradtheta
        call compute_vExBgradr(geom,init_magnet,init_curr, &
          dJ0Phidtheta,dJ0Phidphi,0,Nr,itheta,itheta, &
          iphi,iphi,ivpar,ivpar,vExB_gradr_1D)
        call compute_vExBgradtheta(geom,init_magnet,init_curr, &
          dJ0Phidr,dJ0Phidphi,0,Nr,itheta,itheta, &
          iphi,iphi,ivpar,ivpar,vExB_gradtheta_1D)
        !-> compute vD.gradr and vD.gradtheta
        call compute_vDgradr(geom,init_magnet,init_curr, &
          0,Nr,itheta,itheta,ivpar,ivpar,vD_gradr_1D)
        call compute_vDgradtheta(geom,init_magnet,init_curr, &
          0,Nr,itheta,itheta,ivpar,ivpar,vD_gradtheta_1D)
        do q = igstart_buf,igend_buf
          if (q.lt.0) then
            ir = 0
          else if (q.gt.Nr) then
            ir = Nr
          else
            ir = q
          endif
          !-> compute bstar.gradr and bstar.gradtheta
#ifdef NOPRECOMPUTE 
            call compute_bstar_contravariant(geom, &
              init_magnet,init_curr,ir,itheta,ivpar, &
              Sbstar_gradx1=bstar_gradr_tmp, &
              Sbstar_gradx2=bstar_gradtheta_tmp)
#else
            call precomputed_bstar_gradr(ir,itheta,ivpar, &
              bstar_gradr_tmp)
            call precomputed_bstar_gradtheta(ir,itheta,ivpar, &
              bstar_gradtheta_tmp)
#endif
          !-> initialisation of the velocity projections
          vExB_gradr_tmp       = vExB_gradr_1D(ir)
          vExB_gradtheta_tmp   = vExB_gradtheta_1D(ir)
          vD_gradr_tmp         = vD_gradr_1D(ir)
          vD_gradtheta_tmp     = vD_gradtheta_1D(ir)
          !-> dx/dr, dx/dtheta, dy/dr and dy/dtheta
          ri           = geom%rg(ir)
          dxdr_tmp     = cos_thetaj
          dxdtheta_tmp = -ri*sin_thetaj
          dydr_tmp     = sin_thetaj
          dydtheta_tmp = ri*cos_thetaj
          !-> dr/dt and dtheta/dt
          drdt_tmp     = vparl*bstar_gradr_tmp + &
            vExB_gradr_tmp + vD_gradr_tmp 
          dthetadt_tmp = vparl*bstar_gradtheta_tmp + &
            vExB_gradtheta_tmp + vD_gradtheta_tmp
          !-> computation of ux and uy
          bifi                = 3*mod(iphi,bloc_phi)
          tmprhs2(q,k,bifi)   = f4D_transp(ir,itheta,iphi,ivpar)
          tmprhs2(q,k,bifi+1) = dxdr_tmp*drdt_tmp + &
              dxdtheta_tmp*dthetadt_tmp
          tmprhs2(q,k,bifi+2) = dydr_tmp*drdt_tmp + &
              dydtheta_tmp*dthetadt_tmp
        end do
      end do
    enddo
!$OMP END DO
      
!$OMP DO SCHEDULE(static) 
    do iphi = 0,3*bloc_phi-1
      bifi       = iphi
      zifi       = bifi
      fillde(1)  = rightderive10( &
        tmprhs2(igstart+1:igstart+stencil,jgstart,zifi),d_r) 
      fillde(2)  = rightderive10( &
        tmprhs2(igstart,jgstart+1:jgstart+stencil,zifi), &
        d_theta)
      fillde(3)  = rightderive10( &
        tmprhs2(igstart+1,jgstart+1:jgstart+stencil,zifi), &
        d_theta)
      fillde(4)  = rightderive10( &
        tmprhs2(igend,jgstart+1:jgstart+stencil,zifi),d_theta)
      fillde(5)  = leftderive10( &
        tmprhs2(igend+1-stencil:igend,jgstart,zifi),d_r) 
      fillde(6)  = rightderive10( &
        tmprhs2(igstart+1:igstart+stencil,jgstart+1,zifi),d_r) 
      fillde(7)  = leftderive10( &
        tmprhs2(igend+1-stencil:igend,jgstart+1,zifi),d_r) 
      fillde(8)  = rightderive10( &
        tmprhs2(igstart+1:igstart+stencil,jgend,zifi),d_r) 
      fillde(9)  = leftderive10( &
        tmprhs2(igend+1-stencil:igend,jgend,zifi),d_r) 
      fillde(10) = leftderive10( &
        tmprhs2(igstart,jgend+1-stencil:jgend,zifi),d_theta)
      fillde(11) = leftderive10( &
        tmprhs2(igstart+1,jgend+1-stencil:jgend,zifi),d_theta)
      fillde(12) = leftderive10( &
        tmprhs2(igend,jgend+1-stencil:jgend,zifi),d_theta)
      vect       = 0
      do k = 1,stencil
        vect(k) = rightderive10( &
          tmprhs2(igstart+1:igstart+stencil,jgstart+k,zifi),d_r)
      enddo
      NWcorner = rightderive10(vect(1:stencil),d_theta)
      vect     = 0
      do k = -stencil,-1
        vect(k) = rightderive10( &
          tmprhs2(igstart+1:igstart+stencil,jgend+1+k,zifi),d_r)
      enddo
      NEcorner = leftderive10(vect(-stencil:-1),d_theta)
      vect     = 0
      do k = 1,stencil
        vect(k) = leftderive10( &
          tmprhs2(igend+1-stencil:igend,jgstart+k,zifi),d_r)
      enddo
      SWcorner = rightderive10(vect(1:stencil),d_theta)
      vect     = 0
      do k = -stencil,-1
        vect(k) = leftderive10( &
          tmprhs2(igend+1-stencil:igend,jgend+1+k,zifi),d_r)
      enddo
      SEcorner = leftderive10(vect(-stencil:-1),d_theta)
      
      !--------------------------------------------------------
      
      do ir = igstart,igend
        i = ir - igstart
        sbufWW(i,bifi)          = tmprhs2(ir,jgstart,zifi)
        sbufWW(igsize+i,bifi)   = rightderive10( &
          tmprhs2(ir,jgstart+1:jgstart+stencil,zifi),d_theta)
        sbufWW(2*igsize+i,bifi) = tmprhs2(ir,jgstart+1,zifi) ! DYS  
        sbufEE(i,bifi)          = tmprhs2(ir,jgend,zifi)
        sbufEE(igsize+i,bifi)   = leftderive10( &
          tmprhs2(ir,jgend+1-stencil:jgend,zifi),d_theta)
      enddo
      sbufWW(3*igsize+0,bifi) = fillde(1) !$1
      sbufWW(3*igsize+1,bifi) = fillde(6) !$6
      sbufWW(3*igsize+2,bifi) = fillde(5) !$5
      sbufWW(3*igsize+3,bifi) = fillde(7) !$7
      sbufWW(3*igsize+4,bifi) = NWcorner  !*5
      sbufWW(3*igsize+5,bifi) = SWcorner  !*11
      
      sbufEE(2*igsize+0,bifi) = fillde(8) !$8
      sbufEE(2*igsize+1,bifi) = fillde(9) !$9
      sbufEE(2*igsize+2,bifi) = NEcorner  !*2
      sbufEE(2*igsize+3,bifi) = SEcorner  !*8
      
      !--------------------------------------------------------
      
      do itheta = jgstart,jgend
        j = itheta - jgstart
        sbufNN(j,bifi)          = tmprhs2(igstart,itheta,zifi)
        sbufNN(jgsize+j,bifi)   = rightderive10( &
          tmprhs2(igstart+1:igstart+stencil,itheta,zifi),d_r) 
        sbufNN(2*jgsize+j,bifi) = &
          tmprhs2(igstart+1,itheta,zifi) ! DYS
      
        sbufSS(j,bifi)          = tmprhs2(igend,itheta,zifi)
        sbufSS(jgsize+j,bifi)   = leftderive10( &
          tmprhs2(igend+1-stencil:igend,itheta,zifi),d_r)
      enddo
      
      sbufNN(3*jgsize+0,bifi) = fillde(2)  !$2
      sbufNN(3*jgsize+1,bifi) = fillde(3)  !$3
      sbufNN(3*jgsize+2,bifi) = fillde(10) !$10
      sbufNN(3*jgsize+3,bifi) = fillde(11) !$11
      sbufNN(3*jgsize+4,bifi) = NWcorner   !*7
      sbufNN(3*jgsize+5,bifi) = NEcorner   !*10
      
      sbufSS(2*jgsize+0,bifi) = fillde(12) !$12
      sbufSS(2*jgsize+1,bifi) = fillde(4)  !$4
      sbufSS(2*jgsize+2,bifi) = SWcorner   !*1
      sbufSS(2*jgsize+3,bifi) = SEcorner   !*4
      
      !----------
      
      sbufNW(0,bifi) = tmprhs2(igstart,jgstart,zifi)
      sbufNW(1,bifi) = tmprhs2(igstart+1,jgstart,zifi)   ! DYS
      sbufNW(2,bifi) = tmprhs2(igstart,jgstart+1,zifi)   ! DYS
      sbufNW(3,bifi) = tmprhs2(igstart+1,jgstart+1,zifi) ! DYS
      sbufNW(4,bifi) = NWcorner  !*12
      sbufNW(5,bifi) = fillde(1) !$1
      sbufNW(6,bifi) = fillde(2) !$2
      sbufNW(7,bifi) = fillde(3) !$3
      sbufNW(8,bifi) = fillde(6) !$6
      
      sbufNE(0,bifi) = tmprhs2(igstart,jgend,zifi)
      sbufNE(1,bifi) = tmprhs2(igstart+1,jgend,zifi) ! DYS
      sbufNE(2,bifi) = NEcorner   !*9
      sbufNE(3,bifi) = fillde(8)  !$8
      sbufNE(4,bifi) = fillde(10) !$10
      sbufNE(5,bifi) = fillde(11) !$11
      
      sbufSW(0,bifi) = tmprhs2(igend,jgstart,zifi)
      sbufSW(1,bifi) = tmprhs2(igend,jgstart+1,zifi) ! DYS
      sbufSW(2,bifi) = SWcorner  !*6
      sbufSW(3,bifi) = fillde(4) !$4
      sbufSW(4,bifi) = fillde(5) !$5
      sbufSW(5,bifi) = fillde(7) !$7
      
      sbufSE(0,bifi) = tmprhs2(igend,jgend,zifi)
      sbufSE(1,bifi) = SEcorner   !*3
      sbufSE(2,bifi) = fillde(9)  !$9
      sbufSE(3,bifi) = fillde(12) !$12
    enddo
!$OMP END DO
      
!$OMP BARRIER
!$OMP MASTER
    rbufEE(:,:) = sbufWW(:,:)
    rbufWW(:,:) = sbufEE(:,:)
    rbufSE(:,:) = sbufNW(:,:)
    rbufSW(:,:) = sbufNE(:,:)
    rbufNW(:,:) = sbufSE(:,:)
    rbufNE(:,:) = sbufSW(:,:)
    rbufNN(:,:) = sbufSS(:,:)
    rbufSS(:,:) = sbufNN(:,:)
!$OMP END MASTER
!$OMP BARRIER
      
      lder_dr = leftconst / d_r
      rder_dr = rightconst / d_r      
!$OMP DO SCHEDULE(static) 
      do iphi = 0,3*bloc_phi-1
        bifi = iphi
        zifi = bifi
        !*** f values ***
        do ir = igstart,igend
          i = ir - igstart
          tmprhs2(ir,jgstart-2,zifi) = rbufWW(i,bifi) 
          tmprhs2(ir,jgend+1,zifi)   = rbufEE(i,bifi)  
          tmprhs2(ir,jgend+3,zifi)   = rbufEE(2*igsize+i,bifi)  
          tmprhs2(ir,jgend+2,zifi)   = rbufEE(igsize+i,bifi) + &
            leftderive10( &
            tmprhs2(ir,jgend+1-stencil:jgend,zifi),d_theta) 
          tmprhs2(ir,jgstart-1,zifi) = rbufWW(igsize+i,bifi) + &
            rightderive10( &
            tmprhs2(ir,jgstart+1:jgstart+stencil,zifi),d_theta)
        end do
        !*** partial difference in r and theta of f ***
        !***   in the corners of f                  ***
        do k= 1,stencil
          vect(k) = rightderive10( &
            tmprhs2(igstart+1:igstart+stencil,jgstart+k,zifi),d_r)
        enddo
        NWcorner = rightderive10(vect(1:stencil),d_theta)
        do k = -stencil,-1
          vect(k) = rightderive10( &
            tmprhs2(igstart+1:igstart+stencil,jgend+1+k,zifi),d_r)
        enddo
        NEcorner = leftderive10(vect(-stencil:-1),d_theta)
        do k = 1,stencil
          vect(k) = leftderive10( &
            tmprhs2(igend+1-stencil:igend,jgstart+k,zifi),d_r)
        enddo
        SWcorner = rightderive10(vect(1:stencil),d_theta)
        do k = -stencil,-1
          vect(k) = leftderive10( &
            tmprhs2(igend+1-stencil:igend,jgend+1+k,zifi),d_r)
        enddo
        SEcorner = leftderive10(vect(-stencil:-1),d_theta)
      
        !-> compact domain in r direction at Nr , 
        !->  f is constant outside the domain
        do itheta = jgstart,jgend
          j = itheta - jgstart
          tmprhs2(igend+2,itheta,zifi) = &
            rder_dr*tmprhs2(igend,itheta,zifi) + leftderive10( &
            tmprhs2(igend+1-stencil:igend,itheta,zifi),d_r)
          tmprhs2(igend+1,itheta,zifi) = tmprhs2(igend,itheta,zifi)
          tmprhs2(igend+3,itheta,zifi) = tmprhs2(igend,itheta,zifi)
        end do
        !-> compact domain in r direction, f is constant
        !->   outside the domain
        tmprhs2(igend+1,jgstart-2,zifi) = &
          tmprhs2(igend,jgstart-2,zifi)
        tmprhs2(igend+3,jgstart-2,zifi) = &
          tmprhs2(igend,jgstart-2,zifi)
        tmprhs2(igend+1,jgend+1,zifi)   = &
          tmprhs2(igend,jgend+1,zifi)
        tmprhs2(igend+3,jgend+1,zifi)   = &
          tmprhs2(igend,jgend+1,zifi)
        tmprhs2(igend+1,jgend+3,zifi)   = &
          tmprhs2(igend,jgend+3,zifi)
        tmprhs2(igend+3,jgend+3,zifi)   = &
          tmprhs2(igend,jgend+3,zifi)
        !-> compact domain in r direction, f is constant 
        !->  outside the domain
        tmprhs2(igend+1,jgstart-1,zifi) = &
          tmprhs2(igend,jgstart-1,zifi)
        tmprhs2(igend+3,jgstart-1,zifi) = &
          tmprhs2(igend,jgstart-1,zifi)
        tmprhs2(igend+2,jgstart-2,zifi) = &
          rbufWW(2*igsize+1,bifi) + &
          rder_dr * tmprhs2(igend,jgstart-2,zifi)!$9 8
        tmprhs2(igend+2,jgend+1,zifi)   = &
          tmprhs2(igend,jgend+1,zifi)*rder_dr + &
          rbufEE(3*igsize+2,bifi) !$1 5
        tmprhs2(igend+1,jgend+2,zifi)   = &
          tmprhs2(igend,jgend+2,zifi)
        tmprhs2(igend+3,jgend+2,zifi)   = &
          tmprhs2(igend,jgend+2,zifi)
        tmprhs2(igend+2,jgend+3,zifi)   = &
          rbufEE(3*igsize+3,bifi) + &
          rder_dr * tmprhs2(igend,jgend+3,zifi) !$7 6
        !-> compact domain in r direction, f is constant 
        !->   outside the domain
        derivtmp = rder_dr * tmprhs2(igend,jgstart-1,zifi)
        tmprhs2(igend+2,jgstart-1,zifi) = SWcorner + derivtmp &
          + rbufWW(2*igsize+3,bifi) 
        derivtmp = rder_dr * tmprhs2(igend,jgend+2,zifi)
        tmprhs2(igend+2,jgend+2,zifi) = SEcorner + derivtmp &
          + rbufEE(3*igsize+5,bifi) 
      
        !-> compact domain in r direction at 0, f is constant 
        !->   outside the domain
        do itheta = jgstart,jgend
          j = itheta - jgstart
          tmprhs2(igstart-1,itheta,zifi) = &
            lder_dr*tmprhs2(igstart,itheta,zifi) + rightderive10( &
            tmprhs2(igstart+1:igstart+stencil,itheta ,zifi),d_r) 
          tmprhs2(igstart-2,itheta,zifi) = &
            tmprhs2(igstart,itheta,zifi)
        end do
        !-> compact domain in r direction, f is constant 
        !->   outside the domain
        tmprhs2(igstart-2,jgstart-2,zifi) = &
          tmprhs2(igstart,jgstart-2,zifi)
        tmprhs2(igstart-2,jgend+1,zifi)   = &
          tmprhs2(igstart,jgend+1,zifi)
        tmprhs2(igstart-2,jgend+3,zifi)   = &
          tmprhs2(igstart,jgend+3,zifi)
        !-> compact domain in r direction, f is constant 
        !->   outside the domain
        tmprhs2(igstart-1,jgstart-2,zifi) = &
          lder_dr*tmprhs2(igstart,jgstart-2,zifi) + &
          rbufWW(2*igsize+0,bifi) !$9 8      
        tmprhs2(igstart-2,jgstart-1,zifi) = &
          tmprhs2(igstart,jgstart-1,zifi)
        tmprhs2(igstart-2,jgend+2,zifi)   = &
          tmprhs2(igstart,jgend+2,zifi)
        tmprhs2(igstart-1,jgend+1,zifi)   = &
          rbufEE(3*igsize+0,bifi) + &
          lder_dr * tmprhs2(igstart,jgend+1,zifi) !$1 5
        tmprhs2(igstart-1,jgend+3,zifi)   = &
          lder_dr * tmprhs2(igstart,jgend+3,zifi) + &
          rbufEE(3*igsize+1,bifi) !$ 7 6
        !-> compact domain in r direction, f is constant 
        !->   outside the domain
        derivtmp = lder_dr * tmprhs2(igstart,jgstart-1,zifi)
        tmprhs2(igstart-1,jgstart-1,zifi) = NWcorner + derivtmp &
          + rbufWW(2*igsize+2,bifi) 
        derivtmp = lder_dr * tmprhs2(igstart,jgend+2,zifi) 
        tmprhs2(igstart-1,jgend+2,zifi) = NEcorner + derivtmp &
          + rbufEE(3*igsize+4,bifi) 
      enddo
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
  end subroutine transpose_fillbuffers
      
  !-----------------------------------------------------------------
  !   computation of the derivates of u
  !     by interpolation by cubic splines of u for a knot of 
  !     the mesh (r(ir),theta(itheta))
  !-----------------------------------------------------------------
  subroutine transpose_advecfield(ir,itheta,geom, &
    ux_scoef,uy_scoef,du)
    integer                        , intent(in)  :: ir
    integer                        , intent(in)  :: itheta
    type(geometry)                 , intent(in)  :: geom
    real(RKIND), dimension(-2:,-2:), intent(in)  :: ux_scoef
    real(RKIND), dimension(-2:,-2:), intent(in)  :: uy_scoef
    real(RKIND), dimension(4)      , intent(out) :: du
      
    real(RKIND)                  :: rp, thetap
    real(RKIND)                  :: duxg_dr
    real(RKIND)                  :: duxg_dth
    real(RKIND)                  :: duyg_dr
    real(RKIND)                  :: duyg_dth
    real(RKIND), dimension(-1:2) :: srbase, stbase
    real(RKIND), dimension(-1:2) :: srbase_prime, stbase_prime
    integer                      :: i1, i2
      
    rp     = geom%rg(ir)
    thetap = geom%thetag(itheta)
      
    !*** Calculate dux_dr(rc,thetac), dux_dth(rc,thetac)  ***
    !***   duy_dr(rc,thetac) and duy_dth(rc,thetac) ***
    srbase_prime(-1) = -1._RKIND/(2*geom%dr) 
    srbase_prime( 0) = 0._RKIND                
    srbase_prime( 1) = 1._RKIND/(2*geom%dr)
    srbase(-1)       = 1._RKIND/6._RKIND
    srbase( 0)       = 4._RKIND/6._RKIND
    srbase( 1)       = 1._RKIND/6._RKIND
    stbase_prime(-1) = -1._RKIND/(2*geom%dtheta)
    stbase_prime( 0) = 0._RKIND
    stbase_prime( 1) = 1._RKIND/(2*geom%dtheta)
    stbase(-1)       = 1._RKIND/6._RKIND
    stbase( 0)       = 4._RKIND/6._RKIND
    stbase( 1)       = 1._RKIND/6._RKIND
      
    duxg_dr  = 0._RKIND
    duxg_dth = 0._RKIND
    duyg_dr  = 0._RKIND
    duyg_dth = 0._RKIND
    do i2 = -1,1
      do i1 = -1,1
        duxg_dr = duxg_dr + ux_scoef(ir+i1,itheta+i2)*  & 
          srbase_prime(i1)*stbase(i2)
        duxg_dth = duxg_dth + ux_scoef(ir+i1,itheta+i2)*  & 
          srbase(i1)*stbase_prime(i2)
        duyg_dr = duyg_dr + uy_scoef(ir+i1,itheta+i2)*  & 
          srbase_prime(i1)*stbase(i2)
        duyg_dth = duyg_dth + uy_scoef(ir+i1,itheta+i2)*  & 
          srbase(i1)*stbase_prime(i2)
      end do
    end do
    du(1) = duxg_dr
    du(2) = duxg_dth
    du(3) = duyg_dr
    du(4) = duyg_dth
  end subroutine transpose_advecfield
      
  !-----------------------------------------------------
  ! 2D advection solving (in cartesian coordinates)
  !  - compute the characteristic feet and 
  !  - interpole the distribution function at this point
  ! (using the explicit taylor development)
  ! -> developped by G. Latu and used in the case of
  !    4D transposition of the distribution function 
  !----------------------------------------------------- 
  subroutine transpose_advec2D_xy(geom,dt, &
    nbsubit,init_prof,init_magnet,init_curr, &
    f,Sfmu_eq,activatefilter,maxalp_xy)
    use globals
    use MPIutils_module
#ifdef _OPENMP
    use OMPutils_module, only : omp_get_thread_num
#endif
    use OMPutils_module, only : Romp1_0Nr_0Ntheta, &
      Romp2_0Nr_0Ntheta
    use fdistribu5d_class
    use init_profile_class
    use init_magnetic_class
    use utils_module
    use clock_module
    use filter_module
!R3 #include "r3_info.h" !R3
    type(geometry)       , intent(in)    :: geom
    real(RKIND)          , intent(in)    :: dt
    integer              , intent(in)    :: nbsubit
    type(init_profile)   , intent(in)    :: init_prof
    type(init_magnetic)  , intent(in)    :: init_magnet
    type(init_current)   , intent(in)    :: init_curr
    type(fdistribu5d)    , intent(inout) :: f
    real(RKIND), &       
      dimension(:,:,:)   , pointer       :: Sfmu_eq
    logical              , intent(in)    :: activatefilter
    real(RKIND)          , intent(out)   :: maxalp_xy
      
    integer     :: ir, itheta, iphi, ivpar, zz, k, q
    integer     :: iter_save
    real(RKIND) :: rfeet, tfeet, xfeet, yfeet
    real(RKIND) :: rmin, rmax
    real(RKIND) :: finterpol, outofdom_percent
      
    real(RKIND), dimension(:,:), pointer :: alpha_ij
    real(RKIND), dimension(:,:), pointer :: beta_ij
      
    integer           :: tid
    integer           :: subit
    integer(TIMEPREC) :: tdeb, tfin
    real(RKIND)       :: tcomm2, tcomm3, ttmprhs, tadvect
      
!R3 call r3_info_begin (r3_info_index_0, 'CPU_transp_advec2d') !R3
      
    call clck_time(bclock_advec2D)
      
    ttmprhs = 0._RKIND
    tcomm2  = 0._RKIND
    tcomm3  = 0._RKIND
    tadvect = 0._RKIND
      
    call clck_time(tdeb)
    !*** Transpose the 4D function from f%values to f4D_transp ***
    call pptranspose_forward(f%values)   
    call clck_time(tfin)
    call clck_diff(tdeb,tfin,tcomm2)
      
    maxalp_xy  = 0._RKIND
    rmin       = geom%rg(0)
    rmax       = geom%rg(geom%Nr)
      
    do ivpar = lstart,lend
      do zz = kstart,kend,bloc_phi
        do subit = 1, nbsubit
          call clck_time(tdeb)
          !*** Derivative computation on actual bloc ***
          !***   of 'zz' (phi variable)              ***
          call transpose_fillbuffers(geom%dr,geom%dtheta,zz,geom, &
            init_prof,init_magnet,init_curr,f,ivpar) 
          call clck_time(tfin)
          call clck_diff(tdeb,tfin,ttmprhs)
      
          call clck_time(tdeb)
#ifdef _OPENMP
!$OMP PARALLEL private(alpha_ij,beta_ij,tid,iphi,itheta,ir,rfeet, &
!$OMP tfeet,xfeet,yfeet,finterpol) default(shared)
!$OMP BARRIER
          tid = 1+omp_get_thread_num()
#else
          tid = 1
#endif
          alpha_ij => Romp1_0Nr_0Ntheta(tid)%val
          beta_ij  => Romp2_0Nr_0Ntheta(tid)%val
          
!$OMP DO SCHEDULE(static,1)
          do iphi = zz, zz+bloc_phi-1
            if (iphi .gt. kend) then
              print *, pglobal_id, "erreur ", zz, iphi, kend
              STOP 
            end if
            call hermite(f%hhspline2D_rtheta(tid),&
              tmprhs2(:,:,3*mod(iphi,bloc_phi)))
            call hermite(f%uxspline2D_rtheta(tid),&
              tmprhs2(:,:,3*mod(iphi,bloc_phi)+1))
            call hermite(f%uyspline2D_rtheta(tid),&
              tmprhs2(:,:,3*mod(iphi,bloc_phi)+2))
      
            !*** Computation of the 2D advection in the     *** 
            !***  cartesian coordinates                     ***
            !***   dx/dt = ux and  dy/dt = uy               ***
            !*** where :                                    ***
            !*** . ux = dx/dr*dr/dt + dx/dtheta*dtheta/dt   ***
            !*** . uy = dy/dr*dr/dt + dy/dtheta*dtheta/dt   ***
            !*** with dx/dr     = cos(theta) ;              ***
            !***      dx/dtheta = -r*sin(theta) and         ***
            !***      dy/dr     = sin(theta) ;              ***
            !***      dy/dtheta = r*cos(theta)              ***
            !*** and  . dr/dt     = vpar*bstar_gradr +      ***
            !***            vExB_gradr + vD_gradr           ***
            !***      . dtheta/dt = vpar*bstar_gradtheta +  ***
            !***            vExB_gradtheta + vD_gradtheta   ***
            call transpose_taylor(geom,dt, &
              f%uxspline2d_rtheta(tid)%coef, &
              f%uyspline2d_rtheta(tid)%coef, &
              alpha_ij,beta_ij, &
              tmprhs2(:,:,3*mod(iphi,bloc_phi)+1), &
              tmprhs2(:,:,3*mod(iphi,bloc_phi)+2))
      
            !*** computation of the maximum displacement ***
            !***  normalized by the space step           *** 
            maxalp_xy = max(maxalp_xy,maxval(abs( &
              alpha_ij(f%rstart_modif: &
              f%rend_modif,0:Ntheta-1)))/geom%dr)
            maxalp_xy = max(maxalp_xy,maxval(abs( &
              beta_ij(f%rstart_modif: &
              f%rend_modif,0:Ntheta-1)))/geom%dtheta)
      
            !*** computation of the maximum displacement ***
            !***  normalized by the space step           ***
            do itheta = 0, Ntheta-1
              do ir = f%rstart_modif, f%rend_modif
                !*** computes the feet characteristic       ***
                rfeet = geom%rg(ir)
                tfeet = geom%thetag(itheta)
                !*** transformation of (r,theta) into (x,y) ***
                call cyl2cart(rfeet,tfeet,xfeet,yfeet)
                !*** 2D advection                           ***
                xfeet = xfeet - alpha_ij(ir,itheta)
                yfeet = yfeet - beta_ij(ir,itheta)
                !*** transformation of (x,y) into (r,theta) ***
                call cart2cyl(xfeet,yfeet,rfeet,tfeet)    
                !*** Check that r and theta are on the      ***
                !***   subdomain treated                    ***
                call r_verif(rmin,rmax,rfeet)
                !*** f interpolation ***
                call interpol2d(f%hhspline2D_rtheta(tid), &
                  rfeet,tfeet,finterpol)
                f4D_transp(ir,itheta,iphi,ivpar) = finterpol
              end do
            end do
          end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
          call clck_time(tfin)
          call clck_diff(tdeb,tfin,tadvect)
        end do
      end do
    end do
      
    if (activatefilter) call hffilter_loc(Sfmu_eq)
      
    call clck_time(tdeb)
    !*** Transpose the 4D function from f4D_transp to f%values ***
    call pptranspose_backward(f%values)   
    call clck_time(tfin)
    call clck_diff(tdeb,tfin,tcomm2)
      
#ifdef TIMER    
    write(6,'(I4,A,F13.7,A,F13.7,A,F13.7,A,F13.7,A,F13.7)') &
      pglobal_id, " Temps advection2D,  init rhs ", &
      ttmprhs," comm2 ",tcomm2," comm3 ",tcomm3," advec2D ",tadvect
#endif
      
    !*** counter for particles out of domain ***
    if (plocal_id.eq.0) then
      iter_save                      = max(iter_run,0)
      nbpart_rleft_muid(iter_save)   = nbpart_rleft_iter(mu_id)
      nbpart_rright_muid(iter_save)  = nbpart_rright_iter(mu_id)
      if ((nbpart_rleft_iter(mu_id) + &
        nbpart_rright_iter(mu_id)).ne.0) then
        outofdom_percent = &
          ((nbpart_rright_iter(mu_id)+nbpart_rleft_iter(mu_id))/ &
          float(Nbpoints_muid))*100._RKIND
        write(6,'(A,I4,A,f6.2,A)') &
          '      -> WARNING : particles out ' // &
          'of domain in r for mu_id     (',mu_id,') = ', &
          outofdom_percent,' %'
      end if
    end if
    call clck_time(eclock_advec2D)
    call clck_diff(bclock_advec2D,eclock_advec2D, &
      global_time_advec2D)
      
!R3 call r3_info_end (r3_info_index_0) !R3
  end subroutine transpose_advec2D_xy
      
  !************************************************************
  ! USED FOR PRINTING THE NUMBER OF PARTICLES
  !  OUTSIDE THE DOMAIN
  !************************************************************
  !--------------------------------------------------------
  !  Create the file to list the iterations where
  !   particles go out of domain in r or theta directions
  !--------------------------------------------------------
  subroutine print_particles_outofdomain()
    use globals, only : plocal_id, mu_id, Nbpoints_muid, &
      nbpart_rleft_muid, nbpart_rright_muid, &
      nbpart_thleft_muid, nbpart_thright_muid, &
      nbiter, iter_glob, deltat
      
    character(LEN=22) :: outofdom_filename
    integer           :: iter, iter_glob_num, ifil, Nbpart_muid
    integer           :: max_rright, max_rleft
    integer           :: max_thright, max_thleft, max_tot
    logical           :: outofdomain = .false.
    real(RKIND)       :: outofdom_percent
      
    ifil = 20+mu_id
    write(outofdom_filename,'("outofdomain_mu",i4.4,".txt")') mu_id
    open(ifil,file=outofdom_filename, &
      position='append',form='formatted')
    write(ifil,'(/,A,I4,A)') &
      '---------------- mu_id = ',mu_id,' ----------------'
    do iter = 0,nbiter
      iter_glob_num = iter_glob-nbiter+iter
      max_rright    = max(nbpart_rright_muid(iter),0)
      max_rleft     = max(nbpart_rleft_muid(iter),0)
      max_thright   = max(nbpart_thright_muid(iter),0)
      max_thleft    = max(nbpart_thleft_muid(iter),0)
      max_tot       = max_rright+max_rleft+max_thright+max_thleft
      if (max_tot.ne.0) then
        outofdomain = .true.
        write(ifil,'(/,A,I8,A,f6.2,A)') &
          ' *** iter ',iter_glob_num, &
          ' -> time = ',iter_glob_num*deltat,' ***'
        if (max_rright.ne.0) then
          outofdom_percent = (nbpart_rright_muid(iter)/ &
            float(Nbpoints_muid))*100._RKIND
          write(ifil,'(A,I8,A,f6.2,A)') &
            '   -> number of particles out of domain ' // &
            'at right in r:     ', nbpart_rright_muid(iter), &
            ' so (', outofdom_percent ,' % )' 
        end if
        if (max_rleft.ne.0) then
          outofdom_percent = (nbpart_rleft_muid(iter)/ &
            float(Nbpoints_muid))*100._RKIND
          write(ifil,'(A,I8,A,f6.2,A)') &
            '   -> number of particles out of domain ' // &
            'at left in r:      ',nbpart_rleft_muid(iter), &
            ' so (', outofdom_percent ,' % )'
        end if
        if (max_thright.ne.0) then
          outofdom_percent = (nbpart_thright_muid(iter)/ &
            float(Nbpoints_muid))*100._RKIND
          write(ifil,'(A,I8,A,f6.2,A)') &
            '   -> number of particles out of domain ' // &
            'at right in theta: ', nbpart_thright_muid(iter), &
            ' so (', outofdom_percent ,'% )' 
        end if
        if (max_thleft.ne.0) then
          outofdom_percent = (nbpart_thleft_muid(iter) / &
            float(Nbpoints_muid))*100._RKIND
          write(ifil,'(A,I8,A,f6.2,A)') &
            '   -> number of particles out of domain ' // &
            'at left in theta:  ', nbpart_thleft_muid(iter), &
            ' so (', outofdom_percent ,' % )' 
        end if
      end if
    end do
    if (.not.outofdomain) then
      write(ifil,'(A,f8.2,A,f8.2,A)') &
        ' *** time between  ',(iter_glob-nbiter)*deltat, &
        ' and ',iter_glob*deltat,' ***'
      write(ifil,'(A)') '   -> no particles out of domain '
    end if
    close(ifil)
  end subroutine print_particles_outofdomain
end module advec2D_BSL_module
