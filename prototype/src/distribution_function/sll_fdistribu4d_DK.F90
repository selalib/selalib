!-----------------------------------------------------------
! SELALIB
!-----------------------------------------------------------
!
! MODULE: sll_fdistribu4d_DK
!
!> @author
!> - Virginie Grandgirard
!
! DESCRIPTION: 
!
!> @brief
!> Initialisation of the distribution function for the
!>  drift-kinetic 4D simulation
!>
!>@details
!-----------------------------------------------------------
module sll_fdistribu4d_DK
#include "sll_working_precision.h"
#include "sll_memory.h"

  use sll_constants

  implicit none

  sll_real64, dimension(1) :: whatever  ! dummy params array

  contains

  !---------------------------------------- 
  ! sech = cosh^-1 definition
  !---------------------------------------- 
  function sech(x)
    sll_real64, intent(in) :: x   
    sll_real64             :: sech
    
    sech = 1._f64/cosh(x)
  end function sech  


  !---------------------------------------- -----------------------
  ! Initialization of the radial density profiles 
  !---------------------------------------- -----------------------
  subroutine init_n0_r(r_peak,inv_Ln,deltarn,n0_rmin, &
    r_grid,n0_1d)
    sll_real64, intent(in) :: r_peak
    sll_real64, intent(in) :: inv_Ln
    sll_real64, intent(in) :: deltarn
    sll_real64, intent(in) :: n0_rmin
    sll_real64, dimension(:), intent(in)    :: r_grid
    sll_real64, dimension(:), intent(inout) :: n0_1d

    sll_int32  :: ierr, ir, Nr
    sll_real64 :: dr, Lr
    sll_real64 :: etai, etai_mid
    sll_real64 :: tmp
    sll_real64 :: n0norm_tmp
    sll_real64 :: deltarn_norm

    Nr  = size(r_grid,1)
    Lr  = r_grid(Nr)-r_grid(1)
    dr  = r_grid(2)-r_grid(1)

    deltarn_norm = deltarn*Lr
    
    !*** compute ns0 solution of :                           ***
    !***  2/(n0(r)+n0(r-1))*(n0(r)-n0(r-1))/dr               ***
    !***                  = -(1/Ln0)*cosh^-2(r-rpeak/deltar) ***
    !***  where (1/Ln) = kappa_n0/R                          ***
    n0_1d(1) = n0_rmin 
    do ir = 2,Nr
      etai     = r_grid(ir-1)
      etai_mid = etai + dr*0.5_f64
      tmp      = -inv_Ln * &
        sech((etai_mid-r_peak)/deltarn_norm)**2
      tmp          = 0.5_f64*dr*tmp
      n0_1d(ir) = (1._f64+tmp)/(1._f64-tmp)*n0_1d(ir-1)
    enddo

    !*** normalisation of the density at int(n0(r)rdr)/int(rdr) ***
    ! -> computation of int(n0(r)rdr)
    n0norm_tmp = 0._f64
    do ir = 2,Nr-1
      n0norm_tmp = n0norm_tmp + n0_1d(ir)*r_grid(ir)
    enddo
    n0norm_tmp = n0norm_tmp + 0.5_f64 * &
      (n0_1d(1)*r_grid(1) + n0_1d(Nr)*r_grid(Nr))
    ! -> division by int(rdr)
    n0norm_tmp = n0norm_tmp*2._f64*dr / & 
      (r_grid(Nr)**2-r_grid(1)**2)

    n0_1d(1:Nr) = n0_1d(1:Nr)/n0norm_tmp
  end subroutine init_n0_r


  !---------------------------------------------------------------
  ! Initialization of the radial temperature profiles 
  !---------------------------------------------------------------
  subroutine init_T_r(r_peak,inv_LT, &
    deltarT,T_rmin,T_scal,r_grid,T_1d)
    sll_real64, intent(in) :: r_peak
    sll_real64, intent(in) :: inv_LT
    sll_real64, intent(in) :: deltarT
    sll_real64, intent(in) :: T_rmin
    sll_real64, intent(in) :: T_scal
    sll_real64, dimension(:), intent(in)    :: r_grid
    sll_real64, dimension(:), intent(inout) :: T_1d
    
    sll_int32  :: ierr, ir, Nr
    sll_real64 :: dr, Lr
    sll_real64 :: etai, etai_mid
    sll_real64 :: tmp
    sll_real64 :: w0, w1, Tnorm_tmp
    sll_real64 :: deltarT_norm

    Nr  = size(r_grid,1)
    Lr  = r_grid(Nr)-r_grid(1)
    dr  = r_grid(2)-r_grid(1)
    deltarT_norm = deltarT*Lr

    !*** compute Ts solution of :                           ***
    !***  2/(Ts(r)+Ts(r-1))*(Ts(r)-Ts(r-1))/dr               ***
    !***                  = -(1/LTs)*cosh^-2(r-rpeak/deltar) ***
    !***  where (1/LT) = kappa_Ts/R                          ***
    T_1d(1) = T_rmin
    do ir = 2,Nr
      etai     = r_grid(ir-1)
      etai_mid = etai + dr*0.5_f64
      tmp      = -inv_LT * &
        sech((etai_mid-r_peak)/deltarT_norm)**2
      tmp      = 0.5_f64*dr*tmp
      T_1d(ir) = (1._f64+tmp)/(1._f64-tmp)*T_1d(ir-1)
    enddo

    !*** normalisation of the temperature to 1 at r=rpeak ***
    ir         = int((r_peak-r_grid(1))/dr)
    w1         = (r_peak-r_grid(ir))/dr
    w0         = 1._f64-w1
    Tnorm_tmp  = w0*T_1d(ir)+w1*T_1d(ir+1)
    T_1d(1:Nr) = (T_1d(1:Nr)/Tnorm_tmp)/T_scal
  end subroutine init_T_r


  !---------------------------------------------------------------
  ! Initialisation of the magnetic field B(r,theta)
  !---------------------------------------------------------------
  subroutine init_Brtheta(r_grid,theta_grid,B_rtheta)
    sll_real64, dimension(:)  , intent(in)    :: r_grid
    sll_real64, dimension(:)  , intent(in)    :: theta_grid
    sll_real64, dimension(:,:), intent(inout) :: B_rtheta

    sll_int32 :: ir, itheta
    sll_int32 :: Nr, Ntheta

    Nr     = size(r_grid,1)
    Ntheta = size(theta_grid,1)

    do itheta = 1,Ntheta
      do ir = 1,Nr
        B_rtheta(ir,itheta) = 1._f64
      end do
    end do
  end subroutine init_Brtheta


  !---------------------------------------------------------------
  ! Computation of the equilibrium distribution function for
  !  drift-kinetic 4D simulation
  !   feq(r,vpar) = n0(r)/(2*pi*Ti(r))**(1/2) * 
  !                    exp(-0.5*vpar**2/Ti(r))
  !  where n0 and Ti are respectively the initial density 
  !  and temperature profiles
  !---------------------------------------------------------------
  function compute_feq_val( r, vpar, n0_r, Ti_r) &
    result(val)

    sll_real64             :: val      ! sll_DK_initializer_4d
    sll_real64, intent(in) :: r
    sll_real64, intent(in) :: vpar
    sll_real64, intent(in) :: n0_r
    sll_real64, intent(in) :: Ti_r

    val = n0_r/sqrt(2._f64*sll_pi*Ti_r) * &
      exp(-0.5_f64*vpar**2/Ti_r)
  end function compute_feq_val


  !----------------------------------------------------
  ! Initialization of the 2D array for the equilibrium
  !  distribution function feq(r,vpar)
  !  feq(r,vpar) = n0(r)/(2*pi*Ti(r))**(1/2) * 
  !                    exp(-0.5*vpar**2/Ti(r))
  !----------------------------------------------------
  subroutine init_fequilibrium(Nr,Nvpar, &
    r_grid,vpar_grid,n0_1d,Ti_1d,feq_2d)
    sll_int32, intent(in) :: Nr
    sll_int32, intent(in) :: Nvpar
    sll_real64, dimension(:)  , intent(in)  :: r_grid
    sll_real64, dimension(:)  , intent(in)  :: vpar_grid
    sll_real64, dimension(:)  , intent(in)  :: n0_1d
    sll_real64, dimension(:)  , intent(in)  :: Ti_1d
    sll_real64, dimension(:,:), intent(out) :: feq_2d
    
    sll_int32  :: ir, ivpar
    sll_real64 :: r, vpar, n0_r, Ti_r

    do ivpar = 1,Nvpar
      vpar = vpar_grid(ivpar)
      do ir = 1,Nr
        r    = r_grid(ir)
        n0_r = n0_1d(ir)
        Ti_r = Ti_1d(ir)
        feq_2d(ir,ivpar) = compute_feq_val(r,vpar,n0_r,Ti_r)
      end do
    end do
  end subroutine init_fequilibrium


  !----------------------------------------------------
  ! Compute func(x,y) from func(r)
  !----------------------------------------------------
  subroutine function_xy_from_r(r_grid,func_r, &
    xgrid_2d,ygrid_2d,func_xy)
    use sll_cubic_splines
    use sll_common_coordinate_transformations, only : &
      polar_eta1
    sll_real64, dimension(:)  , intent(in)  :: r_grid
    sll_real64, dimension(:)  , intent(in)  :: func_r
    sll_real64, dimension(:,:), intent(in)  :: xgrid_2d
    sll_real64, dimension(:,:), intent(in)  :: ygrid_2d
    sll_real64, dimension(:,:), intent(out) :: func_xy

    sll_int32  :: Nr, Npt1, Npt2
    sll_int32  :: ix, iy
    sll_real64 :: r, x, y
 
    type(sll_cubic_spline_1D), pointer :: sp1d_r

    Nr   = size(r_grid,1)
    Npt1 = size(xgrid_2d,1)
    Npt2 = size(xgrid_2d,2)

    sp1d_r => new_spline_1d(Npt1, &
      r_grid(1),r_grid(Nr),SLL_HERMITE)
    call compute_spline_1D(func_r,sp1d_r)

    do iy = 1,Npt2
      do ix = 1,Npt1
        x = xgrid_2d(ix,iy)
        y = ygrid_2d(ix,iy)
        r = polar_eta1(x,y,(/0.0_f64/)) ! params doesn't matter for polar_eta1 
        r = min(max(r,r_grid(1)),r_grid(Nr))
        func_xy(ix,iy) = interpolate_value(r,sp1d_r)
      end do
    end do
    call delete(sp1d_r)
  end subroutine function_xy_from_r


  !----------------------------------------------------
  ! Compute func(x,y) from func(r,theta)
  !----------------------------------------------------
  subroutine function_xy_from_rtheta(r_grid,theta_grid, &
    func_rtheta,xgrid_2d,ygrid_2d,func_xy)
    use sll_constants
    use sll_cubic_splines
    use sll_common_coordinate_transformations, only : &
      polar_eta1, polar_eta2
    sll_real64, dimension(:)  , intent(in)  :: r_grid
    sll_real64, dimension(:)  , intent(in)  :: theta_grid
    sll_real64, dimension(:,:), intent(in)  :: func_rtheta
    sll_real64, dimension(:,:), intent(in)  :: xgrid_2d
    sll_real64, dimension(:,:), intent(in)  :: ygrid_2d
    sll_real64, dimension(:,:), intent(out) :: func_xy

    sll_int32  :: Nr, Ntheta
    sll_int32  :: Npt1, Npt2
    sll_int32  :: ix, iy
    sll_real64 :: r, theta, x, y
 
    type(sll_cubic_spline_2D), pointer :: sp2d_rtheta

    Nr     = size(r_grid,1)
    Ntheta = size(theta_grid,1)
    Npt1   = size(xgrid_2d,1)
    Npt2   = size(xgrid_2d,2)

    sp2d_rtheta => new_spline_2d(Npt1,Npt2, &
      r_grid(1),r_grid(Nr), &
      theta_grid(1),theta_grid(Ntheta), &
      SLL_HERMITE,SLL_PERIODIC)
    call compute_spline_2D(func_rtheta,sp2d_rtheta)

    do iy = 1,Npt2
      do ix = 1,Npt1
        x     = xgrid_2d(ix,iy)
        y     = ygrid_2d(ix,iy)
        r     = polar_eta1(x,y,whatever)
        r     = min(max(r,r_grid(1)),r_grid(Nr))
        theta = polar_eta2(x,y,whatever)
        theta = modulo(theta,2._f64*sll_pi)
        func_xy(ix,iy) = interpolate_value_2D(r,theta,sp2d_rtheta)
      end do
    end do
    call delete(sp2d_rtheta)
  end subroutine function_xy_from_rtheta


  !----------------------------------------------------
  ! Initialization of the 3D array for the equilibrium
  !  distribution function feq(x,y,vpar)
  !  feq(x,y,vpar) = n0(x,y)/(2*pi*Ti(x,y))**(1/2) * 
  !                    exp(-0.5*vpar**2/Ti(x,y))
  !----------------------------------------------------
  subroutine init_fequilibrium_xy(xgrid_2d,ygrid_2d, &
    vpar_grid,n0_xy,Ti_xy,feq_xyvpar)
    sll_real64, dimension(:,:)  , intent(in)  :: xgrid_2d
    sll_real64, dimension(:,:)  , intent(in)  :: ygrid_2d
    sll_real64, dimension(:)    , intent(in)  :: vpar_grid
    sll_real64, dimension(:,:)  , intent(in)  :: n0_xy
    sll_real64, dimension(:,:)  , intent(in)  :: Ti_xy
    sll_real64, dimension(:,:,:), intent(out) :: feq_xyvpar

    sll_int32 :: Npt1, Npt2, Nvpar
    sll_int32 :: ix, iy, ivpar

    Npt1  = size(xgrid_2d,1)
    Npt2  = size(xgrid_2d,2)
    Nvpar = size(vpar_grid,1)

    do ivpar = 1,Nvpar
      do iy = 1,Npt2
        do ix = 1,Npt1
           feq_xyvpar(ix,iy,ivpar) = n0_xy(ix,iy) * &
                exp(-0.5_f64*vpar_grid(ivpar)**2/Ti_xy(ix,iy))&
                /sqrt(2._f64*sll_pi*Ti_xy(ix,iy))
        end do
      end do
    end do
  end subroutine init_fequilibrium_xy
end module sll_fdistribu4d_DK
