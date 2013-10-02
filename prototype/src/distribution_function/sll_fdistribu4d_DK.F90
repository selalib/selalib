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
  subroutine init_n0_eta1(eta1_peak,kappan,deltarn,n0_eta1min, &
    eta1_grid,n0_1d)
    sll_real64, intent(in) :: eta1_peak
    sll_real64, intent(in) :: kappan
    sll_real64, intent(in) :: deltarn
    sll_real64, intent(in) :: n0_eta1min
    sll_real64, dimension(:), intent(in)    :: eta1_grid
    sll_real64, dimension(:), intent(inout) :: n0_1d

    sll_int32  :: ierr, ieta1, Neta1
    sll_real64 :: delta_eta1, etai, etai_mid
    sll_real64 :: tmp, inv_Ln0
    sll_real64 :: n0norm_tmp

    Neta1      = size(eta1_grid,1)
    delta_eta1 = eta1_grid(2)-eta1_grid(1)
    
    !*** compute ns0 solution of :                           ***
    !***  2/(n0(r)+n0(r-1))*(n0(r)-n0(r-1))/dr               ***
    !***                  = -(1/Ln0)*cosh^-2(r-rpeak/deltar) ***
    inv_Ln0  = kappan           !??? (1/Ln) = kappa_n0/R
    n0_1d(1) = n0_eta1min 
    do ieta1 = 2,Neta1
      etai     = eta1_grid(ieta1-1)
      etai_mid = etai + delta_eta1*0.5_f64
      tmp      = -inv_Ln0 * &
        sech((etai_mid-eta1_peak)/deltarn)**2
      tmp          = 0.5_f64*delta_eta1*tmp
      n0_1d(ieta1) = (1._f64+tmp)/(1._f64-tmp)*n0_1d(ieta1-1)
    enddo

    !*** normalisation of the density at int(n0(r)rdr)/int(rdr) ***
    ! -> computation of int(n0(r)rdr)
    n0norm_tmp = 0._f64
    do ieta1 = 2,Neta1-1
      n0norm_tmp = n0norm_tmp + n0_1d(ieta1)*eta1_grid(ieta1)
    enddo
    n0norm_tmp = n0norm_tmp + 0.5_f64 * &
      (n0_1d(1)*eta1_grid(1) + n0_1d(Neta1)*eta1_grid(Neta1))
    ! -> division by int(rdr)
    n0norm_tmp = n0norm_tmp*2._f64*delta_eta1 / & 
      (eta1_grid(Neta1)**2-eta1_grid(1)**2)

    n0_1d(1:Neta1) = n0_1d(1:Neta1)/n0norm_tmp
  end subroutine init_n0_eta1


  !---------------------------------------------------------------
  ! Initialization of the radial temperature profiles 
  !---------------------------------------------------------------
  subroutine init_T_eta1(eta1_peak,kappaT, &
    deltarT,T_eta1min,T_scal,eta1_grid,T_1d)
    sll_real64, intent(in) :: eta1_peak
    sll_real64, intent(in) :: kappaT
    sll_real64, intent(in) :: deltarT
    sll_real64, intent(in) :: T_eta1min
    sll_real64, intent(in) :: T_scal
    sll_real64, dimension(:), intent(in)    :: eta1_grid
    sll_real64, dimension(:), intent(inout) :: T_1d
    
    sll_int32  :: ierr, ieta1, Neta1
    sll_real64 :: delta_eta1, etai, etai_mid
    sll_real64 :: tmp, inv_LT
    sll_real64 :: w0, w1, Tnorm_tmp

    Neta1      = size(eta1_grid,1)
    delta_eta1 = eta1_grid(2)-eta1_grid(1)
    
    !*** compute ns0 solution of :                           ***
    !***  2/(n0(r)+n0(r-1))*(n0(r)-n0(r-1))/dr               ***
    !***                  = -(1/Ln0)*cosh^-2(r-rpeak/deltar) ***
    inv_LT  = kappaT           !??? (1/Ln) = kappa_n0/R
    T_1d(1) = T_eta1min
    do ieta1 = 2,Neta1
      etai     = eta1_grid(ieta1-1)
      etai_mid = etai + delta_eta1*0.5_f64
      tmp      = -inv_LT * &
        sech((etai_mid-eta1_peak)/deltarT)**2
      tmp           = 0.5_f64*delta_eta1*tmp
      T_1d(ieta1) = (1._f64+tmp)/(1._f64-tmp)*T_1d(ieta1-1)
    enddo

    !*** normalisation of the temperature to 1 at r=rpeak ***
    ieta1         = int((eta1_peak-eta1_grid(1))/delta_eta1)
    w1            = (eta1_peak-eta1_grid(ieta1))/delta_eta1
    w0            = 1._f64-w1
    Tnorm_tmp     = w0*T_1d(ieta1)+w1*T_1d(ieta1+1)
    T_1d(1:Neta1) = (T_1d(1:Neta1)/Tnorm_tmp)/T_scal
  end subroutine init_T_eta1


  !---------------------------------------------------------------
  ! Computation of the equilibrium distribution function for
  !  drift-kinetic 4D simulation
  !   feq(eta1,vpar) = n0(eta1)/(2*pi*Ti(eta1))**(1/2) * 
  !                    exp(-0.5*eta4**2/Ti(eta1)
  !  where n0 and Ti are respectively the initial density 
  !  and temperature profiles
  !---------------------------------------------------------------
  function compute_feq_val( eta1, vpar, n0_eta1, Ti_eta1) &
    result(val)

    sll_real64             :: val      ! sll_DK_initializer_4d
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: vpar
    sll_real64, intent(in) :: n0_eta1
    sll_real64, intent(in) :: Ti_eta1

    val = n0_eta1*sqrt(2._f64*sll_pi*Ti_eta1) * &
      exp(-0.5_f64*vpar**2/Ti_eta1)
  end function compute_feq_val


  !----------------------------------------------------
  ! Initialization of the 2D array for the equilibrium
  !  distribution function feq(eta1,vpar)
  !----------------------------------------------------
  subroutine init_fequilibrium(Neta1,Nvpar, &
    eta1_grid,vpar_grid,n0_1d,Ti_1d,feq_2d)
    sll_int32, intent(in) :: Neta1
    sll_int32, intent(in) :: Nvpar
    sll_real64, dimension(:)  , intent(in)  :: eta1_grid
    sll_real64, dimension(:)  , intent(in)  :: vpar_grid
    sll_real64, dimension(:)  , intent(in)  :: n0_1d
    sll_real64, dimension(:)  , intent(in)  :: Ti_1d
    sll_real64, dimension(:,:), intent(out) :: feq_2d
    
    sll_int32  :: ieta1, ivpar
    sll_real64 :: eta1, vpar, n0_eta1, Ti_eta1

    do ivpar = 1,Nvpar
      vpar = vpar_grid(ivpar)
      do ieta1 = 1,Neta1
        eta1    = eta1_grid(ieta1)
        n0_eta1 = n0_1d(ieta1)
        Ti_eta1 = Ti_1d(ieta1)
        feq_2d(ieta1,ivpar) = compute_feq_val(eta1,vpar,n0_eta1,Ti_eta1)
      end do
    end do
  end subroutine init_fequilibrium

end module sll_fdistribu4d_DK
