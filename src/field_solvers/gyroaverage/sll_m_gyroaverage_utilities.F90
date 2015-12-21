module sll_m_gyroaverage_utilities
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_gauss_legendre_integration, only: &
    sll_f_gauss_legendre_points, &
    sll_f_gauss_legendre_weights

  use sll_m_gauss_lobatto_integration, only: &
    sll_f_gauss_lobatto_points, &
    sll_f_gauss_lobatto_weights

  implicit none

  public :: &
    sll_s_compute_init_f_polar, &
    sll_s_compute_mu, &
    sll_s_compute_shape_circle, &
    sll_s_zero_bessel_dir_dir

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains

  subroutine sll_s_compute_shape_circle(points,N_points)
    sll_int32,intent(in) :: N_points
    sll_real64,dimension(:,:) ::points
    sll_int32 :: i
    sll_real64 :: x
    do i=1,N_points
      x = 2._f64*sll_p_pi*real(i,f64)/(real(N_points,f64))
      points(1,i) = cos(x)
      points(2,i) = sin(x)
      points(3,i) = 1._f64/real(N_points,f64)
    enddo
   
  end subroutine sll_s_compute_shape_circle


  subroutine sll_s_compute_init_f_polar(f,mode,N,eta_min,eta_max)
    sll_real64,dimension(:,:),intent(out)::f
    sll_int32,intent(in)::N(2),mode(2)
    sll_real64,intent(in)::eta_min(2),eta_max(2)
    sll_int32::i,j
    sll_real64::eta(2),delta_eta(2),kmode,val

    call sll_s_zero_bessel_dir_dir(mode,eta_min(1),eta_max(1),val)
    delta_eta(1)=(eta_max(1)-eta_min(1))/real(N(1),f64)
    delta_eta(2)=(eta_max(2)-eta_min(2))/real(N(2),f64)
    
    kmode=real(mode(2),f64)*(2._f64*sll_p_pi)/(eta_max(2)-eta_min(2))
    
        
    do j=1,N(2)+1
      eta(2)=eta_min(2)+real(j-1,f64)*delta_eta(2)
      do i=1,N(1)+1
        eta(1)=eta_min(1)+real(i-1,f64)*delta_eta(1)
        eta(1)=val*eta(1)/eta_max(1)
        f(i,j)= 0._f64 !temporary, because DBESJ not recognized on helios & curie
        !f(i,j)=cos(kmode*eta(2))
        !f(i,j) = (DBESJN(mode(2),val)*DBESYN(mode(2),eta(1))-DBESYN(mode(2),val)*DBESJN(mode(2),eta(1)))*cos(kmode*eta(2))
      enddo
    enddo   
    
  end subroutine sll_s_compute_init_f_polar
  
  
   subroutine sll_s_zero_bessel_dir_dir(mode,eta_min,eta_max,val)
    sll_real64,intent(in)::eta_min,eta_max
    sll_int32,intent(in)::mode(2)
    sll_real64,intent(out)::val
    sll_real64::alpha,tmp
    sll_int32::mode_max(2),i,j
    logical::is_file
    val=0._f64
    INQUIRE(FILE="zeros_bessel.txt", EXIST=is_file)
    if((is_file).eqv.(.false.))then
      print *,'#file zeros_bessel.txt does not exist'
      return
    endif  
    open(27,file='zeros_bessel.txt',action="read")
      read(27,*) mode_max(1),mode_max(2),alpha
    close(27) 
    if((mode(1)<1).or.(mode(1)>mode_max(1)))then
      print *,'#bad value of mode(1) vs mode_max(1)',mode(1),mode_max(1)
      return
    endif
    if((mode(2)<0).or.(mode(2)>mode_max(2)))then
      print *,'#bad value of mode(2) vs mode_max(2)',mode(2),mode_max(2)
      return
    endif
    if(abs(alpha-eta_min/eta_max)>1.e-12)then
      print *,'#bad value of rmin/rmax w.r.t zeros_bessel.txt',eta_min/eta_max,alpha
      return
    endif
    open(27,file='zeros_bessel.txt',action="read")
      read(27,*) mode_max(1),mode_max(2),alpha
      read(27,*) i,j,tmp
      do while((i.ne.mode(1)).or.(j.ne.mode(2)))
        read(27,*) i,j,tmp
      enddo
    close(27) 
    val = tmp
      
  end subroutine sll_s_zero_bessel_dir_dir

  subroutine sll_s_compute_mu( &
    quadrature_case, &
    mu_points, &
    mu_weights, &
    N_mu, &
    mu_min, &
    mu_max, &
    quadrature_points_per_cell)
    character(len=256), intent(in) :: quadrature_case
    sll_real64, dimension(:), intent(out) :: mu_points
    sll_real64, dimension(:), intent(out) :: mu_weights
    sll_int32, intent(in) :: N_mu
    sll_real64, intent(in) :: mu_min
    sll_real64, intent(in) :: mu_max
    sll_int32, intent(in) :: quadrature_points_per_cell
    sll_int32 :: i
    sll_int32 :: s
    sll_int32 :: num_cells
    
    select case(quadrature_case)
      case ("SLL_RECTANGLE")
        do i=1,N_mu
          mu_points(i) = mu_min+real(i-1,f64)*(mu_max-mu_min)/real(N_mu,f64)
          mu_weights(i) = (mu_max-mu_min)/real(N_mu,f64)*exp(-mu_points(i))
        enddo
        if(N_mu==1)then
          mu_weights(1) = 1._f64
        endif
      case ("SLL_GAUSS_LOBATTO")
        num_cells = N_mu/quadrature_points_per_cell
        mu_points(1:N_mu) = mu_max
        mu_weights(1:N_mu) = 0._f64
        s=1
        do i=1,num_cells
          mu_points(s:s+quadrature_points_per_cell-1) = &
            sll_f_gauss_lobatto_points( &
              quadrature_points_per_cell, &
              mu_min+real(i-1,f64)/real(num_cells,f64)*(mu_max-mu_min), &
              mu_min+real(i,f64)/real(num_cells,f64)*(mu_max-mu_min) )
          mu_weights(s:s+quadrature_points_per_cell-1) = &
            sll_f_gauss_lobatto_weights( &
              quadrature_points_per_cell, &
              mu_min+real(i-1,f64)/real(num_cells,f64)*(mu_max-mu_min), &
              mu_min+real(i,f64)/real(num_cells,f64)*(mu_max-mu_min) )
          s=s+quadrature_points_per_cell        
        enddo
        !mu_points(1:N_mu) = sll_f_gauss_lobatto_points( N_mu, 0._f64, mu_max )
        !mu_weights(1:N_mu) = sll_f_gauss_lobatto_weights( N_mu, 0._f64, mu_max )
        do i=1,N_mu
          mu_weights(i) = mu_weights(i)*exp(-mu_points(i))
        enddo       
      case ("SLL_GAUSS_LEGENDRE")
        num_cells = N_mu/quadrature_points_per_cell
        mu_points(1:N_mu) = mu_max
        mu_weights(1:N_mu) = 0._f64
        s=1
        do i=1,num_cells
          mu_points(s:s+quadrature_points_per_cell-1) = &
            sll_f_gauss_legendre_points( &
              quadrature_points_per_cell, &
              mu_min+real(i-1,f64)/real(num_cells,f64)*(mu_max-mu_min), &
              mu_min+real(i,f64)/real(num_cells,f64)*(mu_max-mu_min) )
          mu_weights(s:s+quadrature_points_per_cell-1) = &
            sll_f_gauss_legendre_weights( &
              quadrature_points_per_cell, &
              mu_min+real(i-1,f64)/real(num_cells,f64)*(mu_max-mu_min), &
              mu_min+real(i,f64)/real(num_cells,f64)*(mu_max-mu_min) )
          s=s+quadrature_points_per_cell        
        enddo
        !mu_points(1:N_mu) = sll_f_gauss_lobatto_points( N_mu, 0._f64, mu_max )
        !mu_weights(1:N_mu) = sll_f_gauss_lobatto_weights( N_mu, 0._f64, mu_max )
        do i=1,N_mu
          mu_weights(i) = mu_weights(i)*exp(-mu_points(i))
        enddo       
      case default
        print *,'#bad quadrature_case',trim(quadrature_case)
        print *,'#not implemented'
        print *,'#in sll_s_compute_mu'
        print*,'#at line and file:',__LINE__,__FILE__
        stop
    end select

    
  end subroutine sll_s_compute_mu


end module sll_m_gyroaverage_utilities
