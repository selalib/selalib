!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************


module sll_gyroaverage_2d_polar
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_fft
  use sll_tridiagonal
  use sll_constants
  use sll_boundary_condition_descriptors

  implicit none

  type sll_plan_gyroaverage_polar
     
     sll_real64          :: eta_min(2)     !< r min et theta min
     sll_real64          :: eta_max(2)     !< r max et theta max
     sll_int32           :: Nc(2)          !< number of cells in r and in theta
     
     sll_int32           :: N_points          !< number of points on the circle
     sll_int32           :: interp_degree(2)  !< interpolation degrees in r,theta

     sll_real64, dimension(:,:), pointer    :: points
     sll_real64, dimension(:,:,:), pointer  :: deriv
     sll_int32, dimension(:), pointer       :: pre_compute_N
     sll_real64, dimension(:,:), pointer    :: pre_compute_coeff
     sll_int32, dimension(:,:), pointer     :: pre_compute_index
     sll_real64, dimension(:), pointer      :: pre_compute_coeff_spl
     sll_real64, dimension(:,:), pointer    :: lunat
     sll_real64, dimension(:), pointer      :: luper
     sll_real64, dimension(:,:,:), pointer  :: A_fft

  end type sll_plan_gyroaverage_polar

contains


  function new_plan_gyroaverage_polar_hermite(eta_min,eta_max,Nc,N_points,interp_degree,deriv_size) result(this)

    implicit none

    sll_real64, intent(in) :: eta_min(2)
    sll_real64, intent(in) :: eta_max(2)
    sll_int32, intent(in)  :: Nc(2)
    sll_int32, intent(in)  :: N_points  
    sll_int32, intent(in)  :: interp_degree(2)
    sll_int32, intent(in)  :: deriv_size
    type(sll_plan_gyroaverage_polar), pointer :: this

    sll_int32 :: err

    SLL_ALLOCATE(this,err)
    SLL_ALLOCATE(this%deriv(deriv_size,Nc(1)+1,Nc(2)+1),err)
    SLL_ALLOCATE(this%points(3,N_points),err)
       
    call compute_shape_circle(this%points,N_points)   
       
    this%eta_min=eta_min
    this%eta_max=eta_max
    this%Nc=Nc
    this%N_points=N_points
    this%interp_degree=interp_degree
    
  end function new_plan_gyroaverage_polar_hermite
  
  
  
  function new_plan_gyroaverage_polar_splines(eta_min,eta_max,Nc,N_points) result(this)

    implicit none

    sll_real64, intent(in) :: eta_min(2)
    sll_real64, intent(in) :: eta_max(2)
    sll_int32, intent(in)  :: Nc(2)
    sll_int32, intent(in)  :: N_points
    type(sll_plan_gyroaverage_polar), pointer :: this

    sll_int32 :: err

    SLL_ALLOCATE(this,err)
    SLL_ALLOCATE(this%points(3,N_points),err)
    
    call compute_shape_circle(this%points,N_points) 
       
    this%eta_min=eta_min
    this%eta_max=eta_max
    this%Nc=Nc
    this%N_points=N_points
    
  end function new_plan_gyroaverage_polar_splines

  
  
  function new_plan_gyroaverage_polar_pade(eta_min,eta_max,Nc) result(this)

    implicit none

    sll_real64, intent(in) :: eta_min(2)
    sll_real64, intent(in) :: eta_max(2)
    sll_int32, intent(in)  :: Nc(2)
    type(sll_plan_gyroaverage_polar), pointer :: this
    sll_int32 :: err
     
    SLL_ALLOCATE(this,err) 
       
    this%eta_min=eta_min
    this%eta_max=eta_max
    this%Nc=Nc

  end function new_plan_gyroaverage_polar_pade

  
  subroutine compute_gyroaverage_points_polar_hermite(gyro,f,rho)
    type(sll_plan_gyroaverage_polar)        :: gyro
    sll_real64,dimension(:,:),intent(inout) :: f
    sll_real64,intent(in)                   :: rho
    
    sll_int32                               :: i,j,k,ii(2),s
    sll_real64                              ::fval,sum_fval,eta_star(2),eta(2),delta_eta(2),x(2)
    
    fval=0._f64
    delta_eta(1)=(gyro%eta_max(1)-gyro%eta_min(1))/real(gyro%Nc(1),f64)
    delta_eta(2)=(gyro%eta_max(2)-gyro%eta_min(2))/real(gyro%Nc(2),f64)
    
    call hermite_coef_nat_per(f(1:gyro%Nc(1)+1,1:gyro%Nc(2)),gyro%deriv,gyro%Nc,gyro%interp_degree)

    
    do j=1,gyro%Nc(2)
     eta(2)=gyro%eta_min(2)+real(j-1,f64)*delta_eta(2)
      do i=1,gyro%Nc(1)+1
        eta(1)=gyro%eta_min(1)+real(i-1,f64)*delta_eta(1)       
        sum_fval = 0._f64
        do k=1,gyro%N_points
          x(1) = eta(1)*cos(eta(2))+rho*gyro%points(1,k)
          x(2) = eta(1)*sin(eta(2))+rho*gyro%points(2,k)
          call localize_polar(x,gyro%eta_min,gyro%eta_max,ii,eta_star,gyro%Nc)
          call interpolate_hermite(gyro%deriv,ii,eta_star,fval,gyro%Nc)
          sum_fval = sum_fval+gyro%points(3,k)*fval
        enddo
        f(i,j) = sum_fval
      enddo
    enddo
    
    f(1:gyro%Nc(1)+1,gyro%Nc(2)+1)=f(1:gyro%Nc(1)+1,1)
    
  end subroutine compute_gyroaverage_points_polar_hermite
  
  
  
  subroutine compute_gyroaverage_points_polar_hermite_c1(gyro,f,rho)
    type(sll_plan_gyroaverage_polar)  :: gyro
    sll_real64,dimension(:,:),intent(inout) :: f
    sll_real64,intent(in)::rho
    sll_int32 ::i,j,k,ii(2),s
    sll_real64::fval,sum_fval,eta_star(2),eta(2),delta_eta(2),x(2)

    fval=0._f64
    delta_eta(1)=(gyro%eta_max(1)-gyro%eta_min(1))/real(gyro%Nc(1),f64)
    delta_eta(2)=(gyro%eta_max(2)-gyro%eta_min(2))/real(gyro%Nc(2),f64)
    
    call hermite_c1_coef_nat_per(f(1:gyro%Nc(1)+1,1:gyro%Nc(2)),gyro%deriv,gyro%Nc,gyro%interp_degree)

    do j=1,gyro%Nc(2)
     eta(2)=gyro%eta_min(2)+real(j-1,f64)*delta_eta(2)
     !eta(2)=gyro%eta_min(2) ! to uncomment for compatibility with precompute 
      do i=1,gyro%Nc(1)+1
        eta(1)=gyro%eta_min(1)+real(i-1,f64)*delta_eta(1)       
        sum_fval = 0._f64
        do k=1,gyro%N_points
          x(1) = eta(1)*cos(eta(2))+rho*gyro%points(1,k)
          x(2) = eta(1)*sin(eta(2))+rho*gyro%points(2,k)
		  eta_star(1)=sqrt(x(1)**2+x(2)**2)
          call localize_nat(ii(1),eta_star(1),gyro%eta_min(1),gyro%eta_max(1),gyro%Nc(1))
    	  eta_star(2)=modulo(datan2(x(2),x(1)),2._f64*sll_pi)
    	  !eta_star(2)=modulo(datan2(x(2),x(1))+real(j-1,f64)*delta_eta(2),2._f64*M_PI) &
          		! to uncomment for compatibility with precompute
   		  call localize_per(ii(2),eta_star(2),gyro%eta_min(2),gyro%eta_max(2),gyro%Nc(2))
          call interpolate_hermite_c1(gyro%deriv,ii,eta_star,fval,gyro%Nc)
          sum_fval = sum_fval+gyro%points(3,k)*fval
        enddo
        f(i,j) = sum_fval
      enddo
    enddo

    !f(1:gyro%Nc(1)+1,gyro%Nc(2)+1)=f(1:gyro%Nc(1)+1,1)
    
  end subroutine compute_gyroaverage_points_polar_hermite_c1
  
  
  
  
   subroutine compute_gyroaverage_points_polar_with_invar_hermite_c1(gyro,f,rho)
    type(sll_plan_gyroaverage_polar)  :: gyro
    sll_real64,dimension(:,:),intent(inout) :: f
    sll_real64,intent(in)::rho
    sll_int32 ::i,j,k,ii(2),s
    sll_real64::fval,eta_star(2),eta(2),delta_eta(2),x(2),angle
	sll_real64,dimension(:,:),allocatable::sum_fval
	sll_int32 ::error

    SLL_ALLOCATE(sum_fval(0:gyro%Nc(1),0:gyro%Nc(2)-1),error)

    sum_fval = 0._f64
    delta_eta(1)=(gyro%eta_max(1)-gyro%eta_min(1))/real(gyro%Nc(1),f64)
    delta_eta(2)=(gyro%eta_max(2)-gyro%eta_min(2))/real(gyro%Nc(2),f64)
    
    call hermite_c1_coef_nat_per(f(1:gyro%Nc(1)+1,1:gyro%Nc(2)),gyro%deriv,gyro%Nc,gyro%interp_degree)

    do i=0,gyro%Nc(1)
      eta(1)=gyro%eta_min(1)+real(i,f64)*delta_eta(1)  ! rayon du centre
      eta(2)=gyro%eta_min(2) ! angle quelconque pour le centre
      do k=1,gyro%N_points
        x(1)=eta(1)*cos(eta(2))+rho*gyro%points(1,k)
        x(2)=eta(1)*sin(eta(2))+rho*gyro%points(2,k)
        eta_star(1)=sqrt(x(1)**2+x(2)**2)
        call localize_nat(ii(1),eta_star(1),gyro%eta_min(1),gyro%eta_max(1),gyro%Nc(1))
        angle=datan2(x(2),x(1))
        do j=0,gyro%Nc(2)-1
          eta_star(2)=modulo(angle+real(j,f64)*delta_eta(2),2._f64*sll_pi)
          call localize_per(ii(2),eta_star(2),gyro%eta_min(2),gyro%eta_max(2),gyro%Nc(2))
	      call interpolate_hermite_c1(gyro%deriv,ii,eta_star,fval,gyro%Nc)
	      sum_fval(i,j)=sum_fval(i,j)+gyro%points(3,k)*fval
	    enddo
	  enddo
	enddo  

    f(1:gyro%Nc(1)+1,1:gyro%Nc(2))=sum_fval(0:gyro%Nc(1),0:gyro%Nc(2)-1)
    f(1:gyro%Nc(1)+1,gyro%Nc(2)+1)=f(1:gyro%Nc(1)+1,1)  
   
    SLL_DEALLOCATE_ARRAY(sum_fval,error)
    
  end subroutine compute_gyroaverage_points_polar_with_invar_hermite_c1
  
  
  
  
  
  subroutine pre_compute_gyroaverage_polar_hermite_c1(gyro,rho)
    type(sll_plan_gyroaverage_polar)  :: gyro
    sll_real64,intent(in)::rho
    sll_int32,dimension(:,:),allocatable :: buf
    sll_int32 ::i,j,k,ell_1,ell_2,ii(2),s,nb,ind(2)
    sll_real64::val(4,0:1,0:1),eta_star(2),eta(2),delta_eta(2),x(2)
    sll_int32 ::error
    
    val=0._f64
    delta_eta(1)=(gyro%eta_max(1)-gyro%eta_min(1))/real(gyro%Nc(1),f64)
    delta_eta(2)=(gyro%eta_max(2)-gyro%eta_min(2))/real(gyro%Nc(2),f64)
    
    SLL_ALLOCATE(buf(0:gyro%Nc(1),0:gyro%Nc(2)),error)    
    
    SLL_ALLOCATE(gyro%pre_compute_N(gyro%Nc(1)+1),error)

    !gyro%pre_compute_N(i)
    !gyro%pre_compute_index(1:2,s)
    !gyro%deriv(ell,ii(1),ii(2))*gyro_pre_compute_coeff(s)



    eta(2)=gyro%eta_min(2)
    buf=0
    nb=0
    do i=1,gyro%Nc(1)+1
      eta(1)=gyro%eta_min(1)+real(i-1,f64)*delta_eta(1)       
      s=0
      do k=1,gyro%N_points
        x(1) = eta(1)*cos(eta(2))+rho*gyro%points(1,k)
        x(2) = eta(1)*sin(eta(2))+rho*gyro%points(2,k)
        call localize_polar(x,gyro%eta_min,gyro%eta_max,ii,eta_star,gyro%Nc)
        do ell_2=0,1
          ind(2)=ii(2)+ell_2
          do ell_1=0,1
            ind(1)=ii(1)+ell_1
            if(buf(ind(1),ind(2)).ne.i)then
              s=s+1
              buf(ind(1),ind(2))=i
            endif
          enddo
        enddo    
      enddo
      gyro%pre_compute_N(i)=s
      nb=nb+s
    enddo
    
    print*,'#N_points pre_compute=',nb,gyro%Nc(1),nb/gyro%Nc(1)

    SLL_ALLOCATE(gyro%pre_compute_index(1:2,nb),error)
    SLL_ALLOCATE(gyro%pre_compute_coeff(4,nb),error)
    
    buf=0
    nb=0
    s=0
    val=0._f64
    eta(2)=gyro%eta_min(2)
    do i=1,gyro%Nc(1)+1
      eta(1)=gyro%eta_min(1)+real(i-1,f64)*delta_eta(1)       
      do k=1,gyro%N_points
        x(1) = eta(1)*cos(eta(2))+rho*gyro%points(1,k)
        x(2) = eta(1)*sin(eta(2))+rho*gyro%points(2,k)
        call localize_polar(x,gyro%eta_min,gyro%eta_max,ii,eta_star,gyro%Nc)
        call contribution_hermite_c1(eta_star,val)
        !print *,eta_star
        !print *,i,k,val(1,0,0),val(1,1,0),ii(1),eta_star(1)
                
        val=val*gyro%points(3,k)
        do ell_2=0,1
          ind(2)=ii(2)+ell_2
          do ell_1=0,1
            ind(1)=ii(1)+ell_1
            j=buf(ind(1),ind(2))
            if(j<=nb)then
              s=s+1
              buf(ind(1),ind(2))=s
              gyro%pre_compute_coeff(1:4,s)=val(1:4,ell_1,ell_2)
              gyro%pre_compute_index(1,s)=ind(1)
              gyro%pre_compute_index(2,s)=ind(2)                     
            else              
              gyro%pre_compute_coeff(1:4,j)=gyro%pre_compute_coeff(1:4,j)+val(1:4,ell_1,ell_2)
            endif
          enddo
        enddo
      enddo      
      nb=s 
    enddo
    
    !print *,gyro%pre_compute_coeff
    !stop
    
    
    print *,'#min/max=',minval(gyro%pre_compute_coeff(1:4,:)),maxval(gyro%pre_compute_coeff(1:4,:))
    SLL_DEALLOCATE_ARRAY(buf,error)
  end subroutine pre_compute_gyroaverage_polar_hermite_c1
  
  
  
  subroutine compute_gyroaverage_pre_compute_polar_hermite_c1(gyro,f)
    type(sll_plan_gyroaverage_polar)  :: gyro
    sll_real64,dimension(:,:),intent(inout) :: f
    sll_int32 ::i,j,k,ell,ii(2),s
    sll_real64::fval,sum_fval,eta_star(2),eta(2),delta_eta(2),x(2)
    
    fval=0._f64
    delta_eta(1)=(gyro%eta_max(1)-gyro%eta_min(1))/real(gyro%Nc(1),f64)
    delta_eta(2)=(gyro%eta_max(2)-gyro%eta_min(2))/real(gyro%Nc(2),f64)
    
    call hermite_c1_coef_nat_per(f(1:gyro%Nc(1)+1,1:gyro%Nc(2)),gyro%deriv,gyro%Nc,gyro%interp_degree)
    
    
    
    do j=1,gyro%Nc(2)
      s=0
      do i=1,gyro%Nc(1)+1
        fval=0._f64
        do k=1,gyro%pre_compute_N(i)
          s=s+1
          do ell=1,4
            ii(1)=gyro%pre_compute_index(1,s)
            ii(2)=modulo(j-1+gyro%pre_compute_index(2,s)+gyro%Nc(2),gyro%Nc(2))
            fval=fval+gyro%deriv(ell,ii(1)+1,ii(2)+1)*gyro%pre_compute_coeff(ell,s)
            !print *,k,ell,s,ii(1),ii(2),gyro%deriv(ell,ii(1)+1,ii(2)+1),gyro%pre_compute_coeff(ell,s)
          enddo  
        enddo
        f(i,j) = fval
        !print *,'#',i,j,fval,gyro%pre_compute_N(i)
        !stop
      enddo
      !stop
    enddo
    
    !print *,f
    !stop
    
    f(1:gyro%Nc(1)+1,gyro%Nc(2)+1)=f(1:gyro%Nc(1)+1,1)
    
  end subroutine compute_gyroaverage_pre_compute_polar_hermite_c1
  
  
  
  
  subroutine compute_gyroaverage_points_polar_spl(gyro,f,rho)
    type(sll_plan_gyroaverage_polar)  :: gyro
    sll_real64,dimension(:,:),intent(inout) :: f
    sll_real64,dimension(:,:),allocatable,target::fbords
    sll_real64,dimension(:,:),pointer::pointer_f_bords
    sll_real64,dimension(:),allocatable,target::buf,dnat,lnat,dper,lper,mper
    sll_real64,dimension(:),pointer::pointer_buf,pointer_dnat,pointer_lnat
    sll_real64,dimension(:),pointer::pointer_dper,pointer_lper,pointer_mper
    sll_real64,intent(in)::rho
    sll_int32 ::i,j,k
    sll_real64::fval,sum_fval,eta(2),delta_eta(2),x(2),xx(2)
    sll_int32 ::error
    
    fval=0._f64
    delta_eta(1)=(gyro%eta_max(1)-gyro%eta_min(1))/real(gyro%Nc(1),f64)
    delta_eta(2)=(gyro%eta_max(2)-gyro%eta_min(2))/real(gyro%Nc(2),f64)
   
    SLL_ALLOCATE(dnat(0:gyro%Nc(1)+2),error)
    SLL_ALLOCATE(lnat(0:gyro%Nc(1)+2),error)
    SLL_ALLOCATE(dper(0:gyro%Nc(2)-1),error)
    SLL_ALLOCATE(lper(0:gyro%Nc(2)-1),error)
    SLL_ALLOCATE(mper(0:gyro%Nc(2)-1),error)
   
    SLL_ALLOCATE(fbords(0:gyro%Nc(1)+2,0:gyro%Nc(2)-1),error)
    SLL_ALLOCATE(buf(0:max(gyro%Nc(1)+2,gyro%Nc(2)-1)),error)
   
    call splcoefnat1d0old(dnat,lnat,gyro%Nc(1))
    call splcoefper1d0old(dper,lper,mper,gyro%Nc(2))
   
    fbords(0,0:gyro%Nc(2)-1)=0._f64
    fbords(1:gyro%Nc(1)+1,0:gyro%Nc(2)-1)=f(1:gyro%Nc(1)+1,1:gyro%Nc(2))
    fbords(gyro%Nc(1)+2,0:gyro%Nc(2)-1)=0._f64
   
	pointer_f_bords => fbords
	pointer_buf => buf
	pointer_dnat => dnat
	pointer_lnat => lnat
	pointer_dper => dper
	pointer_lper => lper
	pointer_mper => mper

    call splcoefnatper2d(pointer_f_bords,pointer_buf,pointer_dnat,pointer_lnat,&
    		pointer_dper,pointer_lper,pointer_mper,gyro%Nc(1),gyro%Nc(2))
    
    
    do j=1,gyro%Nc(2)
     eta(2)=gyro%eta_min(2)+real(j-1,f64)*delta_eta(2)
     !eta(2)=gyro%eta_min(2) ! to uncomment for compatibility with precompute 
      do i=1,gyro%Nc(1)+1
        eta(1)=gyro%eta_min(1)+real(i-1,f64)*delta_eta(1)       
        sum_fval = 0._f64
        do k=1,gyro%N_points
          x(1) = eta(1)*cos(eta(2))+rho*gyro%points(1,k)
          x(2) = eta(1)*sin(eta(2))+rho*gyro%points(2,k)
          xx(1)=sqrt(x(1)**2+x(2)**2)
          xx(2)=modulo(datan2(x(2),x(1)),2._f64*sll_pi)
          !xx(2)=modulo(datan2(x(2),x(1))+real(j-1,f64)*delta_eta(2),2._f64*M_PI) & 
          		! to uncomment for compatibility with precompute
          call splnatper2d(pointer_f_bords,xx(1),gyro%eta_min(1),gyro%eta_max(1),&
				xx(2),gyro%eta_min(2),gyro%eta_max(2),fval,gyro%Nc(1),gyro%Nc(2))
          sum_fval = sum_fval+gyro%points(3,k)*fval
        enddo
        f(i,j) = sum_fval
      enddo
    enddo
    f(1:gyro%Nc(1)+1,gyro%Nc(2)+1)=f(1:gyro%Nc(1)+1,1)
   
   
    SLL_DEALLOCATE_ARRAY(dnat,error)
    SLL_DEALLOCATE_ARRAY(lnat,error)   
    SLL_DEALLOCATE_ARRAY(dper,error)
    SLL_DEALLOCATE_ARRAY(lper,error)
    SLL_DEALLOCATE_ARRAY(mper,error)
    
    SLL_DEALLOCATE_ARRAY(fbords,error)
    SLL_DEALLOCATE_ARRAY(buf,error) 
    
  end subroutine compute_gyroaverage_points_polar_spl
  
  
  
  subroutine compute_gyroaverage_points_polar_with_invar_spl(gyro,f,rho)
   type(sll_plan_gyroaverage_polar)  :: gyro
   sll_real64,dimension(:,:),intent(inout) :: f
   sll_real64,intent(in)::rho
   sll_int32 ::i,j,k
   sll_real64::fval,eta(2),delta_eta(2),x(2),xx(2),angle
   sll_real64,dimension(:,:),allocatable::buf,sum_fval
   sll_real64,dimension(:),allocatable::buf2
   sll_int32 ::error

   SLL_ALLOCATE(buf(0:gyro%Nc(1)+2,0:gyro%Nc(2)),error)
   SLL_ALLOCATE(buf2(0:gyro%Nc(2)-1),error)
   SLL_ALLOCATE(sum_fval(0:gyro%Nc(1),0:gyro%Nc(2)-1),error)
   SLL_ALLOCATE(gyro%lunat(0:1,0:gyro%Nc(1)+2),error)
   SLL_ALLOCATE(gyro%luper(0:3*gyro%Nc(2)-1),error)

   delta_eta(1)=(gyro%eta_max(1)-gyro%eta_min(1))/real(gyro%Nc(1),f64)
   delta_eta(2)=(gyro%eta_max(2)-gyro%eta_min(2))/real(gyro%Nc(2),f64)

   call splcoefnat1d0(gyro%lunat,gyro%Nc(1))
   call splcoefper1d0(gyro%luper,gyro%Nc(2))

   do j=0,gyro%Nc(2)
	  buf(0,j)=0._f64
	  do i=0,gyro%Nc(1)
	    buf(i+1,j)=f(i+1,j+1)
	  enddo
	  buf(gyro%Nc(1)+2,j)=0._f64
	enddo
	
   do j=0,gyro%Nc(2)-1
	  call splcoefnat1d(buf(:,j),gyro%lunat,gyro%Nc(1))
	enddo
	buf(0:gyro%Nc(1)+2,gyro%Nc(2))=buf(0:gyro%Nc(1)+2,0)
	
   !spline decomposition in theta
   do i=0,gyro%Nc(1)+2
     call splcoefper1d(buf(i,:),gyro%luper,gyro%Nc(2))
   enddo

   sum_fval = 0._f64

   do i=0,gyro%Nc(1)
     eta(1)=gyro%eta_min(1)+real(i,f64)*delta_eta(1)  ! rayon du centre
     eta(2)=gyro%eta_min(2) ! angle quelconque pour le centre
     do k=1,gyro%N_points
       x(1)=eta(1)*cos(eta(2))+rho*gyro%points(1,k)
       x(2)=eta(1)*sin(eta(2))+rho*gyro%points(2,k)
       xx(1)=sqrt(x(1)**2+x(2)**2)
       xx(2)=datan2(x(2),x(1))
		do j=0,gyro%Nc(2)-1
	      call splnat1d(buf(:,j),xx(1),gyro%eta_min(1),gyro%eta_max(1),buf2(j),gyro%Nc(1))
	    enddo
       !call splcoefper1d(buf2,gyro%luper,gyro%N(2))
       do j=0,gyro%Nc(2)-1
         angle=modulo(xx(2)+real(j,f64)*delta_eta(2),2._f64*sll_pi)
	      call splper1d(buf2,angle,gyro%eta_min(2),gyro%eta_max(2),fval,gyro%Nc(2))
	      sum_fval(i,j)=sum_fval(i,j)+gyro%points(3,k)*fval
	    enddo
	  enddo
	enddo  

   f(1:gyro%Nc(1)+1,1:gyro%Nc(2))=sum_fval(0:gyro%Nc(1),0:gyro%Nc(2)-1)
   f(1:gyro%Nc(1)+1,gyro%Nc(2)+1)=f(1:gyro%Nc(1)+1,1)  

   SLL_DEALLOCATE_ARRAY(gyro%lunat,error)
   SLL_DEALLOCATE_ARRAY(gyro%luper,error)
   SLL_DEALLOCATE_ARRAY(buf,error)
   SLL_DEALLOCATE_ARRAY(buf2,error)
   SLL_DEALLOCATE_ARRAY(sum_fval,error)

  end subroutine compute_gyroaverage_points_polar_with_invar_spl


 subroutine pre_compute_gyroaverage_polar_spl(gyro,rho)
    type(sll_plan_gyroaverage_polar)  :: gyro
    sll_real64,intent(in)::rho
    sll_int32,dimension(:,:),allocatable :: buf
    sll_int32 ::i,j,k,ell_1,ell_2,ii(2),s,nb,ind(2)
    sll_real64::val(-1:2,-1:2),eta_star(2),eta(2),delta_eta(2),x(2)
    sll_int32 ::error
    
    
    delta_eta(1)=(gyro%eta_max(1)-gyro%eta_min(1))/real(gyro%Nc(1),f64)
    delta_eta(2)=(gyro%eta_max(2)-gyro%eta_min(2))/real(gyro%Nc(2),f64)
    
    SLL_ALLOCATE(buf(0:gyro%Nc(1)+2,0:gyro%Nc(2)-1),error)
       
    SLL_ALLOCATE(gyro%pre_compute_N(gyro%Nc(1)+1),error)

    eta(2)=gyro%eta_min(2)
    buf=0
    nb=0
    do i=1,gyro%Nc(1)+1
      eta(1)=gyro%eta_min(1)+real(i-1,f64)*delta_eta(1)       
      s=0
      do k=1,gyro%N_points
        x(1) = eta(1)*cos(eta(2))+rho*gyro%points(1,k)
        x(2) = eta(1)*sin(eta(2))+rho*gyro%points(2,k)
        call localize_polar(x,gyro%eta_min,gyro%eta_max,ii,eta_star,gyro%Nc)
        do ell_2=-1,2
          ind(2)=modulo(ii(2)+ell_2,gyro%Nc(2))
          do ell_1=-1,2
            ind(1)=ii(1)+1+ell_1
            if(buf(ind(1),ind(2)).ne.i)then
              s=s+1
              buf(ind(1),ind(2))=i
            endif
          enddo
        enddo    
      enddo
      gyro%pre_compute_N(i)=s
      nb=nb+s
    enddo
    
    print*,'#N_points pre_compute=',nb,gyro%Nc(1),nb/gyro%Nc(1)

    SLL_ALLOCATE(gyro%pre_compute_index(1:2,nb),error)
    SLL_ALLOCATE(gyro%pre_compute_coeff_spl(nb),error)
    
    buf=0
    nb=0
    s=0
    val=0._f64
    eta(2)=gyro%eta_min(2)
    do i=1,gyro%Nc(1)+1
      eta(1)=gyro%eta_min(1)+real(i-1,f64)*delta_eta(1)       
      do k=1,gyro%N_points
        x(1) = eta(1)*cos(eta(2))+rho*gyro%points(1,k)
        x(2) = eta(1)*sin(eta(2))+rho*gyro%points(2,k)
        call localize_polar(x,gyro%eta_min,gyro%eta_max,ii,eta_star,gyro%Nc)

        call contribution_spl(eta_star,val)
        
        val=val*gyro%points(3,k)
        
        do ell_2=-1,2
          ind(2)=modulo(ii(2)+ell_2,gyro%Nc(2))
          do ell_1=-1,2
            ind(1)=ii(1)+1+ell_1
            j=buf(ind(1),ind(2))
            if(j<=nb)then
              s=s+1
              buf(ind(1),ind(2))=s
              gyro%pre_compute_coeff_spl(s)=val(ell_1,ell_2)
              gyro%pre_compute_index(1,s)=ind(1)
              gyro%pre_compute_index(2,s)=ind(2)                     
            else              
              gyro%pre_compute_coeff_spl(j)=gyro%pre_compute_coeff_spl(j)+val(ell_1,ell_2)
            endif
          enddo
        enddo
      enddo      
      nb=s 
    enddo
    
    print *,'#min/max=',minval(gyro%pre_compute_coeff_spl(:)),maxval(gyro%pre_compute_coeff_spl(:))
    SLL_DEALLOCATE_ARRAY(buf,error)
  end subroutine pre_compute_gyroaverage_polar_spl



  subroutine compute_gyroaverage_pre_compute_polar_spl(gyro,f)
    type(sll_plan_gyroaverage_polar)  :: gyro
    sll_real64,dimension(:,:),intent(inout) :: f
    sll_real64,dimension(:,:),allocatable,target::fbords
    sll_real64,dimension(:,:),pointer::pointer_f_bords
    sll_real64,dimension(:),allocatable,target::buf,dnat,lnat,dper,lper,mper
    sll_real64,dimension(:),pointer::pointer_buf,pointer_dnat,pointer_lnat
    sll_real64,dimension(:),pointer::pointer_dper,pointer_lper,pointer_mper
    sll_real64,dimension(:,:),allocatable ::tmp_f
    sll_int32 ::i,j,k,ell,ii(2),s
    sll_real64::fval,sum_fval,eta_star(2),eta(2),delta_eta(2),x(2)
    sll_int32 ::error
    
    SLL_ALLOCATE(dnat(0:gyro%Nc(1)+2),error)
    SLL_ALLOCATE(lnat(0:gyro%Nc(1)+2),error)
    SLL_ALLOCATE(dper(0:gyro%Nc(2)-1),error)
    SLL_ALLOCATE(lper(0:gyro%Nc(2)-1),error)
    SLL_ALLOCATE(mper(0:gyro%Nc(2)-1),error)
    SLL_ALLOCATE(fbords(0:gyro%Nc(1)+2,0:gyro%Nc(2)-1),error)
    SLL_ALLOCATE(buf(0:max(gyro%Nc(1)+2,gyro%Nc(2)-1)),error)
    SLL_ALLOCATE(tmp_f(1:gyro%Nc(1)+1,1:gyro%Nc(2)),error)
    
    fval=0._f64
    delta_eta(1)=(gyro%eta_max(1)-gyro%eta_min(1))/real(gyro%Nc(1),f64)
    delta_eta(2)=(gyro%eta_max(2)-gyro%eta_min(2))/real(gyro%Nc(2),f64)    
   
    call splcoefnat1d0old(dnat,lnat,gyro%Nc(1))
    call splcoefper1d0old(dper,lper,mper,gyro%Nc(2))
   
    fbords(0,0:gyro%Nc(2)-1)=0._f64
    fbords(1:gyro%Nc(1)+1,0:gyro%Nc(2)-1)=f(1:gyro%Nc(1)+1,1:gyro%Nc(2))
    fbords(gyro%Nc(1)+2,0:gyro%Nc(2)-1)=0._f64
   
	pointer_f_bords => fbords
	pointer_buf => buf
	pointer_dnat => dnat
	pointer_lnat => lnat
	pointer_dper => dper
	pointer_lper => lper
	pointer_mper => mper

    call splcoefnatper2d(pointer_f_bords,pointer_buf,pointer_dnat,pointer_lnat,&
    		pointer_dper,pointer_lper,pointer_mper,gyro%Nc(1),gyro%Nc(2))    
    
    do j=1,gyro%Nc(2)
      s=0
      do i=1,gyro%Nc(1)+1
        fval=0._f64
        do k=1,gyro%pre_compute_N(i)
          s=s+1
          ii(1)=gyro%pre_compute_index(1,s)
          ii(2)=modulo(j-1+gyro%pre_compute_index(2,s),gyro%Nc(2))
          fval=fval+pointer_f_bords(ii(1),ii(2))*gyro%pre_compute_coeff_spl(s) 
        enddo
        tmp_f(i,j) = fval
      enddo
    enddo 
    
    f(1:gyro%Nc(1)+1,1:gyro%Nc(2))=tmp_f(1:gyro%Nc(1)+1,1:gyro%Nc(2))
    f(1:gyro%Nc(1)+1,gyro%Nc(2)+1)=f(1:gyro%Nc(1)+1,1)
       
    SLL_DEALLOCATE_ARRAY(dnat,error)
    SLL_DEALLOCATE_ARRAY(lnat,error)
    SLL_DEALLOCATE_ARRAY(dper,error)
    SLL_DEALLOCATE_ARRAY(lper,error)
    SLL_DEALLOCATE_ARRAY(mper,error)
    SLL_DEALLOCATE_ARRAY(fbords,error)
    SLL_DEALLOCATE_ARRAY(buf,error)
    SLL_DEALLOCATE_ARRAY(tmp_f,error)
    
  end subroutine compute_gyroaverage_pre_compute_polar_spl


subroutine pre_compute_gyroaverage_polar_spl_FFT(gyro,rho)
    type(sll_plan_gyroaverage_polar)  :: gyro
    sll_real64,intent(in)::rho
    sll_int32,dimension(:,:),allocatable :: buf
    sll_real64,dimension(:),allocatable::buf_fft
    sll_int32 ::i,j,k,ell_1,ell_2,ii(2),s,nb,ind(2),i1,i2
    sll_real64::val(-1:2,-1:2),eta_star(2),eta(2),delta_eta(2),x(2)
    sll_int32 ::error
      
    delta_eta(1)=(gyro%eta_max(1)-gyro%eta_min(1))/real(gyro%Nc(1),f64)
    delta_eta(2)=(gyro%eta_max(2)-gyro%eta_min(2))/real(gyro%Nc(2),f64)
  
    SLL_ALLOCATE(buf(1:gyro%Nc(1)+1,0:gyro%Nc(1)+2),error)
    
    SLL_ALLOCATE(gyro%A_fft(1:gyro%Nc(1)+1,0:gyro%Nc(1)+2,0:gyro%Nc(2)-1),error)
    SLL_ALLOCATE(buf_fft(1:2*gyro%Nc(2)+15),error)


    eta(2)=gyro%eta_min(2)
    buf=0
    nb=0
    do i=1,gyro%Nc(1)+1
      eta(1)=gyro%eta_min(1)+real(i-1,f64)*delta_eta(1)       
      s=0
      do k=1,gyro%N_points
        x(1) = eta(1)*cos(eta(2))+rho*gyro%points(1,k)
        x(2) = eta(1)*sin(eta(2))+rho*gyro%points(2,k)
        call localize_polar(x,gyro%eta_min,gyro%eta_max,ii,eta_star,gyro%Nc)
        do ell_2=-1,2
          ind(2)=modulo(ii(2)+ell_2,gyro%Nc(2))
          do ell_1=-1,2
            ind(1)=ii(1)+1+ell_1
            if(buf(i,ind(1)).ne.1)then
              s=s+1
              buf(i,ind(1))=1
            endif
          enddo
        enddo    
      enddo
      nb=nb+s
    enddo

    SLL_ALLOCATE(gyro%pre_compute_index(1:2,nb),error)
    
    gyro%A_fft=0._f64
    buf=0
    s=0
    val=0._f64
    eta(2)=gyro%eta_min(2)
    do i=1,gyro%Nc(1)+1
      eta(1)=gyro%eta_min(1)+real(i-1,f64)*delta_eta(1)       
      do k=1,gyro%N_points
        x(1) = eta(1)*cos(eta(2))+rho*gyro%points(1,k)
        x(2) = eta(1)*sin(eta(2))+rho*gyro%points(2,k)
        call localize_polar(x,gyro%eta_min,gyro%eta_max,ii,eta_star,gyro%Nc)

        call contribution_spl(eta_star,val)
        
        val=val*gyro%points(3,k)
        
        do ell_2=-1,2
          ind(2)=modulo(ii(2)+ell_2,gyro%Nc(2))
          do ell_1=-1,2
            ind(1)=ii(1)+1+ell_1
            j=buf(i,ind(1))
			if(j==0)then
              s=s+1
              buf(i,ind(1))=1
              gyro%pre_compute_index(1,s)=i
              gyro%pre_compute_index(2,s)=ind(1)    
              gyro%A_fft(i,ind(1),ind(2))=val(ell_1,ell_2)                 
            else
           	  gyro%A_fft(i,ind(1),ind(2))=gyro%A_fft(i,ind(1),ind(2))+val(ell_1,ell_2) 
            endif
          enddo
        enddo
      enddo      
    enddo
    		

	call dffti(gyro%Nc(2),buf_fft)
	do s=1,nb
		call dfftf(gyro%Nc(2),gyro%A_fft(gyro%pre_compute_index(1,s),gyro%pre_compute_index(2,s),:),buf_fft)
	enddo

    SLL_DEALLOCATE_ARRAY(buf,error)
    SLL_DEALLOCATE_ARRAY(buf_fft,error)
    
  end subroutine pre_compute_gyroaverage_polar_spl_FFT


subroutine compute_gyroaverage_pre_compute_polar_spl_FFT(gyro,f)
    type(sll_plan_gyroaverage_polar)  :: gyro
    sll_real64,dimension(:,:),intent(inout) :: f
    sll_real64,dimension(:,:),allocatable,target::fbords
    sll_real64,dimension(:,:),pointer::pointer_f_bords
    sll_real64,dimension(:),allocatable,target::buf,dnat,lnat,dper,lper,mper
    sll_real64,dimension(:),pointer::pointer_buf,pointer_dnat,pointer_lnat
    sll_real64,dimension(:),pointer::pointer_dper,pointer_lper,pointer_mper
    sll_real64,dimension(:,:),allocatable ::tmp_f
    sll_real64,dimension(:),allocatable::buf_fft
    sll_int32 ::i,j,k,ell,ii2,s
    sll_real64::fval,sum_fval,eta_star(2),eta(2),delta_eta(2),x(2)
    sll_int32 ::error
    
    SLL_ALLOCATE(dnat(0:gyro%Nc(1)+2),error)
    SLL_ALLOCATE(lnat(0:gyro%Nc(1)+2),error)
    SLL_ALLOCATE(dper(0:gyro%Nc(2)-1),error)
    SLL_ALLOCATE(lper(0:gyro%Nc(2)-1),error)
    SLL_ALLOCATE(mper(0:gyro%Nc(2)-1),error)
    SLL_ALLOCATE(fbords(0:gyro%Nc(1)+2,0:gyro%Nc(2)-1),error)
    SLL_ALLOCATE(buf(0:max(gyro%Nc(1)+2,gyro%Nc(2)-1)),error)
    SLL_ALLOCATE(tmp_f(1:gyro%Nc(1)+1,0:gyro%Nc(2)-1),error)
    SLL_ALLOCATE(buf_fft(1:2*gyro%Nc(2)+15),error)

    
    fval=0._f64
    tmp_f=0._f64
    delta_eta(1)=(gyro%eta_max(1)-gyro%eta_min(1))/real(gyro%Nc(1),f64)
    delta_eta(2)=(gyro%eta_max(2)-gyro%eta_min(2))/real(gyro%Nc(2),f64)    
   
    call splcoefnat1d0old(dnat,lnat,gyro%Nc(1))
    call splcoefper1d0old(dper,lper,mper,gyro%Nc(2))
   
    fbords(0,0:gyro%Nc(2)-1)=0._f64
    fbords(1:gyro%Nc(1)+1,0:gyro%Nc(2)-1)=f(1:gyro%Nc(1)+1,1:gyro%Nc(2))
    fbords(gyro%Nc(1)+2,0:gyro%Nc(2)-1)=0._f64
   
	pointer_f_bords => fbords
	pointer_buf => buf
	pointer_dnat => dnat
	pointer_lnat => lnat
	pointer_dper => dper
	pointer_lper => lper
	pointer_mper => mper

    call splcoefnatper2d(pointer_f_bords,pointer_buf,pointer_dnat,pointer_lnat,&
    		pointer_dper,pointer_lper,pointer_mper,gyro%Nc(1),gyro%Nc(2))    

 
 	call dffti(gyro%Nc(2),buf_fft)
	do i=0,gyro%Nc(1)+2
		call dfftf(gyro%Nc(2),pointer_f_bords(i,:),buf_fft)
	enddo

	do j=0,gyro%Nc(2)-1
		do s=1,size(gyro%pre_compute_index(1,:))
			tmp_f(gyro%pre_compute_index(1,s),j)=tmp_f(gyro%pre_compute_index(1,s),j)+&
			gyro%A_fft(gyro%pre_compute_index(1,s),&
			gyro%pre_compute_index(2,s),j)*pointer_f_bords(gyro%pre_compute_index(2,s),j)
		enddo
	enddo

	do i=1,gyro%Nc(1)+1
		call dfftb(gyro%Nc(2),tmp_f(i,:),buf_fft)
	enddo
	
	tmp_f=tmp_f/real(gyro%Nc(2),f64)
 
    f(1:gyro%Nc(1)+1,1:gyro%Nc(2))=tmp_f(1:gyro%Nc(1)+1,0:gyro%Nc(2)-1)
    f(1:gyro%Nc(1)+1,gyro%Nc(2)+1)=f(1:gyro%Nc(1)+1,1)
       
    SLL_DEALLOCATE_ARRAY(dnat,error)
    SLL_DEALLOCATE_ARRAY(lnat,error)
    SLL_DEALLOCATE_ARRAY(dper,error)
    SLL_DEALLOCATE_ARRAY(lper,error)
    SLL_DEALLOCATE_ARRAY(mper,error)
    SLL_DEALLOCATE_ARRAY(fbords,error)
    SLL_DEALLOCATE_ARRAY(buf,error)
    SLL_DEALLOCATE_ARRAY(buf_fft,error)
    SLL_DEALLOCATE_ARRAY(tmp_f,error)
    
  end subroutine compute_gyroaverage_pre_compute_polar_spl_FFT


 subroutine compute_gyroaverage_pade_polar(gyro,f,rho)
    type(sll_plan_gyroaverage_polar)  :: gyro
    sll_real64,dimension(:,:),intent(inout) :: f
    sll_real64,dimension(:,:),allocatable :: fcomp
    sll_real64,dimension(:),allocatable :: buf,diagm1,diag,diagp1
    sll_real64,intent(in)::rho
    sll_int32 ::i,j,k
    sll_real64::sum_fval,eta_star(2),eta(2),x(2),dr
    sll_int32 ::error
    
    dr=(gyro%eta_max(1)-gyro%eta_min(1))/real(gyro%Nc(1),f64)

	SLL_ALLOCATE(buf(1:2*gyro%Nc(2)+15),error)
	SLL_ALLOCATE(fcomp(1:gyro%Nc(1)+1,1:gyro%Nc(2)),error)
	SLL_ALLOCATE(diagm1(1:gyro%Nc(1)+1),error)
	SLL_ALLOCATE(diag(1:gyro%Nc(1)+1),error)
	SLL_ALLOCATE(diagp1(1:gyro%Nc(1)+1),error)
	

	fcomp(1:gyro%Nc(1)+1,1:gyro%Nc(2))=f(1:gyro%Nc(1)+1,1:gyro%Nc(2))
	

    !*** Perform FFT 1D in theta direction of ***
    !***   the system solution                ***
	call dffti(gyro%Nc(2),buf)
	do i=1,gyro%Nc(1)+1
		call dfftf(gyro%Nc(2),fcomp(i,:),buf)
	enddo
	fcomp=fcomp/real(gyro%Nc(2),f64)

 	!***POISSON
	do k=1,gyro%Nc(2)
	  do i=1,gyro%Nc(1)
	    diagm1(i+1)=-(rho**2/4)*(1/dr**2-1/(2*dr*(gyro%eta_min(1)+ &
	    	(gyro%eta_max(1)-gyro%eta_min(1))*real(i,f64)/real(gyro%Nc(1),f64))))
	    diag(i)=1-(rho**2/4)*(-(2/dr**2)-((floor(k/2._f64)*1._f64)/ &
	    	(gyro%eta_min(1)+(gyro%eta_max(1)-gyro%eta_min(1))*real(i-1,f64)/real(gyro%Nc(1),f64)))**2)
	    diagp1(i)=-(rho**2/4)*(1/dr**2+1/(2*dr*(gyro%eta_min(1)+ &
	    	(gyro%eta_max(1)-gyro%eta_min(1))*real(i-1,f64)/real(gyro%Nc(1),f64))))
	  enddo
	  diagm1(1)=0._f64
	  diagp1(gyro%Nc(1)+1)=0._f64
	  diag(1)=1._f64
	  diag(gyro%Nc(1)+1)=1._f64
	  !***  Dirichlet boundary conditions ***	  
	  diagp1(1)=0._f64
	  diagm1(gyro%Nc(1)+1)=0._f64
	  !***  Neumann boundary conditions ***
!	  diagp1(1)=-1._f64
!	  diagm1(gyro%Nc(1)+1)=-1._f64
	  call solve_tridiag(diagm1,diag,diagp1,fcomp(1:gyro%Nc(1)+1,k),f(1:gyro%Nc(1)+1,k),gyro%Nc(1)+1)
	enddo

	!*** Perform FFT 1D inverse ***
	do i=1,gyro%Nc(1)+1
	  call dfftb(gyro%Nc(2),f(i+1,1:gyro%Nc(2)),buf)
	enddo
         
    !*** duplicate periodic value ***
    f(1:gyro%Nc(1)+1,gyro%Nc(2)+1)=f(1:gyro%Nc(1)+1,1)
    
  end subroutine compute_gyroaverage_pade_polar











subroutine hermite_coef_nat_per(f,buf3d,N,d)
    sll_int32,intent(in)::N(2),d(2)
    sll_real64,dimension(N(1)+1,N(2)),intent(in)::f
    sll_real64,dimension(9,N(1)+1,N(2)+1),intent(out)::buf3d
    sll_real64 ::w_left_1(-d(1)/2:(d(1)+1)/2),w_right_1((-d(1)+1)/2:d(1)/2+1)
    sll_real64 ::w_left_2(-d(2)/2:(d(2)+1)/2),w_right_2((-d(2)+1)/2:d(2)/2+1)
    sll_real64 ::tmp
    sll_int32  ::i,j,r,s,ii,r_left(2),r_right(2),s_left(2),s_right(2),ind 
    r_left=-d/2
    s_left=(d+1)/2
    r_right=(-d+1)/2
    s_right=d/2+1
    
    
    call compute_w_hermite(w_left_1,r_left(1),s_left(1))
    call compute_w_hermite(w_left_2,r_left(2),s_left(2))    
    if((2*(d(1)/2)-d(1))==0)then
      w_right_1(r_right(1):s_right(1)) = w_left_1(r_left(1):s_left(1))
    else
      w_right_1(r_right(1):s_right(1)) = -w_left_1(s_left(1):r_left(1):-1)
    endif    

    if((2*(d(2)/2)-d(2))==0)then
      w_right_2(r_right(2):s_right(2)) = w_left_2(r_left(2):s_left(2))
    else
      w_right_2(r_right(2):s_right(2)) = -w_left_2(s_left(2):r_left(2):-1)
    endif    
    
    !print *,'w(',r_left(1),':',s_left(1),')=',w_left_1(r_left(1):s_left(1))
    !print *,'w(',r_right(1),':',s_right(1),')=',w_right_1(r_right(1):s_right(1))

    
    do j=1,N(2)
      do i=1,N(1)+1
        buf3d(1,i,j)=f(i,j) !f(0,0)
        tmp=0._f64
        do ii=r_left(1),s_left(1)
          ind=i+ii;if(ind>N(1)+1) ind=2*(N(1)+1)-ind;if(ind<1) ind=2-ind
          tmp=tmp+w_left_1(ii)*f(ind,j)
        enddo
        buf3d(2,i,j)=tmp !fx(0,0)
        tmp=0._f64
        do ii=r_right(1),s_right(1)
          ind=i+ii;if(ind>N(1)+1) ind=2*(N(1)+1)-ind;if(ind<1) ind=2-ind
          tmp=tmp+w_right_1(ii)*f(ind,j)
        enddo
        buf3d(3,i,j)=tmp !fx(1,0)       
      enddo
    enddo
    do i=1,N(1)+1
      do j=1,N(2)
        tmp=0._f64
        do ii=r_left(2),s_left(2)
          ind=modulo(j+ii-1+N(2),N(2))+1
          tmp=tmp+w_left_2(ii)*f(i,ind)
        enddo
        buf3d(4,i,j)=tmp !fy(0,0)
        tmp=0._f64
        do ii=r_right(2),s_right(2)
          ind=modulo(j+ii-1+N(2),N(2))+1
          tmp=tmp+w_right_2(ii)*f(i,ind)
        enddo
        buf3d(5,i,j)=tmp !fy(0,1)               
      enddo
    enddo

    do j=1,N(2)
      do i=1,N(1)+1
        tmp=0._f64
        do ii=r_left(1),s_left(1)
          ind=i+ii;if(ind>N(1)+1) ind=2*(N(1)+1)-ind;if(ind<1) ind=2-ind
          tmp=tmp+w_left_1(ii)*buf3d(4,ind,j)
        enddo
        buf3d(6,i,j)=tmp !fxy(0,0)
        tmp=0._f64
        do ii=r_right(1),s_right(1)
          ind=i+ii;if(ind>N(1)+1) ind=2*(N(1)+1)-ind;if(ind<1) ind=2-ind
          tmp=tmp+w_right_1(ii)*buf3d(4,ind,j)
        enddo
        buf3d(7,i,j)=tmp !fxy(1,0)       
        tmp=0._f64
        do ii=r_left(1),s_left(1)
          ind=i+ii;if(ind>N(1)+1) ind=2*(N(1)+1)-ind;if(ind<1) ind=2-ind
          tmp=tmp+w_left_1(ii)*buf3d(5,ind,j)
        enddo
        buf3d(8,i,j)=tmp  !fxy(0,1)
        tmp=0._f64
        do ii=r_right(1),s_right(1)
          ind=i+ii;if(ind>N(1)+1) ind=2*(N(1)+1)-ind;if(ind<1) ind=2-ind
          tmp=tmp+w_right_1(ii)*buf3d(5,ind,j)
        enddo
        buf3d(9,i,j)=tmp !fxy(1,1)        
      enddo
    enddo

    buf3d(:,:,N(2)+1)=buf3d(:,:,1)
    
  end subroutine hermite_coef_nat_per



  subroutine hermite_c1_coef_nat_per(f,buf3d,N,d)
    integer,intent(in)::N(2),d(2)
    sll_real64,dimension(N(1)+1,N(2)),intent(in)::f
    sll_real64,dimension(4,N(1)+1,N(2)+1),intent(out)::buf3d
    sll_real64,dimension(:),allocatable ::w_left_1,w_left_2
    sll_real64 ::tmp
    sll_int32::i,j,r,s,ii,r_left(2),s_left(2),ind,dd(2) 
    sll_int32::error
    dd(1)=2*((d(1)+1)/2)
    dd(2)=2*((d(2)+1)/2)
    
    r_left=-dd/2
    s_left=(dd+1)/2
    
    SLL_ALLOCATE(w_left_1(-dd(1)/2:dd(1)/2),error)
    SLL_ALLOCATE(w_left_2(-dd(2)/2:dd(2)/2),error)
    
    call compute_w_hermite(w_left_1,r_left(1),s_left(1))
    call compute_w_hermite(w_left_2,r_left(2),s_left(2))    
    
    !print *,'w(',r_left(1),':',s_left(1),')=',w_left_1(r_left(1):s_left(1))
    !print *,'w(',r_right(1),':',s_right(1),')=',w_right_1(r_right(1):s_right(1))

    
    do j=1,N(2)
      do i=1,N(1)+1
        buf3d(1,i,j)=f(i,j) !f(0,0)
        tmp=0._f64
        do ii=r_left(1),s_left(1)
          ind=i+ii;if(ind>N(1)+1) ind=2*(N(1)+1)-ind;if(ind<1) ind=2-ind
          tmp=tmp+w_left_1(ii)*f(ind,j)
        enddo
        buf3d(2,i,j)=tmp !fx(0,0)
      enddo
    enddo
    do i=1,N(1)+1
      do j=1,N(2)
        tmp=0._f64
        do ii=r_left(2),s_left(2)
          ind=modulo(j+ii-1+N(2),N(2))+1
          tmp=tmp+w_left_2(ii)*f(i,ind)
        enddo
        buf3d(3,i,j)=tmp !fy(0,0)
      enddo
    enddo

    do j=1,N(2)
      do i=1,N(1)+1
        tmp=0._f64
        do ii=r_left(1),s_left(1)
          ind=i+ii;if(ind>N(1)+1) ind=2*(N(1)+1)-ind;if(ind<1) ind=2-ind
          tmp=tmp+w_left_1(ii)*buf3d(3,ind,j)
        enddo
        buf3d(4,i,j)=tmp !fxy(0,0)
      enddo
    enddo

    buf3d(:,:,N(2)+1)=buf3d(:,:,1)
    
	SLL_DEALLOCATE_ARRAY(w_left_1,error)
	SLL_DEALLOCATE_ARRAY(w_left_2,error)
    
  end subroutine hermite_c1_coef_nat_per





subroutine compute_w_hermite(w,r,s)
    sll_int32,intent(in)::r,s
    sll_real64,dimension(r:s),intent(out)::w
    sll_int32 ::i,j
    sll_real64::tmp

    !maple code for generation of w
    !for k from r to -1 do
    !  C[k]:=product((k-j),j=r..k-1)*product((k-j),j=k+1..s):
    !  C[k]:=1/C[k]*product((-j),j=r..k-1)*product((-j),j=k+1..-1)*product((-j),j=1..s):
    !od:
    !for k from 1 to s do
    !  C[k]:=product((k-j),j=r..k-1)*product((k-j),j=k+1..s):
    !  C[k]:=1/C[k]*product((-j),j=r..-1)*product((-j),j=1..k-1)*product((-j),j=k+1..s):
    !od:
    !C[0]:=-add(C[k],k=r..-1)-add(C[k],k=1..s):
    
    do i=r,-1
      tmp=1._f64
      do j=r,i-1
        tmp=tmp*real(i-j,f64)
      enddo
      do j=i+1,s
        tmp=tmp*real(i-j,f64)
      enddo
      tmp=1._f64/tmp
      do j=r,i-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=i+1,-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=1,s
        tmp=tmp*real(-j,f64)
      enddo
      w(i)=tmp      
    enddo


    do i=1,s
      tmp=1._f64
      do j=r,i-1
        tmp=tmp*real(i-j,f64)
      enddo
      do j=i+1,s
        tmp=tmp*real(i-j,f64)
      enddo
      tmp=1._f64/tmp
      do j=r,-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=1,i-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=i+1,s
        tmp=tmp*real(-j,f64)
      enddo
      w(i)=tmp      
    enddo
    tmp=0._f64
    do i=r,-1
      tmp=tmp+w(i)
    enddo
    do i=1,s
      tmp=tmp+w(i)
    enddo
    w(0)=-tmp

    !print *,'w',w
    !do ii=r,s
    !  print *,ii,w(r+s-ii)
    !enddo
    !

  
  end subroutine compute_w_hermite


  subroutine localize_polar(x,eta_min,eta_max,ii,eta,N)
    sll_real64,intent(in)::x(2),eta_min(2),eta_max(2)
    sll_int32,intent(out)::ii(2)
    sll_int32,intent(in)::N(2)
    sll_real64,intent(out)::eta(2)
    
    eta(1)=sqrt(x(1)**2+x(2)**2)
    call localize_nat(ii(1),eta(1),eta_min(1),eta_max(1),N(1))
    eta(2)=datan2(x(2),x(1))
    call localize_per(ii(2),eta(2),eta_min(2),eta_max(2),N(2))
  end subroutine localize_polar

  subroutine localize_per(i,x,xmin,xmax,N)
    sll_int32,intent(out)::i
    sll_real64,intent(inout)::x
    sll_real64,intent(in)::xmin,xmax
    sll_int32,intent(in)::N
    x=(x-xmin)/(xmax-xmin)
    x=x-real(floor(x),f64)
    x=x*real(N,f64)
    i=floor(x)
    x=x-real(i,f64)
    if(i==N)then
      i=0
      x=0._f64
    endif
  end subroutine localize_per
  

  subroutine localize_nat(i,x,xmin,xmax,N)
    sll_int32,intent(out)::i
    sll_real64,intent(inout)::x
    sll_real64,intent(in)::xmin,xmax
    sll_int32,intent(in)::N
    x=(x-xmin)/(xmax-xmin)
    x=x*real(N,f64)
    if(x>=real(N,f64))then
      x=real(N,f64)
    endif
    if(x<=0._f64)then
      x=0._f64
    endif    
    i=floor(x)
    x=x-real(i,f64)
    if(i==N)then
      i=N-1
      x=1._f64
    endif
  end subroutine localize_nat


 subroutine interpolate_hermite(f,i,x,fval,N)
    sll_int32,intent(in)::i(2),N(2)
    !real(f64),intent(in)::xx(2),xmin(2),xmax(2)
    sll_real64,intent(in)::x(2)
    sll_real64,intent(out)::fval
    sll_real64,dimension(0:8,0:N(1),0:N(2))::f
    !integer::i(2),i1(2),s
    sll_int32::i1(2),s
    sll_real64::w(2,0:3),tmp(0:3)
    sll_real64::g(0:3,0:3)
    
    !fval =f(0,i(1),i(2))!real(i(1),f64)
    
    !return
    
    do s=1,2
      w(s,0)=(2._f64*x(s)+1)*(1._f64-x(s))*(1._f64-x(s));
      w(s,1)=x(s)*x(s)*(3._f64-2._f64*x(s))
      w(s,2)=x(s)*(1._f64-x(s))*(1._f64-x(s))
      w(s,3)=x(s)*x(s)*(x(s)-1._f64)
      i1(s)=i(s)+1
    enddo

    
    g(0,0)=f(0,i(1),i(2))          !f(0,0)
    g(1,0)=f(0,i1(1),i(2))         !f(1,0)
    g(2,0)=f(1,i(1),i(2))          !fx(0,0)
    g(3,0)=f(2,i(1),i(2))          !fx(1,0)
    g(0,1)=f(0,i(1),i1(2))         !f(0,1) 
    g(1,1)=f(0,i1(1),i1(2))        !f(1,1)
    g(2,1)=f(1,i(1),i1(2))         !fx(0,1)
    g(3,1)=f(2,i(1),i1(2))         !fx(1,1)
    g(0,2)=f(3,i(1),i(2))          !fy(0,0)
    g(1,2)=f(3,i1(1),i(2))         !fy(1,0)
    g(2,2)=f(5,i(1),i(2))          !fxy(0,0)
    g(3,2)=f(6,i(1),i(2))          !fxy(1,0)
    g(0,3)=f(4,i(1),i(2))          !fy(0,1) 
    g(1,3)=f(4,i1(1),i(2))         !fy(1,1)
    g(2,3)=f(7,i(1),i(2))          !fxy(0,1) 
    g(3,3)=f(8,i(1),i(2))          !fxy(1,1)



    do s=0,3
      tmp(s)=w(1,0)*g(0,s)+w(1,1)*g(1,s)+w(1,2)*g(2,s)+w(1,3)*g(3,s)
    enddo  

    fval=w(2,0)*tmp(0)+w(2,1)*tmp(1)+w(2,2)*tmp(2)+w(2,3)*tmp(3)
  
    !print *,fval,' t',f    
  end subroutine interpolate_hermite


  subroutine interpolate_hermite_c1(f,i,x,fval,N)
    sll_int32,intent(in)::i(2),N(2)
    !real(f64),intent(in)::xx(2),xmin(2),xmax(2)
    sll_real64,intent(in)::x(2)
    sll_real64,intent(out)::fval
    sll_real64,dimension(0:3,0:N(1),0:N(2))::f
    !integer::i(2),i1(2),s
    sll_int32::i1(2),s
    sll_real64::w(2,0:3),tmp(0:3)
    sll_real64::g(0:3,0:3)
    
    !fval =f(0,i(1),i(2))!real(i(1),f64)
    
    !return
    
    do s=1,2
      w(s,0)=(2._f64*x(s)+1)*(1._f64-x(s))*(1._f64-x(s));
      w(s,1)=x(s)*x(s)*(3._f64-2._f64*x(s))
      w(s,2)=x(s)*(1._f64-x(s))*(1._f64-x(s))
      w(s,3)=x(s)*x(s)*(x(s)-1._f64)
      i1(s)=i(s)+1
    enddo
    
    g(0,0)=f(0,i(1),i(2))          !f(0,0)
    g(1,0)=f(0,i1(1),i(2))         !f(1,0)
    g(2,0)=f(1,i(1),i(2))          !fx(0,0)
    g(3,0)=f(1,i1(1),i(2))          !fx(1,0)
    g(0,1)=f(0,i(1),i1(2))         !f(0,1) 
    g(1,1)=f(0,i1(1),i1(2))        !f(1,1)
    g(2,1)=f(1,i(1),i1(2))         !fx(0,1)
    g(3,1)=f(1,i1(1),i1(2))         !fx(1,1)
    g(0,2)=f(2,i(1),i(2))          !fy(0,0)
    g(1,2)=f(2,i1(1),i(2))         !fy(1,0)
    g(2,2)=f(3,i(1),i(2))          !fxy(0,0)
    g(3,2)=f(3,i1(1),i(2))          !fxy(1,0)
    g(0,3)=f(2,i(1),i1(2))          !fy(0,1) 
    g(1,3)=f(2,i1(1),i1(2))         !fy(1,1)
    g(2,3)=f(3,i(1),i1(2))          !fxy(0,1) 
    g(3,3)=f(3,i1(1),i1(2))          !fxy(1,1)

    do s=0,3
      tmp(s)=w(1,0)*g(0,s)+w(1,1)*g(1,s)+w(1,2)*g(2,s)+w(1,3)*g(3,s)
    enddo  
    
    fval=w(2,0)*tmp(0)+w(2,1)*tmp(1)+w(2,2)*tmp(2)+w(2,3)*tmp(3)
    !print *,fval,' t',f    
  end subroutine interpolate_hermite_c1
  
  
  subroutine contribution_hermite_c1(x,val)
    sll_real64,intent(in)::x(0:1)
    sll_real64,intent(out)::val(4,0:1,0:1)
    sll_int32::i,ell1,ell2,s
    sll_real64::w(0:3,0:1)!,tmp(0:3)
    do s=0,1
      w(0,s)=(2._f64*x(s)+1)*(1._f64-x(s))*(1._f64-x(s));
      w(1,s)=x(s)*x(s)*(3._f64-2._f64*x(s))
      w(2,s)=x(s)*(1._f64-x(s))*(1._f64-x(s))
      w(3,s)=x(s)*x(s)*(x(s)-1._f64)
    enddo
    
    
    val(1,0,0)=w(0,0)*w(0,1)  
    val(2,0,0)=w(2,0)*w(0,1)  
    val(3,0,0)=w(0,0)*w(2,1)  
    val(4,0,0)=w(2,0)*w(2,1)  

    val(1,1,0)=w(1,0)*w(0,1)  
    val(2,1,0)=w(3,0)*w(0,1)  
    val(3,1,0)=w(1,0)*w(2,1)  
    val(4,1,0)=w(3,0)*w(2,1)  
    
    val(1,0,1)=w(0,0)*w(1,1)  
    val(2,0,1)=w(2,0)*w(1,1)  
    val(3,0,1)=w(0,0)*w(3,1)  
    val(4,0,1)=w(2,0)*w(3,1)  
    
    val(1,1,1)=w(1,0)*w(1,1)  
    val(2,1,1)=w(3,0)*w(1,1)  
    val(3,1,1)=w(1,0)*w(3,1)  
    val(4,1,1)=w(3,0)*w(3,1)  
    
        
    !print *,val
    !stop

  end subroutine contribution_hermite_c1
  
  
  subroutine contribution_spl(x,val)
    sll_real64,intent(in)::x(0:1)
    sll_real64,intent(out)::val(-1:2,-1:2)
    sll_int32::s
    sll_real64::w(-1:2,0:1)
    do s=0,1
      w(-1,s)=(1._f64/6._f64)*(1._f64-x(s))*(1._f64-x(s))*(1._f64-x(s));
      w(0,s)=1._f64/6._f64+0.5_f64*(1._f64-x(s))*(-(1._f64-x(s))*&
    	 (1._f64-x(s))+(1._f64-x(s))+1._f64);
      w(1,s)=1._f64/6._f64+0.5_f64*x(s)*(-x(s)*x(s)+x(s)+1._f64);
      w(2,s)=(1._f64/6._f64)*x(s)*x(s)*x(s);
    enddo
    
    
    val(-1,-1)=w(-1,0)*w(-1,1)  
    val(-1,0)=w(-1,0)*w(0,1)  
    val(-1,1)=w(-1,0)*w(1,1)  
    val(-1,2)=w(-1,0)*w(2,1)  

    val(0,-1)=w(0,0)*w(-1,1)  
    val(0,0)=w(0,0)*w(0,1)  
    val(0,1)=w(0,0)*w(1,1)  
    val(0,2)=w(0,0)*w(2,1)  
    
    val(1,-1)=w(1,0)*w(-1,1)  
    val(1,0)=w(1,0)*w(0,1)  
    val(1,1)=w(1,0)*w(1,1)  
    val(1,2)=w(1,0)*w(2,1)  
    
    val(2,-1)=w(2,0)*w(-1,1)  
    val(2,0)=w(2,0)*w(0,1)  
    val(2,1)=w(2,0)*w(1,1)  
    val(2,2)=w(2,0)*w(2,1)  

  end subroutine contribution_spl
  
  
  subroutine splcoefnat1d0(lunat,N)
    sll_int32,intent(in)::N
    sll_real64,dimension(0:1,0:N+2),intent(out)::lunat
    sll_int32::i
    sll_real64::a(0:2)
    a(0)=-3._f64;a(1)=0._f64;a(2)=3._f64;
    
    lunat(0,0)=a(0);
    lunat(1,0)=1._f64/lunat(0,0);
    lunat(0,1)=4._f64-a(1)*lunat(1,0);
    lunat(1,1)=1._f64/lunat(0,1);
    lunat(0,2)=4._f64-lunat(1,1)*(1._f64-a(2)/a(0));
    do i=2,N
      lunat(1,i)=1._f64/lunat(0,i);
      lunat(0,i+1)=4._f64-lunat(1,i);
    enddo    
    lunat(1,N+2)=a(0)/lunat(0,N);
    lunat(1,N+1)=(a(1)-lunat(1,N+2))/lunat(0,N+1);
    lunat(0,N+2)=a(2)-lunat(1,N+1);
  end subroutine splcoefnat1d0

  subroutine splcoefnat1d(p,lunat,N)
    sll_int32,intent(in)::N
    sll_real64,dimension(0:1,0:N+2),intent(in)::lunat
    sll_real64,dimension(0:N+2),intent(inout)::p
    sll_int32::i
    sll_real64::a(0:2)
    a(0)=-3._f64;a(1)=0._f64;a(2)=3._f64;    
    !p(0)=0._f64;
    !p(N+2)=0._f64;
    do i=0,N+2;
      p(i)=6._f64*p(i);
    enddo;
    do i=1,N+1;p(i)=p(i)-lunat(1,i-1)*p(i-1);enddo
    p(N+2)=p(N+2)-(lunat(1,N+1)*p(N+1)+lunat(1,N+2)*p(N));
    p(N+2)=p(N+2)/lunat(0,N+2);
    do i=N+1,2,-1;p(i)=(p(i)-p(i+1))/lunat(0,i);enddo
    p(1)=(p(1)-(1._f64-a(2)/a(0))*p(2))/lunat(0,1);
    p(0)=(p(0)-a(1)*p(1)-a(2)*p(2))/lunat(0,0);
    !p(i-1)+4*p(i)+p(i+1)=6ftab(i-1,j), i=1..N+1
    !a(0)*p(0)+a(1)*p(1)+a(2)*p(2)=f'(rmin)=0);
    !a(0)*p(Na)+a(1)*p(Na+1)+a(2)*p(Na+2)=f'(rmax)=0);
  end subroutine splcoefnat1d
  
  
   subroutine splcoefnat1d0old(dnat,lnat,N)
    sll_int32,intent(in)::N
    sll_real64,dimension(0:N+2),intent(inout)::dnat,lnat
    sll_int32::i
    sll_real64::a(0:2)
    a(0)=-3._f64;a(1)=0._f64;a(2)=3._f64;
    
    dnat(0)=a(0);
    lnat(0)=1._f64/dnat(0);
    dnat(1)=4._f64-a(1)*lnat(0);
    lnat(1)=1._f64/dnat(1);
    dnat(2)=4._f64-lnat(1)*(1._f64-a(2)/a(0));
    do i=2,N
      lnat(i)=1._f64/dnat(i);
      dnat(i+1)=4._f64-lnat(i);
    enddo    
    lnat(N+2)=a(0)/dnat(N);
    lnat(N+1)=(a(1)-lnat(N+2))/dnat(N+1);
    dnat(N+2)=a(2)-lnat(N+1);
  end subroutine splcoefnat1d0old



subroutine splcoefnatper2d(f,buf,dnatx,lnatx,dpery,lpery,mpery,Nx,Ny)
    sll_int32,intent(in)::Nx,Ny
    sll_real64,dimension(:,:),pointer::f
    sll_real64,dimension(:),pointer::buf,dpery,lpery,mpery,dnatx,lnatx
    sll_int32::i,j    
    !natural spline coefficients in x    
    do j=0,Ny-1
      do i=0,Nx+2;buf(i)=f(i,j);enddo
      call splcoefnat1dold(buf,dnatx,lnatx,Nx)
      do i=0,Nx+2;f(i,j)=buf(i);enddo      
    enddo
    !periodic spline coefficients in y    
    do i=0,Nx+2      
      do j=0,Ny-1;buf(j)=f(i,j);enddo
      call splcoefper1dold(buf,dpery,lpery,mpery,Ny)
      do j=0,Ny-1;f(i,j)=buf(j);enddo      
    enddo
    
  end subroutine splcoefnatper2d

  subroutine splnatper2d(f,xx,xmin,xmax,yy,ymin,ymax,fval,Nx,Ny)
    sll_int32,intent(in)::Nx,Ny
    sll_real64,intent(in)::xx,xmin,xmax,yy,ymin,ymax
    sll_real64,intent(out)::fval
    sll_real64,dimension(:,:),pointer::f
    sll_int32::i,j
    sll_int32::ix(0:3),iy(0:3)
    sll_real64::x,y
    sll_real64::wx(0:3),wy(0:3),tmp(0:3) 

    x=(xx-xmin)/(xmax-xmin)
    !x=x-real(floor(x),f64)
    if(x<0._f64)x=0._f64;
    if(x>=1._f64)x=1._f64-1.e-12_f64;
    x=x*real(Nx,f64)
    i=floor(x)
    x=x-real(i,f64)

    y=(yy-ymin)/(ymax-ymin)
    y=y-real(floor(y),f64)
    y=y*real(Ny,f64)
    j=floor(y)
    y=y-real(j,f64)

    wx(0)=(1._f64/6._f64)*(1._f64-x)*(1._f64-x)*(1._f64-x);
    wx(1)=1._f64/6._f64+0.5_f64*(1._f64-x)*(-(1._f64-x)*&
    	 (1._f64-x)+(1._f64-x)+1._f64);
    wx(2)=1._f64/6._f64+0.5_f64*x*(-x*x+x+1._f64);
    wx(3)=(1._f64/6._f64)*x*x*x;

    wy(0)=(1._f64/6._f64)*(1._f64-y)*(1._f64-y)*(1._f64-y);
    wy(1)=1._f64/6._f64+0.5_f64*(1._f64-y)*(-(1._f64-y)*&
    	 (1._f64-y)+(1._f64-y)+1._f64);
    wy(2)=1._f64/6._f64+0.5_f64*y*(-y*y+y+1._f64);
    wy(3)=(1._f64/6._f64)*y*y*y;

    iy(0)=mod(j+Ny-1,Ny)
    iy(1)=j
    iy(2)=mod(j+1,Ny)
    iy(3)=mod(j+2,Ny)

    ix(0)=i
    ix(1)=i+1
    ix(2)=i+2
    ix(3)=i+3

    
    tmp(0)=wx(0)*f(ix(0),iy(0))+wx(1)*f(ix(1),iy(0))&
    +wx(2)*f(ix(2),iy(0))+wx(3)*f(ix(3),iy(0))
    tmp(1)=wx(0)*f(ix(0),iy(1))+wx(1)*f(ix(1),iy(1))&
    +wx(2)*f(ix(2),iy(1))+wx(3)*f(ix(3),iy(1))
    tmp(2)=wx(0)*f(ix(0),iy(2))+wx(1)*f(ix(1),iy(2))&
    +wx(2)*f(ix(2),iy(2))+wx(3)*f(ix(3),iy(2))
    tmp(3)=wx(0)*f(ix(0),iy(3))+wx(1)*f(ix(1),iy(3))&
    +wx(2)*f(ix(2),iy(3))+wx(3)*f(ix(3),iy(3))
    
    fval=wy(0)*tmp(0)+wy(1)*tmp(1)+wy(2)*tmp(2)+wy(3)*tmp(3)
        
  end subroutine splnatper2d




subroutine splper1d(f,xx,xmin,xmax,fval,N)
    sll_int32,intent(in)::N
    sll_real64,intent(in)::xx,xmin,xmax
    sll_real64,intent(out)::fval
    sll_real64,dimension(0:N-1),intent(inout)::f
    sll_int32::i
    real(f64)::x
    real(f64)::w(0:3) 
    x=(xx-xmin)/(xmax-xmin)
    x=x-real(floor(x),f64)
    x=x*real(N,f64)
    i=floor(x)
    x=x-real(i,f64)
    w(0)=(1._f64/6._f64)*(1._f64-x)*(1._f64-x)*(1._f64-x);
    w(1)=1._f64/6._f64+0.5_f64*(1._f64-x)*(-(1._f64-x)*&
    	 (1._f64-x)+(1._f64-x)+1._f64);
    w(2)=1._f64/6._f64+0.5_f64*x*(-x*x+x+1._f64);
    w(3)=(1._f64/6._f64)*x*x*x;
    fval=w(0)*f(mod(i+N-1,N))+w(1)*f(i)+w(2)*f(mod(i+1,N))+w(3)*f(mod(i+2,N))
  end subroutine splper1d


 subroutine splnat1d(f,xx,xmin,xmax,fval,N)
    sll_int32,intent(in)::N
    sll_real64,intent(in)::xx,xmin,xmax
    sll_real64,intent(out)::fval
    sll_real64,dimension(0:N+2),intent(inout)::f
    sll_int32::i
    sll_real64::x
    sll_real64::w(0:3) 
    x=(xx-xmin)/(xmax-xmin)
    !x=x-real(floor(x),f64)
    if(x<0._f64)x=0._f64;
    if(x>1._f64)x=1._f64-1.e-12_f64;
    x=x*real(N,f64)
    i=floor(x)
    x=x-real(i,f64)

    w(0)=(1._f64/6._f64)*(1._f64-x)*(1._f64-x)*(1._f64-x);
    w(1)=1._f64/6._f64+0.5_f64*(1._f64-x)*(-(1._f64-x)*&
    	 (1._f64-x)+(1._f64-x)+1._f64);
    w(2)=1._f64/6._f64+0.5_f64*x*(-x*x+x+1._f64);
    w(3)=(1._f64/6._f64)*x*x*x;
    fval=w(0)*f(i)+w(1)*f(i+1)+w(2)*f(i+2)+w(3)*f(i+3)
  end subroutine splnat1d
  
  
  subroutine splcoefper1d0old(dper,lper,mper,N)
    sll_int32,intent(in)::N
    sll_real64,dimension(0:N-1),intent(out)::dper,lper,mper
    sll_int32::i
    
    dper(0)=4._f64
    mper(0)=0.25_f64
    do i=0,N-2
      lper(i)=1._f64/dper(i)
      dper(i+1)=4._f64-lper(i)
      mper(i+1)=-mper(i)/dper(i+1)
    enddo
    dper(N-1)=dper(N-1)-(lper(N-2)+2._f64*mper(N-2))  
    do i=0,N-1
      dper(i)=1._f64/dper(i)
    enddo
  end subroutine splcoefper1d0old
  
  
  
  subroutine splcoefper1d0(luper,N)
    sll_int32,intent(in)::N
    sll_real64,dimension(0:3*N-1),intent(out)::luper
    sll_int32::i
    
    luper(0+3*0)=4._f64
    luper(2+3*0)=0.25_f64
    do i=0,N-2
      luper(1+3*i)=1._f64/luper(0+3*i)
      luper(0+3*(i+1))=4._f64-luper(1+3*i)
      luper(2+3*(i+1))=-luper(2+3*i)/luper(0+3*(i+1))
    enddo
    luper(0+3*(N-1))=luper(0+3*(N-1))-(luper(1+3*(N-2))+2._f64*luper(2+3*(N-2)))  
    do i=0,N-1
      luper(0+3*i)=1._f64/luper(0+3*i)
    enddo
  end subroutine splcoefper1d0
  
  
  subroutine splcoefper1d(f,luper,N)
    sll_int32,intent(in)::N
    sll_real64,dimension(0:3*N-1),intent(in)::luper
    sll_real64,dimension(0:N-1),intent(inout)::f
    sll_int32::i
    do i=0,N-1;f(i)=6._f64*f(i);enddo;
    do i=1,N-1
      f(i)=f(i)-f(i-1)*luper(1+3*(i-1))
    enddo
    do i=0,N-2
      f(N-1)=f(N-1)-luper(2+3*i)*f(i)
    enddo
    f(N-1)=f(N-1)*luper(0+3*(N-1));f(N-2)=luper(0+3*(N-2))*(f(N-2)-(1._f64-luper(2+3*(N-3)))*f(N-1))
    do i=N-3,1,-1
      f(i)=luper(0+3*i)*(f(i)-f(i+1)+luper(2+3*(i-1))*f(N-1))
    enddo
    f(0)=luper(0+3*0)*(f(0)-f(1)-f(N-1));
  end subroutine splcoefper1d
  


subroutine splcoefnat1dold(p,dnat,lnat,N)
    sll_int32,intent(in)::N
    sll_real64,dimension(0:N+2),intent(in)::dnat,lnat
    sll_real64,dimension(0:N+2),intent(inout)::p
    sll_int32::i
    sll_real64::a(0:2)
    a(0)=-3._f64;a(1)=0._f64;a(2)=3._f64;
    
    !p(0)=0._f64;
    !p(N+2)=0._f64;
    do i=0,N+2;
      p(i)=6._f64*p(i);
    enddo;
    do i=1,N+1;p(i)=p(i)-lnat(i-1)*p(i-1);enddo
    p(N+2)=p(N+2)-(lnat(N+1)*p(N+1)+lnat(N+2)*p(N));
    p(N+2)=p(N+2)/dnat(N+2);
    do i=N+1,2,-1;p(i)=(p(i)-p(i+1))/dnat(i);enddo
    p(1)=(p(1)-(1._f64-a(2)/a(0))*p(2))/dnat(1);
    p(0)=(p(0)-a(1)*p(1)-a(2)*p(2))/dnat(0);
    !p(i-1)+4*p(i)+p(i+1)=6ftab(i-1,j), i=1..N+1
    !a(0)*p(0)+a(1)*p(1)+a(2)*p(2)=f'(rmin)=0);
    !a(0)*p(Na)+a(1)*p(Na+1)+a(2)*p(Na+2)=f'(rmax)=0);
  end subroutine splcoefnat1dold
  
  subroutine splcoefper1dold(f,dper,lper,mper,N)
    sll_int32,intent(in)::N
    sll_real64,dimension(0:N-1),intent(in)::dper,lper,mper
    sll_real64,dimension(0:N-1),intent(inout)::f
    sll_int32::i
    do i=0,N-1;f(i)=6._f64*f(i);enddo;
    do i=1,N-1
      f(i)=f(i)-f(i-1)*lper(i-1)
    enddo
    do i=0,N-2
      f(N-1)=f(N-1)-mper(i)*f(i)
    enddo
    f(N-1)=f(N-1)*dper(N-1);f(N-2)=dper(N-2)*(f(N-2)-(1._f64-mper(N-3))*f(N-1))
    do i=N-3,1,-1
      f(i)=dper(i)*(f(i)-f(i+1)+mper(i-1)*f(N-1))
    enddo
    f(0)=dper(0)*(f(0)-f(1)-f(N-1));
  end subroutine splcoefper1dold
  
  
  subroutine solve_tridiag(a,b,c,v,x,n)
! Using Thomas' Algorithm
	implicit none
!        a - sub-diagonal (means it is the diagonal below the main diagonal)
!        b - the main diagonal
!        c - sup-diagonal (means it is the diagonal above the main diagonal)
!        v - right part
!        x - the answer
!        n - number of equations
 
	sll_int32,intent(in) :: n
	sll_real64,dimension(n),intent(in) :: a,b,c,v
	sll_real64,dimension(n),intent(out) :: x
	sll_real64,dimension(n) :: bp,vp
	sll_real64 :: m
	sll_int32 i
 
	! Make copies of the b and v variables so that they are unaltered by this sub
	bp(1) = b(1)
	vp(1) = v(1)
 
	!The first pass (setting coefficients):
    firstpass: do i = 2,n
	m = a(i)/bp(i-1)
	bp(i) = b(i) - m*c(i-1)
	vp(i) = v(i) - m*vp(i-1)
    end do firstpass
 
	x(n) = vp(n)/bp(n)
	!The second pass (back-substition)
    backsub:do i = n-1, 1, -1
	x(i) = (vp(i) - c(i)*x(i+1))/bp(i)
    end do backsub

  end subroutine solve_tridiag



  subroutine compute_shape_circle(points,N_points)
    sll_int32,intent(in) :: N_points
    sll_real64,dimension(:,:) ::points
    sll_int32 :: i
    sll_real64 :: x
    do i=1,N_points
      x = 2._f64*sll_pi*real(i,f64)/(real(N_points,f64))
      points(1,i) = cos(x)
      points(2,i) = sin(x)
      points(3,i) = 1._f64/real(N_points,f64)
    enddo
   
  end subroutine compute_shape_circle




end module sll_gyroaverage_2d_polar
