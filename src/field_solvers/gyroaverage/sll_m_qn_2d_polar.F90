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



module sll_m_qn_2d_polar
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

! use F77_fftpack, only: &
!   zfftb, &
!   zfftf, &
!   zffti

! use F77_lapack, only: &
!   zgetrf, &
!   zgetri

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_gnuplot, only: &
    sll_o_gnuplot_2d

  use sll_m_gyroaverage_utilities, only: &
    sll_s_compute_shape_circle, &
    sll_s_zero_bessel_dir_dir

  implicit none

  public :: &
    sll_s_compute_error, &
    sll_s_compute_gamma0, &
    sll_f_compute_gamma0_quadrature, &
    sll_s_compute_splines_coefs_matrix_nat_1d, &
    sll_s_compute_splines_coefs_matrix_per_1d, &
    sll_s_contribution_spl, &
    sll_s_initialize_mu_quadr_for_phi, &
    sll_s_localize_polar, &
    sll_s_matrix_product_compf, &
    sll_f_new_plan_qn_polar_splines, &
    sll_s_precompute_gyroaverage_index, &
    sll_s_precompute_inverse_qn_matrix_polar_splines, &
    sll_t_plan_qn_polar, &
    sll_s_solve_qn_polar_splines, &
    sll_s_splcoefnat1d0old, &
    sll_s_splcoefper1d0old, &
    sll_s_test_solve_qn_polar_splines

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type sll_t_plan_qn_polar
     
     sll_real64          :: eta_min(2)     !< r min et theta min
     sll_real64          :: eta_max(2)     !< r max et theta max
     sll_int32           :: Nc(2)          !< number of cells in r and in theta
     
     sll_int32           :: N_points          !< number of points on the circle
     
     sll_real64, dimension(:,:), pointer    :: points
     sll_int32, dimension(:,:), pointer       :: pre_compute_N
     sll_real64, dimension(:,:,:), pointer    :: pre_compute_coeff
     sll_int32, dimension(:,:,:), pointer     :: pre_compute_index
     sll_real64, dimension(:,:), pointer      :: pre_compute_coeff_spl
     sll_int32 :: size_pre_compute

     sll_real64, dimension(:,:,:), pointer    :: mat_gyro
     sll_real64, dimension(:,:,:), pointer    :: mat_double_gyro 
     sll_real64, dimension(:,:,:,:), pointer    :: mat_gyro_circ
     sll_real64, dimension(:,:,:,:), pointer    :: mat_double_gyro_circ

     sll_comp64, dimension(:,:,:), pointer    :: mat_qn_inverse
     
     sll_real64, dimension(:), allocatable    :: lambda
     sll_real64, dimension(:), allocatable    :: T_i
     ! solve \lambda(r)\phi-\tilde{\phi} = second member
     
     sll_real64, dimension(:), pointer        :: mu_points_for_phi
     sll_real64, dimension(:), pointer        :: mu_weights_for_phi
     sll_int32                                :: N_mu_for_phi

  end type sll_t_plan_qn_polar

contains


  function sll_f_new_plan_qn_polar_splines(eta_min,eta_max,Nc,N_points,lambda,T_i) result(this)

    implicit none

    sll_real64, intent(in) :: eta_min(2)
    sll_real64, intent(in) :: eta_max(2)
    sll_int32, intent(in)  :: Nc(2)
    sll_int32, intent(in)  :: N_points
    sll_real64, dimension(:), intent(in)    :: lambda
    sll_real64, dimension(:), intent(in)    :: T_i
    type(sll_t_plan_qn_polar), pointer :: this

    sll_int32 :: err

    SLL_ALLOCATE(this,err)
    SLL_ALLOCATE(this%points(3,N_points),err)
    SLL_ALLOCATE(this%lambda(1:Nc(1)+1),err)
    SLL_ALLOCATE(this%T_i(1:Nc(1)+1),err)
    
    call sll_s_compute_shape_circle(this%points,N_points) 
       
    this%eta_min=eta_min
    this%eta_max=eta_max
    this%Nc=Nc
    this%N_points=N_points
    
  
    this%lambda=lambda
    this%T_i=T_i
    
  end function sll_f_new_plan_qn_polar_splines

  
  
  subroutine sll_s_compute_splines_coefs_matrix_nat_1d(mat,dnat,lnat,N)
    sll_int32,intent(in)::N
    sll_real64,dimension(0:N+2,1:N+1),intent(inout) :: mat
    sll_real64,dimension(0:N+2),intent(in)::dnat,lnat
    sll_real64,dimension(:,:),allocatable :: matb
    sll_int32 :: err,i,j
    
    ! Conditions aux bords nulles
    
    SLL_ALLOCATE(matb(0:N+2,0:N+2),err)
    
    matb = 0._f64
    
    do j = 0,N+2
      matb(j,j) = 6._f64
    enddo
    
    do i = 1,N+1
      do j = 0,N+2
        matb(i,j) = matb(i,j)-lnat(i-1)*matb(i-1,j)
      enddo
    enddo
    
    do j = 0,N+2
      matb(N+2,j) = (matb(N+2,j)-lnat(N+1)*matb(N+1,j)-lnat(N+2)*matb(N,j))/dnat(N+2)
    enddo

    do i = N+1,2,-1
      do j = 0,N+2
        matb(i,j) = (matb(i,j)-matb(i+1,j))/dnat(i)
      enddo
    enddo
    
    do j = 0,N+2
      matb(1,j) = (matb(1,j) - 2._f64*matb(2,j))/dnat(1) 
    enddo
    
    do j = 0,N+2
      matb(0,j) = (matb(0,j) - 3._f64*matb(2,j))/dnat(0) 
    enddo 
    
    mat(0:N+2,1:N+1) = matb(0:N+2,1:N+1)
  
  end subroutine sll_s_compute_splines_coefs_matrix_nat_1d
  
  
  
  subroutine sll_s_compute_splines_coefs_matrix_per_1d(mat,dper,lper,mper,N)
    sll_int32,intent(in)::N
    sll_real64,dimension(0:N-1,0:N-1),intent(inout) :: mat
    sll_real64,dimension(0:N+2),intent(in)::dper,lper,mper
    sll_int32 :: i,j

    mat = 0._f64
    
    do j = 0,N-1
      mat(j,j) = 6._f64
    enddo

    do i = 1,N-1
      do j = 0,N-1
      mat(i,j) = mat(i,j)-lper(i-1)*mat(i-1,j)
      enddo
    enddo

    do i = 0,N-2
      do j = 0,N-1
      mat(N-1,j) = mat(N-1,j)-mper(i)*mat(i,j)
      enddo
    enddo
    
    do j = 0,N-1
      mat(N-1,j) = mat(N-1,j)*dper(N-1) 
    enddo
    
    do j = 0,N-1
      mat(N-2,j) = dper(N-2)*(mat(N-2,j)-(1._f64-mper(N-3))*mat(N-1,j))
    enddo

    do i=N-3,1,-1
      do j = 0,N-1
        mat(i,j) = dper(i)*(mat(i,j)-mat(i+1,j)+mper(i-1)*mat(N-1,j))
      enddo
    enddo
        
    do j = 0,N-1
      mat(0,j) = dper(0)*(mat(0,j)-mat(1,j)-mat(N-1,j))
    enddo
  
  end subroutine sll_s_compute_splines_coefs_matrix_per_1d
  
  
  subroutine kronecker_product(A,NxA,NyA,B,NxB,NyB,kronecker)
    sll_int32,intent(in)::NxA,NyA,NxB,NyB
    sll_real64,dimension(0:NxA-1,0:NyA-1),intent(in)::A
    sll_real64,dimension(0:NxB-1,0:NyB-1),intent(in)::B
    sll_real64,dimension(0:NxA*NxB-1,0:NyA*NyB-1),intent(inout) :: kronecker
    sll_int32 :: iA,iB,jA,jB
    
    do iA = 0,NxA-1
      do iB = 0,NxB-1
        do jA = 0,NyA-1
          do jB = 0,NyB-1
            kronecker(iB+NxB*iA,jB+NyB*jA)=A(iA,jA)*B(iB,jB)
          enddo
        enddo
      enddo
    enddo
    
  end subroutine kronecker_product
  
  subroutine matrix_product(A,NxA,NyA,B,NxB,NyB,prod)
    sll_int32,intent(in)::NxA,NyA,NxB,NyB
    sll_real64,dimension(0:NxA-1,0:NyA-1),intent(in)::A
    sll_real64,dimension(0:NxB-1,0:NyB-1),intent(in)::B
    sll_real64,dimension(0:NxA-1,0:NyB-1),intent(inout) :: prod
    sll_int32 :: i,j,k
    sll_real64 :: result
    
    if (NyA/=NxB) then
      print *,'#incompatible sizes in matrix_product'
      stop
    else
      do i = 0,NxA-1
        do j = 0,NyB-1
          result = 0._f64
          do k = 0,NyA-1
            result = result + A(i,k)*B(k,j)
          enddo
          prod(i,j)=result 
        enddo
      enddo
    endif
    
  end subroutine matrix_product
 
   subroutine sll_s_matrix_product_compf(A,NxA,NyA,B,NxB,NyB,prod)
    sll_int32,intent(in)::NxA,NyA,NxB,NyB
    sll_comp64,dimension(:,:),intent(in)::A
    sll_comp64,dimension(:,:),intent(in)::B
    sll_comp64,dimension(:,:),intent(inout) :: prod
    sll_int32 :: i,j,k
    sll_comp64 :: result
    
    if (NyA/=NxB) then
      print *,'#incompatible sizes in matrix_product'
      stop
    else
      do i = 1,NxA
        do j = 1,NyB
          result = (0.0_f64,0.0_f64)
          do k = 1,NyA
            result = result + A(i,k)*B(k,j)
          enddo
          prod(i,j)=result 
        enddo
      enddo
    endif
    
  end subroutine sll_s_matrix_product_compf
 
 
  
  
    subroutine matrix_product_comp(A,NxA,NyA,B,NxB,NyB,prod)
    sll_int32,intent(in)::NxA,NyA,NxB,NyB
    sll_comp64,dimension(0:NxA-1,0:NyA-1),intent(in)::A
    sll_comp64,dimension(0:NxB-1,0:NyB-1),intent(in)::B
    sll_comp64,dimension(0:NxA-1,0:NyB-1),intent(inout) :: prod
    sll_int32 :: i,j,k
    sll_comp64 :: result
    
    if (NyA/=NxB) then
      print *,'#incompatible sizes in matrix_product'
      stop
    else
      do i = 0,NxA-1
        do j = 0,NyB-1
          result = (0._f64,0.0_f64)
          do k = 0,NyA-1
            result = result + A(i,k)*B(k,j)
          enddo
          prod(i,j)=result 
        enddo
      enddo
    endif
    
  end subroutine matrix_product_comp
  
  
    subroutine matrix_product_circ(A,NxA,NyA,B,NxB,NyB,prod,Ntheta)
    sll_int32,intent(in)::Ntheta,NxA,NyA,NxB,NyB
    sll_real64,dimension(0:Ntheta-1,0:NxA-1,0:NyA-1),intent(in)::A
    sll_real64,dimension(0:Ntheta-1,0:NxB-1,0:NyB-1),intent(in)::B
    sll_real64,dimension(0:Ntheta-1,0:NxA-1,0:NyB-1),intent(inout) :: prod
    sll_real64,dimension(:,:),allocatable :: mat_stock,mat_stock_sum 
    sll_int32 :: i,j,error
    !sll_real64 :: result

    SLL_ALLOCATE(mat_stock(0:NxA-1,0:NyB-1),error)
    SLL_ALLOCATE(mat_stock_sum(0:NxA-1,0:NyB-1),error)
   
    if (NyA/=NxB) then
      print *,'#incompatible sizes in matrix_product_circ'
      stop
    else
      do i = 0,Ntheta-1
      mat_stock_sum = 0._f64
        do j = 0,Ntheta-1
          call matrix_product(A(modulo(i+j,Ntheta),:,:),NxA,NyA,B(j,:,:),NxB,NyB,mat_stock)
          mat_stock_sum = mat_stock_sum + mat_stock
        enddo
      prod(i,:,:) = mat_stock_sum 
      enddo
    endif
    
  end subroutine matrix_product_circ


  
   subroutine sll_s_precompute_gyroaverage_index(quasineutral,rho,N_rho)
    type(sll_t_plan_qn_polar)  :: quasineutral
    sll_int32 :: N_rho
    sll_real64,dimension(1:N_rho),intent(in) :: rho
    sll_int32,dimension(:,:),allocatable :: buf
    sll_int32 ::i,j,k,p,ell_1,ell_2,ii(2),s,nb,ind(2)
    sll_real64::val(-1:2,-1:2),eta_star(2),eta(2),delta_eta(2),x(2)
    sll_int32 ::error,max_nb
    
    
    delta_eta(1)=(quasineutral%eta_max(1)-quasineutral%eta_min(1))/real(quasineutral%Nc(1),f64)
    delta_eta(2)=(quasineutral%eta_max(2)-quasineutral%eta_min(2))/real(quasineutral%Nc(2),f64)
    
    SLL_ALLOCATE(buf(0:quasineutral%Nc(1)+2,0:quasineutral%Nc(2)-1),error)      
    SLL_ALLOCATE(quasineutral%pre_compute_N(1:N_rho,quasineutral%Nc(1)+1),error)

    eta(2)=quasineutral%eta_min(2)
    
    max_nb=0
 
    do p=1,N_rho
    
    buf=0
    nb=0
    do i=1,quasineutral%Nc(1)+1
      eta(1)=quasineutral%eta_min(1)+real(i-1,f64)*delta_eta(1)       
      s=0
      do k=1,quasineutral%N_points
        x(1) = eta(1)*cos(eta(2))+rho(p)*quasineutral%points(1,k)
        x(2) = eta(1)*sin(eta(2))+rho(p)*quasineutral%points(2,k)
        call sll_s_localize_polar(x,quasineutral%eta_min,quasineutral%eta_max,ii,eta_star,quasineutral%Nc)
        do ell_2=-1,2
          ind(2)=modulo(ii(2)+ell_2,quasineutral%Nc(2))
          do ell_1=-1,2
            ind(1)=ii(1)+1+ell_1
            if(buf(ind(1),ind(2)).ne.i)then
              s=s+1
              buf(ind(1),ind(2))=i
            endif
          enddo
        enddo    
      enddo
      quasineutral%pre_compute_N(p,i)=s
      nb=nb+s
    enddo
    max_nb = max(max_nb,nb)
    
    quasineutral%size_pre_compute = max_nb    
   ! print*,'#N_points pre_compute=',nb,quasineutral%Nc(1),nb/quasineutral%Nc(1)

    enddo ! p=1,N_rho


    SLL_ALLOCATE(quasineutral%pre_compute_index(1:N_rho,1:2,1:max_nb),error)
    SLL_ALLOCATE(quasineutral%pre_compute_coeff_spl(1:N_rho,1:max_nb),error)

    do p=1,N_rho
    
    buf=0
    nb=0
    s=0
    val=0._f64
    eta(2)=quasineutral%eta_min(2)
    do i=1,quasineutral%Nc(1)+1
      eta(1)=quasineutral%eta_min(1)+real(i-1,f64)*delta_eta(1)       
      do k=1,quasineutral%N_points
        x(1) = eta(1)*cos(eta(2))+rho(p)*quasineutral%points(1,k)
        x(2) = eta(1)*sin(eta(2))+rho(p)*quasineutral%points(2,k)
        call sll_s_localize_polar(x,quasineutral%eta_min,quasineutral%eta_max,ii,eta_star,quasineutral%Nc)

        call sll_s_contribution_spl(eta_star,val)
        
        val=val*quasineutral%points(3,k)
        
        do ell_2=-1,2
          ind(2)=modulo(ii(2)+ell_2,quasineutral%Nc(2))
          do ell_1=-1,2
            ind(1)=ii(1)+1+ell_1
            j=buf(ind(1),ind(2))
            if(j<=nb)then
              s=s+1
              buf(ind(1),ind(2))=s
              quasineutral%pre_compute_coeff_spl(p,s)=val(ell_1,ell_2)
              quasineutral%pre_compute_index(p,1,s)=ind(1)
              quasineutral%pre_compute_index(p,2,s)=ind(2)                     
            else              
              quasineutral%pre_compute_coeff_spl(p,j)=quasineutral%pre_compute_coeff_spl(p,j)+val(ell_1,ell_2)
            endif
          enddo
        enddo
      enddo      
      nb=s 
    enddo
    
    enddo ! p=1,N_rho
    
 !   print *,'#min/max=',minval(quasineutral%pre_compute_coeff_spl(:)),maxval(quasineutral%pre_compute_coeff_spl(:))
  
    SLL_DEALLOCATE_ARRAY(buf,error)
  
  end subroutine sll_s_precompute_gyroaverage_index  



   subroutine sll_s_precompute_inverse_qn_matrix_polar_splines(quasineutral,mu_points,mu_weights,N_mu)
    type(sll_t_plan_qn_polar) :: quasineutral
    sll_int32,intent(in) :: N_mu
    sll_real64,dimension(1:N_mu),intent(in) :: mu_points
    sll_real64,dimension(1:N_mu),intent(in) :: mu_weights
    sll_comp64,dimension(:,:),allocatable :: mat_stock1,mat_stock2
    sll_real64,dimension(:),allocatable,target::dnat,lnat,dper,lper,mper
    sll_real64,dimension(:),pointer::pointer_dnat,pointer_lnat
    sll_real64,dimension(:),pointer::pointer_dper,pointer_lper,pointer_mper
    sll_real64,dimension(:,:),allocatable,target::mat_nat,mat_per
    sll_real64,dimension(:,:),pointer::pointer_mat_nat,pointer_mat_per
    sll_real64,dimension(:,:,:),allocatable,target::mat_spl2D_circ
    sll_real64,dimension(:,:,:),pointer::pointer_mat_spl2D_circ
    sll_real64,dimension(:,:,:,:),allocatable,target::mat_contribution_circ
    sll_real64,dimension(:,:,:,:),pointer::pointer_mat_contribution_circ
    sll_comp64,dimension(:,:,:),allocatable::D_spl2D
    sll_comp64,dimension(:,:,:,:),allocatable::D_contr
    sll_int32,dimension(:),allocatable :: IPIV
    sll_comp64,dimension(:),allocatable :: WORK
    sll_real64 :: mode
    sll_comp64 :: exp_comp
    sll_int32 :: INFO
    sll_int32 :: Nr,Ntheta
    sll_int32 :: error
    sll_int32 :: i,j,k,m,p,s
    sll_int32 :: ii(2)
  
 ! Taille du maillage  
    Nr = quasineutral%Nc(1)
    Ntheta = quasineutral%Nc(2)
    
 ! Allocate
    SLL_ALLOCATE(mat_stock1(0:Nr,0:Nr),error)
    SLL_ALLOCATE(mat_stock2(0:Nr,0:Nr),error)
    SLL_ALLOCATE(quasineutral%mat_qn_inverse(0:Ntheta-1,0:Nr,0:Nr),error)    
    SLL_ALLOCATE(dnat(0:Nr+2),error)
    SLL_ALLOCATE(lnat(0:Nr+2),error)
    SLL_ALLOCATE(dper(0:Ntheta-1),error)
    SLL_ALLOCATE(lper(0:Ntheta-1),error)
    SLL_ALLOCATE(mper(0:Ntheta-1),error)
    SLL_ALLOCATE(mat_nat(0:Nr+2,0:Nr),error)
    SLL_ALLOCATE(mat_per(0:Ntheta-1,0:Ntheta-1),error)
    SLL_ALLOCATE(mat_spl2D_circ(0:Ntheta-1,0:Nr+2,0:Nr),error)
    SLL_ALLOCATE(D_spl2D(0:Ntheta-1,0:Nr+2,0:Nr),error)
    SLL_ALLOCATE(mat_contribution_circ(0:N_mu,0:Ntheta-1,0:Nr,0:Nr+2),error)
    SLL_ALLOCATE(D_contr(0:N_mu,0:Ntheta-1,0:Nr,0:Nr+2),error)
    SLL_ALLOCATE(IPIV(Nr+1),error)
    SLL_ALLOCATE(WORK((Nr+1)**2),error) 
    
    SLL_ALLOCATE(quasineutral%mat_double_gyro_circ(1:N_mu,0:Ntheta-1,0:Nr,0:Nr),error)
    SLL_ALLOCATE(quasineutral%mat_gyro_circ(1:N_mu,0:Ntheta-1,0:Nr,0:Nr),error)

 
 ! Initialise les coefficients de splines  
    call sll_s_splcoefnat1d0old(dnat,lnat,Nr)
    call sll_s_splcoefper1d0old(dper,lper,mper,Ntheta)
    
 ! Pointeurs    
    pointer_dnat => dnat
  pointer_lnat => lnat
  pointer_dper => dper
  pointer_lper => lper
  pointer_mper => mper  
  pointer_mat_nat => mat_nat
    pointer_mat_per => mat_per
    pointer_mat_spl2D_circ => mat_spl2D_circ
    pointer_mat_contribution_circ => mat_contribution_circ
   
 ! Construction de la matrice des coefficients de splines    
    call sll_s_compute_splines_coefs_matrix_nat_1d(pointer_mat_nat,pointer_dnat,pointer_lnat,Nr)    
    call sll_s_compute_splines_coefs_matrix_per_1d(pointer_mat_per,pointer_dper,pointer_lper,pointer_mper,Ntheta) 
    do j=0,Ntheta-1  
      pointer_mat_spl2D_circ(j,:,:)=pointer_mat_per(0,j)*pointer_mat_nat
    enddo
    
 ! Construction de la matrice de contribution   
    pointer_mat_contribution_circ = 0._f64
    do p=1,N_mu   
      s=0
      do i=1,Nr+1
        do k=1,quasineutral%pre_compute_N(p,i)
          s=s+1
          ii(1)=quasineutral%pre_compute_index(p,1,s)
          ii(2)=modulo(quasineutral%pre_compute_index(p,2,s),Ntheta)
          pointer_mat_contribution_circ(p,ii(2),i-1,ii(1))=pointer_mat_contribution_circ(p,ii(2),i-1,ii(1))+quasineutral%pre_compute_coeff_spl(p,s) 
        enddo
      enddo         
    enddo  

 ! Matrices D^{spl} et D^{contr}
    D_spl2D = (0.0_f64,0.0_f64)
    D_contr = (0.0_f64,0.0_f64)
    do m=0,Ntheta-1
      do j=0,Ntheta-1
        mode=real(-2._f64*sll_p_pi*real(j,f64)*real(m,f64)/real(Ntheta,f64),f64)
        exp_comp = cmplx( cos(mode), sin(mode), kind=f64 )
        D_spl2D(m,:,:) = D_spl2D(m,:,:) + pointer_mat_spl2D_circ(j,:,:)*exp_comp
        do p=1,N_mu
          D_contr(p,m,:,:) = D_contr(p,m,:,:) + pointer_mat_contribution_circ(p,j,:,:)*exp_comp
        enddo
      enddo
    enddo
 
 ! Construction de la matrice à inverser
    quasineutral%mat_qn_inverse = (0._f64,0._f64)
    do m=0,Ntheta-1
      do p=1,N_mu
        mat_stock1 = (0._f64,0._f64)
        mat_stock2 = (0._f64,0._f64)
        call matrix_product_comp(D_contr(p,m,:,:),Nr+1,Nr+3,D_spl2D(m,:,:),Nr+3,Nr+1,mat_stock1)
        call matrix_product_comp(mat_stock1(:,:),Nr+1,Nr+1,mat_stock1(:,:),Nr+1,Nr+1,mat_stock2)

        !call matrix_product_comp(transpose(mat_stock1(:,:)),Nr+1,Nr+1,mat_stock1(:,:),Nr+1,Nr+1,mat_stock2)

        do i=0,Nr
          mat_stock2(i,i) =   mat_stock2(i,i) - quasineutral%lambda(i+1)*(1._f64,0._f64)
        enddo
        do i=0,Nr
          quasineutral%mat_qn_inverse(m,i,:) = &
          quasineutral%mat_qn_inverse(m,i,:)   &
          - mu_weights(p)*mat_stock2(i,:)*  &
          cmplx(exp(-mu_points(p)/quasineutral%T_i(i+1)),0._f64,f64)
        enddo
      enddo 
    enddo     
 ! Inversion des blocs
    do m=0,Ntheta-1
      call ZGETRF(Nr+1,Nr+1,quasineutral%mat_qn_inverse(m,:,:),Nr+1,IPIV,INFO)
      call ZGETRI(Nr+1,quasineutral%mat_qn_inverse(m,:,:),Nr+1,IPIV,WORK,(Nr+1)**2,INFO)
    enddo
    
 ! Deallocate    
    SLL_DEALLOCATE_ARRAY(mat_stock1,error)
    SLL_DEALLOCATE_ARRAY(mat_stock2,error) 
    SLL_DEALLOCATE_ARRAY(dnat,error)
    SLL_DEALLOCATE_ARRAY(lnat,error)
    SLL_DEALLOCATE_ARRAY(dper,error)
    SLL_DEALLOCATE_ARRAY(lper,error)
    SLL_DEALLOCATE_ARRAY(mper,error)    
    SLL_DEALLOCATE_ARRAY(mat_nat,error)
    SLL_DEALLOCATE_ARRAY(mat_per,error)
    SLL_DEALLOCATE_ARRAY(mat_spl2D_circ,error)
    SLL_DEALLOCATE_ARRAY(D_spl2D,error)
    SLL_DEALLOCATE_ARRAY(mat_contribution_circ,error)
    SLL_DEALLOCATE_ARRAY(D_contr,error)
    SLL_DEALLOCATE_ARRAY(IPIV,error)
    SLL_DEALLOCATE_ARRAY(WORK,error) 
    
  end subroutine sll_s_precompute_inverse_qn_matrix_polar_splines






 
 
   subroutine sll_s_solve_qn_polar_splines(quasineutral,phi)
    type(sll_t_plan_qn_polar) :: quasineutral
    sll_real64,dimension(1:quasineutral%Nc(1)+1,1:quasineutral%Nc(2)),intent(inout) :: phi
    sll_comp64,dimension(:,:),allocatable :: phi_comp,phi_old
    sll_real64,dimension(:),allocatable::buf_fft
    !sll_int32,dimension(:),allocatable :: IPIV
    !sll_comp64,dimension(:),allocatable :: WORK
    !sll_real64 :: mode
    sll_comp64 :: result
    !sll_int32 :: INFO
    sll_int32 :: Nr,Ntheta
    sll_int32 :: error
    sll_int32 :: i,j,m
    !sll_int32 :: ii(2)
  
 ! Taille du maillage  
    Nr = quasineutral%Nc(1)
    Ntheta = quasineutral%Nc(2)
    
 ! Allocate
    SLL_ALLOCATE(phi_comp(1:Nr+1,1:Ntheta),error)
    SLL_ALLOCATE(phi_old(1:Nr+1,1:Ntheta),error)
    SLL_ALLOCATE(buf_fft(1:4*Ntheta+15),error)
 
 ! FFT(PHI)
    phi_comp=phi*(1._f64,0._f64)
    call zffti(Ntheta,buf_fft)
    do i=1,Nr+1
    call zfftf(Ntheta,phi_comp(i,:),buf_fft)
    enddo   

 ! Produit matrice/vecteur 
    phi_old=phi_comp
    do m = 1,Ntheta
      do i = 1,Nr+1
        result = (0._f64,0.0_f64)
        do j = 1,Nr+1
          result = result + quasineutral%mat_qn_inverse(m-1,i-1,j-1)*phi_old(j,m)
        enddo
        phi_comp(i,m) = result
      enddo 
    enddo  

    
 ! FFT^-1
    do i=1,Nr+1
    call zfftb(Ntheta,phi_comp(i,:),buf_fft)
  enddo
  phi=real(phi_comp,f64)/real(Ntheta,f64)
    
 ! Sorties     
 !   print *,"phi_min : ",minval(phi)
 !   print *,"phi_max : ",maxval(phi)
    
 ! Deallocate    
    SLL_DEALLOCATE_ARRAY(phi_comp,error)
    SLL_DEALLOCATE_ARRAY(phi_old,error)
    SLL_DEALLOCATE_ARRAY(buf_fft,error)
    
  end subroutine sll_s_solve_qn_polar_splines



 subroutine sll_s_test_solve_qn_polar_splines(Nc,eta_min,eta_max,mu_points,mu_weights,N_mu,mode,lambda,T_i,phi_init,phi_qn)
  sll_int32,intent(in)  :: Nc(2)
  sll_real64,intent(in) :: eta_min(2)
  sll_real64,intent(in) :: eta_max(2)
  sll_int32,intent(in) :: N_mu
  sll_int32,intent(in)  :: mode(2)
  sll_real64,dimension(1:N_mu),intent(in) :: mu_points
  sll_real64,dimension(1:N_mu),intent(in) :: mu_weights
  sll_real64,dimension(1:Nc(1)+1,1:Nc(2)+1),intent(in) :: phi_init
  sll_real64,dimension(1:Nc(1)+1,1:Nc(2)+1),intent(in) :: phi_qn
  sll_real64,dimension(1:Nc(1)+1),intent(in) :: lambda,T_i
  sll_int32  :: N_min(2),N_max(2)
  sll_real64,dimension(:),allocatable :: gamma0
  sll_real64  :: tmp1
  sll_real64 :: eps,rho2d(2),mu2dmax(2),error(3)
  sll_int32 :: i,ierr,p
  
  SLL_ALLOCATE(gamma0(1:Nc(1)+1),ierr)
  
  N_min = 1
  N_max(1) = Nc(1)
  N_max(2) = Nc(2)
  
  mu2dmax(1) = maxval(mu_points)
  mu2dmax(2) = maxval(mu_points)
  
  eps=1.d-10
  gamma0 = 0._f64
  
  call compute_N_bounds_polar_circle(N_min(1),N_max(1),Nc(1),2._f64*sqrt(2._f64*mu2dmax),eta_min(1),eta_max(1))

  do p = 1, N_mu
    rho2d(1) = sqrt(2._f64*mu_points(p))
    rho2d(2) = sqrt(2._f64*mu_points(p))
    call solution_polar_circle(rho2d,mode,eta_min,eta_max,tmp1)
    do i = 1, Nc(1)+1
      gamma0(i) = gamma0(i) + mu_weights(p)*exp(-mu_points(p)/T_i(i))*(lambda(i)-tmp1**2)
    enddo
  enddo
  do i = 1, Nc(1)+1
    gamma0(i) = 1._f64/gamma0(i)
  enddo
  
  ! print *,'#gamma0val=',gamma0

  print *,'#N_min(1:2) / N_max(1:2) = ',N_min,' / ',N_max
  call compute_error_1D(phi_qn,phi_init,gamma0,error,N_min,N_max)
  print *,'#error subdomain=',error
  call compute_error_1D(phi_qn,phi_init,gamma0,error,(/1,1/),Nc)
  print *,'#error whole domain=',error


             call sll_o_gnuplot_2d( &
                  eta_min(1), &
                  eta_max(1), &
                  Nc(1)+1, &
                  eta_min(2), &
                  eta_max(2), &
                  Nc(2)+1, &
                  phi_qn, &
                  'fdiff', &
                  1, &
                  ierr)


  end subroutine sll_s_test_solve_qn_polar_splines




  function sll_f_compute_gamma0_quadrature( &
    Nc, &
    eta_min, &
    eta_max, &
    mu_points, &
    mu_weights, &
    N_mu, &
    mode ) &
    result(gamma0) 
    sll_int32,intent(in)  :: Nc(2)
    sll_real64,intent(in) :: eta_min(2)
    sll_real64,intent(in) :: eta_max(2)
    sll_int32,intent(in) :: N_mu
    sll_int32,intent(in)  :: mode(2)
    sll_real64,dimension(1:N_mu),intent(in) :: mu_points
    sll_real64,dimension(1:N_mu),intent(in) :: mu_weights
    sll_real64 :: gamma0
    sll_real64 :: tmp1
    sll_real64 :: rho2d(2)
    sll_int32 :: p
    
    gamma0 = 0._f64
    do p = 1, N_mu
      rho2d(1) = sqrt(2._f64*mu_points(p))
      rho2d(2) = sqrt(2._f64*mu_points(p))
      call solution_polar_circle(rho2d,mode,eta_min,eta_max,tmp1)
      !gamma0 = gamma0 + mu_weights(p)*exp(-mu_points(p))*(1._f64-tmp1**2)
      gamma0 = gamma0 + mu_weights(p)*(1._f64-tmp1**2)
      !gamma0 = gamma0 + mu_weights(p)*(tmp1**2)
     enddo
    
    !gamma0 = 1._f64-gamma0

    return
    print*, nc
    
  end function sll_f_compute_gamma0_quadrature 
  

subroutine solve_circulant_system(Ntheta,Nr,mat_circ,sol)
  ! Solve mat_circ*X=sol where mat_circ(0:Ntheta-1,0:Nr,0:Nr) 
  ! is a circulent matrix of size Ntheta*(Nr+1) Ntheta*(Nr+1)
  ! and sol(1:Nr+1,1:Ntheta)
  sll_int32 :: Ntheta,Nr,i,j,m,error
  sll_comp64,dimension(:,:,:),allocatable :: Dm
  sll_comp64,dimension(:,:),allocatable :: sol_comp,sol_old
  sll_real64,dimension(:),allocatable::buf_fft
  sll_real64 :: mode
  sll_comp64 :: exp_comp,result
  sll_int32,dimension(:),allocatable :: IPIV
  sll_comp64,dimension(:),allocatable :: WORK
  sll_int32 :: INFO
  sll_real64,dimension(0:Ntheta-1,0:Nr,0:Nr),intent(in) :: mat_circ
  sll_real64,dimension(1:Nr+1,1:Ntheta),intent(inout) :: sol

 SLL_ALLOCATE(Dm(0:Ntheta-1,0:Nr,0:Nr),error)
 SLL_ALLOCATE(sol_comp(1:Nr+1,1:Ntheta),error)
 SLL_ALLOCATE(sol_old(1:Nr+1,1:Ntheta),error)
 SLL_ALLOCATE(buf_fft(1:4*Ntheta+15),error)
 SLL_ALLOCATE(IPIV(Nr+1),error)
 SLL_ALLOCATE(WORK((Nr+1)**2),error)

 sol_comp=sol*(1._f64,0._f64)
 
 ! FFT(PHI)
  call zffti(Ntheta,buf_fft)
  do i=1,Nr+1
  call zfftf(Ntheta,sol_comp(i,:),buf_fft)
  enddo   

 ! Matrices Dm  
    Dm = (0._f64  ,0.0_f64)
    do m=0,Ntheta-1
      do j=0,Ntheta-1
        mode=real(-2._f64*sll_p_pi*real(j,f64)*real(m,f64)/real(Ntheta,f64),f64)
        exp_comp = cmplx( cos(mode), sin(mode), kind=f64 )
        Dm(m,:,:) = Dm(m,:,:) + mat_circ(j,:,:)*exp_comp
      enddo
      
      call ZGETRF(Nr+1,Nr+1,Dm(m,:,:),Nr+1,IPIV,INFO)
      call ZGETRI(Nr+1,Dm(m,:,:),Nr+1,IPIV,WORK,(Nr+1)**2,INFO)
              
    enddo
 
   ! Produit Dm 
    sol_old=sol_comp
    do m = 1,Ntheta
      do i = 1,Nr+1
        result = (0._f64,0.0_f64)
        do j = 1,Nr+1
          result = result + Dm(m-1,i-1,j-1)*sol_old(j,m)
        enddo
        sol_comp(i,m) = result
      enddo 
    enddo
    
    ! FFT^-1
    do i=1,Nr+1
    call zfftb(Ntheta,sol_comp(i,:),buf_fft)
  enddo
  
  sol=real(sol_comp/cmplx(Ntheta,0._f64,f64),f64)

end subroutine solve_circulant_system




  subroutine test_solve_circulant_system(mat,Nr,Ntheta)
    sll_int32,intent(in) :: Nr,Ntheta
    sll_real64,dimension(0:Ntheta-1,0:Nr,0:Nr),intent(in) :: mat
    sll_real64,dimension(:,:,:),allocatable :: inv,test_id
    sll_real64,dimension(:,:),allocatable :: phi
    sll_int32 :: i,j,m, error


    SLL_ALLOCATE(inv(0:Ntheta-1,0:Nr,0:Nr),error)
    SLL_ALLOCATE(test_id(0:Ntheta-1,0:Nr,0:Nr),error)
    SLL_ALLOCATE(phi(1:Nr+1,1:Ntheta),error)

  
  do m = 1,Nr+1
    phi = 0._f64
    phi(m,1) = 1._f64
  
    call solve_circulant_system(Ntheta,Nr,mat,phi)
  
    do i=0,Nr
      do j=0,Ntheta-1
        inv(modulo(Ntheta-j,Ntheta),i,m-1) = phi(i+1,j+1)
      enddo
    enddo    
  enddo
  
  call matrix_product_circ(mat,Nr+1,Nr+1,inv,Nr+1,Nr+1,test_id,Ntheta)
  
  do i=0,Nr
    test_id(0,i,i)=test_id(0,i,i)-1._f64
  enddo
  
  print *,"Inversion error = ",maxval(test_id)
  
  end subroutine test_solve_circulant_system









 !----------------------------------------




subroutine hermite_coef_nat_per(f,buf3d,N,d)
    sll_int32,intent(in)::N(2),d(2)
    sll_real64,dimension(N(1)+1,N(2)),intent(in)::f
    sll_real64,dimension(9,N(1)+1,N(2)+1),intent(out)::buf3d
    sll_real64 ::w_left_1(-d(1)/2:(d(1)+1)/2),w_right_1((-d(1)+1)/2:d(1)/2+1)
    sll_real64 ::w_left_2(-d(2)/2:(d(2)+1)/2),w_right_2((-d(2)+1)/2:d(2)/2+1)
    sll_real64 ::tmp
    sll_int32  ::i,j,ii,r_left(2),r_right(2),s_left(2),s_right(2),ind 
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
    sll_int32::i,j,ii,r_left(2),s_left(2),ind,dd(2) 
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


  subroutine sll_s_localize_polar(x,eta_min,eta_max,ii,eta,N)
    sll_real64,intent(in)::x(2),eta_min(2),eta_max(2)
    sll_int32,intent(out)::ii(2)
    sll_int32,intent(in)::N(2)
    sll_real64,intent(out)::eta(2)
    
    eta(1)=sqrt(x(1)**2+x(2)**2)
    call localize_nat(ii(1),eta(1),eta_min(1),eta_max(1),N(1))
    eta(2)=atan2(x(2),x(1))
    call localize_per(ii(2),eta(2),eta_min(2),eta_max(2),N(2))
  end subroutine sll_s_localize_polar

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
    sll_int32::s
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
  
  
  subroutine sll_s_contribution_spl(x,val)
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

  end subroutine sll_s_contribution_spl
  
  
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
  
  
   subroutine sll_s_splcoefnat1d0old(dnat,lnat,N)
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
  end subroutine sll_s_splcoefnat1d0old



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
  
  
  subroutine sll_s_splcoefper1d0old(dper,lper,mper,N)
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
  end subroutine sll_s_splcoefper1d0old
  
  
  
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





  subroutine compute_N_bounds_polar_circle(N_min,N_max,N,rho,eta_min,eta_max)
    ! version modifiée : N -> N+1
    sll_int32,intent(out)::N_min,N_max
    sll_int32,intent(in)::N
    sll_real64,intent(in)::rho(2),eta_min,eta_max
    sll_real64::delta_eta
    sll_int32::N_rho 
    delta_eta=(eta_max-eta_min)/real(N,f64)
    N_rho=floor(max(rho(1),rho(2))/delta_eta)+1
    N_min=1+N_rho
    N_max=N+1-N_rho
    if((N_min>N+1).or.(N_max<1).or.(N_min>N_max))then
      print *,'#Warning: rho is too big'
      print *,'#Bad computation of N_min and N_max'
      N_min=1
      N_max=N+1
    endif
  end subroutine compute_N_bounds_polar_circle



  subroutine compute_factor_bounds(f,f_init,eps,bounds,N_min,N_max)
    sll_real64,dimension(:,:),intent(in)::f,f_init
    sll_int32,intent(in)::N_min(2),N_max(2)
    sll_real64,intent(in)::eps
    sll_real64,intent(out)::bounds(2)
    sll_int32::i,j
    sll_real64::tmp
    
    bounds(1)=1.e10_f64
    bounds(2)=-bounds(1)
    do j=N_min(2),N_max(2)+1
      do i=N_min(1),N_max(1)+1
        if(abs(f_init(i,j))>eps)then
          tmp=f(i,j)/f_init(i,j)
          if(tmp>bounds(2))bounds(2)=tmp;
          if(tmp<bounds(1))bounds(1)=tmp;
        endif  
      enddo
    enddo   
  end subroutine compute_factor_bounds

  subroutine sll_s_compute_error(f,f_init,J_factor,err,N_min,N_max)
    sll_real64,dimension(:,:),intent(in)::f,f_init
    sll_int32,intent(in)::N_min(2),N_max(2)
    sll_real64,intent(in)::J_factor
    sll_real64,intent(out)::err(3)
    sll_int32::i,j
    sll_real64::tmp,delta
    
    err=0._f64
    
    do j=N_min(2),N_max(2)+1
      do i=N_min(1),N_max(1)+1
        tmp=f(i,j)-J_factor*f_init(i,j)
        err(1)=err(1)+abs(tmp)
        err(2)=err(2)+abs(tmp)**2
        err(3)=max(err(3),abs(tmp))        
      enddo
    enddo
    delta=1._f64/(real((N_max(2)-N_min(2)),f64)* real((N_max(1)-N_min(1)),f64))
    err(1)=err(1)*delta
    err(2)=sqrt(err(2)*delta)
      
  end subroutine sll_s_compute_error


  subroutine compute_error_1D(f,f_init,J_factor,err,N_min,N_max)
    sll_real64,dimension(:,:),intent(in)::f,f_init
    sll_int32,intent(in)::N_min(2),N_max(2)
    sll_real64,dimension(:),intent(in) :: J_factor
    sll_real64,intent(out)::err(3)
    sll_int32::i,j
    sll_real64::tmp,delta
    
    err=0._f64
    
    do j=N_min(2),N_max(2)+1
      do i=N_min(1),N_max(1)+1
        tmp=f(i,j)-J_factor(i)*f_init(i,j)
        err(1)=err(1)+abs(tmp)
        err(2)=err(2)+abs(tmp)**2
        err(3)=max(err(3),abs(tmp))        
      enddo
    enddo
    delta=1._f64/(real((N_max(2)-N_min(2)),f64)* real((N_max(1)-N_min(1)),f64))
    err(1)=err(1)*delta
    err(2)=sqrt(err(2)*delta)
      
  end subroutine compute_error_1D


  subroutine solution_polar_circle(rho,mode,eta_min,eta_max,val)
    sll_real64,intent(in)::rho(2),eta_min(2),eta_max(2)
    sll_int32,intent(in)::mode(2)
    sll_real64,intent(out)::val
    sll_real64::tmp
    !integer::i
    !logical::is_file
    val=0._f64
    if(abs(rho(1)-rho(2))>1.e-12)then
      print *,'#for the moment rho(1)=rho(2) is needed'
      return
    endif
    call sll_s_zero_bessel_dir_dir(mode,eta_min(1),eta_max(1),tmp)
    val = 0._f64 !temporary because DBESJ not recognized on helios and curie
    !val = DBESJN(0,tmp*rho(1)/eta_max(1))
    !print *,i,j,mode_max,alpha,tmp      
  end subroutine solution_polar_circle

  subroutine sll_s_compute_gamma0(mode,eta_min,eta_max,val)
    sll_real64,intent(in)::eta_min(2),eta_max(2)
    sll_int32,intent(in)::mode(2)
    sll_real64,intent(out)::val
    sll_real64::tmp!,mu,mu_max
    !integer::N_approx
    !logical::is_file
!    N_approx = 2**15
!    mu_max = 50._f64
!    sum1=0._f64
!    sum2=0._f64
    call sll_s_zero_bessel_dir_dir(mode,eta_min(1),eta_max(1),tmp)
    !print *,"tmp=",tmp
    !print *,"eta_max(1)=",eta_max(1)
    if(abs(tmp-3.9651944700162098_f64)>1.e-12)then
      print *,'#error : gamma0 not computed for this mode : ',mode
      return
    endif
    if(abs(eta_max(1)-18._f64)>1.e-12)then
      print *,'#error : gamma0 not computed for this eta_max : ',eta_max(1)
      return
    endif
    val = 0.95319247622058476357_f64
!    delta_x=real(real(mu_max,f64)/real(N_approx,f64),f64)
!    do i=1,N_approx/2-1
!      mu = real(2._f64*real(i,f64)*delta_x,f64)
!      sum1 = sum1 + DBESJN(0,tmp*sqrt(2._f64*mu)/eta_max(1))**2*exp(-mu)
!    enddo
!    do i=1,N_approx/2
!      mu = real((2._f64*real(i,f64)-1._f64)*delta_x,f64)
!      sum2 = sum2 + DBESJN(0,tmp*sqrt(2._f64*mu)/eta_max(1))**2*exp(-mu)
!    enddo
!    val = 2._f64*sum1 + 4._f64*sum2 + DBESJN(0,0._f64)**2 + DBESJN(0,tmp*sqrt(2._f64*mu_max)/eta_max(1))**2*exp(-mu_max)
!    val = val*real(delta_x/3._f64,f64)
  end subroutine sll_s_compute_gamma0

  
  
  subroutine sll_s_initialize_mu_quadr_for_phi( &
    quasineutral, &
    mu_quadr_for_phi_case, &
    N_mu_for_phi, &
    mu_max_for_phi, &    
    mu_points_user_defined, &
    mu_weights_user_defined, &
    N_mu_user_defined)
    
    type(sll_t_plan_qn_polar)                        :: quasineutral
    character(len=256), intent(in)                 :: mu_quadr_for_phi_case 
    sll_int32, intent(in)                          :: N_mu_for_phi
    sll_real64, intent(in)                         :: mu_max_for_phi
    sll_real64, dimension(:), intent(in), optional :: mu_points_user_defined
    sll_real64, dimension(:), intent(in), optional :: mu_weights_user_defined
    sll_int32, intent(in)                          :: N_mu_user_defined
    sll_real64 :: h
    sll_int32 :: ierr,i 



    if(mu_quadr_for_phi_case=="SLL_USER_DEFINED") then
      quasineutral%N_mu_for_phi = N_mu_user_defined
      SLL_ALLOCATE(quasineutral%mu_points_for_phi(1:N_mu_user_defined),ierr)
      SLL_ALLOCATE(quasineutral%mu_weights_for_phi(1:N_mu_user_defined),ierr)
      if( .not.( (present(mu_points_user_defined)) &
        .and.(present(mu_weights_user_defined)) )) then
        print *,'#provide mu_points_user_defined, mu_weights_user_defined'
        print *,'#in sll_s_initialize_mu_quadr_for_phi'
        stop
      endif  
!    else
!      if(&
!        (present(mu_points_user_defined)) &
!        .or.(present(mu_weights_user_defined)) &
!        )then
!        print *,'# do not provide mu_points_user_defined, mu_weights_user_defined'
!        print *,'#in sll_s_initialize_mu_quadr_for_phi'
!        stop       
!      endif
    endif  

    select case (mu_quadr_for_phi_case)    
      case ("SLL_USER_DEFINED")
        if(size(mu_points_user_defined) /= N_mu_user_defined) then
          print *,'#incompatible sizes for mu_points_user_defined :',size(mu_points_user_defined)
          print *,'#and N_mu_user_defined :',N_mu_user_defined
          print *,'#in sll_s_initialize_mu_quadr_for_phi'
          stop
        elseif(size(mu_weights_user_defined) /= N_mu_user_defined) then
          print *,'#incompatible sizes for mu_weights_user_defined :',size(mu_weights_user_defined)
          print *,'#and N_mu_user_defined :',N_mu_user_defined
          print *,'#in sll_s_initialize_mu_quadr_for_phi'
          stop
        else
          quasineutral%mu_points_for_phi(1:N_mu_user_defined)=mu_points_user_defined
          quasineutral%mu_weights_for_phi(1:N_mu_user_defined)=mu_weights_user_defined
        endif
      case ("SLL_LIN_LEE") 
        quasineutral%N_mu_for_phi = N_mu_for_phi
        SLL_ALLOCATE(quasineutral%mu_points_for_phi(1:N_mu_for_phi),ierr)
        SLL_ALLOCATE(quasineutral%mu_weights_for_phi(1:N_mu_for_phi),ierr)
        select case (N_mu_for_phi)
          case (1)
            quasineutral%mu_points_for_phi(1) = 1._f64
            quasineutral%mu_weights_for_phi(1) = 1._f64
          case (2)
            quasineutral%mu_points_for_phi(1) = 0.4167845_f64
            quasineutral%mu_points_for_phi(2) = 2.495154605_f64
            quasineutral%mu_weights_for_phi(1) = 0.7194_f64
            quasineutral%mu_weights_for_phi(2) = 0.2806_f64
          case (3)
            quasineutral%mu_points_for_phi(1) = 0.148131245_f64
            quasineutral%mu_points_for_phi(2) = 1._f64
            quasineutral%mu_points_for_phi(3) = 3.15959522_f64
            quasineutral%mu_weights_for_phi(1) = 0.3583_f64
            quasineutral%mu_weights_for_phi(2) = 0.5004_f64 
            quasineutral%mu_weights_for_phi(3) = 0.1413_f64 
          case default  
            print *,'# bad N_mu_user_defined in SLL_LIN_LEE (must be 1, 2 or 3) : ', N_mu_for_phi
            print *,'#in sll_s_initialize_mu_quadr_for_phi'
            stop  
        end select
      case ("SLL_RECTANGLES") 
        quasineutral%N_mu_for_phi = N_mu_for_phi
        SLL_ALLOCATE(quasineutral%mu_points_for_phi(1:N_mu_for_phi),ierr)
        SLL_ALLOCATE(quasineutral%mu_weights_for_phi(1:N_mu_for_phi),ierr)
        do i=1,N_mu_for_phi
          quasineutral%mu_points_for_phi(i) = mu_max_for_phi*real(i-1,f64)/real(N_mu_for_phi,f64)
          quasineutral%mu_weights_for_phi(i) = mu_max_for_phi/real(N_mu_for_phi,f64)
        enddo  
      case ("SLL_SIMPSON") 
        quasineutral%N_mu_for_phi = N_mu_for_phi
        SLL_ALLOCATE(quasineutral%mu_points_for_phi(1:N_mu_for_phi),ierr)
        SLL_ALLOCATE(quasineutral%mu_weights_for_phi(1:N_mu_for_phi),ierr)
        do i=1,N_mu_for_phi
          quasineutral%mu_points_for_phi(i) = mu_max_for_phi*real(i-1,f64)/real(N_mu_for_phi-1,f64)
        enddo     
        h = mu_max_for_phi/(3._f64*real(N_mu_for_phi,f64))
        quasineutral%mu_weights_for_phi(1) = h
        do i=1,(N_mu_for_phi-1)/2-1
          quasineutral%mu_weights_for_phi(2*i) = 4._f64*h
          quasineutral%mu_weights_for_phi(2*i+1) = 2._f64*h
        enddo 
        quasineutral%mu_weights_for_phi(N_mu_for_phi-1) = 4._f64*h  
        quasineutral%mu_weights_for_phi(N_mu_for_phi) = h
      case default
        print *,'#mu_quadr_for_phi_case not defined'
        print *,'#in sll_s_initialize_mu_quadr_for_phi'
        stop    
    end select
  end subroutine sll_s_initialize_mu_quadr_for_phi
  
  


end module sll_m_qn_2d_polar
