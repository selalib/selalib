module sll_m_qn_2d_polar_precompute
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

! use F77_blas, only: &
!   dgemm, &
!   zgemm

! use F77_fftpack, only: &
!   zfftb, &
!   zfftf, &
!   zffti

! use F77_lapack, only: &
!   zgetrf, &
!   zgetri

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_qn_2d_polar, only: &
    sll_s_compute_splines_coefs_matrix_nat_1d, &
    sll_s_compute_splines_coefs_matrix_per_1d, &
    sll_s_contribution_spl, &
    sll_s_localize_polar, &
    sll_s_matrix_product_compf, &
    sll_s_splcoefnat1d0old, &
    sll_s_splcoefper1d0old

  implicit none

  public :: &
    sll_s_compute_qns_inverse_polar_splines, &
    sll_s_compute_qns_matrix_polar_splines, &
    sll_s_solve_qns_polar_splines

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains

  subroutine pre_precompute_double_gyroaverage_coeff_polar_splines( &
    r_min, &
    r_max, &
    num_cells_r, &
    num_cells_theta, &
    rho, &
    points, &
    N_points, &
    pre_compute_N, &
    size_pre_compute )
    sll_real64, intent(in) :: r_min
    sll_real64, intent(in) :: r_max
    sll_int32, intent(in) :: num_cells_r
    sll_int32, intent(in) :: num_cells_theta
    sll_real64, dimension(:), intent(in) :: rho
    sll_real64, dimension(:,:), intent(in) :: points
    sll_int32, intent(in) :: N_points
    sll_int32, dimension(:), intent(out) :: pre_compute_N
    sll_int32, intent(out) :: size_pre_compute
    sll_int32,dimension(:,:),allocatable :: buf
    sll_int32 ::i,k,ell_1,ell_2,ii(2),s,nb,ind(2)
    sll_real64::eta_star(2),eta(2),delta_eta(2),x(2)
    sll_int32 ::error,max_nb
    sll_real64 :: eta_min(2), eta_max(2)
    sll_int32 :: Nc(2)
    
    delta_eta(1)=(r_max-r_min)/real(num_cells_r,f64)
    delta_eta(2)=2._f64*sll_p_pi/real(num_cells_theta,f64)
    
    Nc(1) = num_cells_r
    Nc(2) = num_cells_theta
    
    eta_min(1) = r_min
    eta_min(2) = 0._f64
    eta_max(1) = r_max
    eta_max(2) = 2._f64*sll_p_pi
    
    SLL_ALLOCATE(buf(0:num_cells_r+2,0:num_cells_theta-1),error)      
    
    eta(2)=0._f64
    
    max_nb=0
    buf=0
    nb=0
    do i=1,Nc(1)+1
      eta(1)=eta_min(1)+real(i-1,f64)*delta_eta(1)       
      s=0
      do k=1,N_points
        x(1) = eta(1)*cos(eta(2))+rho(i)*points(1,k)
        x(2) = eta(1)*sin(eta(2))+rho(i)*points(2,k)
        call sll_s_localize_polar( &
          x, &
          eta_min, &
          eta_max, &
          ii, &
          eta_star, &
          Nc)
        do ell_2=-1,2
          ind(2)=modulo(ii(2)+ell_2,Nc(2))
          do ell_1=-1,2
            ind(1)=ii(1)+1+ell_1
            if(buf(ind(1),ind(2)).ne.i)then
              s=s+1
              buf(ind(1),ind(2))=i
            endif
          enddo
        enddo    
      enddo
      pre_compute_N(i)=s
      nb=nb+s
    enddo
    max_nb = max(max_nb,nb)
    
    size_pre_compute = max_nb    
  
    SLL_DEALLOCATE_ARRAY(buf,error)
  
  end subroutine pre_precompute_double_gyroaverage_coeff_polar_splines  



   subroutine precompute_double_gyroaverage_coeff_polar_splines( &
    r_min, &
    r_max, &
    num_cells_r, &
    num_cells_theta, &
    rho, &
    points, &
    N_points, &
    size_pre_compute, &
    pre_compute_index, &
    pre_compute_coeff_spl)
    sll_real64, intent(in) :: r_min
    sll_real64, intent(in) :: r_max
    sll_int32, intent(in) :: num_cells_r
    sll_int32, intent(in) :: num_cells_theta
    sll_real64, dimension(:), intent(in) :: rho
    sll_real64, dimension(:,:), intent(in) :: points
    sll_int32, intent(in) :: N_points
    sll_int32, intent(in) :: size_pre_compute
    sll_int32, dimension(:,:), intent(out) :: pre_compute_index
    sll_real, dimension(:), intent(out) :: pre_compute_coeff_spl
    sll_int32,dimension(:,:),allocatable :: buf
    sll_int32 ::i,j,k,ell_1,ell_2,ii(2),s,nb,ind(2)
    sll_real64::val(-1:2,-1:2),eta_star(2),eta(2),delta_eta(2),x(2)
    sll_int32 ::error
    sll_real64 :: eta_min(2), eta_max(2)
    sll_int32 :: Nc(2)
    
    delta_eta(1)=(r_max-r_min)/real(num_cells_r,f64)
    delta_eta(2)=2._f64*sll_p_pi/real(num_cells_theta,f64)
    
    Nc(1) = num_cells_r
    Nc(2) = num_cells_theta
    
    eta_min(1) = r_min
    eta_min(2) = 0._f64
    eta_max(1) = r_max
    eta_max(2) = 2._f64*sll_p_pi
    
    SLL_ALLOCATE(buf(0:num_cells_r+2,0:num_cells_theta-1),error)      
    
    eta(2)=0._f64
    
 
 
    
    buf=0
    nb=0
    s=0
    val=0._f64
    eta(2)=0._f64
    do i=1,Nc(1)+1
      eta(1)=eta_min(1)+real(i-1,f64)*delta_eta(1)       
      do k=1,N_points
        x(1) = eta(1)*cos(eta(2))+rho(i)*points(1,k)
        x(2) = eta(1)*sin(eta(2))+rho(i)*points(2,k)
        call sll_s_localize_polar(x,eta_min,eta_max,ii,eta_star,Nc)

        call sll_s_contribution_spl(eta_star,val)
        
        val=val*points(3,k)
        
        do ell_2=-1,2
          ind(2)=modulo(ii(2)+ell_2,Nc(2))
          do ell_1=-1,2
            ind(1)=ii(1)+1+ell_1
            j=buf(ind(1),ind(2))
            if(j<=nb)then
              s=s+1
              buf(ind(1),ind(2))=s
              pre_compute_coeff_spl(s)=val(ell_1,ell_2)
              pre_compute_index(1,s)=ind(1)
              pre_compute_index(2,s)=ind(2)                     
            else              
              pre_compute_coeff_spl(j)=pre_compute_coeff_spl(j)+val(ell_1,ell_2)
            endif
          enddo
        enddo
      enddo      
      nb=s 
    enddo




      
    SLL_DEALLOCATE_ARRAY(buf,error)
    return
    print*, size_pre_compute
  
  end subroutine precompute_double_gyroaverage_coeff_polar_splines  

  subroutine compute_contribution_matrix( &
    num_cells_r, &
    num_cells_theta, &
    pre_compute_N, &
    pre_compute_index, &
    pre_compute_coeff_spl, &
    size_pre_compute, &
    pointer_mat_contribution_circ)
    sll_int32, intent(in) :: num_cells_r
    sll_int32, intent(in) :: num_cells_theta
    sll_int32, dimension(:), intent(in) :: pre_compute_N
    sll_int32, dimension(:,:), intent(in) :: pre_compute_index
    sll_real64, dimension(:), intent(in) :: pre_compute_coeff_spl
    sll_int32, intent(in) :: size_pre_compute
    sll_real64, dimension(:,:,:), intent(out) :: pointer_mat_contribution_circ
    sll_int32 :: i
    sll_int32 :: s
    sll_int32 :: k
    sll_int32 :: ii(2)
    pointer_mat_contribution_circ = 0._f64
    s=0
    do i=1,num_cells_r+1
      do k=1,pre_compute_N(i)
        s=s+1
        ii(1)=pre_compute_index(1,s)+1
        ii(2)=modulo(pre_compute_index(2,s),num_cells_theta)+1
        pointer_mat_contribution_circ(ii(2),i,ii(1))= &
          pointer_mat_contribution_circ(ii(2),i,ii(1))+ &
          pre_compute_coeff_spl(s) 
      enddo
    enddo         
    return
    print*, size_pre_compute
  end subroutine compute_contribution_matrix

  subroutine compute_D_spl2D( &
    num_cells_r, &
    num_cells_theta, &
    D_spl2D)
    sll_int32, intent(in) :: num_cells_r
    sll_int32, intent(in) :: num_cells_theta
    sll_int32 :: ierr
    sll_int32 :: Nr
    sll_int32 :: Ntheta
    sll_comp64, dimension(:,:,:), intent(out) :: D_spl2D
    sll_real64,dimension(:),allocatable,target::dnat,lnat,dper,lper,mper
    sll_real64,dimension(:),pointer::pointer_dnat,pointer_lnat
    sll_real64,dimension(:),pointer::pointer_dper,pointer_lper,pointer_mper
    sll_real64,dimension(:,:),allocatable,target::mat_nat,mat_per
    sll_real64,dimension(:,:),pointer::pointer_mat_nat,pointer_mat_per
    sll_real64,dimension(:,:,:),allocatable,target::mat_spl2D_circ
    sll_real64,dimension(:,:,:),pointer::pointer_mat_spl2D_circ
    sll_int32 :: j
    !sll_int32 :: m
    !sll_comp64 :: exp_comp
    !sll_real64 :: mode
    sll_comp64, dimension(:), allocatable :: fft_array
    sll_real64, dimension(:), allocatable :: buf_fft
    sll_int32 :: k
    
    
    Nr = num_cells_r
    Ntheta = num_cells_theta

    SLL_ALLOCATE(dnat(0:Nr+2),ierr)
    SLL_ALLOCATE(lnat(0:Nr+2),ierr)
    SLL_ALLOCATE(dper(0:Ntheta-1),ierr)
    SLL_ALLOCATE(lper(0:Ntheta-1),ierr)
    SLL_ALLOCATE(mper(0:Ntheta-1),ierr)
    SLL_ALLOCATE(mat_nat(0:Nr+2,0:Nr),ierr)
    SLL_ALLOCATE(mat_per(0:Ntheta-1,0:Ntheta-1),ierr)
    SLL_ALLOCATE(mat_spl2D_circ(0:Ntheta-1,0:Nr+2,0:Nr),ierr)
    !SLL_ALLOCATE(D_spl2D(0:Ntheta-1,0:Nr+2,0:Nr),error)
    SLL_ALLOCATE(buf_fft(4*Ntheta+15),ierr)
    SLL_ALLOCATE(fft_array(Ntheta),ierr)



    pointer_dnat => dnat
    pointer_lnat => lnat
    pointer_dper => dper
    pointer_lper => lper
    pointer_mper => mper    
    pointer_mat_nat => mat_nat
    pointer_mat_per => mat_per
    pointer_mat_spl2D_circ => mat_spl2D_circ

    call sll_s_splcoefnat1d0old(dnat,lnat,Nr)
    call sll_s_splcoefper1d0old(dper,lper,mper,Ntheta)
 
    call sll_s_compute_splines_coefs_matrix_nat_1d(pointer_mat_nat,pointer_dnat,pointer_lnat,Nr)    
    call sll_s_compute_splines_coefs_matrix_per_1d(pointer_mat_per,pointer_dper,pointer_lper,pointer_mper,Ntheta) 
    do j=0,Ntheta-1  
      pointer_mat_spl2D_circ(j,:,:)=pointer_mat_per(0,j)*pointer_mat_nat
    enddo

    call zffti(Ntheta,buf_fft)

    do k=0,Nr
      do j=0,Nr+2
        fft_array(1:Ntheta)=pointer_mat_spl2D_circ(0:Ntheta-1,j,k)*(1._f64,0._f64)
    call zfftf(Ntheta,fft_array(1:Ntheta),buf_fft)
    D_spl2D(1:Ntheta,j+1,k+1) = fft_array(1:Ntheta)
    enddo  
    enddo   


    
    
!    D_spl2D = 0._f64
!    do m=0,Ntheta-1
!      do j=0,Ntheta-1
!        mode=real(-2._f64*sll_p_pi*real(j,f64)*real(m,f64)/real(Ntheta,f64),f64)
!        exp_comp = cmplx( cos(mode), sin(mode), kind=f64 )
!        D_spl2D(m,:,:) = D_spl2D(m,:,:) + pointer_mat_spl2D_circ(j,:,:)*exp_comp
!      enddo
!    enddo


  end subroutine compute_D_spl2D


  subroutine compute_D_contr( &
    num_cells_r, &
    num_cells_theta, &
    pointer_mat_contribution_circ, &
    D_contr)
    sll_int32, intent(in) :: num_cells_r
    sll_int32, intent(in) :: num_cells_theta
    sll_real64, dimension(:,:,:), intent(in) :: pointer_mat_contribution_circ
    sll_int32 :: ierr
    sll_int32 :: Nr
    sll_int32 :: Ntheta
    sll_comp64, dimension(:,:,:), intent(out) :: D_contr
    !sll_int32 :: m
    sll_int32 :: j
    !sll_real64 :: mode
    !sll_comp64 :: exp_comp
    sll_comp64, dimension(:), allocatable :: fft_array
    sll_real64, dimension(:), allocatable :: buf_fft
    sll_int32 :: k

    
    Nr = num_cells_r
    Ntheta = num_cells_theta

    SLL_ALLOCATE(buf_fft(4*Ntheta+15),ierr)
    SLL_ALLOCATE(fft_array(Ntheta),ierr)


    call zffti(Ntheta,buf_fft)

    do k=0,Nr+2
      do j=0,Nr
        fft_array(1:Ntheta)=pointer_mat_contribution_circ(1:Ntheta,j+1,k+1)*(1._f64,0._f64)
    call zfftf(Ntheta,fft_array(1:Ntheta),buf_fft)
    D_contr(1:Ntheta,j+1,k+1) = fft_array(1:Ntheta)
    enddo  
    enddo   


  
!    D_contr = 0._f64  
!    do m=0,Ntheta-1
!      do j=0,Ntheta-1
!        mode=real(-2._f64*sll_p_pi*real(j,f64)*real(m,f64)/real(Ntheta,f64),f64)
!        exp_comp = cmplx( cos(mode), sin(mode), kind=f64 )
!          D_contr(m,:,:) = D_contr(m,:,:) + pointer_mat_contribution_circ(j,:,:)*exp_comp
!      enddo
!    enddo

  
  
  end subroutine compute_D_contr

  subroutine compute_double_gyroaverage_matrix( &
    D_spl2D, &
    D_contr, &
    num_cells_r, &
    num_cells_theta, &
    mat)
    sll_comp64, dimension(:,:,:), intent(in) :: D_spl2D
    sll_comp64, dimension(:,:,:), intent(in) :: D_contr
    sll_int32, intent(in) :: num_cells_r
    sll_int32, intent(in) :: num_cells_theta
    sll_comp64, dimension(:,:,:), intent(out) :: mat
    sll_int32 :: Nr
    sll_int32 :: Ntheta
    sll_int32 :: m
    sll_int32 :: i
    sll_comp64, dimension(:,:), allocatable :: mat_stock1
    sll_comp64, dimension(:,:), allocatable :: mat_stock2
    sll_comp64, dimension(:,:), allocatable :: mat_stock3
    sll_int32 :: ierr
    
    
    Nr = num_cells_r
    Ntheta = num_cells_theta
    
    SLL_ALLOCATE(mat_stock1(Nr+1,Nr+1), ierr)
    SLL_ALLOCATE(mat_stock2(Nr+1,Nr+1), ierr)
    SLL_ALLOCATE(mat_stock3(Nr+1,Nr+1), ierr)
    
     
    mat = (0._f64,0._f64)
    do m=1,Ntheta
      mat_stock1 = (0._f64,0._f64)
      mat_stock2 = (0._f64,0._f64)
      call sll_s_matrix_product_compf( &
        D_contr(m,:,:), &
        Nr+1, &
        Nr+3, &
        D_spl2D(m,:,:), &
        Nr+3, &
        Nr+1, &
        mat_stock1)
 !     call matrix_product_comp( &
!        mat_stock1(:,:), &
!        Nr+1, &
!        Nr+1, &
!        mat_stock1(:,:), &
!        Nr+1, &
!        Nr+1, &
!        mat_stock2)
!
     !print *,'mat_stock1 l1',m,mat_stock1(1,1),mat_stock1(1,2)
     !print *,'mat_stock1 l1',m,mat_stock1(2,1),mat_stock1(2,2)

      !mat_stock3 = transpose(mat_stock1)
      mat_stock3 = mat_stock1
      call sll_s_matrix_product_compf( &
        mat_stock3, &
        Nr+1, &
        Nr+1, &
        mat_stock1, &
        Nr+1, &
        Nr+1, &
        mat_stock2)
     ! mat_stock2 = mat_stock1
        
!  if(m==1)then

        do i=1,Nr+1
          mat_stock2(i,i) =  mat_stock2(i,i) - (1._f64,0._f64)
        enddo


!      endif  
      mat(m,:,:) = -mat_stock2
      !mat(m,:,:) = mat_stock2


     !print *,'mat_stock2 l1',m,mat_stock2(1,1),mat_stock2(1,2)
     !print *,'mat_stock2 l1',m,mat_stock2(2,1),mat_stock2(2,2)


    enddo     

  end subroutine compute_double_gyroaverage_matrix


  subroutine compute_double_gyroaverage_matrix_blas( &
    D_spl2D, &
    D_contr, &
    num_cells_r, &
    num_cells_theta, &
    mat)
    sll_comp64, dimension(:,:,:), intent(in) :: D_spl2D
    sll_comp64, dimension(:,:,:), intent(in) :: D_contr
    sll_int32, intent(in) :: num_cells_r
    sll_int32, intent(in) :: num_cells_theta
    sll_comp64, dimension(:,:,:), intent(out) :: mat
    sll_int32 :: Nr
    sll_int32 :: Ntheta
    sll_int32 :: m
    sll_int32 :: i
    sll_comp64, dimension(:,:), allocatable :: mat_stock1
    sll_comp64, dimension(:,:), allocatable :: mat_stock2


    sll_comp64, dimension(:,:), allocatable :: mat_A
    sll_comp64, dimension(:,:), allocatable :: mat_B
    sll_comp64, dimension(:,:), allocatable :: mat_C

    sll_int32 :: ierr
    
    
    Nr = num_cells_r
    Ntheta = num_cells_theta
    
    SLL_ALLOCATE(mat_stock1(Nr+1,Nr+1), ierr)
    SLL_ALLOCATE(mat_stock2(Nr+1,Nr+1), ierr)
    SLL_ALLOCATE(mat_A(Nr+1,Nr+3), ierr)
    SLL_ALLOCATE(mat_B(Nr+3,Nr+1), ierr)
    SLL_ALLOCATE(mat_C(Nr+1,Nr+1), ierr)
    
     
    mat = (0._f64,0._f64)
    do m=1,Ntheta
      mat_stock1 = (0._f64,0._f64)
      mat_A(1:Nr+1,1:Nr+3) = D_contr(m,:,:)
      mat_B(1:Nr+3,1:Nr+1) = D_spl2D(m,:,:)
      !A MxK
      !B KxN
      !C MxN
      !M = Nr+1, K = Nr+3, N = Nr+1
      !ALPHA=1 BETA=0
      !ZGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M)       
      CALL ZGEMM( &
        'N', &
        'N', &
        Nr+1, &
        Nr+1, &
        Nr+3, &
        1._f64, &
        mat_A, &
        Nr+1, &
        mat_B, &
        Nr+3, &
        0._f64, &
        mat_C, &
        Nr+1)
      
!      call matrix_product_comp( &
!        D_contr(m,:,:), &
!        Nr+1, &
!        Nr+3, &
!        D_spl2D(m,:,:), &
!        Nr+3, &
!        Nr+1, &
!        mat_stock1)
      
      mat_stock2 = (0._f64,0._f64)
      mat_stock1(:,:) = mat_C(1:Nr+1,1:Nr+1)
      CALL ZGEMM( &
        'T', &
        'N', &
        Nr+1, &
        Nr+1, &
        Nr+1, &
        1._f64, &
        mat_C, &
        Nr+1, &
        mat_stock1(:,:), &
        Nr+1, &
        0._f64, &
        mat_stock2, &
        Nr+1)


!      call matrix_product_comp( &
!        transpose(mat_stock1(:,:)), &
!        Nr+1, &
!        Nr+1, &
!        mat_stock1(:,:), &
!        Nr+1, &
!        Nr+1, &
!        mat_stock2)
      do i=1,Nr+1
        mat_stock2(i,i) =  mat_stock2(i,i) - (1._f64,0._f64)
      enddo
      mat(m,:,:) = -mat_stock2
    enddo     

  end subroutine compute_double_gyroaverage_matrix_blas

  subroutine compute_double_gyroaverage_matrix_blas_real( &
    D_spl2D, &
    D_contr, &
    num_cells_r, &
    num_cells_theta, &
    mat)
    sll_real64, dimension(:,:,:), intent(in) :: D_spl2D
    sll_real64, dimension(:,:,:), intent(in) :: D_contr
    sll_int32, intent(in) :: num_cells_r
    sll_int32, intent(in) :: num_cells_theta
    sll_real64, dimension(:,:,:), intent(out) :: mat
    sll_int32 :: Nr
    sll_int32 :: Ntheta
    sll_int32 :: m
    sll_int32 :: i
    sll_real64, dimension(:,:), allocatable :: mat_stock1
    sll_real64, dimension(:,:), allocatable :: mat_stock2


    sll_real64, dimension(:,:), allocatable :: mat_A
    sll_real64, dimension(:,:), allocatable :: mat_B
    sll_real64, dimension(:,:), allocatable :: mat_C

    sll_int32 :: ierr
    
    
    Nr = num_cells_r
    Ntheta = num_cells_theta
    
    SLL_ALLOCATE(mat_stock1(Nr+1,Nr+1), ierr)
    SLL_ALLOCATE(mat_stock2(Nr+1,Nr+1), ierr)
    SLL_ALLOCATE(mat_A(Nr+1,Nr+3), ierr)
    SLL_ALLOCATE(mat_B(Nr+3,Nr+1), ierr)
    SLL_ALLOCATE(mat_C(Nr+1,Nr+1), ierr)
    
     
    mat = 0._f64
    do m=1,Ntheta
      mat_stock1 = 0._f64
      mat_A(1:Nr+1,1:Nr+3) = D_contr(m,:,:)
      mat_B(1:Nr+3,1:Nr+1) = D_spl2D(m,:,:)
      !A MxK
      !B KxN
      !C MxN
      !M = Nr+1, K = Nr+3, N = Nr+1
      !ALPHA=1 BETA=0
      !ZGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M)       
      CALL DGEMM( &
        'N', &
        'N', &
        Nr+1, &
        Nr+1, &
        Nr+3, &
        1._f64, &
        mat_A, &
        Nr+1, &
        mat_B, &
        Nr+3, &
        0._f64, &
        mat_C, &
        Nr+1)
      
!      call matrix_product_comp( &
!        D_contr(m,:,:), &
!        Nr+1, &
!        Nr+3, &
!        D_spl2D(m,:,:), &
!        Nr+3, &
!        Nr+1, &
!        mat_stock1)
      
      mat_stock2 = 0._f64
      mat_stock1(:,:) = mat_C(1:Nr+1,1:Nr+1)
      CALL DGEMM( &
        'T', &
        'N', &
        Nr+1, &
        Nr+1, &
        Nr+1, &
        1._f64, &
        mat_C, &
        Nr+1, &
        mat_stock1(:,:), &
        Nr+1, &
        0._f64, &
        mat_stock2, &
        Nr+1)


!      call matrix_product_comp( &
!        transpose(mat_stock1(:,:)), &
!        Nr+1, &
!        Nr+1, &
!        mat_stock1(:,:), &
!        Nr+1, &
!        Nr+1, &
!        mat_stock2)
      do i=1,Nr+1
        mat_stock2(i,i) =  mat_stock2(i,i) - 1._f64
      enddo
      mat(m,:,:) = -mat_stock2
    enddo     

  end subroutine compute_double_gyroaverage_matrix_blas_real




  subroutine sll_s_compute_qns_matrix_polar_splines( &
    r_min, &
    r_max, &
    num_cells_r, &
    num_cells_theta, &
    rho, &
    points, &
    N_points, &
    mat)
    sll_real64, intent(in) :: r_min
    sll_real64, intent(in) :: r_max
    sll_int32, intent(in) :: num_cells_r
    sll_int32, intent(in) :: num_cells_theta
    sll_real64, dimension(:), intent(in) :: rho
    sll_real64, dimension(:,:), intent(in) :: points
    sll_int32, intent(in) :: N_points
    sll_comp64, dimension(:,:,:), intent(out) :: mat
    sll_int32, dimension(:), allocatable :: pre_compute_N
    sll_int32 :: size_pre_compute
    sll_int32, dimension(:,:), allocatable :: pre_compute_index
    sll_real64, dimension(:), allocatable :: pre_compute_coeff_spl
    sll_real64, dimension(:,:,:), allocatable :: pointer_mat_contribution_circ
    sll_comp64, dimension(:,:,:), allocatable :: D_spl2D 
    sll_comp64, dimension(:,:,:), allocatable :: D_contr
    sll_int32 :: ierr 
    !sll_real64, dimension(:,:,:), allocatable :: real_mat
    !sll_real64, dimension(:,:,:), allocatable :: real_D_spl2D 
    !sll_real64, dimension(:,:,:), allocatable :: real_D_contr




    SLL_ALLOCATE(pre_compute_N(num_cells_r+1),ierr)
  
  
  
    call pre_precompute_double_gyroaverage_coeff_polar_splines( &
      r_min, &
      r_max, &
      num_cells_r, &
      num_cells_theta, &
      rho, &
      points, &
      N_points, &
      pre_compute_N, &
      size_pre_compute )
  
    !print *,'#size_pre_compute=',size_pre_compute
  
    SLL_ALLOCATE(pre_compute_index(2,size_pre_compute),ierr)
    SLL_ALLOCATE(pre_compute_coeff_spl(size_pre_compute),ierr)
  
    call precompute_double_gyroaverage_coeff_polar_splines( &
      r_min, &
      r_max, &
      num_cells_r, &
      num_cells_theta, &
      rho, &
      points, &
      N_points, &
      size_pre_compute, &
      pre_compute_index, &
      pre_compute_coeff_spl)

    !print *,'#precompute_double_gyroaverage_coeff_polar_splines done'

  
    SLL_ALLOCATE(pointer_mat_contribution_circ(num_cells_theta, num_cells_r+1, num_cells_r+3),ierr)

    call compute_contribution_matrix( &
      num_cells_r, &
      num_cells_theta, &
      pre_compute_N, &
      pre_compute_index, &
      pre_compute_coeff_spl, &
      size_pre_compute, &
      pointer_mat_contribution_circ)

    !print *,'#compute_contribution_matrix done'


  
    SLL_ALLOCATE(D_spl2D(num_cells_theta,num_cells_r+3,num_cells_r+1), ierr)

    call compute_D_spl2D( &
      num_cells_r, &
      num_cells_theta, &
      D_spl2D)

    !print *,'#compute_D_spl2D done'


    SLL_ALLOCATE(D_contr(num_cells_theta,num_cells_r+1,num_cells_r+3), ierr)

    call compute_D_contr( &
      num_cells_r, &
      num_cells_theta, &
      pointer_mat_contribution_circ, &
      D_contr)

    !print *,'#compute_D_contr done'


    !print *,'#complex values of D_spl2D',minval(aimag(D_spl2D)),maxval(aimag(D_spl2D))
    !print *,'#complex values of D_contr',minval(aimag(D_contr)),maxval(aimag(D_contr))
    !print *,'#real values of D_spl2D',minval(real(D_spl2D)),maxval(real(D_spl2D))
    !print *,'#real values of D_contr',minval(real(D_contr)),maxval(real(D_contr))



   !print *,'#begin compute_double_gyroaverage_matrix'

    call compute_double_gyroaverage_matrix( &
      D_spl2D, &
      D_contr, &
      num_cells_r, &
      num_cells_theta, &
      mat)
   !print *,'#end compute_double_gyroaverage_matrix'


!    call compute_double_gyroaverage_matrix_blas( &
!      D_spl2D, &
!      D_contr, &
!      num_cells_r, &
!      num_cells_theta, &
!      mat)

!    SLL_ALLOCATE(real_D_spl2D(num_cells_theta,num_cells_r+3,num_cells_r+1), ierr)
!    SLL_ALLOCATE(real_D_contr(num_cells_theta,num_cells_r+1,num_cells_r+3), ierr)
!    SLL_ALLOCATE(real_mat(num_cells_theta,num_cells_r+1,num_cells_r+1), ierr)
!    
!    real_D_spl2D = real(D_spl2D,f64)
!    real_D_contr = real(D_contr,f64)
!    
!    real_mat = 0._f64
!
!
!    
!    call compute_double_gyroaverage_matrix_blas_real( &
!      real_D_spl2D, &
!      real_D_contr, &
!      num_cells_r, &
!      num_cells_theta, &
!      real_mat)
!    
!    mat = real_mat*(1._f64,0._f64)
   

    !print *,'#compute_double_gyroaverage_matrix'



    
  end subroutine sll_s_compute_qns_matrix_polar_splines

  subroutine sll_s_compute_qns_inverse_polar_splines( &
    mat, &
    lambda, &
    num_cells_r, &
    num_cells_theta)
    sll_comp64, dimension(:,:,:), intent(inout) :: mat
    sll_real64, dimension(:), intent(in) :: lambda
    sll_int32, intent(in) :: num_cells_r
    sll_int32, intent(in) :: num_cells_theta
    sll_int32 :: ierr
    sll_int32, dimension(:), allocatable :: IPIV
    sll_comp64, dimension(:), allocatable :: WORK
    sll_int32 :: INFO
    sll_int32 :: i
    sll_int32 :: m
    
    
    
    SLL_ALLOCATE(IPIV(num_cells_r+1),ierr)
    SLL_ALLOCATE(WORK((num_cells_r+1)**2),ierr) 
    do i=1,num_cells_r+1
      do m = 1,num_cells_theta
        mat(m,i,i) = mat(m,i,i)+lambda(i)
      enddo  
   enddo

    do m=1,num_cells_theta
    call ZGETRF( &
    num_cells_r+1, &
    num_cells_r+1, &
    mat(m,:,:), &
    num_cells_r+1, &
    IPIV, &
    INFO)
    call ZGETRI( &
    num_cells_r+1, &
    mat(m,:,:), &
    num_cells_r+1, &
    IPIV, &
    WORK, &
    (num_cells_r+1)**2, &
    INFO)
    enddo

    print *,'#here is INFO'
    print *,INFO


  end subroutine sll_s_compute_qns_inverse_polar_splines

  subroutine sll_s_solve_qns_polar_splines( &
    mat, &
    phi, &
    num_cells_r, &
    num_cells_theta)
    sll_comp64, dimension(:,:,:), intent(in) :: mat
    sll_real64, dimension(:,:), intent(inout) :: phi
    sll_comp64,dimension(:,:),allocatable :: phi_comp
    sll_comp64,dimension(:,:),allocatable :: phi_old
    sll_real64,dimension(:),allocatable::buf_fft
    sll_int32, intent(in) :: num_cells_r
    sll_int32, intent(in) :: num_cells_theta
    sll_int32 :: Nr
    sll_int32 :: Ntheta
    sll_int32 :: m
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: ierr
    sll_comp64 :: result
    
   
    Nr = num_cells_r
    Ntheta = num_cells_theta
    
    SLL_ALLOCATE(phi_comp(1:Nr+1,1:Ntheta),ierr)
    SLL_ALLOCATE(phi_old(1:Nr+1,1:Ntheta),ierr)
    SLL_ALLOCATE(buf_fft(1:4*Ntheta+15),ierr)
 
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
        result = (0._f64,0._f64)
        do j = 1,Nr+1
          result = result + mat(m,i,j)*phi_old(j,m)
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
    SLL_DEALLOCATE_ARRAY(phi_comp,ierr)
    SLL_DEALLOCATE_ARRAY(phi_old,ierr)
    SLL_DEALLOCATE_ARRAY(buf_fft,ierr)
        
  end subroutine sll_s_solve_qns_polar_splines
  
  
!PN  subroutine precompute_matrix( &
!PN    Ti, &
!PN    r_min, &
!PN    r_max, &
!PN    num_cells_r, &
!PN    num_cells_theta, &
!PN    mu_points, &
!PN    mu_weights, &
!PN    N_mu, &
!PN    N_points, &
!PN    mat)
!PN    sll_real64,dimension(:),intent(in) :: Ti
!PN    sll_real64, intent(in) :: r_min
!PN    sll_real64, intent(in) :: r_max
!PN    sll_int32, intent(in) :: num_cells_r
!PN    sll_int32, intent(in) :: num_cells_theta
!PN    sll_real64, dimension(:), intent(in) :: mu_points
!PN    sll_real64, dimension(:), intent(in) :: mu_weights
!PN    sll_int32, intent(in) :: N_mu
!PN    sll_int32, intent(in) :: N_points
!PN    sll_real64, dimension(:,:,:), intent(out) :: mat
!PN    
!PN    
!PN  end subroutine precompute_matrix

end module sll_m_qn_2d_polar_precompute
