program contrib_rho_tester

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use cubic_nonuniform_splines
  use numeric_constants
  use contrib_rho_module

  !use init1d_module
  
  implicit none
  
  sll_int ::N_int1,N_int2,i,j,N_cells_eta1,N_cells_eta2,N_size,N_store,s,jj!,Nr,Nth
  
  sll_real64,dimension(:),allocatable::store_contrib_rho
  integer,dimension(:),allocatable::store_index_contrib_rho
  integer,dimension(:,:),allocatable::size_contrib_rho
  sll_real64,dimension(:,:,:),allocatable::integration_points
  
  
  sll_real64::dx1,dx2,x1_min,x1_max,x2_min,x2_max,x1,x2,eta1_min,eta1_max,eta2_min,eta2_max
  
  sll_int :: test_case
  sll_real64 :: tmp
  
  sll_real64,dimension(:), pointer :: node_positions_eta1,node_positions_eta2,f_eta1,f_eta2
  sll_real64,dimension(:,:), pointer :: f,coef
  type(cubic_nonunif_spline_1D), pointer :: spl_eta1, spl_eta2 
  sll_int :: err,f_case
  sll_real64,dimension(:), pointer :: rho, rho_exact
  sll_real64 :: eta_1,eta_2
  
  test_case=3
  f_case = 4
  N_int1 = 32
  N_int2 = 2000
  
  N_cells_eta1 = 160
  N_cells_eta2 = 160
  
  N_store = N_cells_eta2 * N_int1 *5
  
  
  
  x1_min=0._f64
  x1_max=2._f64
  
  x2_min=-1._f64
  x2_max=3._f64
  
  x1_min=0._f64
  x1_max=1._f64
  
  x2_min=0._f64
  x2_max=1._f64
  
  
  dx1=(x1_max-x1_min)/real(N_int1,f64)
  dx2=(x2_max-x2_min)/real(N_int2,f64)
  
  SLL_ALLOCATE(integration_points(2,N_int2+1,N_int1+1),err)
  SLL_ALLOCATE(size_contrib_rho(-1:N_cells_eta2+1,N_int1+1),err)
  SLL_ALLOCATE(store_index_contrib_rho(N_store),err)
  SLL_ALLOCATE(store_contrib_rho(N_store),err)


  spl_eta1 =>  new_cubic_nonunif_spline_1D( N_cells_eta1, HERMITE_SPLINE)
  spl_eta2 =>  new_cubic_nonunif_spline_1D( N_cells_eta2, PERIODIC_SPLINE)

  SLL_ALLOCATE(node_positions_eta1(N_cells_eta1+1), err)
  SLL_ALLOCATE(node_positions_eta2(N_cells_eta2+1), err)
  SLL_ALLOCATE(f_eta1(N_cells_eta1+1), err)
  SLL_ALLOCATE(f_eta2(N_cells_eta2+1), err)
  SLL_ALLOCATE(f(N_cells_eta1+1,N_cells_eta2+1), err)
  SLL_ALLOCATE(coef(-1:N_cells_eta1+1,-1:N_cells_eta2+1), err)
  SLL_ALLOCATE(rho(N_int1+1), err)
  SLL_ALLOCATE(rho_exact(N_int1+1), err)

  do i=1,N_cells_eta1+1
    node_positions_eta1(i)=real(i-1,f64)/real(N_cells_eta1,f64)
  enddo

  do i=1,N_cells_eta2+1
    node_positions_eta2(i)=real(i-1,f64)/real(N_cells_eta2,f64)
  enddo

  
  ! compute initial function and its spline coefficients
  
  f_eta1 = 0._f64
  f_eta2 = 0._f64
  call compute_spline_nonunif( f_eta1, spl_eta1, node_positions_eta1)
  call compute_spline_nonunif( f_eta2, spl_eta2, node_positions_eta2)


  
  
  
  if(test_case==1)then
    do i=1,N_int1+1
      do j=1,N_int2+1
        integration_points(1,j,i) = x1_min+real(i-1,f64)*dx1
        integration_points(2,j,i) = x2_min+real(j-1,f64)*dx2
      enddo
    enddo
    eta1_min = x1_min
    eta1_max = x1_max
    eta2_min = x2_min
    eta2_max = x2_max
    
    integration_points(1,:,:) = (integration_points(1,:,:)-eta1_min)/(eta1_max-eta1_min)
    integration_points(2,:,:) = (integration_points(2,:,:)-eta2_min)/(eta2_max-eta2_min)
    
  endif
  
  if(test_case==2)then
    call random_number(integration_points(1:2,1:N_int2+1,1:N_int1+1))
  endif
  
  if(test_case==3)then
    x1_min=-5._f64
    x1_max=5._f64
  
    x2_min=-3._f64
    x2_max=3._f64
  
    dx1=(x1_max-x1_min)/real(N_int1,f64)
    dx2=(x2_max-x2_min)/real(N_int2,f64)
    !polar coordinates
    do i=1,N_int1+1
      do j=1,N_int2+1
        x1 = x1_min+real(i-1,f64)*dx1
        x2 = x2_min+real(j-1,f64)*dx2
        integration_points(1,j,i) = sqrt(x1*x1+x2*x2)
        if(abs(x1)>1.e-12)then
          integration_points(2,j,i) = datan(x2/x1)
        else
          integration_points(2,j,i) = 0._f64
        endif  
      enddo
    enddo
    
    eta1_min = integration_points(1,1,1) 
    eta1_max = eta1_min
    eta2_min = integration_points(2,1,1) 
    eta2_max = eta2_min
    do i=1,N_int1+1
      do j=1,N_int2+1
        x1=integration_points(1,j,i)
        if(x1>eta1_max)then
          eta1_max=x1
        endif
        if(x1<eta1_min)then
          eta1_min=x1
        endif
        x2=integration_points(2,j,i)
        if(x2>eta2_max)then
          eta2_max=x2
        endif
        if(x2<eta2_min)then
          eta2_min=x2
        endif        
      enddo
    enddo
    
    !print *,'#',eta1_min,eta1_max,eta2_min,eta2_max  
    
    integration_points(1,:,:) = (integration_points(1,:,:)-eta1_min)/(eta1_max-eta1_min)
    integration_points(2,:,:) = (integration_points(2,:,:)-eta2_min)/(eta2_max-eta2_min)
    
  endif

  do j=1,N_cells_eta2+1
    do i=1,N_cells_eta1+1
      eta_1 = eta1_min+(eta1_max-eta1_min)*node_positions_eta1(i)
      eta_2 = eta2_min+(eta2_max-eta2_min)*node_positions_eta2(i)
      if(f_case==1) then      
        f(i,j) = 1._f64        
      endif
      if(f_case==2) then
      f(i,j) = eta_1
      endif
      if(f_case==3) then
        f(i,j) = eta_1*cos(eta_2)
      endif
      if(f_case==4) then
        f(i,j) = exp(-10*eta_1**2)
      endif  
    enddo
  enddo  
  
  do j=1,N_cells_eta2+1
     f_eta1 = f(1:N_cells_eta1+1,j)     
     call compute_spline_nonunif( f_eta1, spl_eta1)
     coef(-1:N_cells_eta1+1,j) = spl_eta1%coeffs(-1:N_cells_eta1+1)
  enddo



  do i=-1,N_cells_eta1+1
     f_eta2 = coef(i,1:N_cells_eta2+1)
     call compute_spline_nonunif( f_eta2, spl_eta2)
     coef(i,-1:N_cells_eta2+1) = spl_eta2%coeffs(-1:N_cells_eta2+1)
  enddo

  
  do i=1,N_int1+1
    rho_exact(i)=0._f64
    do j=1,N_int2+1
      eta_1 = integration_points(1,j,i)
      eta_2 = integration_points(2,j,i)
      eta_1 = eta1_min+(eta1_max-eta1_min)*eta_1
      eta_2 = eta2_min+(eta2_max-eta2_min)*eta_2
      if(f_case==1) then      
        tmp = 1._f64        
      endif
      if(f_case==2) then
        tmp = eta_1
      endif
      if(f_case==3) then
        tmp = eta_1*cos(eta_2)
      endif
      if(f_case==4) then
        tmp = exp(-10*eta_1**2)
      endif  
      !x1 = x1_min+real(i-1,f64)*dx1
      !x2 = x2_min+real(j-1,f64)*dx2
      !tmp = exp(-10._f64*eta_1**2)
      !tmp = exp(-10._f64*(x1**2+x2**2))
      rho_exact(i) = rho_exact(i)+tmp
    enddo
    rho_exact(i) = rho_exact(i)*(x2_max-x2_min)/real(N_int2+1,f64)
  enddo
  
  
  !compute splines
  
  size_contrib_rho =0
  
  
  !N_size =  
  N_size =  compute_contrib_rho(integration_points,N_int1,N_int2,size_contrib_rho,&
     & store_index_contrib_rho,store_contrib_rho,N_store,N_cells_eta1,N_cells_eta2)  
  
  print *,'#',N_size,N_store
  
  if(N_size>N_store)then
    deallocate(store_index_contrib_rho)
    deallocate(store_contrib_rho)
    N_store=N_size
    allocate(store_index_contrib_rho(N_store))
    allocate(store_contrib_rho(N_store))    
    N_size =  compute_contrib_rho(integration_points,N_int1,N_int2,size_contrib_rho,&
     & store_index_contrib_rho,store_contrib_rho,N_store,N_cells_eta1,N_cells_eta2)  
  endif
  
  if(N_size>N_store)then
    print *,'N_size=',N_size,'N_store',N_store
    stop
  endif
  
  print *,'#',N_size,N_store
  
  
  !stop
  
  s=0
  do i=1,N_int1+1
    tmp=0._f64
    do j=-1,N_cells_eta2+1
      !print *,i,j,size_contrib_rho(j,i)
      do jj=1,size_contrib_rho(j,i)
        s=s+1
        !print *,j,store_index_contrib_rho(s),store_contrib_rho(s)
        tmp = tmp+coef(store_index_contrib_rho(s),j)*store_contrib_rho(s)
      enddo
    enddo
    rho(i) = tmp*(x2_max-x2_min)/real(N_int2+1,f64)
    print *,x1_min+real(i-1,f64)*dx1,rho(i),rho_exact(i)!tmp*dx2,tmp*(x2_max-x2_min)/real(N_int2+1,f64)   
  enddo  
  
  

end program