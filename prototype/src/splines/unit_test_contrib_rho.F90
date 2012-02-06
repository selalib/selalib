program contrib_rho_tester

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"


!function
!f(H)=mu/(sqrt(2._f64*sll_pi))*(2._f64-2._f64*xi)/(3._f64-2._f64*xi)*(1._f64+H/(1._f64-xi))*exp(-H)
!mu=0.92_f64
!xi=0.90_f64
!L=14.71_f64
!Hmin 
!Hmax
!integration
!1._f64-mu*(3._f64-2._f64*xi+2*eta_1)/(3._f64-2._f64*xi)*exp(-eta_1)

  use cubic_nonuniform_splines
  use numeric_constants
  use contrib_rho_module

  !use init1d_module
  
  implicit none
  
  sll_int ::N_int1,N_int2,i,j,N_cells_eta1,N_cells_eta2,N_size,N_store,s,jj!,Nr,Nth
  
  sll_int :: N_int1_quarter,N_int2_quarter
  
  sll_real64,dimension(:),allocatable::store_contrib_rho
  integer,dimension(:),allocatable::store_index_contrib_rho
  integer,dimension(:,:),pointer::size_contrib_rho
  sll_real64,dimension(:,:,:),pointer::integration_points
  
  
  sll_real64::dx1,dx2,x1_min,x1_max,x2_min,x2_max,x1,x2,eta1_min,eta1_max,eta2_min,eta2_max
  
  sll_real64,dimension(:),pointer::x2_min_tab,x2_max_tab
  
  sll_int :: test_case
  sll_real64 :: tmp
  
  sll_real64,dimension(:), pointer :: node_positions_eta1,node_positions_eta2,f_eta1,f_eta2
  sll_real64,dimension(:), pointer :: node_positions_eta1_non_unif
  sll_real64,dimension(:,:), pointer :: f,coef
  type(cubic_nonunif_spline_1D), pointer :: spl_eta1, spl_eta2 
  sll_int :: err,f_case,N_eta1_non_unif,ii
  sll_real64,dimension(:), pointer :: rho, rho_exact
  sll_real64 :: eta_1,eta_2
  sll_real64 :: h_min,h_max
  
  sll_real64 :: L,mu,xi
  
  sll_real64 :: M
  
  

  test_case=-1
  mu=0.92_f64
  xi=0.90_f64
  L=14.71_f64
  
  M = 1._f64
  
  f_case = 5
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
  
  
  SLL_ALLOCATE(x2_min_tab(1:N_int1+1),err)
  SLL_ALLOCATE(x2_max_tab(1:N_int1+1),err)
  
  x2_min_tab = x2_min
  x2_max_tab = x2_max

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


  if(test_case==0)then
    
    open(unit=900,file='integration_points.dat')
    read(900,*) x1,x2
    x1_min = x1
    x1_max = x2
    do s=1,N_int1+1
      read(900,*) i,x1
      x2_max_tab(i) =x1
      x2_min_tab(i) =-x1
    enddo
    do s=1,(N_int1+1)*(N_int2+1)    
      read(900,* ) i,j,x1,x2
      integration_points(1,j,i) = x1
      integration_points(2,j,i) = x2
      !print *,i,j,x1,x2
    enddo
    close(900)
  endif

  if(test_case==-1)then
    
    open(unit=900,file='integration_points_store.dat')
    !read(900,*) x1,x2
    !x1_min = x1
    !x1_max = x2
    !do s=1,N_int1+1
    !  read(900,*) i,x1
    !  x2_max_tab(i) =x1
    !  x2_min_tab(i) =-x1
    !enddo
    read(900,*) N_int1_quarter,N_int2_quarter,h_min,h_max,N_eta1_non_unif
    
    SLL_DEALLOCATE(integration_points,err)
    SLL_DEALLOCATE(x2_max_tab,err)
    SLL_DEALLOCATE(x2_min_tab,err)
    SLL_DEALLOCATE(size_contrib_rho,err)
    SLL_DEALLOCATE(rho, err)
    SLL_DEALLOCATE(rho_exact, err)
    
    N_int1 = 2*N_int1_quarter
    N_int2 = 2*N_int2_quarter
    
    SLL_ALLOCATE(integration_points(2,N_int2+1,N_int1+1),err)
    SLL_ALLOCATE(x2_max_tab(N_int1+1),err)
    SLL_ALLOCATE(x2_min_tab(N_int1+1),err)
    SLL_ALLOCATE(size_contrib_rho(-1:N_cells_eta2+1,N_int1+1),err)
    SLL_ALLOCATE(rho(N_int1+1), err)
    SLL_ALLOCATE(rho_exact(N_int1+1), err)
  
 
    !read(900,*) x1,x2
    !x1_min = x1
    !x1_max = x2
    x1_min = 0._f64
    x1_max= L
    do s=1,N_int1_quarter+1
      read(900,*) i,x1
      !print *,i,x1
      x2_max_tab(i) =x1
      x2_min_tab(i) =-x1
    enddo
    
    if(N_eta1_non_unif>=0)then
      SLL_ALLOCATE(node_positions_eta1_non_unif(N_eta1_non_unif),err)    
    endif
    
    do s=1,N_eta1_non_unif+1
      read(900,*) i, node_positions_eta1_non_unif(i)
    enddo
    
    !stop
    
    do s=1,(N_int1_quarter+1)*(N_int2_quarter+1)    
      read(900,* ) i,j,x1,x2
      !print *,i,j,x1,x2
      
      integration_points(1,j+N_int2_quarter,i) = x1
      integration_points(2,j+N_int2_quarter,i) = x2
      !print *,i,j,x1,x2
    enddo    
                    
    close(900)
    
    !print *,N_int1_quarter,N_int2_quarter
    
    !stop
    
    do i = 1,N_int1/2
      integration_points(1,N_int2+1,i)=integration_points(1,N_int2,i)
      integration_points(2,N_int2+1,i)=integration_points(2,N_int2,i)
    enddo
    
    
    do j = 1,N_int2/2+1
      integration_points(1,j,N_int1/2+1)=integration_points(1,j,N_int1/2)
      integration_points(2,j,N_int1/2+1)=integration_points(2,j,N_int1/2)
    enddo
    
    integration_points(1,N_int2+1,N_int1/2+1)=integration_points(1,N_int2+1,N_int1/2)
    integration_points(2,N_int2+1,N_int1/2+1)=integration_points(2,N_int2+1,N_int1/2)
    
    
    !stop
    

    do s=1,N_int1_quarter
      x2_max_tab(N_int1-s+2) =x2_max_tab(s)
      x2_min_tab(N_int1-s+2) =x2_min_tab(s)
      !print *,s,N_int1-s+1,x2_max_tab(s),x2_max_tab(N_int1-s+1)
    enddo

    ! A  ! C !  vmax 
    ! B  ! D !   0
    !0   L/2 L  vmin
    
    ! A is filled 
    
    
    ! fill B
    do i=1,N_int1_quarter+1      
      do j=1,N_int2_quarter
        integration_points(1,j,i)=integration_points(1,N_int2+2-j,i)
        integration_points(2,j,i)=-integration_points(2,N_int2+2-j,i)+2._f64*sll_pi
      enddo
    enddo
    
    ! fill C

    do i=1,N_int1_quarter      
      do j=1,N_int2_quarter+1
        integration_points(1,j+N_int2_quarter,i+N_int1_quarter+1)=integration_points(1,j+N_int2_quarter,-i+N_int1_quarter+1)
        integration_points(2,j+N_int2_quarter,i+N_int1_quarter+1)=-integration_points(2,j+N_int2_quarter,-i+N_int1_quarter+1)+sll_pi
      enddo
    enddo

    ! fill D

    do i=1,N_int1_quarter      
      do j=1,N_int2_quarter
        integration_points(1,j,i+N_int1_quarter+1)=integration_points(1,N_int2+2-j,i+N_int1_quarter+1)
        integration_points(2,j,i+N_int1_quarter+1)=-integration_points(2,N_int2+2-j,i+N_int1_quarter+1)+2._f64*sll_pi
      enddo
    enddo

    
    !eta1_min = integration_points(1,1,1)
    !eta1_max = eta1_min
    !eta2_min = integration_points(2,1,1)
    !eta2_max = eta2_min
    ! normalization of integration points
    integration_points(1,:,:) = (integration_points(1,:,:)-h_min)/(h_max-h_min)
    integration_points(2,:,:) = (integration_points(2,:,:))/(2.0_f64*sll_pi)
    
    eta1_min = h_min
    eta1_max = h_max
    eta2_min = 0._f64
    eta2_max = 2.0_f64*sll_pi
    
    
    do i=1,N_int1+1
      do j=1,N_int2+1
        x1=integration_points(1,j,i)
        x2=integration_points(2,j,i)
        !if(x1<eta1_min)then
        !  eta1_min = x1
        !endif
        !if(x2<eta2_min)then
        !  eta2_min = x2
        !endif
        !if(x1>eta1_max)then
        !  eta1_max = x1
        !endif
        !if(x2>eta2_max)then
        !  eta2_max = x2
        !endif
        !print *,i,j,x1,x2
      enddo
    enddo
    
    
    
    print *,'#',eta1_min,eta1_max,eta2_min,eta2_max
    
    x2_max_tab(N_int1_quarter+1)=x2_max_tab(N_int1_quarter)
    x2_min_tab(N_int1_quarter+1)=x2_min_tab(N_int1_quarter)
    !do i=1,N_int1+1
    !  print *,i,x2_max_tab(i),x2_min_tab(i)
    !enddo
    
    !print *,h_min,h_max
    
  endif

    

  
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

    open(unit=900,file='integration_points.dat')
    write(900,*) x1_min,x1_max
    do i=1,N_int1+1
      write(900,*) i,x2_max_tab(i)
    enddo
    
    do i=1,N_int1+1
      do j=1,N_int2+1
        write(900,* ) i,j,integration_points(1,j,i),integration_points(2,j,i)
      enddo  
    enddo
    close(900)

    
  endif
  
  if(test_case==2)then
    call random_number(integration_points(1:2,1:N_int2+1,1:N_int1+1))
  endif
  
  if(test_case==3)then
    x1_min=-5._f64
    x1_max=5._f64
  
    x2_min=-3._f64
    x2_max=3._f64

    x2_min_tab = x2_min
    x2_max_tab = x2_max

  
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
      if(N_eta1_non_unif>=0)then      
        eta_1 = real(N_eta1_non_unif,f64)*node_positions_eta1(i)
        ii=floor(eta_1)
        eta_1 = eta_1-real(ii,f64)
        eta_1 = (1._f64-eta_1)*node_positions_eta1_non_unif(ii+1)+eta_1*node_positions_eta1_non_unif(ii+2)
      endif  
      eta_2 = eta2_min+(eta2_max-eta2_min)*node_positions_eta2(j)
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
      if(f_case==5) then
        f(i,j) = mu/(sqrt(2._f64*sll_pi))*(2._f64-2._f64*xi)/(3._f64-2._f64*xi)*(1._f64+eta_1/(1._f64-xi))*exp(-M*eta_1)
      endif
    enddo
  enddo  
  
  do j=1,N_cells_eta2+1
     f_eta1 = f(1:N_cells_eta1+1,j)     
     call compute_spline_nonunif( f_eta1, spl_eta1)
     coef(-1:N_cells_eta1+1,j) = spl_eta1%coeffs(-1:N_cells_eta1+1)
  enddo

  print *,'#',eta1_min,eta1_max,h_min,h_max,exp(-h_max*M)
  
  !stop

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
      !if(N_eta1_non_unif<0)then
      !  eta_1 = eta1_min+(eta1_max-eta1_min)*eta_1
      !else
      !  eta_1 = real(N_eta1_non_unif,f64)*eta_1
      !  ii=floor(eta_1)
      !  eta_1 = eta_1-real(ii,f64)
      !  eta_1 = (1._f64-eta_1)*node_positions_eta1_non_unif(ii+1)+eta_1*node_positions_eta1_non_unif(ii+2)      
      !endif
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
      if(f_case==5) then
        tmp = mu/(sqrt(2._f64*sll_pi))*(2._f64-2._f64*xi)/(3._f64-2._f64*xi)*(1._f64+eta_1/(1._f64-xi))*exp(-eta_1)
      endif      
      !x1 = x1_min+real(i-1,f64)*dx1
      !x2 = x2_min+real(j-1,f64)*dx2
      !tmp = exp(-10._f64*eta_1**2)
      !tmp = exp(-10._f64*(x1**2+x2**2))
      rho_exact(i) = rho_exact(i)+tmp
      !print *,i,j,tmp,eta_1
    enddo
    !print *,''
    rho_exact(i) = rho_exact(i)*(x2_max_tab(i)-x2_min_tab(i))/real(N_int2+1,f64)
  enddo


  !stop
  do i=1,N_int1+1
    eta_1 = integration_points(1,N_int2/2+1,i)
    !if(N_eta1_non_unif<0)then
    !  eta_1 = eta1_min+(eta1_max-eta1_min)*eta_1
    !else
    !  eta_1 = real(N_eta1_non_unif,f64)*eta_1
    !  ii=floor(eta_1)
    !  eta_1 = eta_1-real(ii,f64)
    !  eta_1 = (1._f64-eta_1)*node_positions_eta1_non_unif(ii+1)+eta_1*node_positions_eta1_non_unif(ii+2)      
    !endif
    eta_1 = eta1_min+(eta1_max-eta1_min)*eta_1
    if(f_case==5)then
      !rho_exact(i)=1._f64-mu*(3._f64-2._f64*xi+2*eta_1)/(3._f64-2._f64*xi)*exp(-eta_1)
      !rho_exact(i)=mu*(3._f64-2._f64*xi+2*eta_1)/(3._f64-2._f64*xi)*exp(-eta_1)      
      rho_exact(i)=mu/(3._f64-2._f64*xi)*exp(-M*eta_1)*(2._f64*(1-xi+eta_1)/sqrt(M)+1._f64/(M*sqrt(M)))      
    endif  
  enddo
  
  !print *,node_positions_eta1_non_unif
  
  !print *,eta1_min,eta1_max
  !stop
  
  !we redefine the integration_points on [0,1], bu using the composition of the two transforms:
  ! uniform mesh-> non uniform mesh -> mapping
  if(test_case==-1)then
    integration_points(1,:,:) = eta1_min+(eta1_max-eta1_min)*integration_points(1,:,:)
    do i=1,N_int1+1
      do j=1,N_int2+1
        x1=integration_points(1,j,i)
        ii=1
        do while(node_positions_eta1_non_unif(ii)<=x1)
          ii = ii+1
        enddo
        ii=ii-1
        x1=(x1-node_positions_eta1_non_unif(ii))/(node_positions_eta1_non_unif(ii+1)-node_positions_eta1_non_unif(ii))
        integration_points(1,j,i) = (real(ii,f64)+x1)/real(N_eta1_non_unif)
      enddo
    enddo
  endif
  
  
  
  !compute splines
  
  size_contrib_rho =0
  
  
  
  !N_size =  
  N_size =  compute_contrib_rho(integration_points,N_int1,N_int2,size_contrib_rho,&
     & store_index_contrib_rho,store_contrib_rho,N_store,N_cells_eta1,N_cells_eta2)  
  
  print *,'#',N_size,N_store

  !stop
  
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
      do jj=1,size_contrib_rho(j,i)
        s=s+1
        tmp = tmp+coef(store_index_contrib_rho(s),j)*store_contrib_rho(s)
      enddo
    enddo
    rho(i) = tmp*(x2_max_tab(i)-x2_min_tab(i))/real(N_int2+1,f64)
  enddo  

  !change for rho(L/2)
  rho(N_int1/2+1)=0.5_f64*(rho(N_int1/2)+rho(N_int1/2+2))
  

  
  s=0
  do i=1,N_int1+1
    tmp=0._f64
    do j=-1,N_cells_eta2+1
      !print *,i,j,size_contrib_rho(j,i)
      do jj=1,size_contrib_rho(j,i)
        s=s+1
        !print *,jj,store_index_contrib_rho(s),store_contrib_rho(s)
        tmp = tmp+coef(store_index_contrib_rho(s),j)*store_contrib_rho(s)
      enddo
    enddo
    !rho(i) = tmp*(x2_max_tab(i)-x2_min_tab(i))/real(N_int2+1,f64)
    print *,x1_min+real(i-1,f64)*(x1_max-x1_min)/real(N_int1,f64),rho(i),rho_exact(i),&
    &2._f64*x2_max_tab(i),&
    &tmp/real(N_int2+1,f64)!tmp*dx2,tmp*(x2_max-x2_min)/real(N_int2+1,f64)   
  enddo  
  
  !print *,'rho0'
  
  !do j=1,N_int2/2+1
  !  x2 = real(j-1,f64)*(x2_max_tab(1)-x2_min_tab(1))/real(N_int2,f64)    
  !  print *,j,integration_points(1,j+N_int2/2,1),integration_points(2,j+N_int2/2,1),(x2**2/2._f64)/h_max
  !enddo


end program