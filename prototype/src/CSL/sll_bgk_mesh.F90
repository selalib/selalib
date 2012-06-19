module bgk_mesh_construction
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use numeric_constants
  use cubic_nonuniform_splines
  !use utils
  implicit none

contains 
  subroutine construct_bgk_mesh()
    use numeric_constants
    implicit none
  
  
  
  
  end subroutine construct_bgk_mesh  
  
  subroutine compute_bgk_phi(L,Nx,mu,xi,tab_phi,tab_dphi)
    use numeric_constants
    implicit none
    sll_real64,intent(in) ::L,mu,xi
    sll_int,intent(in) :: Nx
    sll_real64 :: xmin,xmax,x,dx,k2,phi(2)
    sll_int:: i,j
    sll_real64,dimension(:), pointer :: tab_phi,tab_dphi 
    xmin = 0._f64
    xmax = L/2._f64
    phi = 0._f64
    x   = xmin
    dx  = (xmax-xmin)/real(Nx,f64)
    k2  = mu*(1._f64-2._f64*xi+2._f64*phi(1))/(3._f64-2._f64*xi)*exp(-phi(1))
    tab_phi(1) = phi(1)
    tab_dphi(1)= phi(2)
    !write(10,* ) x, phi(1), phi(2),k2,v1,v2,1.d0/phi(2)
    do j=2, Nx+1
      x =xmin+ (j-1)*dx
      call Solve_Runge_K4(mu,xi,x, phi, dx, 2, Diff_Eq,0._f64,1._f64)
      k2=mu*(1._f64-2._f64*xi+2._f64*phi(1))/(3._f64-2._f64*xi)*exp(-phi(1))
      tab_phi(j) = phi(1)
      tab_dphi(j)= phi(2)
      ! write(10,* ) x, phi(1), phi(2),k2,v1,v2,1.d0/phi(2)
    enddo!j 
  
    !open(unit=10,file='phi.dat')
    !  write(10,*) Nx,L
    !  do j=1,Nx+1
    !    write(10,*) j,tab_phi(j),tab_dphi(j)
    !  enddo
    !close(10)
    
  end subroutine compute_bgk_phi
  
  subroutine construct_mesh_bgk(tab_phi,h_positions,tab_phi_positions,h_theta_positions,&
  Nx,Nen,Nen2,h_max,L,N_theta)
    sll_real64,dimension(:), pointer :: tab_phi,h_positions 
    sll_real64,dimension(:,:), pointer :: tab_phi_positions 
    sll_real64,dimension(:,:,:), pointer :: h_theta_positions 
    sll_int,intent(in)::Nx,Nen,Nen2,N_theta
    sll_int::jj,i,j,N_h,k
    sll_real64,intent(in) :: h_max,L
    sll_real64 :: delta_v,v,h,dx,xmin,xmax,x_factor,dtheta
    sll_real64 :: y_factor,total_length,tmp_loc,alpha,length
    
    
    N_h = Nen+Nen2
    
    xmin = 0._f64
    xmax= L/2._f64
    
    dx  = (xmax-xmin)/real(Nx,f64)
    
    !check compatibility between Nx and Nen
    jj= Nx/Nen
    if(jj==0)then
      print *,'bad compatibility between Nx=',Nx,' and Nen=',Nen
      print *,'Nx/Nen should be a non zero integer'
      stop
    endif
    if(jj*Nen/=Nx)then
      print *,'bad compatibility between Nx=',Nx,' and Nen=',Nen
      print *,'Nx/Nen should be a non zero integer'
      stop
    endif
     
    ! first define the h_positions for the grid
    jj=(Nx/Nen)
    do j=1,Nen+1
      h_positions(j)=tab_phi(1+(j-1)*jj)
    enddo


  
    delta_v = sqrt(2._f64*(h_positions(Nen+1)))/real(Nen,f64)
    do j=2,Nen
      v = real(j-1,f64)*delta_v
      !h_positions(j)=0.5_f64*v*v
    enddo
  
  
    do j=1,Nen2
      h_positions(Nen+1+j) = (sqrt(h_positions(Nen+1+j-1))+(sqrt(h_max)-sqrt(h_positions(Nen+1)))/real(Nen2,8))**2
    enddo

  
    print *,'#h_min error=',h_positions(1)
    print *,'#h_max error=',h_positions(N_h+1)-h_max
    h_positions(1)=0._f64
    h_positions(N_h+1)=h_max
  
  
    ! compute length for a given h
    do i=2,N_h+1
      h=h_positions(i)
      j=1
      do while((h-tab_phi(j)>=0._f64).and.(j<=Nx))
        tab_phi_positions(2,j) = sqrt(2._f64*(h-tab_phi(j)))
        tab_phi_positions(1,j) = xmin+real(j-1,f64)*dx
        j=j+1
      enddo
      if((j==Nx+1).and.(h-tab_phi(j)>=0._f64))then
        tab_phi_positions(1,j) = xmin+real(j-1,f64)*dx
        tab_phi_positions(2,j) = sqrt(2._f64*(h-tab_phi(j)))    
        j=j+1
      endif    
      jj=j-1
      !compute length
      y_factor = tab_phi_positions(1,jj)**2
      x_factor = tab_phi_positions(2,1)**2
      total_length = 0._f64
      do j=1,jj-1
        tmp_loc = y_factor*(tab_phi_positions(2,j+1)-tab_phi_positions(2,j))**2
        tmp_loc = tmp_loc+x_factor*(tab_phi_positions(1,j+1)-tab_phi_positions(1,j))**2
        total_length = total_length+sqrt(tmp_loc)    
      enddo
  
  
      if(total_length<=0._f64)then
        print *,'length is zero',total_length
        stop
      endif
      print *,'#length=',total_length
      !dtheta(i) =total_length/real(N_theta,f64)
      dtheta =total_length/real(N_theta,f64)
     !index_theta_positions(1,i) = 1
      length=0._f64
      j=1
      do k=2,N_theta
        do while(length<=real(k-1,8)*dtheta)
          tmp_loc = y_factor*(tab_phi_positions(2,j+1)-tab_phi_positions(2,j))**2
          tmp_loc = tmp_loc+x_factor*(tab_phi_positions(1,j+1)-tab_phi_positions(1,j))**2
          tmp_loc= sqrt(tmp_loc)
          length = length+tmp_loc
          j=j+1
        enddo  
        length = length-tmp_loc
        j=j-1
        !index_theta_positions(k,1) = j
        if((length>real(k-1,f64)*dtheta).or.(real(k-1,f64)*dtheta>=length+tmp_loc))then
          print *,'Problem of localization',j,length,real(k-1,f64)*dtheta,length+tmp_loc
          stop
        endif
        alpha= (real(k-1,f64)*dtheta-length)/tmp_loc
        if((alpha>=1._f64).or.(alpha<0._f64))then
          print *,'bad value of alpha=',alpha
        endif
      
        if(j>Nx)then
          print *,'j=',j,' Nx=',Nx
          print *,'length=',length,'Ntheta*dtheta=',real(N_theta,f64)*dtheta
          stop
        endif
        h_theta_positions(1,k,i) =  alpha*tab_phi_positions(1,j+1)+(1._f64-alpha)*tab_phi_positions(1,j)
        h_theta_positions(2,k,i) =  alpha*tab_phi_positions(2,j+1)+(1._f64-alpha)*tab_phi_positions(2,j)
        !print *,j,tmp,real(k-1,8)*dtheta(i),tmp+tmp_loc
      enddo
      !k=1
      h_theta_positions(1,1,i) =  tab_phi_positions(1,1)
      h_theta_positions(2,1,i) =  tab_phi_positions(2,1)
      !k=N_theta+1
      h_theta_positions(1,N_theta+1,i) =  tab_phi_positions(1,jj)
      h_theta_positions(2,N_theta+1,i) =  tab_phi_positions(2,jj)
    
    enddo

    open(unit=10,file='meshi.dat')
    !write(10,*) Nx,L
    do i= 1, N_h+1
      do k=1,N_theta+1
        write(10,* ) i,k, h_theta_positions(1,k,i), h_theta_positions(2,k,i)
      enddo
    enddo
    close(10)

    
    
  end subroutine construct_mesh_bgk
  
    

  

!*****************************************************************************************
!  Subroutine: To solve  the differential Equation with Runge Kutta 4 method             *
!*****************************************************************************************
  subroutine Solve_Runge_K4(mu,xi,X,Phi,dx,Neq,Diff_Eq,Y,C)
    use numeric_constants
    implicit none
    !External Diff_Eq 
    sll_int   :: Neq
    sll_real64   :: x, phi(Neq), dx,mu,xi,y,c
    sll_real64     :: Step1(Neq),Step2(Neq),Step3(Neq), Step4(Neq), phi_temp(Neq)
    sll_int   :: i

    call Diff_Eq(mu,xi,x,phi,Step1,Neq,y,c)
    do i=1, Neq
      phi_temp(i) = phi(i) + Step1(i)*(dx/2._f64)
    enddo
    call Diff_Eq(mu,xi,x+dx/2._f64,phi_temp,Step2,Neq,y,c)
    do i=1, Neq
      phi_temp(i) = phi(i) + Step2(i)*(dx/2._f64)
    enddo
    call Diff_Eq(mu,xi,x+dx/2._f64,phi_temp,Step3,Neq,y,c)
    do i=1, Neq
      phi_temp(i) = phi(i) + Step3(i)*dx
    enddo
    call Diff_Eq(mu,xi,x+dx,phi_temp,Step4,Neq,y,c)
    do i=1, Neq
      phi(i)=phi(i)+(dx/6._f64)*(Step1(i)+2._f64*Step2(i)+2._f64*Step3(i) + Step4(i))
    enddo
    return
  end subroutine Solve_Runge_K4
!*****************************************************************************************
!  Subroutine: To define  the differential Equation                                      *
!*****************************************************************************************
  subroutine Diff_Eq(mu,xi,X,Phi,Dphi, Neq,Y,C)
    use numeric_constants
    implicit none
    sll_real64   :: mu,xi,y,c
    sll_int  :: Neq
    sll_real64   :: x, phi(Neq), dphi(Neq)
    ! Differential Equations
    if((mu.ne.0).and.(xi.ne.0)) then
      dphi(1) =phi(2)
      dphi(2) = -mu*(3._f64-2._f64*xi+2._f64*phi(1))/(3._f64-2._f64*xi)*exp(-phi(1))+1._f64
    else
      dphi(1) =-C*phi(1)/Y
      dphi(2) =0._f64
    endif
    return
  end subroutine Diff_Eq




end module bgk_mesh_construction