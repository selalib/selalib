!gfortran -O3 maillage.f90
!./a.out <in.param
Program Curvilinear_mesh2D

  Implicit None
  External Solve_Runge_K4, Diff_Eq, Bary
  !External Diff_Eq
  !External Bary
  Integer                 :: i,j,k,n,m, Nx, Neq,index,ii,s
  Integer                 :: Nen,Nen2,Nth
  Parameter                  (Neq = 2)
  Real*8                  :: L,mu,xi, Pi,k2,v1,v2,v3,alpha,y
  Real*8                  :: phi(Neq), x,x2, dx,r,dr,theta,dth,&
       Xmax,Xmin,fmin,fmax,dE, th_max, th_min,stheta
  Real*8                  :: Enmin, Enmax,xmin1,xmax1, Enl2
  real*8, allocatable     :: tab_phi(:), tab_dphi(:),xx(:,:),vv(:,:), xx_node(:,:),vv_node(:,:),&
       jacobian_cell(:,:), xx_cell(:,:),vv_cell(:,:),integration_points(:,:,:),vrho_max_tab(:)
  Real*8                  :: x_index,eps,En,eps_nw,eps_dich,Et, Air_T1, Air_T2
  Real*8                  :: c,f,fp,g,gp,a,b,gc,ga,K0,dx1,dif1,de2, dth2
  Character(6)            :: ichar
  Logical                 :: test,loc
  Integer                 :: Nth_in,  Nh1_in,  Nh2_in
  Real*8                  :: xrho_min, xrho_max, x0_min,x0_max,x0,hrho,xrho,vrho, h_ap
  Real*8                  :: dxrho,dvrho,  vrho_min, vrho_max,xrho_tmp,vrho_tmp
  Real*8                  :: baryc,CoefA, CoefB, CoefC,xA,vA,xB,vB,xD,vD, xM,vM,xC,vC
  Integer                 :: indxh, Nx_rho,Nv_rho,vrho_max_index, icount
  Namelist /param/ Nth_in, Nh1_in , Nh2_in, mu,&
       xi,L,Enmin, eps_dich,eps, eps_nw, Nx, loc,&
      Nv_rho,Nx_rho    
  read(*,NML=param);open(unit=900,file="in.param");write(900,NML=param);close(900)
  Pi   = 2.d0*asin(1.d0)
  ENL2 = 0.48708847265908323_8
  !Nx0  =0.d0
  if( mod(Nth_in,4).ne.0.d0) then
     Stop" Nth_in must be divisible by 4"
  endif
  vrho_min=0;  vrho_max=0    
  xrho_min=0;  xrho_max=0  
  baryc  =99.99
  Xmin   = 0.d0     ; Xmax   = L/2.d0
  Nen    = 2*Nh1_in ; Nen2   = 2*Nh2_in
  Nth    = 2*(Nth_in/4)
  !C= 0.21739130434782611 ! 1/4.6
  C= 1/4.6d0    
  !Enmin   = 0.001d0!0.01d0!1.0e-5
  !Enl2   = 0.48708847265908323
  test=.false.
  th_min =0.d0 ; th_max = Pi/2.d0
  dth=(th_max-th_min)/Nth
  index=0
  allocate(tab_phi(Nx+1))
  allocate(tab_dphi(Nx+1))
  allocate(xx(4*(Nth+1),Nen+Nen2+2))
  allocate(vv(4*(Nth+1),Nen+Nen2+2))
  allocate(xx_node(2*(Nth+2),(Nh1_in+Nh2_in)+2))
  allocate(vv_node(2*(Nth+2),(Nh1_in+Nh2_in)+2))
  allocate(xx_cell(2*Nth+2,(Nh1_in+Nh2_in)+2))
  allocate(vv_cell(2*Nth+2,(Nh1_in+Nh2_in)+2))
  allocate(jacobian_cell(2*Nth,(Nh1_in+Nh2_in)+1))
  allocate(integration_points(2,1:Nv_rho+1,1:Nx_rho+1))
  
  allocate(vrho_max_tab(1:Nx))
  
  tab_phi     = 0.d0 ; tab_dphi    = 0.d0
  phi(1)      = 0.d0 ; phi(2)      = 0.d0
  xx_cell     = 0.d0 ; xx_node     = 0.d0
  vv_cell     = 0.d0 ; vv_node     = 0.d0
  
  ! open(unit=10,file='rg.dat')
  x   = xmin
  dx  = (xmax-xmin)/Nx
  k2  = mu*(1.d0-2.d0*xi+2.d0*phi(1))/(3.d0-2.d0*xi)*exp(-phi(1))
  tab_phi(1) = phi(1)
  tab_dphi(1)= phi(2)
  !write(10,* ) x, phi(1), phi(2),k2,v1,v2,1.d0/phi(2)
  do j=2, Nx+1
     x =xmin+ (j-1)*dx
     call Solve_Runge_K4(mu,xi,x, phi, dx, 2, Diff_Eq,y,c)
     k2=mu*(1.d0-2.d0*xi+2.d0*phi(1))/(3.d0-2.d0*xi)*exp(-phi(1))
     tab_phi(j) = phi(1)
     tab_dphi(j)= phi(2)
     ! write(10,* ) x, phi(1), phi(2),k2,v1,v2,1.d0/phi(2)
  enddo!j 
  
  ! close(10)
  
  
  
  
  Enl2=tab_phi(Nx+1)
  dE=(Enl2-Enmin)/Nen
  Enmax=Enl2+Nen2*de
  de=(Enmax-Enmin)/(Nen+Nen2)

  print *,'#',Enl2,Enmax,exp(-Enl2),exp(-Enmax)
  
  !stop
  
  !************************************
  !   TEST
  !************************************
  if(test) then
     stop
     do i=1,Nen+Nen2
        En=Enmin+(i-1)*dE
        write(ichar,'(I3.3)'),i
        open(unit=11,file='ligne'//ichar//'.dat') 
        s=0
        do j=1, Nx+1
           x =xmin+ (j-1)*dx
           phi(1)=tab_phi(j); phi(2)=tab_dphi(j)
           if(En.ge. phi(1)) then
              v1 = sqrt(2.d0*(En-phi(1)));  v2 =-v1 ; v3=x
           else
              s=s+1
              v1 =0.d0; v2 =0.d0; index=i
              if(s.eq.1) then
                 x_index=x                
              endif
           endif
           write(11,* ) x,v1,v2     
        enddo!j
        close(11)
     enddo !i
  endif
  !************************************
  !   END TEST
  !************************************
  
  !************************************
  !    Axis points 
  !************************************
  print *,Enmin
  if(Enmin.ge. 1.d-9) then
     n=0
     r=L/4.d0  
     do m=1, Nen
        r=L/4.d0
        En= Enmin+(m-1)*de
        theta=0.d0     
        g =1.d0
        do while (abs(g).gt.eps_nw) 
           x=r*cos(theta)
           x2=(x-xmin)/dx
           ii=abs(floor(x2))
           alpha = x2-ii
           !g  = (alpha*tab_phi(ii+1)+(1.d0-alpha)*tab_phi(ii))-En
           !gp = (alpha*tab_dphi(ii+1)+(1.d0-alpha)*tab_dphi(ii))         
           !warning new definition
           g  = (alpha*tab_phi(ii+2)+(1.d0-alpha)*tab_phi(ii+1))-En
           gp = (alpha*tab_dphi(ii+2)+(1.d0-alpha)*tab_dphi(ii+1))         
           r=r-g/gp 
           n=n+1
        enddo !n  
        xx(1,m)=r
        vv(1,m)=0.d0
        xx(2*Nth+2,m)=L-r
        vv(2*Nth+2,m)=0.d0
        if(m==1)then
          print *,m,ii,alpha,g,r
        endif
     enddo !m 
  endif
  xrho_min=xx(1,1)     
  xx(1,Nen+1)=L/2.d0 ; vv(1,Nen+1)=0.d0
  do m=Nen+2, Nen+Nen2+1
     En= Enmin+(m-1)*de
     xx(1,m)=L/2.d0       ; vv(1,m)=sqrt(2.d0*(En-tab_phi(Nx+1)))   
     xx(2*Nth+1,m)=L/2.d0 ;vv(2*Nth+1,m)=sqrt(2.d0*(En-tab_phi(Nx+1))) 
  enddo
  !************************************
  !   End:  Axis points 
  !************************************
  En= Enmin
  r = (th_min+th_max)/2.d0 
  !open(unit=12,file='alpha.dat')
  do k=1,Nth+1
     theta= th_min+ (k-1)*dth
     g =1.d0
     gp=0.d0  
     if(k.ne. Nth+1) then 
        n=0
        do while (abs(g).gt.eps_nw) 
           x=r*cos(theta)
           x2=(x-xmin)/dx
           ii=abs(floor(x2))
           alpha = x2-ii
           !g  = r**2*sin(theta)**2/2.d0+(alpha*tab_phi(ii+1)+(1.d0-alpha)*tab_phi(ii))-En
           !gp = r*sin(theta)**2 +cos(theta)*(alpha*tab_dphi(ii+1)+(1.d0-alpha)*tab_dphi(ii))    
           !warning new definition
           g  = r**2*sin(theta)**2/2.d0+(alpha*tab_phi(ii+2)+(1.d0-alpha)*tab_phi(ii+1))-En
           gp = r*sin(theta)**2 +cos(theta)*(alpha*tab_dphi(ii+2)+(1.d0-alpha)*tab_dphi(ii+1))    
           r=r-g/gp 
           n=n+1
        enddo !n
     else 
        r=sqrt(2.d0*Enmin/sin(theta))
     endif
     xx(k,1)=r*cos(theta)
     if(mod(k,2).eq. 0) then 
        xx_cell(k/2,1)=r*cos(theta)/2.d0; vv_cell(k/2,1)=r*sin(theta)/2.d0
     endif
     vv(k,1)=r*sin(theta)
  enddo !k
  xx(Nth+1,1)=0.d0; vv(Nth+1,1)=sqrt(2*Enmin) 
  Enmax=Enl2+Nen2*de
  !close(12) 
  phi(1)=1.d0/(sqrt(2*(Enmin-tab_phi(1))))**(1.d0/c)
  do k=1,Nth+1
     !write(ichar,'(I3.3)'),k
     !open(unit=13,file='rgg'//ichar//'.dat')
     !open(unit=13,file='rg2.dat')
     xmin1 = xx(k,1); xmax1 = L/2.d0
     x     = xmin1
     K0    = vv(k,1)**(c-1.d0/c)   
     dx1   = (xmax1-xmin1)/real(Nx)
     ii    = floor((xmin1-xmin)/dx)
     phi(1)= 1.d0/vv(k,1)**(1.d0/c)
     v1 =1.d0/(phi(1))**(c)
     xx(4*Nth+5-k,1)=x  ; vv(4*Nth+5-k,1)=-v1
     xx(2*Nth+3-k,1)=L-x; vv(2*Nth+3-k,1)=v1
     xx(3*Nth+4-k,1)=L-x; vv(3*Nth+4-k,1)=-v1
     j=2
     x =xmin1+ (j-1)*dx
     m=2
     En=Enmin+1.d0*de
     v1=0.d0
     if(k.ne. Nth+1) then
        v2=v1- sqrt(2*(Enmax-tab_phi(ii)))
     else
        v2=v1- sqrt(2*(Enmax))
     endif
     do while((x.lt.L/2.d0).and.( v2.lt.0))
        x =xmin1+ (j-1)*dx1
        x2=(x-xmin)/(dx)
        ii=abs(floor(x2))
        alpha = x2-ii
        y=(alpha*tab_dphi(ii+1)+(1.d0-alpha)*tab_dphi(ii))      
        call Solve_Runge_K4(0.d0,0.d0,x, phi, dx1, 1, Diff_Eq,y,c)
        v1 =1.d0/(phi(1))**(c)  
        v2=v1- sqrt(2*(Enmax-tab_phi(ii)))
        if(En.ge. tab_phi(ii)) then
           v3=  sqrt(2.d0*(En-tab_phi(ii)))
        else
           v3=0.d0
        endif
        dif1=v1-v3
        !write(13,* ) x  , phi(1), v1,dif1
        if(dif1.ge. 0.d0) then
           xx(k,m)= x
           x2=(x-dx1-xmin)/(dx)
           ii=abs(floor(x2))
           alpha = x2-ii     
           vv(k,m)= v1
           xx(2*Nth-k+3,m)= L-x       ;  vv(2*Nth-k+3,m)= vv(k,m)
           xx(3*Nth-k+4,m)= L-x       ; vv(3*Nth-k+4,m)= -v1
           xx(4*Nth-k+5,m)= x         ; vv(4*Nth-k+5,m)= -v1
           xx(Nth+1,m)=0.d0           ; vv(Nth+1,m)=sqrt(2*En)
           xx(Nth+2,m)  =L            ;vv(Nth+2,m)  =vv(Nth+1,m)
           xx(2*Nth+2,m)=L-xx(1,m)    ; vv(2*Nth+2,m)=vv(1,m)
           xx(2*Nth+3,m)=L-xx(1,m)    ; vv(2*Nth+3,m)=-vv(1,m)
           xx(3*Nth+3,m)=L-xx(Nth+1,m); vv(3*Nth+3,m)=-vv(Nth+1,m)       
           xx(4*Nth+4,m)=xx(1,m)      ; vv(4*Nth+4,m)=-vv(1,m)
           xx(3*Nth+4,m)=xx(Nth+1,m)  ; vv(3*Nth+4,m)=-vv(Nth+1,m)
           m=m+1
           En=Enmin+(m-1)*de       
        endif
        j=j+1
     enddo!j   
     !close(13)  
  enddo !k  
  
  
  !stop
  !loc=.true.
  if(loc) then
     integration_points=-1.d0
     xrho_max=L/2.d0
     !**************************************************************************************
     !                   Rho: Localisation  using Barycentric coordinates 
     !**************************************************************************************
     !                    A=(k+1,m)                     
     !           C=(k+1,m)**********************C=(k+1,m+1)
     !                    *  *        T2       *
     ! eta2=theta     *       *            *
     !                    *    T1       *      *
     !                    *                  * *
     !            A=(k,m) **********************B=(k,m+1)
     !                       eta1=H
     Xrho  = 3.d0
     Vrho  = 1.d0/2.0d0
     !Nx_rho=1000
     !Nv_rho=1000
     CoefA = 0.d0
     CoefB = 0.d0
     CoefC = 0.d0
     vrho_max=sqrt(2.d0*Enmax)
     dxrho =(xrho_max-xrho_min)/Nx_rho
     dvrho =(vrho_max-vrho_min)/Nv_rho
     !print*,"xrho_max",xrho_max, dxrho, dvrho
     !stop
     !print*,xx(1,1)/dx,680*dx,xx(1,1),681*dx
     
     !stop
     
     !we redefine x(1,1) using piecewise linear representation of phi
     i=1
     do while(tab_phi(i)<=Enmin)
       i=i+1
     enddo
     i=i-1
     !print *,i,tab_phi(i),Enmin,tab_phi(i+1)
     alpha=(Enmin-tab_phi(i))/(tab_phi(i+1)-tab_phi(i))
     print *,xx(1,1),xmin+(real(i-1,8)+alpha)*dx,xx(1,1)-(xmin+(real(i-1,8)+alpha)*dx)
     print *,alpha,i
     !xx(1,1)=xmin+(real(i-1,8)+alpha)*dx
     !stop
     
     vrho_max_index=Nth+1
     
     !integration_points = 0._8
     
     do i=1,Nx_rho       
       xrho=xrho_min+(i-1)*dxrho
       x2=(xrho-xmin)/dx
       ii=abs(floor(x2))
       alpha = abs(x2-ii)
      ! if(xrho.lt.1.d-7) then
       if(ii.eq.0) then
        y=0.d0
       else
       !y=(alpha*tab_phi(ii+1)+(1.d0-alpha)*tab_phi(ii))
       ! warning change here
       y=(alpha*tab_phi(ii+2)+(1.d0-alpha)*tab_phi(ii+1))
       endif
       vrho_max=sqrt(2.d0*(Enmax-y))
       !warning change of the value of vrho_max
       !xx(Nth+1,Nen+Nen2+1),vv(Nth+1,Nen+Nen2+1)
       if(xx(vrho_max_index,Nen+Nen2+1)>xrho)then
         vrho_max_index=Nth+1
       endif
       do while(xx(vrho_max_index,Nen+Nen2+1)<=xrho)
         vrho_max_index=vrho_max_index-1
       enddo
       vrho_max_index=vrho_max_index+1
       k=vrho_max_index-1
       !print *,xx(k+1,Nen+Nen2+1),xrho,xx(k,Nen+Nen2+1)
       !print *,vrho_max , vv(k+1,Nen+Nen2+1)+(xrho-xx(k+1,Nen+Nen2+1))/(xx(k,Nen+Nen2+1)-xx(k+1,Nen+Nen2+1))&
       !   &*(vv(k,Nen+Nen2+1)-vv(k+1,Nen+Nen2+1))
       vrho_max = vv(k+1,Nen+Nen2+1)+(xrho-xx(k+1,Nen+Nen2+1))/(xx(k,Nen+Nen2+1)-xx(k+1,Nen+Nen2+1))&
          &*(vv(k,Nen+Nen2+1)-vv(k+1,Nen+Nen2+1))
       vrho_max_tab(i) = vrho_max   
       dvrho =(vrho_max-vrho_min)/Nv_rho
       !stop
      do j=1,Nv_rho+1
        vrho=vrho_min+(j-1)*dvrho
        !print *,vrho,sqrt(2.d0*(Enmin-y))
        !print *,tab_phi(ii),tab_phi(ii+1)
        !print *,xrho,xx(1,1)
        !print *,y,Enmin,Enmin-y
        !stop
        if(Enmin>=y)then
          if(vrho.lt.sqrt(2.d0*(Enmin-y))) then
            vrho = sqrt(2.d0*(Enmin-y))
          endif
        endif
        !stop
        hrho=(vrho**2)/2.d0+y    
        indxh=int((hrho-Enmin)/de)+1
        
        if((i==1).and.j==1)then
          xrho_tmp=xrho
          vrho_tmp=vrho
          print *,'xrho= ',xrho,xx(1,1),xrho-xx(1,1),L/2._8,'vrho= ',vrho,vrho_min+(j-1)*dvrho,vrho_max-vrho
          print *,'vrho_max= ',vrho_max,'Enmax',Enmax,'vrho_min',vrho_min
          print *,'hrho= ',hrho,'hmin= ',Enmin
          print *,'indxh= ',indxh
          print *,'h_inf',Enmin+real(indxh-1,8)*de
          print *,'h_next',Enmin+real(indxh,8)*de
          print *,'y=',y,y-Enmin,'Nth',Nth
          print *,'x(Nth+1)',xx(Nth+1,Nen+Nen2+1),vv(Nth+1,Nen+Nen2+1)
          print *,'x(Nth)',xx(Nth,Nen+Nen2+1),vv(Nth,Nen+Nen2+1)          
          print *,'val',vv(Nth+1,Nen+Nen2+1)+(xrho-xx(Nth+1,Nen+Nen2+1))/(xx(Nth,Nen+Nen2+1)-xx(Nth+1,Nen+Nen2+1))&
          &*(vv(Nth,Nen+Nen2+1)-vv(Nth+1,Nen+Nen2+1))
        !stop
        endif
        
       k=1
       do while(k.le. Nth)    
          xA=xx(k,indxh)    ; vA=vv(k,indxh)
          xB=xx(k,indxh+1)  ; vB=vv(k,indxh+1)
          xC=xx(k+1,indxh)  ; vC=vv(k+1,indxh)
          xD=xx(k+1,indxh+1); vD=vv(k+1,indxh+1)
          xM=xrho           ; vM=vrho
          call Bary(xA,vA,xB,vB,xC,vC, xM,vM,CoefA, CoefB, CoefC)
          baryc=CoefA+CoefB+CoefC
          
          if(abs(baryc-1).le. 1.d-10) then
             h_ap  =(CoefA*(enmin+(indxh-1)*de)+CoefB*(enmin+(indxh)*de)+CoefC*(enmin+(indxh-1)*de))
             theta =(CoefA*(th_min+(k-1)*dth)+CoefB*(th_min+(k-1)*dth)+CoefC*(th_min+(k)*dth))
             !write(*,*)"Tr1",baryc, k,theta,indxh,hrho,abs(hrho-h_ap)
             integration_points(1,j,i) = hrho
             integration_points(2,j,i) = theta             
             k=Nth+1
          else
             xA=xC; vA=VC
             xC=xD; vC=vD
             call Bary(xA,vA,xB,vB,xC,vC, xM,vM,CoefA, CoefB, CoefC)
             baryc=CoefA+CoefB+CoefC 
             if(abs(baryc-1).le. 1.d-10) then
                h_ap  =(CoefA*(enmin+(indxh-1)*de)+CoefB*(enmin+(indxh)*de)+CoefC*(enmin+(indxh)*de))
                theta =(CoefA*(th_min+(k)*dth)+CoefB*(th_min+(k-1)*dth)+CoefC*(th_min+(k)*dth))
                !write(*,*)"Tr2",baryc, k,theta,indxh,hrho, abs(hrho-h_ap)
                integration_points(1,j,i) = hrho
                integration_points(2,j,i) = theta                             
                k=Nth+1
               
             else
                k=k+1
                if((k==Nth+1).and.(integration_points(1,j,i)==-1))then
                  k=1
                  if(indxh<Nen+Nen2)then 
                    indxh=indxh+1
                    integration_points(1,j,i)=-2
                  else
                   integration_points(1,j,i)=-4                    
                  endif
                  
                endif
                if((k==Nth+1).and.(integration_points(1,j,i)==-2))then
                  k=1
                  if(indxh>2)then
                    indxh=indxh-2
                  endif  
                  integration_points(1,j,i)=-3
                endif
                if((k==Nth+1).and.(integration_points(1,j,i)==-4))then
                  k=1
                  if(indxh>2)then
                    indxh=indxh-1
                  endif  
                  integration_points(1,j,i)=-3
                endif                
                if((k==Nth+1).and.(integration_points(1,j,i)==-3))then
                  print *,'BIG PROBLEM ','i=',i,'j=',j
                  !stop
                endif
                
                !if(k>=Nth+1)then
                !  print *,'point not found ','i=',i,'j',j
                !  stop
                !endif
             endif
          endif
          
       enddo
    enddo
 enddo
 endif
 k=0
     do i=1,Nx_rho
      do j=1,Nv_rho+1
       xrho=integration_points(1,j,i)
       vrho=integration_points(2,j,i)
       if(xrho==-3)then
         print *,i,j,xrho,vrho
         k=1     
       endif
      enddo  
    enddo
 
 if(k==1)then
   print *,'localization problem'
   !stop
 endif
 
 
 
 open(unit=900,file='integration_points.dat')
     write(900,*),Nx_rho,Nv_rho,Enmin,Enmax
     
     do i=1,Nx_rho+1
       write(900,*) i,vrho_max_tab(i)  
     enddo
     
     do i=1,Nx_rho+1
      do j=1,Nv_rho+1
        write(900,* ) i,j,integration_points(1,j,i),integration_points(2,j,i)
      enddo  
    enddo

 close(900)
 
 !stop
 
 deallocate(tab_phi)
 deallocate(tab_dphi)
 open(unit=900,file='test_mesh.dat') 
 do i=1,Nth+1
   do j=1,Nen+Nen2+1
        write(900,* ) i,j,xx(i,j),vv(i,j)     
   enddo
 enddo
 close(900)

 open(unit=900,file='test_point.dat') 
        write(900,* ) xrho_tmp,vrho_tmp   
 close(900)


 !stop

 !*****************************************
 !   First part
 !*****************************************
 do k=1,Nth+1
    do m=1, Nen+Nen2+1
       ! write(*,* ),k/2,m/2, xx(k,m), vv(k,m)
       if ( (mod(k,2) .eq. 0 ).and.(mod(m,2) .eq. 0 )) then
          xx_cell(int(k/2.d0),int(m/2)+1)= xx(k,m)
          vv_cell(int(k/2.d0),int(m/2)+1)= vv(k,m)
          !write(*,* ),k,m, xx(k,m), vv(k,m)
          if(vv_cell(k/2,int(m/2)+1).le.1.d-10) then
             write(*,*),"Insufficient precision: You  must  increase the value of Nx"
             stop
          endif
       endif
       if ( (mod(k,2) .ne. 0).and.(mod(m,2) .ne. 0 )) then   
          xx_node(int(k/2.d0)+1,int(m/2)+1+1)= xx(k,m)
          vv_node(int(k/2.0) +1 ,int(m/2)+1+1)=  vv(k,m)
          !write(*,* ),k/2+1,m/2+1, xx(k,m), vv(k,m)
          !write(*,* ),k,m, xx(k,m), vv(k,m)
       endif
    enddo
 enddo

 !*****************************************
 !   Second part
 !*****************************************
 do k=Nth+2,2*Nth+2
    if(mod(k,2).ne. 0) then    
       xx_cell(k/2,1)=L-xx_cell(Nth+1-k/2,1)
       vv_cell(k/2,1)=  vv_cell(Nth+1-k/2,1)       
    endif
    do m=1, Nen+Nen2+1
       if ( (mod(k,2) .ne. 0 ).and.(mod(m,2) .eq. 0 )) then
          xx_cell(int(k/2.d0),int(m/2.d0)+1) = xx(k,m)
          vv_cell(int(k/2.d0),int(m/2.d0)+1) = vv(k,m)
          ! write(*,* ),k/2,m/2, xx(k,m), vv(k,m)
       endif
       
       if ( (mod(k,2) .eq. 0).and.(mod(m,2) .ne. 0 )) then
          xx_node(int(k/2.d0)+1,1)=L
          xx_node(int(k/2.d0)+1,int(m/2.d0)+1+1)= xx(k,m)
          vv_node(int(k/2.0) +1,int(m/2.d0)+1+1)= vv(k,m)
          ! write(*,* ),int(k/2.d0)+1,int(m/2.d0)+1, xx(k,m), vv(k,m)
       endif
    enddo
 enddo
 
 !*****************************************
 !   Third part
 !*****************************************
 do k=2*Nth+3,3*Nth+3 
    if(mod(k,2).eq. 0) then 
       xx_cell(k/2-1,1)=L-xx_cell(k/2-Nth-1,1)
       vv_cell(k/2-1,1)=  -vv_cell(k/2-Nth-1,1)   
       !print*,k,k/2-1,k/2-Nth-1         
    endif
    do m=1, Nen+Nen2+1 
       if ( (mod(k,2) .eq. 0 ).and.(mod(m,2) .eq. 0 )) then
          xx_cell(int(k/2.d0)-1,int(m/2.d0)+1)= xx(k,m)
          vv_cell(int(k/2.d0)-1,int(m/2.d0)+1)= vv(k,m)
          ! write(*,* ),k/2-1,m/2-1, xx(k,m), vv(k,m)
       endif
       
       if ( (mod(k,2) .ne. 0).and.(mod(m,2) .ne. 0 )) then
          xx_node(int(k/2.d0)+2,1)= L
          xx_node(int(k/2.d0)+2,int(m/2.d0)+1+1)= xx(k,m)
          vv_node(int(k/2.0) +2,int(m/2.d0)+1+1)= vv(k,m)
          ! write(*,* )"par3",int(k/2.0) +2,int(m/2.d0)+1, xx(k,m), vv(k,m)
       endif
    enddo
 enddo
 
 !*****************************************
 !  Fourth part
 !*****************************************
 do k=3*Nth+4,4*Nth+4
    if(mod(k,2).ne. 0) then 
       xx_cell(k/2-1,1) =xx_cell(2*Nth+2-k/2,1)
       vv_cell(k/2-1,1) =  -vv_cell(2*Nth+2-k/2,1)   
    endif
    do m=1, Nen+Nen2+1
       if ( (mod(k,2) .ne. 0 ).and.(mod(m,2) .eq. 0 )) then
          xx_cell(int(k/2.d0)-1,m/2+1) = xx(k,m)
          vv_cell(int(k/2.d0)-1,m/2+1) = vv(k,m)
          ! write(*,* ),k/2-1,m/2, xx(k,m), vv(k,m)
       endif
       if ( (mod(k,2) .eq. 0).and.(mod(m,2) .ne. 0 )) then
          xx_node(int(k/2.d0)+2,int(m/2.d0)+1+1) =  xx(k,m)
          vv_node(int(k/2.0)+2, int(m/2.d0)+1+1)  = vv(k,m)
          ! write(*,* )"part4",int(k/2.0)+2,int(m/2.d0)+1, xx(k,m), vv(k,m)
       endif
    enddo
 enddo
   
 !************************************* 
 ! Sqrt(g) At the center of mesh      *
 !*************************************

 ! (k,m+1)******************(k+1,k+1)
 !        *  *     T2      *
 !        *       *        *
 !        *  T1       *    *
 !        *              * *
 !  (k,m) ******************(k,m+1)
 
 de=2*de
 dth=2*dth
 do k=1,Nth_in/4 
    do m=1,(Nh1_in+Nh2_in)+1
       Air_T1=(1.d0/2.d0)*abs((xx_node(k+1,m)-xx_node(k,m))*(vv_node(k,m+1)-vv_node(k,m))&
            -(xx_node(k,m+1)-xx_node(k,m))*(vv_node(k+1,m)-vv_node(k,m))     )
       Air_T2=(1.d0/2.d0)*abs((xx_node(k+1,m+1)-xx_node(k+1,m))*(vv_node(k,m+1)-vv_node(k+1,m))&
            -(xx_node(k,m+1)-xx_node(k+1,m))*(vv_node(k+1,m+1)-vv_node(k+1,m)) )
       jacobian_cell(k,m)=(Air_T1+Air_T2)/(dth*de)
       jacobian_cell(Nth+1-k,m)=jacobian_cell(k,m)
       jacobian_cell(Nth+k,m)=jacobian_cell(k,m)
       jacobian_cell(2*Nth+1-k,m)=jacobian_cell(k,m)     
    enddo
 enddo

 deallocate(xx)
 deallocate(vv)
 !************************************* 
 !            Diagnostic              *
 !*************************************
 open(unit=12,file='xxvv_cell.dat')
 open(unit=13,file='jac_cell.dat')
 
!!$ do k=1,Nth_in
!!$    do m=1,Nh1_in+Nh2_in+1
!!$       write(12,* ),k,m, xx_cell(k,m), vv_cell(k,m)
!!$       write(13,* ),k,m, jacobian_cell(k,m)
!!$    enddo
!!$ enddo
 icount = 1
 do k= 1, Nth/2 
    do m=1,(Nh1_in+Nh2_in)+1
       write(12,* ),k,m, xx_cell(k,m), vv_cell(k,m)
       write(13,* ),k,m, jacobian_cell(k,m)
    enddo
    icount=icount+1
 enddo
 do k= Nth/2 + 1, Nth 
    do m=1,(Nh1_in+Nh2_in)+1
       write(12,* ),k,m, xx_cell(k,m)-L, vv_cell(k,m)
       write(13,* ),k,m, jacobian_cell(k,m)
    enddo
    icount=icount+1
 enddo
 k = Nth+1
 do m=1,(Nh1_in+Nh2_in)+1
    write(12,* ),k,m, xx_node(k+1,m)-L, -vv_node(k+1,m)
    write(13,* ),k,m, 0.0
 enddo
 icount=Nth+2
 do k= 3*Nth/2, Nth + 1,-1 
    do m=1,Nh1_in+Nh2_in+1
       write(12,* ), icount,m, xx_cell(k,m)-L, vv_cell(k,m)
       write(13,* ), icount,m, jacobian_cell(k,m)
    enddo
    icount=icount+1
 enddo
 do k= 3*Nth/2 + 1, 2*Nth
    do m=1,(Nh1_in+Nh2_in)+1
       write(12,* ),icount,m, xx_cell(k,m), vv_cell(k,m)
       write(13,* ), icount,m, jacobian_cell(k,m)
    enddo
    icount=icount+1
 enddo
 k = 2*Nth+2
 do m=1,(Nh1_in+Nh2_in)+1
    write(12,* ),k,m, xx_node(k+2,m), vv_node(k+2,m)
    write(13,* ),k,m, 0.0
 enddo

 close(12)
 close(13)
 
 open(unit=14,file='xxvv_node.dat')
 icount = 1
 do k= 1, Nth/2 + 1
    do m=1,(Nh1_in+Nh2_in)+1+1
       write(14,* ),k,m, xx_node(k,m), vv_node(k,m)
    enddo
    icount=icount+1
 enddo
 do k= Nth/2 + 3, Nth + 3
    do m=1,(Nh1_in+Nh2_in)+1+1
       write(14,* ),k-1,m, xx_node(k,m)-L, vv_node(k,m)
    enddo
    icount=icount+1
 enddo
 do k= 3*Nth/2 + 2, Nth + 4,-1 
    do m=1,Nh1_in+Nh2_in+2
       write(14,* ), icount ,m, xx_node(k,m)-L, vv_node(k,m)
    enddo
    icount=icount+1
 enddo
 do k= 3*Nth/2 + 4, 2*Nth + 4
    do m=1,(Nh1_in+Nh2_in)+1+1
       write(14,* ),icount,m, xx_node(k,m), vv_node(k,m)
    enddo
    icount=icount+1
 enddo
 close(14)
 write(*,*)'The curvilinear  mesh has been created'
 return
End program Curvilinear_mesh2D

!*****************************************************************************************
!  Subroutine: To solve  the differential Equation with Runge Kutta 4 method             *
!*****************************************************************************************
Subroutine Solve_Runge_K4(mu,xi,X,Phi,dx,Neq,Diff_Eq,Y,C)
  
  Implicit none
  External Diff_Eq 
  Integer   :: Neq
  Real*8    :: x, phi(Neq), dx,mu,xi,y,c
  Real*8    :: Step1(Neq),Step2(Neq),Step3(Neq), Step4(Neq), phi_temp(Neq)
  Integer   :: i

  call Diff_Eq(mu,xi,x,phi,Step1,Neq,y,c)
  do i=1, Neq
     phi_temp(i) = phi(i) + Step1(i)*(dx/2.d0)
  enddo
  call Diff_Eq(mu,xi,x+dx/2.d0,phi_temp,Step2,Neq,y,c)
  do i=1, Neq
     phi_temp(i) = phi(i) + Step2(i)*(dx/2.d0)
  enddo
  call Diff_Eq(mu,xi,x+dx/2.d0,phi_temp,Step3,Neq,y,c)
  do i=1, Neq
     phi_temp(i) = phi(i) + Step3(i)*dx
  enddo
  call Diff_Eq(mu,xi,x+dx,phi_temp,Step4,Neq,y,c)
  do i=1, Neq
     phi(i)=phi(i)+(dx/6.d0)*(Step1(i)+2.d0*Step2(i)+2.d0*Step3(i) + Step4(i))
  enddo
  return
End Subroutine Solve_Runge_K4
!*****************************************************************************************
!  Subroutine: To define  the differential Equation                                      *
!*****************************************************************************************
Subroutine Diff_Eq(mu,xi,X,Phi,Dphi, Neq,Y,C)
  Implicit None
  Real*8   :: mu,xi,y,c
  Integer  :: Neq
  Real*8   :: x, phi(Neq), dphi(Neq)
  ! Differential Equations
  if((mu.ne.0).and.(xi.ne.0)) then
     dphi(1) =phi(2)
     dphi(2) = -mu*(3.d0-2.d0*xi+2.d0*phi(1))/(3.d0-2.d0*xi)*exp(-phi(1))+1.d0
  else
     dphi(1) =-C*phi(1)/Y
     dphi(2) =0.d0
  endif
  return
End subroutine Diff_Eq
!*****************************************************************************************
!  Subroutine: To find  the barycentric coordinates                                      *
!*****************************************************************************************
Subroutine Bary(xA,vA,xB,vB,xC,vC, xM,vM, CoefA, CoefB, CoefC)
  Implicit None
  Real*8   :: xA,vA,xB,vB,xC,vC, xM,vM,xD,Vd
  Real*8   :: CoefA,CoefB,CoefC,Baryc
  Real*8   :: AirT_ABC, AirT_ABM,AirT_BCM,AirT_AMC
  
  AirT_ABC = (1.d0/2.d0)*abs((xC-xA)*(vB-vA) -(xB-xA)*(vC-vA))
  AirT_ABM = (1.d0/2.d0)*abs((xM-xA)*(vB-vA) -(xB-xA)*(vM-vA))
  AirT_BCM = (1.d0/2.d0)*abs((xM-xB)*(vC-vB) -(xC-xB)*(vM-vB))
  AirT_AMC = (1.d0/2.d0)*abs((xC-xA)*(vM-vA) -(xM-xA)*(vC-vA))  
  CoefA    = AirT_BCM/airT_ABC
  CoefB    = AirT_AMC/airT_ABC
  CoefC    = AirT_ABM/airT_ABC
  return    
End subroutine Bary
