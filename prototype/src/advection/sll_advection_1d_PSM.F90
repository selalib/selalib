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

! in development; should be at least cubic splines
! attached with computation of characteristics


module sll_module_advection_1d_PSM
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_boundary_condition_descriptors
use sll_module_advection_1d_base
use sll_module_characteristics_1d_base
use sll_module_interpolators_1d_base
use sll_constants
implicit none

  type,extends(sll_advection_1d_base) :: PSM_1d_advector
    sll_real64, dimension(:), pointer :: buf1d    
    sll_real64, dimension(:), pointer :: buf1d_out
    sll_real64, dimension(:), pointer :: dtab     
    sll_real64, dimension(:), pointer :: ltab     
    sll_real64, dimension(:), pointer :: mtab     
    sll_real64, dimension(:), pointer :: alphax   
    sll_real64, dimension(:), pointer :: prim     
    sll_real64, dimension(:), pointer :: bufout0  
    sll_real64, dimension(:), pointer :: bufout   
    sll_real64, dimension(:), pointer :: p        
    sll_real64 :: eta_min
    sll_real64 :: eta_max
    sll_int32 :: Npts
    sll_int32 :: Nbdr
  contains
    procedure, pass(adv) :: initialize => &
       initialize_PSM_1d_advector
    procedure, pass(adv) :: advect_1d => &
      PSM_advect_1d
    procedure, pass(adv) :: advect_1d_constant => &
      PSM_advect_1d_constant
    procedure, pass(adv) :: delete => delete_PSM_1d_advector
  end type PSM_1d_advector
   




contains
  function new_PSM_1d_advector( &
    Npts, &
    eta_min, &
    eta_max, &
    Nbdr) &  
    result(adv)      
    type(PSM_1d_advector), pointer :: adv
    sll_int32, intent(in) :: Npts
    sll_real64, intent(in) :: eta_min
    sll_real64, intent(in) :: eta_max
    sll_int32, intent(in), optional :: Nbdr
    sll_int32 :: ierr
    
    SLL_ALLOCATE(adv,ierr)
        
    call initialize_PSM_1d_advector(&
      adv, &
      Npts, &
      eta_min, &
      eta_max, &
      Nbdr)    
    
  end function  new_PSM_1d_advector


  subroutine initialize_PSM_1d_advector(&
    adv, &
    Npts, &
    eta_min, &
    eta_max, &
    Nbdr)    
    class(PSM_1d_advector), intent(inout) :: adv
    sll_int32, intent(in) :: Npts
    sll_real64, intent(in) :: eta_min
    sll_real64, intent(in) :: eta_max
    sll_int32, intent(in), optional :: Nbdr
    sll_int32 :: ierr
    sll_int32 :: i
    sll_real64 :: delta_eta
    sll_int32 :: Nx
    
    
    
    adv%Npts = Npts
    adv%eta_min = eta_min
    adv%eta_max = eta_max
    
    SLL_ALLOCATE(adv%buf1d(Npts),ierr)
    SLL_ALLOCATE(adv%buf1d_out(Npts),ierr)
    Nx = Npts-1
 
    
    if(present(Nbdr))then
      adv%Nbdr = Nbdr
    else
      adv%Nbdr = 100  
    endif
    
    
    !print *,Nx,adv%Nbdr
    
    
    
    SLL_ALLOCATE(adv%dtab(0:Nx-1),ierr)
    SLL_ALLOCATE(adv%ltab(0:Nx-2),ierr)
    SLL_ALLOCATE(adv%mtab(0:Nx-1),ierr)
    SLL_ALLOCATE(adv%alphax(0:Nx-1),ierr)
    SLL_ALLOCATE(adv%prim(-adv%Nbdr:Nx-1+adv%Nbdr),ierr)
    SLL_ALLOCATE(adv%bufout0(0:Nx-1+2*adv%Nbdr),ierr)
    SLL_ALLOCATE(adv%bufout(0:Nx),ierr)
    SLL_ALLOCATE(adv%p(0:Nx-1),ierr)


    adv%dtab(0) = 4._f64
    adv%mtab(0) = 1._f64
    do i=0,Nx-3
      adv%ltab(i)=1._f64/adv%dtab(i)
      adv%dtab(i+1)=4._f64-adv%ltab(i)
      adv%mtab(i+1)=-adv%mtab(i)/adv%dtab(i)      
    enddo
    adv%ltab(Nx-2)=(1._f64+adv%mtab(Nx-2))/adv%dtab(Nx-2)
    adv%dtab(Nx-1)=4._f64-adv%ltab(Nx-2)*(adv%mtab(Nx-2)+1._f64)
    do i=0,Nx-3
      adv%dtab(Nx-1)=adv%dtab(Nx-1)+adv%mtab(i)*adv%mtab(i+1)
    enddo  
    do i=0,Nx-1
      adv%dtab(i)=1._f64/adv%dtab(i)
    enddo
     
!  dtab= (double*) malloc(sizeof(double)*(Nx));ltab= (double*) malloc(sizeof(double)*(Nx-1));
!mtab= (double*) malloc(sizeof(double)*(Nx));  
!  dtab[0]=4.;mtab[0]=1.;for(i=0;i<Nx-2 ;i++){ltab[i]=1./dtab[i];dtab[i+1]=4.-ltab[i];mtab[i+1]=-mtab[i]/dtab[i];}
!  ltab[Nx-2]=(1.+mtab[Nx-2])/dtab[Nx-2];dtab[Nx-1]=4-ltab[Nx-2]*(mtab[Nx-2]+1);
!  for(i=0;i<Nx-2;i++)dtab[Nx-1]=dtab[Nx-1]+mtab[i]*mtab[i+1];for(i=0;i<Nx;i++)dtab[i]=1./dtab[i];

          
  end subroutine initialize_PSM_1d_advector

  subroutine PSM_advect_1d(&
    adv, &
    A, &
    dt, &
    input, &
    output)
    class(PSM_1d_advector) :: adv
    sll_real64, dimension(:), intent(in) :: A
    sll_real64, intent(in) :: dt 
    sll_real64, dimension(:), intent(in) :: input
    sll_real64, dimension(:), intent(out) :: output      
    sll_int32 :: Npts
    sll_int32 :: Nbdr
    sll_int32 :: i
    sll_int32 :: Nx
    sll_int32 :: s
    sll_real64 :: dtmp2
    sll_real64 :: x
    sll_real64 :: dx
    sll_int32 :: ix
    sll_int32 :: ix1
    sll_real64 :: result
    sll_real64 :: M
    sll_real64 :: w(0:2)
    sll_real64 :: tmp
    sll_real64 :: df0
    sll_real64 :: df1
    sll_real64 :: dt_loc
    
    dt_loc = 0.5_f64*dt
    
    Npts = adv%Npts
    
    Nbdr = adv%Nbdr
    
    Nx = Npts-1
    
    dtmp2 = 1._f64/(adv%eta_max-adv%eta_min)
    
    dx = 1._f64/real(Nx,f64)
    
    adv%alphax(0:Nx-1) = 0._f64

    
    !print *,size(A),Nx,Npts,size(adv%alphax)
    
    do i=0,Nx-1
      do s=0,9
        x = real(i,f64)*dx-adv%alphax(i)
        if(x>=1._f64)then
          x = x-1._f64
        endif
        if(x<0._f64)then
          x =x+1._f64
        endif
        if(x==1)then
          x = 0._f64
        endif
        ix = floor(x*real(Nx,f64))
        if((x>=1._f64).or.(x<0._f64))then
          print *,'#x too big or too small'
          print *,'#in  PSM_advect_1d'
          print *,x,Nx,ix,i,s,adv%alphax(i)
          stop
        endif
        ix1=ix+1
        if(ix1==Nx)then
          ix1=0
        endif
        
        !print *,'ix=',i,s,ix,ix1,x,dtmp2,adv%alphax(i),dt
        
        result=(1._f64-x)*A(ix+1)+x*A(ix1+1)
        adv%alphax(i) = result*dtmp2*dt_loc
      enddo
    enddo
        
    dtmp2=adv%alphax(Nx-1)
    do i=Nx-1,1,-1
      adv%alphax(i)=adv%alphax(i)+adv%alphax(i-1)
    enddo
    adv%alphax(0)=adv%alphax(0)+dtmp2
    


!  dtmp2=1./(xmax-xmin);dtmp3=1./(ymax-ymin);
!  #ifdef PTFX1D
!  for(i=0;i<Nx*Ny;i++)alphax[i]=0.;
!  //for(s=0;s<10;s++)for(i=0;i<Nx;i++)for(j=1;j<Ny;j++){
!  for(j=0;j<Ny;j++){p=&alphax[Nx*j];for(i=0;i<Nx;i++)for(s=0;s<10;s++){
!    //localization of x,y ix and iy
!    x=i*dx-p[i];if(x>=1.)x-=1.;if(x<0.)x+=1.;if(x>=1.)x-=1.;ix=(int)(x*Nx);
!    if(x>=1. || x<0.){fprintf(stderr,"x too big/small %1.1000lg i=%d j=%d\n",x,i,j);exit(1);}
!    assert(x>=0.&&x<1.);assert(ix>=0 &&ix<Nx);x=x*Nx-ix;assert(x>=0 &&x<1.);
!    ix1=ix+1;if(ix1==Nx)ix1=0;
!    //result=(1.-x)*Ey[j+(Ny+1)*ix]+x*Ey[j+(Ny+1)*ix1];alphax[(i+Nx*j)]=result*dtmp2*dt;
!    result=(1.-x)*Ey[ix+Nx*j]+x*Ey[ix1+Nx*j];p[i]=result*dtmp2*dt;
!  }}  

!  if(interp==5 || interp==11 || interp==6 || interp<0){
!    for(j=0;j<Ny;j++){p=&alphax[Nx*j];dtmp2=p[Nx-1];for(i=Nx-1;i>0;i--)p[i]=(p[i]+p[i-1]);p[0]=(p[0]+dtmp2);        
!  }

    
    adv%p(0:Nx-1) = input(1:Nx)
    M = sum(adv%p(0:Nx-1))
    adv%prim(0)=adv%p(0)
    do i=1,Nx-1
      adv%prim(i)=adv%prim(i-1)+adv%p(i)
    enddo
    
    do i=-Nbdr,-1
      adv%prim(i)=adv%prim(i+Nx)-M
    enddo
    do  i=Nx,Nx+Nbdr-1
      adv%prim(i)=adv%prim(i-Nx)+M
    enddo
    
    do i=0,Nx-1
      adv%bufout0(i+Nbdr)=3._f64*(adv%p(i)+adv%p(mod(i+Nx-1,Nx)))
    enddo
    do i=1,Nx-1
      adv%bufout0(i+Nbdr)= adv%bufout0(i+Nbdr)-adv%bufout0(i-1+Nbdr)*adv%ltab(i-1)
    enddo  
    do i=0,Nx-3
      adv%bufout0(Nx-1+Nbdr) = adv%bufout0(Nx-1+Nbdr)+adv%mtab(i+1)*adv%bufout0(i+Nbdr)
    enddo
    adv%bufout0(Nx-1+Nbdr)=adv%bufout0(Nx-1+Nbdr)*adv%dtab(Nx-1)
    adv%bufout0(Nx-2+Nbdr)=adv%dtab(Nx-2)*(adv%bufout0(Nx-2+Nbdr) &
      -(1._f64+adv%mtab(Nx-2))*adv%bufout0(Nx-1+Nbdr))
    do i=Nx-3,0,-1 
      adv%bufout0(i+Nbdr)=adv%dtab(i) &
        *(adv%bufout0(i+Nbdr)-adv%bufout0(i+1+Nbdr)-adv%mtab(i)*adv%bufout0(Nx-1+Nbdr))
    enddo    
    do i=0,Nbdr-1
      adv%bufout0(i)=adv%bufout0(Nx+i)
    enddo  
    do i=0,Nbdr-1
      adv%bufout0(Nx+Nbdr+i)=adv%bufout0(Nbdr+i)
    enddo


    

!    M=0.;for(i=0;i<Nx;i++)M+=p[i];prim[0]=p[0];
!    for(i=1;i<Nx;i++)prim[i]=prim[i-1]+p[i];
!    for(i=-Nbdr;i<0;i++)prim[i]=prim[i+Nx]-M;for(i=Nx;i<Nx+Nbdr;i++)prim[i]=prim[i-Nx]+M;
!    for(i=0;i<Nx;i++)bufout0[i+Nbdr]=3*(p[i]+p[(i+Nx-1)%Nx]); 
!     /*resolution of L phi = ftab*/    
!    for(i=1;i<Nx;i++) bufout0[i+Nbdr]-=bufout0[i-1+Nbdr]*ltab[i-1];
!    for(i=0;i<Nx-2;i++)bufout0[Nx-1+Nbdr]+=mtab[i+1]*bufout0[i+Nbdr];
!    /*resolution of U eta =phi*/
!    bufout0[Nx-1+Nbdr]*=dtab[Nx-1];
!    bufout0[Nx-2+Nbdr]=dtab[Nx-2]*(bufout0[Nx-2+Nbdr]-(1.+mtab[Nx-2])*bufout0[Nx-1+Nbdr]);
!    for(i=Nx-3;i>=0;i--) bufout0[i+Nbdr]=dtab[i]*(bufout0[i+Nbdr]-bufout0[i+1+Nbdr]-mtab[i]*bufout0[Nx-1+Nbdr]);
!    for(i=0;i<Nbdr;i++)bufout0[i]=bufout0[Nx+i];
!    for(i=0;i<Nbdr;i++)bufout0[Nx+Nbdr+i]=bufout0[Nbdr+i];


     do i=0,Nx
       x=real(i,f64)*dx-adv%alphax(mod(i,Nx))
       ix=floor(x*real(Nx,f64))
       x=x*real(Nx,f64)-real(ix,f64)
       if(x==1)then
         x = 0._f64
         ix = ix-1
       endif
       if((x<0._f64).or.(x>=1._f64))then
         print *,'#Localization problem in PSM_advect_1d'
         print *,x,ix,Nx,adv%alphax(mod(i,Nx))
         stop
       endif
       df0=adv%bufout0(ix+Nbdr)
       df1=adv%bufout0(ix+1+Nbdr)
       tmp=adv%p(mod(ix+Nx,Nx))
       w(0)=x*(x-1._f64)*(x-1._f64)
       w(1)=x*x*(x-1._f64)
       w(2)=x*x*(3._f64-2._f64*x)
       result=w(0)*df0+w(1)*df1+w(2)*tmp
       adv%bufout(i)=result+adv%prim(ix-1)       
     enddo
     
     do i=0,Nx-1
       output(i+1) = adv%bufout(i+1)-adv%bufout(i)
     enddo
     
     output(Nx+1) = output(1)

!    for(i=0;i<Nx+1;i++){
!      //i1=(i-1+Nx)%Nx;
!      x=(i)*dx-alphax[i%Nx];
!      ix=floor(x*Nx);
!      assert(x>-1.&&x<2.);assert(ix-1>=-Nbdr &&ix+1<Nx+Nbdr);x=x*Nx-ix;assert(x>=0 &&x<1.);
!      if(x>=2. || x<=-1.){fprintf(stderr,"x too big/small %1.1000lg i=%d j=%d\n",x,i,j);exit(1);}
!      df0=bufout0[ix+Nbdr];df1=bufout0[ix+1+Nbdr];tmp=p[(ix+Nx)%Nx];      
!      lim(tmp,df0,df1);
!      w[0]=x*(x-1)*(x-1);
!      w[1]=x*x*(x-1);
!      w[2]=x*x*(3-2*x);
!      result=w[0]*df0+w[1]*df1+w[2]*tmp;
!      bufout[i]=result+prim[ix-1];
!    }
!for(i=0;i<Nx;i++)p[i]=bufout[i+1]-bufout[i];

    
          
  end subroutine PSM_advect_1d


  subroutine PSM_advect_1d_constant(&
    adv, &
    A, &
    dt, &
    input, &
    output)
    class(PSM_1d_advector) :: adv
    sll_real64, intent(in) :: A
    sll_real64, intent(in) :: dt 
    sll_real64, dimension(:), intent(in) :: input
    sll_real64, dimension(:), intent(out) :: output      
    sll_int32 :: ierr
    
    !this version is not optimized
    

          
  end subroutine PSM_advect_1d_constant

  subroutine delete_PSM_1d_advector( adv )
    class(PSM_1d_advector), intent(inout) :: adv
    sll_int32 :: ierr
    SLL_DEALLOCATE(adv%buf1d,ierr)    
    SLL_DEALLOCATE(adv%buf1d_out,ierr)    
    SLL_DEALLOCATE(adv%dtab,ierr)         
    SLL_DEALLOCATE(adv%ltab,ierr)         
    SLL_DEALLOCATE(adv%mtab,ierr)         
    SLL_DEALLOCATE(adv%alphax,ierr)       
    SLL_DEALLOCATE(adv%prim,ierr)         
    SLL_DEALLOCATE(adv%bufout0,ierr)      
    SLL_DEALLOCATE(adv%bufout,ierr)       
    SLL_DEALLOCATE(adv%p,ierr)    
  end subroutine delete_PSM_1d_advector



end module sll_module_advection_1d_PSM
