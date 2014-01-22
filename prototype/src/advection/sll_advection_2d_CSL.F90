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

!in development

module sll_module_advection_2d_CSL
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_boundary_condition_descriptors
use sll_module_advection_2d_base
use sll_module_characteristics_2d_base
use sll_module_interpolators_2d_base
implicit none

  type,extends(sll_advection_2d_base) :: CSL_2d_advector
  
    class(sll_characteristics_2d_base), pointer  :: charac
    sll_real64, dimension(:), pointer :: eta1_coords
    sll_real64, dimension(:), pointer :: eta2_coords
    sll_real64, dimension(:,:), pointer :: charac_feet1
    sll_real64, dimension(:,:), pointer :: charac_feet2
    sll_int32 :: Npts1
    sll_int32 :: Npts2  
    sll_int32 :: nbmax
    sll_real64, dimension(:,:), pointer :: tt
    sll_real64, dimension(:,:), pointer :: tcell
    sll_real64, dimension(:), pointer :: intx
    sll_real64, dimension(:), pointer :: inty
    sll_real64, dimension(:,:), pointer :: dir
    sll_int32, dimension(:,:,:), pointer :: cell
     

  contains
     procedure, pass(adv) :: initialize => &
       initialize_CSL_2d_advector
    procedure, pass(adv) :: advect_2d => &
      CSL_advect_2d
  
  end type CSL_2d_advector
   




contains
  function new_CSL_2d_advector( &
    charac, &
    Npts1, &
    Npts2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    eta1_coords, &
    eta2_coords) &  
    result(adv)      
    type(CSL_2d_advector), pointer :: adv
    class(sll_characteristics_2d_base), pointer  :: charac
    sll_int32, intent(in) :: Npts1
    sll_int32, intent(in) :: Npts2
    sll_real64, intent(in), optional :: eta1_min
    sll_real64, intent(in), optional :: eta1_max
    sll_real64, intent(in), optional :: eta2_min
    sll_real64, intent(in), optional :: eta2_max
    sll_real64, dimension(:), pointer, optional :: eta1_coords
    sll_real64, dimension(:), pointer, optional :: eta2_coords
    sll_int32 :: ierr
    
    SLL_ALLOCATE(adv,ierr)
        
    call initialize_CSL_2d_advector(&
      adv, &
      charac, &
      Npts1, &
      Npts2, &
      eta1_min, &
      eta1_max, &
      eta2_min, &
      eta2_max, &
      eta1_coords, &
      eta2_coords)    
    
  end function  new_CSL_2d_advector


  subroutine initialize_CSL_2d_advector(&
    adv, &
    charac, &
    Npts1, &
    Npts2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    eta1_coords, &
    eta2_coords)    
    class(CSL_2d_advector), intent(inout) :: adv
    class(sll_characteristics_2d_base), pointer  :: charac
    sll_int32, intent(in) :: Npts1
    sll_int32, intent(in) :: Npts2
    sll_real64, intent(in), optional :: eta1_min
    sll_real64, intent(in), optional :: eta1_max
    sll_real64, intent(in), optional :: eta2_min
    sll_real64, intent(in), optional :: eta2_max
    sll_real64, dimension(:), pointer, optional :: eta1_coords
    sll_real64, dimension(:), pointer, optional :: eta2_coords
    sll_int32 :: ierr
    sll_int32 :: i
    sll_real64 :: delta_eta1
    sll_real64 :: delta_eta2
    
    
    adv%Npts1 = Npts1
    adv%Npts2 = Npts2
    adv%charac => charac
    !SLL_ALLOCATE(adv%x1_mesh(Npts1),ierr)
    !SLL_ALLOCATE(adv%x2_mesh(Npts2),ierr)
    SLL_ALLOCATE(adv%eta1_coords(Npts1),ierr)
    SLL_ALLOCATE(adv%eta2_coords(Npts2),ierr)

    SLL_ALLOCATE(adv%charac_feet1(Npts1,Npts2),ierr)
    SLL_ALLOCATE(adv%charac_feet2(Npts1,Npts2),ierr)

    if(present(eta1_min).and.present(eta1_max))then
      if(present(eta1_coords))then
        print *,'#provide either eta1_coords or eta1_min and eta1_max'
        print *,'#and not both in subroutine initialize_CSL_2d_advector'
        stop
      else
        delta_eta1 = (eta1_max-eta1_min)/real(Npts1-1,f64)
        do i=1,Npts1
          adv%eta1_coords(i) = eta1_min+real(i-1,f64)*delta_eta1
        enddo
      endif
    else if(present(eta1_coords))then
      if(size(eta1_coords,1)<Npts1)then
        print *,'#bad size for eta1_coords in initialize_CSL_2d_advector'
        stop
      else
        adv%eta1_coords(1:Npts1) = eta1_coords(1:Npts1)
      endif     
    else
      print *,'#Warning, we assume eta1_min = 0._f64 eta1_max = 1._f64'
      delta_eta1 = 1._f64/real(Npts1-1,f64)
      do i=1,Npts1
          adv%eta1_coords(i) = real(i-1,f64)*delta_eta1
      enddo                      
    endif


    if(present(eta2_min).and.present(eta2_max))then
      if(present(eta2_coords))then
        print *,'#provide either eta2_coords or eta2_min and eta2_max'
        print *,'#and not both in subroutine initialize_CSL_2d_advector'
        stop
      else
        delta_eta2 = (eta2_max-eta2_min)/real(Npts2-1,f64)
        do i=1,Npts2
          adv%eta2_coords(i) = eta2_min+real(i-1,f64)*delta_eta2
        enddo
      endif
    else if(present(eta2_coords))then
      if(size(eta2_coords,1)<Npts2)then
        print *,'#bad size for eta2_coords in initialize_CSL_2d_advector'
        stop
      else
        adv%eta2_coords(1:Npts2) = eta2_coords(1:Npts2)
      endif     
    else
      print *,'#Warning, we assume eta2_min = 0._f64 eta2_max = 1._f64'
      delta_eta2 = 1._f64/real(Npts2-1,f64)
      do i=1,Npts2
          adv%eta2_coords(i) = real(i-1,f64)*delta_eta2
      enddo                      
    endif

    adv%nbmax=30
    !print*,f(1,33)
    !return
    SLL_ALLOCATE(adv%tt(adv%nbmax,2),ierr)
    SLL_ALLOCATE(adv%cell(2,adv%nbmax,4),ierr)
    SLL_ALLOCATE(adv%tcell(adv%nbmax,4),ierr)
    SLL_ALLOCATE(adv%intx(0:adv%nbmax),ierr)
    SLL_ALLOCATE(adv%inty(0:adv%nbmax),ierr)
    SLL_ALLOCATE(adv%dir(adv%nbmax,2),ierr)
      
  end subroutine initialize_CSL_2d_advector

  subroutine CSL_advect_2d(&
    adv, &
    A1, &
    A2, &
    dt, &
    input, &
    output)
    class(CSL_2d_advector) :: adv
    sll_real64, dimension(:,:), intent(in) :: A1
    sll_real64, dimension(:,:), intent(in) :: A2
    sll_real64, intent(in) :: dt 
    sll_real64, dimension(:,:), intent(in) :: input
    sll_real64, dimension(:,:), intent(out) :: output      
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: Nc_x1
    sll_int32 :: Nc_x2
    sll_real64 :: x1_min
    sll_real64 :: x1_max
    sll_real64 :: x2_min
    sll_real64 :: x2_max
    sll_real64 :: xx(4)
    sll_real64 :: yy(4)
    sll_int32 :: ii(4)
    sll_int32 :: jj(4)
    sll_int32 :: imin
    sll_int32 :: jmin
    sll_int32 :: imax
    sll_int32 :: jmax
    sll_int32 :: ell
    sll_int32 :: ell1
    sll_real64 :: xA
    sll_real64 :: xB
    sll_real64 :: yA
    sll_real64 :: yB
    sll_int32 :: i0
    sll_int32 :: j0
    sll_int32 :: i1
    sll_int32 :: j1
    sll_int32 :: s
    sll_int32 :: nbx(4)
    sll_int32 :: nby(4)
    sll_int32 :: dirx
    sll_int32 :: diry
    sll_int32 :: k
    sll_int32 :: sx
    sll_int32 :: sy
    sll_real64 :: xA_loc
    sll_real64 :: yA_loc
    sll_real64 :: xB_loc
    sll_real64 :: yB_loc
    sll_int32 :: i0_loc
    sll_int32 :: j0_loc
    sll_real64 :: res
    sll_real64 :: minfl
    sll_real64 :: maxfl
    
    
    
    Nc_x1 = adv%Npts1-1 
    Nc_x2 = adv%Npts2-1 
    
    x1_min = adv%eta1_coords(1)
    x1_max = adv%eta1_coords(adv%Npts1)
    x2_min = adv%eta2_coords(1)
    x2_max = adv%eta2_coords(adv%Npts2)
    
    
    call adv%charac%compute_characteristics( &
      A1, &
      A2, &
      dt, &
      adv%eta1_coords, &
      adv%eta2_coords, &
      adv%charac_feet1, &
      adv%charac_feet2)
   
    do j=0,Nc_x2
      do i=0,Nc_x1
        
        output(i+1,j+1)=0._f64
        xx(1)=(xx(1)-x1_min)/(x1_max-x1_min)*real(Nc_x1,f64)
        xx(2)=(xx(2)-x1_min)/(x1_max-x1_min)*real(Nc_x1,f64)
        xx(3)=(xx(3)-x1_min)/(x1_max-x1_min)*real(Nc_x1,f64)
        xx(4)=(xx(4)-x1_min)/(x1_max-x1_min)*real(Nc_x1,f64)

        yy(1)=(yy(1)-x2_min)/(x2_max-x2_min)*real(Nc_x2,f64)
        yy(2)=(yy(2)-x2_min)/(x2_max-x2_min)*real(Nc_x2,f64)
        yy(3)=(yy(3)-x2_min)/(x2_max-x2_min)*real(Nc_x2,f64)
        yy(4)=(yy(4)-x2_min)/(x2_max-x2_min)*real(Nc_x2,f64)
    

        xx=xx+0.5_f64
        yy=yy+0.5_f64 
	        
        ii(1)=floor(xx(1))
        ii(2)=floor(xx(2))
        ii(3)=floor(xx(3))
        ii(4)=floor(xx(4))
	
        jj(1)=floor(yy(1))
        jj(2)=floor(yy(2))
        jj(3)=floor(yy(3))
        jj(4)=floor(yy(4))
	
        imin=min(ii(1),ii(2),ii(3),ii(4))
        jmin=min(jj(1),jj(2),jj(3),jj(4))
        imax=max(ii(1),ii(2),ii(3),ii(4))
        jmax=max(jj(1),jj(2),jj(3),jj(4))



        do ell=1,4

          ell1=ell+1
          if(ell1==5)then
            ell1=1  
          endif
          xA=xx(ell)
          yA=yy(ell)
          xB=xx(ell1)
          yB=yy(ell1)
          i0=ii(ell)
          j0=jj(ell)
          i1=ii(ell1)
          j1=jj(ell1)

          s=1
          if(i0<i1)then
            do k=i0+1,i1
              adv%tt(s,1)=(real(k,f64)-xA)/(xB-xA)
              s=s+1
            enddo
            dirx=1
          endif
          if(i0>i1)then
            do k=i0,i1+1,-1
              adv%tt(s,1)=(real(k,f64)-xA)/(xB-xA)
              s=s+1
            enddo
            dirx=-1
          endif
          nbx(ell)=s-1
          s=1
          if(j0<j1)then
            do k=j0+1,j1
              adv%tt(s,2)=(real(k,f64)-yA)/(yB-yA)
              s=s+1
            enddo
            diry=1
          endif
          if(j0>j1)then
            do k=j0,j1+1,-1
              adv%tt(s,2)=(real(k,f64)-yA)/(yB-yA)
              s=s+1
            enddo
            diry=-1
          endif
          nby(ell)=s-1
    
          adv%cell(1,1,ell)=i0
          adv%cell(2,1,ell)=j0
          adv%tcell(1,ell)=0._f64
          sx=1
          sy=1
          s=1

          do while((sx<=nbx(ell)).and.(sy<=nby(ell)))
            if(adv%tt(sx,1)<adv%tt(sy,2))then
              s=s+1
              adv%cell(1,s,ell)=adv%cell(1,s-1,ell)+dirx
              adv%cell(2,s,ell)=adv%cell(2,s-1,ell)
              adv%tcell(s,ell)=adv%tt(sx,1)
              adv%intx(2*adv%cell(1,s-1,ell)+(dirx-1)/2-2*imin)=yA+adv%tt(sx,1)*(yB-yA)
              adv%dir(s,1)=1
              adv%dir(s,2)=dirx  
              sx=sx+1
            else
              s=s+1
              adv%cell(1,s,ell)=adv%cell(1,s-1,ell)
              adv%cell(2,s,ell)=adv%cell(2,s-1,ell)+diry
              adv%tcell(s,ell)=adv%tt(sy,2)
              adv%inty(2*adv%cell(2,s-1,ell)+(diry-1)/2-2*jmin)=xA+adv%tt(sy,2)*(xB-xA)
              adv%dir(s,1)=2
              adv%dir(s,2)=diry
              sy=sy+1
            endif
          enddo
          do while(sx<=nbx(ell))
            s=s+1
            adv%cell(1,s,ell)=adv%cell(1,s-1,ell)+dirx
            adv%cell(2,s,ell)=adv%cell(2,s-1,ell)
            adv%tcell(s,ell)=adv%tt(sx,1)
            adv%intx(2*adv%cell(1,s-1,ell)+(dirx-1)/2-2*imin)=yA+adv%tt(sx,1)*(yB-yA)
            adv%dir(s,1)=1
            adv%dir(s,2)=dirx
            sx=sx+1
          enddo  
          do while(sy<=nby(ell))
            s=s+1
            adv%cell(1,s,ell)=adv%cell(1,s-1,ell)
            adv%cell(2,s,ell)=adv%cell(2,s-1,ell)+diry
            adv%tcell(s,ell)=adv%tt(sy,2)
            adv%inty(2*adv%cell(2,s-1,ell)+(diry-1)/2-2*jmin)=xA+adv%tt(sy,2)*(xB-xA)
            adv%dir(s,1)=2
            adv%dir(s,2)=diry
            sy=sy+1
          enddo        

!Computation of extern edges


        if ((ell==1) .or. (ell==4) .or. &
          ((ell==2) .and. (i==Nc_x1)) .or. &
          ((ell==3) .and. (j==Nc_x2))) then
          xB_loc=xx(ell)
          yB_loc=yy(ell)
          do k=2,nbx(ell)+nby(ell)+2
            xA_loc=xB_loc
            yA_loc=yB_loc
            i0_loc=adv%cell(1,k-1,ell)
            j0_loc=adv%cell(2,k-1,ell)
            if(k==nbx(ell)+nby(ell)+2)then
              xB_loc=xx(ell1)
              yB_loc=yy(ell1)
            else
              if(adv%dir(k,1)==1)then
                xB_loc=real(adv%cell(1,k,ell)+(1-adv%dir(k,2))/2,f64)
                yB_loc=yy(ell)+adv%tcell(k,ell)*(yy(ell1)-yy(ell))
              else
                xB_loc=xx(ell)+adv%tcell(k,ell)*(xx(ell1)-xx(ell))
                yB_loc=real(adv%cell(2,k,ell)+(1-adv%dir(k,2))/2,f64)
              endif
            endif

            res = 0._f64
            !call calcule_coeff(N0,N1,buf2d,i0_loc,j0_loc,xA_loc,yA_loc,xB_loc,yB_loc,res,aretesh,aretesv,sommets,areteshb,areteshh,&
            !  aretesvg,aretesvd,sommetsbg,sommetsbd,sommetshg,sommetshd,interp_case)

            output(i+1,j+1)=output(i+1,j+1)+res

            if ((i-1>=0) .and. (ell==4)) then
              output(i,j+1)=output(i,j+1)-res
            endif
            if ((j-1>=0) .and. (ell==1)) then
              output(i+1,j)=output(i+1,j)-res
            endif
          enddo
        endif
      end do

!#endif
!
!Computation of vertical intern edges
!
!#ifdef COMPUTE_V
    do ell=0,imax-imin-1
    minfl=min(floor(adv%intx(2*ell)),floor(adv%intx(2*ell+1)))
    maxfl=max(floor(adv%intx(2*ell)),floor(adv%intx(2*ell+1)))

        i0_loc=imin+ell
        yB_loc=min(adv%intx(2*ell),adv%intx(2*ell+1))
        do k=0,maxfl-minfl
          yA_loc=yB_loc
          j0_loc=minfl+k  
          if(k==maxfl-minfl)then
            yB_loc=max(adv%intx(2*ell),adv%intx(2*ell+1))
          else
            yB_loc=real(minfl+k+1,f64)
          endif

          !call calcule_coeffv(N0,N1,buf2d,i0_loc,j0_loc,yA_loc,yB_loc,res,aretesh,aretesv,sommets,areteshb,areteshh,&
        !aretesvg,aretesvd,sommetsbg,sommetsbd,sommetshg,sommetshd,interp_case)
          output(i+1,j+1)=output(i+1,j+1)+res
      enddo
    enddo
!
!#endif
!
!Computation of horizontal intern edges
!
!#ifdef COMPUTE_X
    do ell=0,jmax-jmin-1
      minfl=min(floor(adv%inty(2*ell)),floor(adv%inty(2*ell+1)))
      maxfl=max(floor(adv%inty(2*ell)),floor(adv%inty(2*ell+1)))

      j0_loc=jmin+ell
      xA_loc=min(adv%inty(2*ell),adv%inty(2*ell+1))
       do k=0,maxfl-minfl
         xB_loc=xA_loc
         i0_loc=minfl+k  
         if(k==maxfl-minfl)then
           xA_loc=max(adv%inty(2*ell),adv%inty(2*ell+1))
         else
           xA_loc=real(minfl+k+1,f64)
         endif
          
         !call calcule_coeffh(N0,N1,buf2d,i0_loc,j0_loc,xA_loc,xB_loc,res,aretesh,aretesv,sommets,areteshb,areteshh,&
      !aretesvg,aretesvd,sommetsbg,sommetsbd,sommetshg,sommetshd,interp_case)
      output(i+1,j+1)=output(i+1,j+1)+res

      enddo
    enddo

!#endif
!        deallocate(tnbr)
!        deallocate(tpts)

!      enddo
!
!    enddo
   
!#ifdef DIAG_TIME
!      f=buf2d
!#endif
!    deallocate(tt,cell,tcell,intx,inty,dir)
!
!    if ((interp_case==3) .or. (interp_case==6) .or. (interp_case==7)) then
!      deallocate(aretesh,aretesv,sommets)
!      deallocate(areteshb,areteshh,aretesvg,aretesvd,sommetsbg,sommetsbd,sommetshg,sommetshd)
!    endif
!    if (interp_case==4) then
!      deallocate(areteshb,areteshh,aretesvg,aretesvd,sommetsbg,sommetsbd,sommetshg,sommetshd)
!    endif















        
        
      enddo
    enddo
   
        

          
  end subroutine CSL_advect_2d



  subroutine aux(N0,N1,f,aretesh,aretesv,sommets,ordre,dom) !used in PPM case

    sll_int32,intent(in)::N0,N1,ordre
    sll_real64,dimension(0:1,0:1),intent(in)::dom
    real(f64),dimension(0:N0,0:N1-1),intent(in)::f
    !    sll_real64,dimension(2,-1:N0+1,-1:N1+1),intent(in)::carac
    sll_int32::i,j,im3,im2,im1,ib,ip1,ip2,jm3,jm2,jm1,jb,jp1,jp2
    real(f64),dimension(0:N0,0:N1-1),intent(inout)::aretesh,aretesv,sommets
!		print*,N0,N1,ordre
!stop
    
    if(dom(1,1)==1e-23)then
      print *,'#is this possible?'
    endif


    do i=0,N0
      do j=0,N1-1
        im3=i-3
        im2=i-2
        im1=i-1
        ip1=i+1
        ip2=i+2
        if(i-3<=0)im3=0;
        if(i-2<=0)im2=0;
        if(i-1<=0)im1=0;
        ib=i;
        if(i+1>=N0)ip1=N0;
        if(i+2>=N0)ip2=N0;
        jm3=modulo(j-3,N1)
        jm2=modulo(j-2,N1)
        jm1=modulo(j-1,N1)
        jb=modulo(j,N1)
        jp1=modulo(j+1,N1)
        jp2=modulo(j+2,N1)
        if (ordre==1) then !PPM1
          aretesv(i,j)=7._f64/12._f64*(f(im1,jb)+f(ib,jb)) &
            -1._f64/12._f64*(f(im2,jb)+f(ip1,jb))
          aretesh(i,j)=7._f64/12._f64*(f(ib,jm1)+f(ib,jb)) &	
            -1._f64/12._f64*(f(ib,jm2)+f(ib,jp1))
        else if (ordre==2) then !PPM2
          aretesv(i,j)=1._f64/60._f64*(f(ip2,jb)+f(im3,jb)) &
            -8._f64/60._f64*(f(ip1,jb)+f(im2,jb)) &
            +37._f64/60._f64*(f(ib,jb)+f(im1,jb))
          aretesh(i,j)=1._f64/60._f64*(f(ib,jp2)+f(ib,jm3)) &
            -8._f64/60._f64*(f(ib,jp1)+f(ib,jm2)) &
            +37._f64/60._f64*(f(ib,jb)+f(ib,jm1))
        else if (ordre==0) then !PPM0
          aretesv(i,j)=1._f64/2._f64*(f(ib,jb)+f(im1,jb))
          aretesh(i,j)=1._f64/2._f64*(f(ib,jb)+f(ib,jm1))
        end if
      end do
    end do
    
    do i=0,N0
      do j=0,N1-1
        im3=i-3
        im2=i-2
        im1=i-1
        ip1=i+1
        ip2=i+2
        if(i-3<=0)im3=0;
        if(i-2<=0)im2=0;
        if(i-1<=0)im1=0;
        ib=i;
        if(i+1>=N0)ip1=N0;
        if(i+2>=N0)ip2=N0;
				
				
				
        jm3=modulo(j-3,N1)
        jm2=modulo(j-2,N1)
        jm1=modulo(j-1,N1)
        jb=modulo(j,N1)
        jp1=modulo(j+1,N1)
        jp2=modulo(j+2,N1)

        if (ordre==1) then !PPM1
          sommets(i,j)=7._f64/12._f64*(aretesv(ib,jm1)+aretesv(ib,jb)) &
            -1._f64/12._f64*(aretesv(ib,jm2)+aretesv(ib,jp1))
        else if (ordre==2) then !PPM2
          sommets(i,j)=1._f64/60._f64*(aretesv(ib,jp2)+aretesv(ib,jm3)) &
            -8._f64/60._f64*(aretesv(ib,jp1)+aretesv(ib,jm2)) &
            +37._f64/60._f64*(aretesv(ib,jb)+aretesv(ib,jm1))
        else if (ordre==0) then !PPM0
          sommets(i,j)=1._f64/2._f64*(aretesv(ib,jb)+aretesv(ib,jm1))
        end if
      end do
    end do

  end subroutine aux

  subroutine compute_aux2_new( &
    N_x1, &
    N_x2, &
    r_x1, &
    s_x1, &
    r_x2, &
    s_x2, &
    input, &
    output)
    sll_int32, intent(in) :: N_x1
    sll_int32, intent(in) :: N_x2
    sll_int32, intent(in) :: r_x1
    sll_int32, intent(in) :: s_x1
    sll_int32, intent(in) :: r_x2
    sll_int32, intent(in) :: s_x2
    sll_real64, dimension(:,:), intent(in) :: input
    sll_real64, dimension(:,:,:,:), intent(out) :: output
    !sll_real64, dimension(:), allocatable :: w_x1
    sll_real64, dimension(:), allocatable :: ww_x1
    !sll_real64, dimension(:), allocatable :: w_x2
    sll_real64, dimension(:), allocatable :: ww_x2
    sll_int32 :: ierr
    
    !SLL_ALLOCATE(w_x1(r_x1,s_x1),ierr)
    !SLL_ALLOCATE(w_x2(r_x2,s_x2),ierr)
    SLL_ALLOCATE(ww_x1(r_x1:s_x1-1),ierr)
    SLL_ALLOCATE(ww_x2(r_x2:s_x2-1),ierr)
    
    call compute_ww(ww_x1,r_x1,s_x1)
    call compute_ww(ww_x2,r_x2,s_x2)
    
  end subroutine compute_aux2_new
  
  
  subroutine compute_ww(ww,r,s)
    sll_real64, dimension(:), intent(out) :: ww
    sll_int32, intent(in) :: r
    sll_int32, intent(in) :: s
    sll_real64, dimension(:), allocatable :: w
    sll_real64 :: tmp
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: ierr

    SLL_ALLOCATE(w(r:s),ierr)

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

    tmp=0._f64
    do i=r,-1
      tmp=tmp+w(i)
      ww(i)=-tmp
    enddo
    tmp=0._f64
    do i=s,1,-1
      tmp=tmp+w(i)
      ww(i-1)=tmp
    enddo

    SLL_DEALLOCATE_ARRAY(w,ierr)
    
  end subroutine compute_ww

  subroutine aux2(N0,N1,f,areteshb,areteshh,aretesvg,aretesvd,sommetsbg,sommetsbd, &
    sommetshg,sommetshd,ordre,carac,dom) !used in PPM case
    integer,intent(in)::N0,N1,ordre
    real(f64),dimension(0:1,0:1),intent(in)::dom
    real(f64),dimension(0:N0,0:N1-1),intent(in)::f
    real(f64),dimension(2,-1:N0+1,-1:N1+1),intent(in)::carac
    integer::i,j,i3
    real(f64),dimension(0:N0,0:N1-1),intent(inout)::areteshb,areteshh,aretesvg,aretesvd,&
      sommetsbg,sommetsbd,sommetshg,sommetshd
    real(f64) ::w(-ordre:ordre+1),tmp,ww(-ordre:ordre)
    integer::r,s,ii,d
    !integer :: ib,im3,im2,im1,ip1,ip2,jm3,jm2,jm1,jb,jp1,jp2
       
    d=ordre
    r=-d
    s=d+1
    if(carac(1,1,1)==1e-23)then
      print *,'#is this possible?'
    endif
    if(dom(1,1)==1e-23)then
      print *,'#is this possible?'
    endif
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

!    do i=r,-1
!      tmp=1._f64
!      !do j=r,i-1
!      !  tmp=tmp*real(i-j,f64)
!      !enddo
!      !do j=i+1,s
!      !  tmp=tmp*real(i-j,f64)
!      !enddo
!      !tmp=1._f64/tmp
!      do j=r,i-1 !-j/(i-j)=j/(j-i)=1/(1-i/j)
!        tmp=tmp*(1._f64-real(i,f64)/real(j,f64))
!      enddo
!      do j=i+1,-1
!        tmp=tmp*(1._f64-real(i,f64)/real(j,f64))
!      enddo
!      do j=1,s
!        tmp=tmp*(1._f64-real(i,f64)/real(j,f64))
!      enddo
!      tmp=tmp*real(i,f64)
!      w(i)=1._f64/tmp      
!    enddo


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
    
    !compute now ww
    !maple code
    !#for conservative formulation
    !tmp:=0:
    !for k from r to -1 do
    !tmp:=tmp+C[k]:
    !CC[k]:=-tmp:
    !od:
    !tmp:=0:
    !for k from s to 1 by -1 do
    !  tmp:=tmp+C[k]:
    !  CC[k-1]:=tmp:
    !od:
    !seq(CC[k],k=r..s-1);
    !evalf(%);

    tmp=0._f64
    do i=r,-1
      tmp=tmp+w(i)
      ww(i)=-tmp
    enddo
    tmp=0._f64
    do i=s,1,-1
      tmp=tmp+w(i)
      ww(i-1)=tmp
    enddo

    !print *,'ww',ww
    !do ii=r,s-1
    !  print *,ii,ww(r+s-1-ii)
    !enddo
    !stop
    do j=0,N1-1
      do i=0,N0
        tmp=0._f64
        do ii=r,s-1
          i3=i+ii-1;if(i3<=0)i3=0;if(i3>=N0)i3=N0          
          tmp=tmp+ww(r+s-1-ii)*f(i3,j)
        enddo
        aretesvg(i,j)=tmp
        tmp=0._f64
        do ii=r,s-1
          i3=i+ii;if(i3<=0)i3=0;if(i3>=N0)i3=N0          
          tmp=tmp+ww(ii)*f(i3,j)
        enddo
        aretesvd(i,j)=tmp
        tmp=0._f64
        do ii=r,s-1
          tmp=tmp+ww(ii)*f(i,modulo(j+ii,N1))
        enddo
        areteshh(i,j)=tmp
        tmp=0._f64
        do ii=r,s-1
          tmp=tmp+ww(r+s-1-ii)*f(i,modulo(j+ii-1,N1))
        enddo
        areteshb(i,j)=tmp
      enddo
    enddo

    do j=0,N1-1
      do i=0,N0
        tmp=0._f64
        do ii=r,s-1
          i3=i+ii;if(i3<=0)i3=0;if(i3>=N0)i3=N0          
          tmp=tmp+ww(ii)*areteshh(i3,j)
        enddo
        sommetshd(i,j)=tmp
        tmp=0._f64
        do ii=r,s-1
          i3=i+ii-1;if(i3<=0)i3=0;if(i3>=N0)i3=N0          
          tmp=tmp+ww(r+s-1-ii)*areteshh(i3,j)
        enddo
        sommetshg(i,j)=tmp
        tmp=0._f64
        do ii=r,s-1
          i3=i+ii;if(i3<=0)i3=0;if(i3>=N0)i3=N0          
          tmp=tmp+ww(ii)*areteshb(i3,j)
        enddo
        sommetsbd(i,j)=tmp
        tmp=0._f64
        do ii=r,s-1
          i3=i+ii-1;if(i3<=0)i3=0;if(i3>=N0)i3=N0          
          tmp=tmp+ww(r+s-1-ii)*areteshb(i3,j)
        enddo
        sommetsbg(i,j)=tmp
      enddo
    enddo

  end subroutine aux2










end module sll_module_advection_2d_CSL
